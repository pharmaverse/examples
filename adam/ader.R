## ----r echo=TRUE, message=FALSE-----------------------------------------------
# Load Packages
library(admiral)
library(admiralonco)
# pharmaverseadam contains example datasets generated from the CDISC pilot
# project SDTM ran through admiral templates
library(pharmaverseadam)
library(dplyr)
library(lubridate)
library(stringr)
library(metacore)
library(metatools)
library(xportr)

## ----r echo=TRUE, message=FALSE-----------------------------------------------
# ---- Load Specs for Metacore ----
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx") %>%
  select_dataset("ADER")

## ----r------------------------------------------------------------------------
# ---- Load source datasets ----
# Load ADRS, ADTTE, ADSL, ADLB, ADVS and ADEX
adrs <- pharmaverseadam::adrs_onco
adtte <- pharmaverseadam::adtte_onco

adsl <- pharmaverseadam::adsl
adlb <- pharmaverseadam::adlb
advs <- pharmaverseadam::advs
adex <- pharmaverseadam::adex %>%
  filter(PARCAT1 == "INDIVIDUAL")

## ----r------------------------------------------------------------------------
# ---- Derivations ----
# For ADTTE censor variables add "IND" to PARAMCD
adttei <- adtte %>%
  mutate(PARAMCD = paste0(PARAMCD, "IND"))

ader_tte <- adsl %>%
  select(!!!get_admiral_option("subject_keys")) %>%
  # Create OS and PFS variables from ADTTE
  derive_vars_transposed(
    dataset_merge = adtte,
    by_vars = get_admiral_option("subject_keys"),
    key_var = PARAMCD,
    value_var = AVAL
  ) %>%
  # Create OS and PFS censor variables
  derive_vars_transposed(
    dataset_merge = adttei,
    by_vars = get_admiral_option("subject_keys"),
    key_var = PARAMCD,
    value_var = CNSR
  )

## ----r------------------------------------------------------------------------
# ---- Add ADRS data ----
# Add response date to ADSL for duration of response calculation
ader_bor <- ader_tte %>%
  derive_vars_merged(
    dataset_add = adrs,
    filter_add = PARAMCD == "BOR" & ANL01FL == "Y",
    by_vars = get_admiral_option("subject_keys"),
    new_vars = exprs(BOR = AVAL, BORC = AVALC)
  )

## ----r------------------------------------------------------------------------
# ---- Add Sequence Number ----
ader_aseq <- ader_bor %>%
  derive_var_obs_number(
    by_vars = get_admiral_option("subject_keys"),
    check_type = "error"
  )

## ----r------------------------------------------------------------------------
# ---- Derive Covariates ----
# Include numeric values for STUDYIDN, USUBJIDN, SEXN, RACEN etc.

covar <- adsl %>%
  create_var_from_codelist(metacore, input_var = STUDYID, out_var = STUDYIDN) %>%
  create_var_from_codelist(metacore, input_var = SEX, out_var = SEXN) %>%
  create_var_from_codelist(metacore, input_var = RACE, out_var = RACEN) %>%
  create_var_from_codelist(metacore, input_var = ETHNIC, out_var = ETHNICN) %>%
  create_var_from_codelist(metacore, input_var = ARMCD, out_var = COHORT) %>%
  create_var_from_codelist(metacore, input_var = ARMCD, out_var = COHORTC) %>%
  create_var_from_codelist(metacore, input_var = ARM, out_var = ARMN) %>%
  create_var_from_codelist(metacore, input_var = ACTARM, out_var = ACTARMN) %>%
  create_var_from_codelist(metacore, input_var = COUNTRY, out_var = COUNTRYN) %>%
  create_var_from_codelist(metacore, input_var = COUNTRY, out_var = COUNTRYL) %>%
  mutate(
    STUDYIDN = as.numeric(word(USUBJID, 1, sep = fixed("-"))),
    SITEIDN = as.numeric(word(USUBJID, 2, sep = fixed("-"))),
    USUBJIDN = as.numeric(word(USUBJID, 3, sep = fixed("-"))),
    SUBJIDN = as.numeric(SUBJID),
    ROUTE = unique(adex$EXROUTE)[1],
    FORM = unique(adex$EXDOSFRM)[1],
    REGION1 = COUNTRY,
    REGION1N = COUNTRYN,
  ) %>%
  create_var_from_codelist(metacore, input_var = FORM, out_var = FORMN) %>%
  create_var_from_codelist(metacore, input_var = ROUTE, out_var = ROUTEN)

## ----r------------------------------------------------------------------------
# ---- Derive additional baselines from ADVS and ADLB ----

labsbl <- adlb %>%
  filter(ABLFL == "Y" & PARAMCD %in% c("CREAT", "ALT", "AST", "BILI")) %>%
  mutate(PARAMCDB = paste0(PARAMCD, "BL")) %>%
  select(STUDYID, USUBJID, PARAMCDB, AVAL)

covar_vslb <- covar %>%
  derive_vars_merged(
    dataset_add = advs,
    filter_add = PARAMCD == "HEIGHT" & ABLFL == "Y",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(HTBL = AVAL)
  ) %>%
  derive_vars_merged(
    dataset_add = advs,
    filter_add = PARAMCD == "WEIGHT" & ABLFL == "Y",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WTBL = AVAL)
  ) %>%
  derive_vars_transposed(
    dataset_merge = labsbl,
    by_vars = exprs(STUDYID, USUBJID),
    key_var = PARAMCDB,
    value_var = AVAL
  ) %>%
  mutate(
    BMIBL = compute_bmi(height = HTBL, weight = WTBL),
    BSABL = compute_bsa(
      height = HTBL,
      weight = WTBL,
      method = "Mosteller"
    ),
    CRCLBL = compute_egfr(
      creat = CREATBL, creatu = "SI", age = AGE, weight = WTBL, sex = SEX,
      method = "CRCL"
    ),
    EGFRBL = compute_egfr(
      creat = CREATBL, creatu = "SI", age = AGE, weight = WTBL, sex = SEX,
      method = "CKD-EPI"
    )
  ) %>%
  rename(TBILBL = BILIBL)

## ----r------------------------------------------------------------------------
# Combine covariates with APPPK data
# Combine covariates with ADER data

ader_prefinal <- ader_aseq %>%
  derive_vars_merged(
    dataset_add = covar_vslb,
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r------------------------------------------------------------------------
ader <- ader_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

ader_xpt <- ader %>%
  xportr_type(metacore, domain = "ADER") %>% # Coerce variable type to match spec
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "ader.xpt")) # Write xpt v5 transport file
