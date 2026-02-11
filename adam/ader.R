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
# Load ADRS, ADTTE, ADSL, ADLB, ADVS, ADEX, ADPP and ADAE
adrs <- pharmaverseadam::adrs_onco
adtte <- pharmaverseadam::adtte_onco
adsl <- pharmaverseadam::adsl
adlb <- pharmaverseadam::adlb
advs <- pharmaverseadam::advs
adex <- pharmaverseadam::adex %>%
  filter(PARCAT1 == "INDIVIDUAL")
adpp <- pharmaverseadam::adpp
adae <- pharmaverseadam::adae

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
# ---- Add Exposure Metrics  ----
covar_auc <- covar_vslb %>%
  derive_vars_transposed(
    dataset_merge = adpp,
    filter = PARAMCD %in% c("AUCLST", "CMAX"),
    by_vars = get_admiral_option("subject_keys"),
    key_var = PARAMCD,
    value_var = AVAL
  ) %>%
  rename(AUCSS = AUCLST, CMAXSS = CMAX)

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
# Combine covariates with ADER data

ader_prefinal <- ader_aseq %>%
  derive_vars_merged(
    dataset_add = covar_auc,
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

## ----r echo=TRUE, message=FALSE-----------------------------------------------
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx") %>%
  select_dataset("ADEE")

## ----r------------------------------------------------------------------------
# ---- Create adee base dataset

# Get variable names from both datasets
adsl_vars <- names(adsl)
adtte_vars <- names(adtte)

# Find common variables
common_vars <- intersect(adsl_vars, adtte_vars)

# Remove key variables to get variables to drop
vars_to_drop <- setdiff(common_vars, c("STUDYID", "USUBJID"))

# Ensure PARAMN exists in ADTTE
if (!"PARAMN" %in% names(adtte)) {
  adtte <- adtte %>%
    mutate(
      PARAMN = case_when(
        PARAMCD == "PFS" ~ 1,
        PARAMCD == "OS" ~ 2,
        PARAMCD == "TTP" ~ 3,
        PARAMCD == "TTNT" ~ 4,
        TRUE ~ 99
      )
    )
}

# ---- Create ADEE Base

adee_base <- adtte %>%
  # Filter to efficacy endpoints
  filter(PARAMCD %in% c("OS", "PFS", "TTP", "TTNT")) %>%
  # Add derived variables
  mutate(
    EVENT = 1 - CNSR,
    AVALU = if_else(!is.na(AVAL), "DAYS", NA_character_),
  ) %>%
  # Remove overlapping variables (use clean method)
  select(-any_of(vars_to_drop))


## ----r------------------------------------------------------------------------
# ---- Add Analysis variables

adee_aseq <- adee_base %>%
  # Analysis flags
  mutate(
    ANL01FL = if_else(PARAMCD == "PFS", "Y", ""),
    ANL02FL = if_else(PARAMCD == "OS", "Y", ""),
    ANL03FL = if_else(PARAMCD == "TTP", "Y", ""),
    ANL04FL = if_else(PARAMCD == "TTNT", "Y", "")
  ) %>%
  # Parameter categories
  mutate(
    PARCAT1 = "EFFICACY",
    PARCAT2 = "TIME TO EVENT"
  ) %>%
  # Sequence number
  derive_var_obs_number(
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD),
    new_var = ASEQ,
    check_type = "error"
  )


## ----r------------------------------------------------------------------------
# Combine covariates with ADER data

adee_prefinal <- adee_aseq %>%
  derive_vars_merged(
    dataset_add = covar_auc,
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r------------------------------------------------------------------------
adee <- adee_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

adee_xpt <- adee %>%
  xportr_type(metacore, domain = "ADEE") %>% # Coerce variable type to match spec
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adee.xpt")) # Write xpt v5 transport file

## ----r------------------------------------------------------------------------
# [Similar structure...]

## ----r------------------------------------------------------------------------
# [Similar structure...]
