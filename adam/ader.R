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
library(purrr)
library(rlang)
library(metacore)
library(metatools)
library(xportr)

## ----r echo=TRUE, message=FALSE-----------------------------------------------
# ---- Load Specs for Metacore ----
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx", quiet = TRUE) %>%
  select_dataset("ADER")

## ----r------------------------------------------------------------------------
# ---- Load source datasets ----
# Load ADRS, ADTTE, ADSL, ADLB, ADVS, ADEX, ADPP, ADAE and ADTR
adrs <- pharmaverseadam::adrs_onco
adtte <- pharmaverseadam::adtte_onco
adsl <- pharmaverseadam::adsl %>%
  filter(TRT01A != "Screen Failure")
adlb <- pharmaverseadam::adlb
advs <- pharmaverseadam::advs
adex <- pharmaverseadam::adex %>%
  filter(PARCAT1 == "INDIVIDUAL")
adpp <- pharmaverseadam::adpp
adae <- pharmaverseadam::adae
adtr <- pharmaverseadam::adtr_onco

## ----r------------------------------------------------------------------------
# ---- Derive Covariates ----
# Include numeric values for STUDYIDN, USUBJIDN, SEXN, RACEN etc.

covar <- purrr::reduce(
  list(
    c("STUDYID", "STUDYIDN"),
    c("SEX", "SEXN"),
    c("RACE", "RACEN"),
    c("ETHNIC", "ETHNICN"),
    c("ARMCD", "COHORT"),
    c("ARMCD", "COHORTC"),
    c("ARM", "ARMN"),
    c("ACTARM", "ACTARMN"),
    c("TRT01A", "TRT01AN"),
    c("TRT01P", "TRT01PN"),
    c("COUNTRY", "COUNTRYN"),
    c("COUNTRY", "COUNTRYL")
  ),
  ~ create_var_from_codelist(.x, metacore,
    input_var = !!rlang::sym(.y[1]),
    out_var   = !!rlang::sym(.y[2])
  ),
  .init = adsl
) %>%
  mutate(
    STUDYIDN = as.numeric(word(USUBJID, 1, sep = fixed("-"))),
    SITEIDN = as.numeric(word(USUBJID, 2, sep = fixed("-"))),
    USUBJIDN = as.numeric(word(USUBJID, 3, sep = fixed("-"))),
    SUBJIDN = as.numeric(SUBJID),
    ROUTE = unique(adex$EXROUTE)[1],
    FORM = unique(adex$EXDOSFRM)[1],
    REGION1 = COUNTRY,
    REGION1N = COUNTRYN
  ) %>%
  create_var_from_codelist(metacore, input_var = FORM, out_var = FORMN) %>%
  create_var_from_codelist(metacore, input_var = ROUTE, out_var = ROUTEN)

## ----r------------------------------------------------------------------------
# ---- Derive additional baselines from ADVS and ADLB ----

labsbl <- adlb %>%
  filter(ABLFL == "Y" & PARAMCD %in% c("CREAT", "ALT", "AST", "BILI")) %>%
  mutate(PARAMCDB = paste0(PARAMCD, "BL")) %>%
  select(!!!get_admiral_option("subject_keys"), PARAMCDB, AVAL)

covar_vslb <- covar %>%
  derive_vars_merged(
    dataset_add = advs,
    filter_add = PARAMCD == "HEIGHT" & ABLFL == "Y",
    by_vars = get_admiral_option("subject_keys"),
    new_vars = exprs(HTBL = AVAL)
  ) %>%
  derive_vars_merged(
    dataset_add = advs,
    filter_add = PARAMCD == "WEIGHT" & ABLFL == "Y",
    by_vars = get_admiral_option("subject_keys"),
    new_vars = exprs(WTBL = AVAL)
  ) %>%
  derive_vars_transposed(
    dataset_merge = labsbl,
    by_vars = get_admiral_option("subject_keys"),
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
# Add response date to ADER for duration of response calculation
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
    by_vars = get_admiral_option("subject_keys")
  )

## ----r------------------------------------------------------------------------
ader <- ader_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the specs
  sort_by_key(metacore) # Sorts the rows by the sort keys

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

ader_xpt <- ader %>%
  xportr_type(metacore, domain = "ADER") %>% # Coerce variable type to match specs
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "ader.xpt")) # Write xpt v5 transport file

## ----r echo=TRUE, message=FALSE-----------------------------------------------
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx", quiet = TRUE) %>%
  select_dataset("ADEE")

## ----r------------------------------------------------------------------------
# ---- Create ADEE base dataset

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
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(PARAMCD),
    new_var = ASEQ,
    check_type = "error"
  )

## ----r------------------------------------------------------------------------
# Combine covariates with ADEE data

adee_prefinal <- adee_aseq %>%
  derive_vars_merged(
    dataset_add = covar_auc,
    by_vars = get_admiral_option("subject_keys")
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
  xportr_type(metacore, domain = "ADEE") %>% # Coerce variable type to match specs
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adee.xpt")) # Write xpt v5 transport file

## ----r echo=TRUE, message=FALSE-----------------------------------------------
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx", quiet = TRUE) %>%
  select_dataset("ADES")

## ----r------------------------------------------------------------------------
# ---- Create ades base dataset

# Derive subject-level summary parameters from ADAE using derive_param_exist_flag()
# NOTE: pharmaverseadam uses ASEV/ASEVN (severity) not AETOXGR/AETOXGRN

# Get all subjects from ADSL
adsl_sub <- adsl %>%
  select(!!!get_admiral_option("subject_keys"))

# Derive binary subject-level AE flag parameters
subject_params <- derive_param_exist_flag(
  dataset_ref = adsl_sub,
  dataset_add = adae,
  condition = TRTEMFL == "Y",
  false_value = "N",
  missing_value = "N",
  set_values_to = exprs(
    AVAL = yn_to_numeric(AVALC),
    PARAMCD = "TEAE",
    PARAM = "Any Treatment-Emergent Adverse Event",
    PARAMN = 1
  )
) %>%
  derive_param_exist_flag(
    dataset_ref = adsl_sub,
    dataset_add = adae,
    condition = TRTEMFL == "Y" & ASEVN == 3,
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "TEAESEV",
      PARAM = "Any Severe Treatment-Emergent Adverse Event",
      PARAMN = 2
    )
  ) %>%
  derive_param_exist_flag(
    dataset_ref = adsl_sub,
    dataset_add = adae,
    condition = AESER == "Y",
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "SAE",
      PARAM = "Any Serious Adverse Event",
      PARAMN = 3
    )
  ) %>%
  derive_param_exist_flag(
    dataset_ref = adsl_sub,
    dataset_add = adae,
    condition = AREL %in% c("POSSIBLE", "PROBABLE", "RELATED"),
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "TRAE",
      PARAM = "Any Treatment-Related Adverse Event",
      PARAMN = 4
    )
  ) %>%
  derive_param_exist_flag(
    dataset_ref = adsl_sub,
    dataset_add = adae,
    condition = AEACN == "DRUG WITHDRAWN",
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "AEDCN",
      PARAM = "AE Leading to Treatment Discontinuation",
      PARAMN = 5
    )
  ) %>%
  arrange(USUBJID, PARAMCD)

## ----r------------------------------------------------------------------------
# ---- Create event level parameters

# Get variable names for clean dropping
adsl_vars <- names(covar_auc)
adae_vars <- names(adae)
common_vars <- intersect(adsl_vars, adae_vars)
vars_to_drop <- setdiff(common_vars, c("STUDYID", "USUBJID"))

# Create event-level records from ADAE
# NOTE: Using actual pharmaverseadam variables (ASEV/ASEVN, AREL)
event_params <- adae %>%
  filter(TRTEMFL == "Y") %>% # Treatment-emergent only
  mutate(
    PARAMCD = "AETERM",
    PARAM = "Adverse Event Term",
    PARAMN = 10,
    AVAL = 1, # Event occurred
    AVALC = "Y",

    # Keep AE-specific variables (8-char names)
    # Using actual pharmaverseadam ADAE variables
    AEDECOD = AEDECOD, # Preferred term
    AEBODSYS = AEBODSYS, # System organ class
    AESEV = ASEV, # Severity (char): MILD, MODERATE, SEVERE
    AESEVN = ASEVN, # Severity (num): 1, 2, 3
    AESER = AESER, # Serious flag: Y/N
    AEREL = AREL, # Relationship (char): NOT RELATED, POSSIBLE, etc.

    # Create numeric relationship for analysis
    AERELN = case_when(
      AREL == "NOT RELATED" ~ 0,
      AREL == "UNLIKELY RELATED" ~ 1,
      AREL == "POSSIBLE" ~ 2,
      AREL == "PROBABLE" ~ 3,
      AREL == "RELATED" ~ 4,
      TRUE ~ NA_real_
    ),
    AESTDT = ASTDT, # AE start date (8 chars)
    AEENDT = AENDT # AE end date (8 chars)
  ) %>%
  select(-any_of(vars_to_drop))

# ---- Combine subject and event levels

# Ensure all AE-specific variables exist in subject_params (as NA)
# This prevents issues when binding with event_params
subject_params_complete <- subject_params %>%
  derive_vars_merged(
    dataset_add = adsl_sub,
    by_vars = get_admiral_option("subject_keys"),
  ) %>%
  mutate(
    # Add event-level variables as NA for subject-level records
    AESTDT = as.Date(NA),
    AEENDT = as.Date(NA),
    AEDECOD = NA_character_,
    AEBODSYS = NA_character_,
    AESEV = NA_character_,
    AESEVN = NA_integer_,
    AESER = NA_character_,
    AEREL = NA_character_,
    AERELN = NA_real_
  )

event_params_complete <- event_params %>%
  derive_vars_merged(
    dataset_add = adsl_sub,
    by_vars = get_admiral_option("subject_keys")
  )

# Combine both levels
ades_base <- bind_rows(
  subject_params_complete,
  event_params_complete
) %>%
  arrange(USUBJID, PARAMCD)

## ----r------------------------------------------------------------------------
# ---- Add analysis variables

ades_flags <- ades_base %>%
  # Analysis flags
  mutate(
    ANL01FL = if_else(PARAMCD == "TEAE", "Y", ""),
    ANL02FL = if_else(PARAMCD == "TEAESEV", "Y", ""),
    ANL03FL = if_else(PARAMCD == "SAE", "Y", ""),
    ANL04FL = if_else(PARAMCD == "TRAE", "Y", ""),
    ANL05FL = if_else(PARAMCD == "AETERM", "Y", "")
  ) %>%
  # Parameter categories
  mutate(
    PARCAT1 = "SAFETY",
    PARCAT2 = case_when(
      PARAMN <= 5 ~ "SUBJECT-LEVEL",
      PARAMN >= 10 ~ "EVENT-LEVEL",
      TRUE ~ NA_character_
    )
  ) %>%
  # Analysis timepoint
  mutate(
    AVISIT = if_else(PARAMN <= 5, "OVERALL", "AT EVENT"),
    AVISITN = if_else(PARAMN <= 5, 99, 0)
  ) %>%
  # Sort and create sequence number
  # Use coalesce to handle NA dates (puts them first in sort)
  arrange(USUBJID, PARAMN, coalesce(AESTDT, as.Date("1900-01-01"))) %>%
  group_by(!!!get_admiral_option("subject_keys")) %>%
  mutate(ASEQ = row_number()) %>%
  ungroup()

## ----r------------------------------------------------------------------------
# Combine covariates with ADES data

ades_prefinal <- ades_flags %>%
  derive_vars_merged(
    dataset_add = covar_auc,
    by_vars = get_admiral_option("subject_keys")
  )

## ----r------------------------------------------------------------------------
## Check Data With metacore and metatools

ades <- ades_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

ades_xpt <- ades %>%
  xportr_type(metacore, domain = "ADES") %>% # Coerce variable type to match specs
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore, domain = "ADES") %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "ades.xpt")) # Write xpt v5 transport file

## ----r echo=TRUE, message=FALSE-----------------------------------------------
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx", quiet = TRUE) %>%
  select_dataset("ADTRR")

## ----r------------------------------------------------------------------------
# ---- Create base tumor size parameter dataset

# Get variable names for clean dropping
adsl_vars <- names(covar_auc)
adtr_vars <- names(adtr)
common_vars <- intersect(adsl_vars, adtr_vars)
vars_to_drop <- setdiff(common_vars, c("STUDYID", "USUBJID"))

tsize_final <- adtr %>%
  filter(PARAMCD == "SDIAM") %>%
  mutate(
    PARAMCD = "TSIZE",
    PARAM = "Target Lesions Sum of Diameters",
    PARAMN = 1
  ) %>%
  # Derive Nominal Relative Time from First Dose (NFRLT)
  derive_var_nfrlt(
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "DAYS",
    visit_day = ADY
  ) %>%
  # Derive Actual Relative Time from First Dose (AFRLT)
  derive_vars_duration(
    new_var = AFRLT,
    start_date = TRTSDT,
    end_date = ADT,
    out_unit = "DAYS",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  select(-any_of(vars_to_drop))

## ----r------------------------------------------------------------------------
# ---- Add BOR from ADRS

adrs_vars <- names(adrs)
common_vars_adrs <- intersect(adsl_vars, adrs_vars)
vars_to_drop_adrs <- setdiff(common_vars_adrs, c("STUDYID", "USUBJID"))

bor <- adrs %>%
  filter(PARAMCD == "BOR" & SAFFL == "Y") %>%
  mutate(
    PARAMN = 2,
    # Create BORN from AVALC if AVAL doesn't exist
    BORN = if ("AVAL" %in% names(.)) {
      AVAL
    } else {
      case_when(
        AVALC == "CR" ~ 4,
        AVALC == "PR" ~ 3,
        AVALC == "SD" ~ 2,
        AVALC == "PD" ~ 1,
        AVALC == "NE" ~ 0,
        TRUE ~ NA_real_
      )
    }
  ) %>%
  select(-any_of(vars_to_drop_adrs))

## ----r------------------------------------------------------------------------
# ---- Derive Nadir

# Calculate nadir from TSIZE records
# Keep BASE, CHG, PCHG from the nadir timepoint
nadir <- tsize_final %>%
  filter(AVISITN > 0 & !is.na(AVAL)) %>%
  group_by(!!!get_admiral_option("subject_keys")) %>%
  filter(AVAL == min(AVAL, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    PARAMCD = "NADIR",
    PARAM = "Nadir Tumor Size",
    PARAMN = 3,
    NADIR = AVAL,
    NADPCHG = PCHG, # Keep PCHG at nadir
    NADVST = AVISIT # Keep visit of nadir
  )

## ----r------------------------------------------------------------------------
# ---- Combine parameters

adtrr_base <- bind_rows(
  tsize_final,
  bor,
  nadir
) %>%
  arrange(USUBJID, PARAMN, AVISITN)

# ---- Add analysis variables

# Ensure AVALU exists before mutating
if (!"AVALU" %in% names(adtrr_base)) {
  adtrr_base <- adtrr_base %>%
    mutate(AVALU = NA_character_)
}

adtrr_seq <- adtrr_base %>%
  # Analysis flags
  mutate(
    # Baseline flag
    ABLFL = case_when(
      !is.na(ABLFL) ~ ABLFL,
      !is.na(AVISITN) & AVISITN == 0 ~ "Y",
      TRUE ~ ""
    ),

    # Post-baseline flag
    ANL01FL = if_else(!is.na(AVISITN) & AVISITN > 0, "Y", ""),

    # Responders (CR or PR)
    ANL02FL = if_else(!is.na(AVALC) & AVALC %in% c("CR", "PR"), "Y", ""),

    # Has change from baseline
    ANL03FL = if_else(!is.na(PCHG), "Y", "")
  ) %>%
  # Parameter categories
  mutate(
    PARCAT1 = "TUMOR RESPONSE",
    PARCAT2 = case_when(
      PARAMCD == "TSIZE" ~ "MEASUREMENT",
      PARAMCD == "BOR" ~ "OVERALL RESPONSE",
      PARAMCD == "NADIR" ~ "BEST RESPONSE",
      TRUE ~ NA_character_
    )
  ) %>%
  # Set AVALU (now guaranteed to exist)
  mutate(
    AVALU = case_when(
      !is.na(AVALU) & AVALU != "" ~ AVALU, # Keep existing non-empty
      PARAMCD == "TSIZE" ~ "mm",
      PARAMCD == "NADIR" ~ "mm",
      TRUE ~ NA_character_
    )
  ) %>%
  # Sequence number
  derive_var_obs_number(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(PARAMN, AVISITN),
    new_var = ASEQ,
    check_type = "error"
  ) %>%
  arrange(USUBJID, PARAMN, AVISITN)

## ----r------------------------------------------------------------------------
# ---- Combine with covariates

adtrr_prefinal <- adtrr_seq %>%
  derive_vars_merged(
    dataset_add = covar_auc,
    by_vars = get_admiral_option("subject_keys")
  )

## ----r------------------------------------------------------------------------
## Check Data With metacore and metatools

adtrr <- adtrr_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

adtrr_xpt <- adtrr %>%
  xportr_type(metacore, domain = "ADTRR") %>% # Coerce variable type to match specs
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore, domain = "ADTRR") %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adtrr.xpt")) # Write xpt v5 transport file
