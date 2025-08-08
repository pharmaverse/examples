## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(metacore)
library(metatools)
library(pharmaversesdtm)
library(pharmaverseadam)
library(admiral)
library(xportr)
library(dplyr)
library(lubridate)
library(stringr)
library(reactable)

# Read in input data
adsl <- pharmaverseadam::adsl
ae <- pharmaversesdtm::ae
ex <- pharmaversesdtm::ex

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values
ae <- convert_blanks_to_na(ae)
ex <- convert_blanks_to_na(ex)

## ----r echo=TRUE--------------------------------------------------------------
# ---- Load Specs for Metacore ----
metacore <- spec_to_metacore(
  path = "./metadata/safety_specs.xlsx",
  # All datasets are described in the same sheet
  where_sep_sheet = FALSE
) %>%
  select_dataset("ADAE")

## ----r------------------------------------------------------------------------
# Select required ADSL variables
adsl_vars <- exprs(TRTSDT, TRTEDT, DTHDT)

# Join ADSL variables with VS
adae <- ae %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r------------------------------------------------------------------------
# Derive ASTDT/ASTDTF/ASTDY and AENDT/AENDTF/AENDY
adae <- adae %>%
  derive_vars_dt(
    new_vars_prefix = "AEN",
    dtc = AEENDTC,
    date_imputation = "last",
    highest_imputation = "M", # imputation is performed on missing days or months
    flag_imputation = "auto" # to automatically create AENDTF variable
  ) %>%
  derive_vars_dt(
    new_vars_prefix = "AST",
    dtc = AESTDTC,
    highest_imputation = "M", # imputation is performed on missing days or months
    flag_imputation = "auto", # to automatically create ASTDTF variable
    min_dates = exprs(TRTSDT), # apply a minimum date for the imputation
    max_dates = exprs(AENDT) # apply a maximum date for the imputation
  ) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ASTDT, AENDT)
  )

## ----r------------------------------------------------------------------------
# Derive ADURN/ADURU
adae <- adae %>%
  derive_vars_duration(
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT
  )

## ----r------------------------------------------------------------------------
# Derive LDOSEDT
# In our ex data the EXDOSFRQ (frequency) is "QD" which stands for once daily
# If this was not the case then we would need to use the admiral::create_single_dose_dataset() function
# to generate single doses from aggregate dose information
# Refer to https://pharmaverse.github.io/admiral/reference/create_single_dose_dataset.html
ex <- ex %>%
  derive_vars_dt(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN"
  )

adae <- adae %>%
  derive_vars_joined(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(EXENDT),
    new_vars = exprs(LDOSEDT = EXENDT),
    join_vars = exprs(EXENDT),
    join_type = "all",
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDT),
    filter_join = EXENDT <= ASTDT,
    mode = "last"
  )

## ----r------------------------------------------------------------------------
# Derive TRTEMFL and ONTRTFL
adae <- adae %>%
  derive_var_trtemfl(
    start_date = ASTDT,
    end_date = AENDT,
    trt_start_date = TRTSDT,
    trt_end_date = TRTEDT
  ) %>%
  derive_var_ontrtfl(
    start_date = ASTDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 30
  )

## ----r------------------------------------------------------------------------
# Derive AOCCIFL
adae <- adae %>%
  # create temporary numeric ASEVN for sorting purpose
  mutate(TEMP_AESEVN = as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))) %>%
  derive_var_extreme_flag(
    new_var = AOCCIFL,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(TEMP_AESEVN, ASTDT, AESEQ),
    mode = "first"
  )

## ----r------------------------------------------------------------------------
queries <- admiral::queries %>%
  filter(PREFIX %in% c("CQ01", "SMQ02"))

## ----r------------------------------------------------------------------------
# Derive CQ01NAM and SMQ02NAM
adae <- adae %>%
  derive_vars_query(dataset_queries = queries)

## ----r eval=TRUE--------------------------------------------------------------
adae <- adae %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r checks, warning=FALSE, message=FALSE-----------------------------------
dir <- tempdir() # Specify the directory for saving the XPT file

adae %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore, na_acceptable = TRUE) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) %>% # Sorts the rows by the sort keys
  xportr_type(metacore, domain = "ADAE") %>% # Coerce variable type to match spec
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adae.xpt"), metadata = metacore, domain = "ADAE")

