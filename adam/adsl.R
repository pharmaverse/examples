## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(metacore)
library(metatools)
library(pharmaversesdtm)
library(admiral)
library(xportr)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# Read in input SDTM data
data("dm")
data("ex")

## ----r metacore, warning=FALSE, results='hold'--------------------------------
# Read in metacore object
load(metacore_example("pilot_ADaM.rda"))
metacore <- metacore %>%
   select_dataset("ADSL")

## ----r------------------------------------------------------------------------
metacore$ds_vars

## ----r demographics-----------------------------------------------------------
adsl_preds <- build_from_derived(metacore,
                                 ds_list = list("dm" = dm),
                                 predecessor_only = FALSE, keep = TRUE)
head(adsl_preds, n=10)

## ----r------------------------------------------------------------------------
get_control_term(metacore, variable = AGEGR1)

## ----r ct---------------------------------------------------------------------
adsl_ct <- adsl_preds %>%
   create_cat_var(metacore, ref_var = AGE,
                  grp_var = AGEGR1, num_grp_var = AGEGR1N) %>%
   create_var_from_codelist(metacore = metacore,
                            input_var = RACE,
                            out_var = RACEN) %>%
   # Removing screen failures from ARM and TRT01P to match the define and FDA guidance
   mutate(ARM = if_else(ARM == "Screen Failure", NA_character_, ARM),
          TRT01P = if_else(TRT01P == "Screen Failure", NA_character_, TRT01P)
   )

head(adsl_ct, n=10)

## ----r exposure---------------------------------------------------------------
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

adsl_raw <- adsl_ct %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & nchar(EXSTDTC) >= 10,
    new_vars = exprs(TRTSDTM = EXSTDTM),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10,
    new_vars = exprs(TRTEDTM = EXENDTM),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
   derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM)) %>%  # Convert Datetime variables to date
   derive_var_trtdurd() %>%
   derive_var_merged_exist_flag(
     dataset_add = ex,
     by_vars = exprs(STUDYID, USUBJID),
     new_var = SAFFL,
     condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
   ) %>%
   drop_unspec_vars(metacore) # This will drop any columns that aren't specified in the metacore object

head(adsl_raw, n=10)

## ----r, warning=FALSE, message=FALSE, include=FALSE---------------------------
# Create dummy variables to match metacore specs to avoid later errors
# In practice these would be mainly created using derivation functions from admiral
adsl_raw <- adsl_raw %>%
  mutate(
    SITEGR1 = NA,
    TRT01PN = NA,
    TRT01A = NA,
    TRT01AN = NA,
    AVGDD = NA,
    CUMDOSE = NA,
    ITTFL = NA,
    EFFFL = NA,
    COMP8FL = NA,
    COMP16FL = NA,
    COMP24FL = NA,
    DISCONFL = NA,
    DSRAEFL = NA,
    BMIBL = NA,
    BMIBLGR1 = NA,
    HEIGHTBL = NA,
    WEIGHTBL = NA,
    EDUCLVL = NA,
    DISONSDT = NA,
    DURDIS = NA,
    DURDSGR1 = NA,
    VISIT1DT = NA,
    VISNUMEN = NA,
    RFENDT = NA,
    DCDECOD = NA,
    EOSSTT = NA,
    DCSREAS = NA,
    MMSETOT = NA
  )

## ----r checks, warning=FALSE, message=FALSE-----------------------------------

adsl_raw %>%
   check_variables(metacore) %>% # Check all variables specified are present and no more
   check_ct_data(metacore, na_acceptable = TRUE) %>% # Checks all variables with CT only contain values within the CT
   order_cols(metacore) %>% # Orders the columns according to the spec
   sort_by_key(metacore) %>% # Sorts the rows by the sort keys
   xportr_type(metacore, domain = "ADSL") %>% # Coerce variable type to match spec
   xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
   xportr_label(metacore) %>% # Assigns variable label from metacore specifications
   xportr_df_label(metacore) # Assigns dataset label from metacore specifications

