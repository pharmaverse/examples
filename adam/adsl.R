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
dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
vs <- pharmaversesdtm::vs
suppdm <- pharmaversesdtm::suppdm

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values
dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
vs <- convert_blanks_to_na(vs)
suppdm <- convert_blanks_to_na(suppdm)

# Combine Parent and Supp - very handy! ----
dm_suppdm <- combine_supp(dm, suppdm)

## ----r metacore, warning=FALSE, results='hold'--------------------------------
# Read in metacore object
metacore <- spec_to_metacore(
  path = "./metadata/safety_specs.xlsx",
  where_sep_sheet = FALSE,
  quiet = TRUE
) %>%
  select_dataset("ADSL")

## ----r demographics-----------------------------------------------------------
adsl_preds <- build_from_derived(metacore,
                                 ds_list = list("dm" = dm_suppdm, "suppdm" = dm_suppdm),
                                 predecessor_only = FALSE, keep = FALSE)

## ----r------------------------------------------------------------------------
get_control_term(metacore, variable = AGEGR1)

## ----r grouping_option_1------------------------------------------------------
adsl_ct <- adsl_preds %>%
  create_cat_var(metacore, ref_var = AGE,
                 grp_var = AGEGR1, num_grp_var = AGEGR1N)


## ----r grouping_option_2------------------------------------------------------
agegr1_lookup <- exprs(
  ~condition,            ~AGEGR1, ~AGEGR1N,
  AGE < 18,                "<18",        1,
  between(AGE, 18, 64),  "18-64",        2,
  !is.na(AGE),             ">64",        3,
  is.na(AGE),          "Missing",        4
)

adsl_cat <- derive_vars_cat(
  dataset = adsl_preds,
  definition = agegr1_lookup
)

## ----r grouping_option_3------------------------------------------------------
format_agegr1 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 64) ~ "18-64",
    var_input > 64 ~ ">64",
    TRUE ~ "Missing"
  )
}

format_agegr1n <- function(var_input) {
  case_when(
    var_input < 18 ~ 1,
    between(var_input, 18, 64) ~ 2,
    var_input > 64 ~ 3,
    TRUE ~ 4
  )
}

adsl_cust <- adsl_preds %>%
  mutate(
    AGEGR1 = format_agegr1(AGE),
    AGEGR1N = format_agegr1n(AGE)
  )

## ----r codelist---------------------------------------------------------------
adsl_ct <- adsl_ct %>%
  create_var_from_codelist(metacore = metacore,
                           input_var = RACE,
                           out_var = RACEN)

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
  # Treatment Start Datetime
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Treatment End Datetime
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Treatment Start and End Date
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM)) %>%  # Convert Datetime variables to date
  # Treatment Start Time
  derive_vars_dtm_to_tm(source_vars = exprs(TRTSDTM)) %>% 
  # Treatment Duration
  derive_var_trtdurd() %>%
  # Safety Population Flag
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = SAFFL,
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  ) %>%
  drop_unspec_vars(metacore) # This will drop any columns that are not specified in the metacore object

## ----r treatment_char, eval=TRUE----------------------------------------------
adsl <- adsl_raw %>%
  mutate(TRT01P = if_else(ARM %in% c("Screen Failure", "Not Assigned", "Not Treated"), "No Treatment", ARM),
         TRT01A = if_else(ACTARM %in% c("Screen Failure", "Not Assigned", "Not Treated"), "No Treatment", ACTARM)
  )

## ----r treatment_num, eval=TRUE-----------------------------------------------
adsl <- adsl %>%
  create_var_from_codelist(metacore, input_var = TRT01P, out_var = TRT01PN) %>%
  create_var_from_codelist(metacore, input_var = TRT01A, out_var = TRT01AN)

## ----r disposition, eval=TRUE-------------------------------------------------
# Convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )

## ----r eval=TRUE--------------------------------------------------------------
format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "DISCONTINUED"
  )
}

## ----r eval=TRUE--------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )

## ----r eval=TRUE--------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC,
    highest_imputation = "M",
    date_imputation = "first"
  ) 

## ----r eval=TRUE--------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(RANDDT = DSSTDT),
    filter_add = DSDECOD == "RANDOMIZED",
  ) %>% 
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(SCRFDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE"
  ) %>% 
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(FRVDT = DSSTDT),
    filter_add = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
  )

## ----r eval=TRUE--------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  ) %>% 
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  )

## ----r eval=TRUE--------------------------------------------------------------
assign_randfl <- function(x) {
  if_else(!is.na(x), "Y", NA_character_)
}

adsl <- adsl %>% 
  mutate(
    RANDFL = assign_randfl(RANDDT)
  )

## ----r death, eval=TRUE-------------------------------------------------------
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD, DTHDOM = "AE"),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM, DTHDOM = "DS"),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS, DTHDOM)
  )

## ----r grouping, eval=TRUE----------------------------------------------------
region1_lookup <- exprs(
  ~condition,                              ~REGION1, ~REGION1N,
  COUNTRY %in% c("CAN", "USA"),     "North America",         1,
  !is.na(COUNTRY),              "Rest of the World",         2,
  is.na(COUNTRY),                         "Missing",         3
)

racegr1_lookup <- exprs(
  ~condition,              ~RACEGR1, ~RACEGR1N,
  RACE %in% c("WHITE"),     "White",        1,
  RACE != "WHITE",      "Non-white",        2,
  is.na(RACE),            "Missing",        3
)

dthcgr1_lookup <- exprs(
  ~condition,                                                                                 ~DTHCGR1, ~DTHCGR1N,
  DTHDOM == "AE",                                                                      "ADVERSE EVENT",         1,
  !is.na(DTHDOM) & str_detect(DTHCAUS, "(PROGRESSIVE DISEASE|DISEASE RELAPSE)"), "PROGRESSIVE DISEASE",         2,
  !is.na(DTHDOM) & !is.na(DTHCAUS),                                                            "OTHER",         3,
  is.na(DTHDOM),                                                                         NA_character_,        NA
)

adsl <- adsl %>% 
  derive_vars_cat(
    definition = region1_lookup
  ) %>% 
  derive_vars_cat(
    definition = racegr1_lookup
  ) %>% 
  derive_vars_cat(
    definition = dthcgr1_lookup
  ) 

## ----r checks, warning=FALSE, message=FALSE-----------------------------------
dir <- tempdir() # Specify the directory for saving the XPT file

adsl %>%
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore, na_acceptable = TRUE) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) %>% # Sorts the rows by the sort keys
  xportr_type(metacore, domain = "ADSL") %>% # Coerce variable type to match spec
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_df_label(metacore) %>%  # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adsl.xpt"), metadata = metacore, domain = "ADSL")
