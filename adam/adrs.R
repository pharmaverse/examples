## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(admiral)
library(admiralonco)
library(metacore)
library(metatools)
library(xportr)
library(pharmaversesdtm)
library(pharmaverseadam)
library(dplyr)
library(lubridate)
library(stringr)

## ----r read-specs-------------------------------------------------------------
metacore <- spec_to_metacore("./metadata/onco_spec.xlsx") %>%
  select_dataset("ADRS")

## ----r load-data--------------------------------------------------------------
adsl <- pharmaverseadam::adsl
rs <- pharmaversesdtm::rs_onco_recist
tu <- pharmaversesdtm::tu_onco_recist

# Convert blanks to NA
rs <- convert_blanks_to_na(rs)
tu <- convert_blanks_to_na(tu)

## ----r merge-adsl-rs----------------------------------------------------------
adsl_vars <- exprs(RANDDT)
adrs_merged <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = exprs(STUDYID, USUBJID)
)

## ----r set-param-details------------------------------------------------------
adrs_ovr <- adrs_merged %>%
  filter(RSEVAL == "INVESTIGATOR" & RSTESTCD == "OVRLRESP") %>%
  mutate(
    PARAMCD = "OVR",
    PARAM = "Overall Response by Investigator",
    PARCAT1 = "Tumor Response",
    PARCAT2 = "Investigator",
    PARCAT3 = "RECIST 1.1"
  )

## ----r impute-dates-----------------------------------------------------------
adrs_imputed <- adrs_ovr %>%
  derive_vars_dt(
    dtc = RSDTC,
    new_vars_prefix = "A",
    highest_imputation = "D",
    date_imputation = "last"
  ) %>%
  mutate(AVISIT = VISIT)

## ----r derive-aval------------------------------------------------------------
adrs_aval <- adrs_imputed %>%
  mutate(
    AVALC = RSSTRESC,
    AVAL = aval_resp(AVALC)
  )

## ----r flag-worst-assessment--------------------------------------------------
worst_resp <- function(arg) {
  case_when(
    arg == "NE" ~ 1,
    arg == "CR" ~ 2,
    arg == "PR" ~ 3,
    arg == "SD" ~ 4,
    arg == "NON-CR/NON-PD" ~ 5,
    arg == "PD" ~ 6,
    TRUE ~ 0
  )
}

adrs_anl01fl <- adrs_aval %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, ADT),
      order = exprs(worst_resp(AVALC), RSSEQ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ADT >= RANDDT
  )

## ----r derive-pd--------------------------------------------------------------
adrs_pd <- adrs_anl01fl %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs_anl01fl,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = PARAMCD == "OVR" & AVALC == "PD" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "first",
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "PD",
      PARAM = "Disease Progression by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----r derive-death-----------------------------------------------------------
adsldth <- adsl %>%
  select(STUDYID, USUBJID, DTHDT, !!!adsl_vars)

adrs_death <- adrs_pd %>%
  derive_extreme_records(
    dataset_ref = adsldth,
    dataset_add = adsldth,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = !is.na(DTHDT),
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "DEATH",
      PARAM = "Death",
      PARCAT1 = "Reference Event",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ADT = DTHDT
    )
  ) %>%
  select(-DTHDT)

## ----r derive-lsta------------------------------------------------------------
adrs_lsta <- adrs_death %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs_death,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = PARAMCD == "OVR" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "LSTA",
      PARAM = "Last Disease Assessment by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      ANL01FL = "Y"
    )
  )

## ----r apply-metadata-check, message=FALSE, warning=FALSE---------------------
adrs_checked <- adrs_lsta %>%
  add_variables(metacore) %>% # Add variables specified in the metadata
  drop_unspec_vars(metacore) %>% # Drop variables not specified in metadata
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Check controlled terminology
  order_cols(metacore) %>% # Order columns according to metadata
  sort_by_key(metacore) # Sort rows by sort keys

## ----r xportr-----------------------------------------------------------------
dir <- tempdir() # Specify the directory for saving the XPT file

adrs_xpt <- adrs_checked %>%
  xportr_type(metacore, domain = "ADRS") %>% # Coerce variable types to match metadata
  xportr_length(metacore) %>% # Assign variable lengths from metadata
  xportr_label(metacore) %>% # Assign variable labels from metadata
  xportr_format(metacore) %>% # Assign variable formats from metadata
  xportr_df_label(metacore) %>% # Assign dataset labels from metadata
  xportr_write(file.path(dir, "adrs.xpt")) # Write the XPT file
