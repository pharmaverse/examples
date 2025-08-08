## ----r message=FALSE, warning=FALSE-------------------------------------------
library(admiral)
library(admiralonco)
library(dplyr)
library(lubridate)
library(metacore)
library(metatools)
library(xportr)
library(pharmaverseadam)

## ----r read-specs-------------------------------------------------------------
# Load metacore specifications
metacore <- spec_to_metacore("./metadata/onco_spec.xlsx") %>%
  select_dataset("ADTTE")

# Load source datasets
adsl <- pharmaverseadam::adsl
adrs <- pharmaverseadam::adrs_onco

adrs <- adrs_onco

## ----r------------------------------------------------------------------------
# Define event and censoring sources
death_event <- event_source(
  dataset_name = "adrs",
  filter = PARAMCD == "DEATH" & AVALC == "Y" & ANL01FL == "Y",
  date = ADT,
  set_values_to = exprs(
    EVNTDESC = "Death",
    SRCDOM = "ADRS",
    SRCVAR = "ADT"
  )
)

pd_event <- event_source(
  dataset_name = "adrs",
  filter = PARAMCD == "PD" & ANL01FL == "Y",
  date = ADT,
  set_values_to = exprs(
    EVNTDESC = "Progressive Disease",
    SRCDOM = "ADRS",
    SRCVAR = "ADT"
  )
)

lastalive_censor <- censor_source(
  dataset_name = "adsl",
  date = LSTALVDT,
  set_values_to = exprs(
    EVNTDESC = "Last Known Alive",
    CNSDTDSC = "Last Known Alive Date",
    SRCDOM = "ADSL",
    SRCVAR = "LSTALVDT"
  )
)

lasta_censor <- censor_source(
  dataset_name = "adrs",
  filter = PARAMCD == "LSTA" & ANL01FL == "Y",
  date = ADT,
  set_values_to = exprs(
    EVNTDESC = "Progression Free Alive",
    CNSDTDSC = "Last Tumor Assessment",
    SRCDOM = "ADRS",
    SRCVAR = "ADT"
  )
)

rand_censor <- censor_source(
  dataset_name = "adsl",
  date = RANDDT,
  set_values_to = exprs(
    EVNTDESC = "Randomization Date",
    CNSDTDSC = "Randomization Date",
    SRCDOM = "ADSL",
    SRCVAR = "RANDDT"
  )
)

## ----r------------------------------------------------------------------------
# Derive Overall Survival (OS)
adtte <- derive_param_tte(
  dataset_adsl = adsl,
  start_date = RANDDT,
  event_conditions = list(death_event),
  censor_conditions = list(lastalive_censor, rand_censor),
  source_datasets = list(adsl = adsl, adrs = adrs),
  set_values_to = exprs(PARAMCD = "OS", PARAM = "Overall Survival")
)

# Derive Progression-Free Survival (PFS)
adtte_pfs <- adtte %>%
  derive_param_tte(
    dataset_adsl = adsl,
    start_date = RANDDT,
    event_conditions = list(pd_event, death_event),
    censor_conditions = list(lasta_censor, rand_censor),
    source_datasets = list(adsl = adsl, adrs = adrs),
    set_values_to = exprs(PARAMCD = "PFS", PARAM = "Progression-Free Survival")
  )

## ----r------------------------------------------------------------------------
# Derive analysis value
adtte_aval <- adtte_pfs %>%
  derive_vars_duration(
    new_var = AVAL,
    start_date = STARTDT,
    end_date = ADT
  )

## ----r------------------------------------------------------------------------
# Derive analysis sequence number
adtte_aseq <- adtte_aval %>%
  derive_var_obs_number(
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD),
    check_type = "error"
  )

## ----r------------------------------------------------------------------------
# Add ADSL variables
adtte_adsl <- adtte_aseq %>%
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r, message=FALSE, warning=FALSE------------------------------------------
# Apply metadata and perform checks
adtte_adsl_checked <- adtte_adsl %>%
  add_variables(metacore) %>%
  drop_unspec_vars(metacore) %>%
  check_variables(metacore) %>%
  check_ct_data(metacore) %>%
  order_cols(metacore) %>%
  sort_by_key(metacore)

# Apply apply labels, formats, and export the dataset to an XPT file.
adtte_final <- adtte_adsl_checked %>%
  xportr_type(metacore, domain = "ADTTE") %>%
  xportr_length(metacore) %>%
  xportr_label(metacore) %>%
  xportr_df_label(metacore)

# Write dataset to XPT file (optional)
dir <- tempdir()
xportr_write(adtte_final, file.path(dir, "adtte.xpt"))
