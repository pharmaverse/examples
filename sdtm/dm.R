## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(sdtm.oak)
library(pharmaverseraw)
library(dplyr)

dm_raw <- pharmaverseraw::dm_raw
ds_raw <- pharmaverseraw::ds_raw
ec_raw <- pharmaverseraw::ec_raw

## ----r------------------------------------------------------------------------
dm_raw <- dm_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "dm_raw"
  )

ds_raw <- ds_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ds_raw"
  )

ec_raw <- ec_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ec_raw"
  )

## ----r, echo = TRUE-----------------------------------------------------------
study_ct <- read.csv("metadata/sdtm_ct.csv")

## ----r------------------------------------------------------------------------
ref_date_conf_df <- tibble::tribble(
  ~raw_dataset_name, ~date_var,     ~time_var,      ~dformat,      ~tformat, ~sdtm_var_name,
  "ec_raw",       "IT.ECSTDAT", NA_character_, "dd-mmm-yyyy", NA_character_,     "RFXSTDTC",
  "ec_raw",       "IT.ECENDAT", NA_character_, "dd-mmm-yyyy", NA_character_,     "RFXENDTC",
  "ec_raw",       "IT.ECSTDAT", NA_character_, "dd-mmm-yyyy", NA_character_,      "RFSTDTC",
  "ec_raw",       "IT.ECENDAT", NA_character_, "dd-mmm-yyyy", NA_character_,      "RFENDTC",
  "dm_raw",            "IC_DT", NA_character_,  "mm/dd/yyyy", NA_character_,      "RFICDTC",
  "ds_raw",          "DSDTCOL",     "DSTMCOL",  "mm-dd-yyyy",         "H:M",     "RFPENDTC",
  "ds_raw",          "DEATHDT", NA_character_,  "mm/dd/yyyy", NA_character_,       "DTHDTC"
)

## ----r------------------------------------------------------------------------
dm <- dm_raw %>%
  mutate(
    SUBJID = substr(PATNUM, 5, 8)
  ) %>%
  select(oak_id, raw_source, patient_number, SUBJID)

## ----r------------------------------------------------------------------------
dm <- dm %>%
  # Map AGE using assign_no_ct
  assign_no_ct(
    raw_dat = dm_raw,
    raw_var = "IT.AGE",
    tgt_var = "AGE",
    id_vars = oak_id_vars()
  ) %>%
  # Map AGEU using hardcode_ct
  hardcode_ct(
    raw_dat = dm_raw,
    raw_var = "IT.AGE",
    tgt_var = "AGEU",
    tgt_val = "Year",
    ct_spec = study_ct,
    ct_clst = "C66781",
    id_vars = oak_id_vars()
  ) %>%
  # Map SEX using assign_ct
  assign_ct(
    raw_dat = dm_raw,
    raw_var = "IT.SEX",
    tgt_var = "SEX",
    ct_spec = study_ct,
    ct_clst = "C66731",
    id_vars = oak_id_vars()
  ) %>%
  # Map ETHNIC using assign_ct
  assign_ct(
    raw_dat = dm_raw,
    raw_var = "IT.ETHNIC",
    tgt_var = "ETHNIC",
    ct_spec = study_ct,
    ct_clst = "C66790",
    id_vars = oak_id_vars()
  ) %>%
  # Map RACE using assign_ct
  assign_ct(
    raw_dat = dm_raw,
    raw_var = "IT.RACE",
    tgt_var = "RACE",
    ct_spec = study_ct,
    ct_clst = "C74457",
    id_vars = oak_id_vars()
  ) %>%
  # Map ARM using assign_ct
  assign_ct(
    raw_dat = dm_raw,
    raw_var = "PLANNED_ARM",
    tgt_var = "ARM",
    ct_spec = study_ct,
    ct_clst = "ARM",
    id_vars = oak_id_vars()
  ) %>%
  # Map ARMCD using assign_no_ct
  assign_no_ct(
    raw_dat = dm_raw,
    raw_var = "PLANNED_ARMCD",
    tgt_var = "ARMCD",
    id_vars = oak_id_vars()
  ) %>%
  # Map ACTARM using assign_ct
  assign_ct(
    raw_dat = dm_raw,
    raw_var = "ACTUAL_ARM",
    tgt_var = "ACTARM",
    ct_spec = study_ct,
    ct_clst = "ARM",
    id_vars = oak_id_vars()
  ) %>%
  # Map ACTARMCD using assign_no_ct
  assign_no_ct(
    raw_dat = dm_raw,
    raw_var = "ACTUAL_ARMCD",
    tgt_var = "ACTARMCD",
    id_vars = oak_id_vars()
  )

## ----r eval=TRUE--------------------------------------------------------------
dm <- dm %>%
  # Derive RFSTDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "RFSTDTC",
    min_max = "min",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  )

## ----r------------------------------------------------------------------------
dm <- dm %>%
  # Derive RFENDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "RFENDTC",
    min_max = "max",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  )

## ----r------------------------------------------------------------------------
dm <- dm %>%
  # Derive RFXSTDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "RFXSTDTC",
    min_max = "min",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  ) %>%
  # Derive RFXENDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "RFXENDTC",
    min_max = "max",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  ) %>%
  # Derive RFICDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "RFICDTC",
    min_max = "min",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  ) %>%
  # Derive RFPENDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "RFPENDTC",
    min_max = "max",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  ) %>%
  # Map DTHDTC using oak_cal_ref_dates
  oak_cal_ref_dates(
    ds_in = .,
    der_var = "DTHDTC",
    min_max = "min",
    ref_date_config_df = ref_date_conf_df,
    raw_source = list(
      ec_raw = ec_raw,
      ds_raw = ds_raw,
      dm_raw = dm_raw
    )
  )

## ----r------------------------------------------------------------------------
dm <- dm %>%
  mutate(
    STUDYID = dm_raw$STUDY,
    DOMAIN = "DM",
    USUBJID = paste0("01-", dm_raw$PATNUM),
    COUNTRY = dm_raw$COUNTRY,
    DTHFL = dplyr::if_else(is.na(DTHDTC), NA_character_, "Y")
  ) %>%
  # Map DMDTC using assign_datetime
  assign_datetime(
    raw_dat = dm_raw,
    raw_var = "COL_DT",
    tgt_var = "DMDTC",
    raw_fmt = c("m/d/y"),
    id_vars = oak_id_vars()
  ) %>%
  # Derive study day
  derive_study_day(
    sdtm_in = .,
    dm_domain = .,
    tgdt = "DMDTC",
    refdt = "RFXSTDTC",
    study_day_var = "DMDY"
  )
