## ----r echo=TRUE, message=FALSE-----------------------------------------------
# Load Packages
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)
library(metacore)
library(metatools)
library(xportr)
library(readr)
library(pharmaversesdtm)
library(pharmaverseadam)

## ----r echo=TRUE, message=FALSE-----------------------------------------------
# ---- Load Specs for Metacore ----
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx") %>%
  select_dataset("ADPPK")

## ----r------------------------------------------------------------------------
# ---- Load source datasets ----
# Load PC, EX, VS, LB and ADSL

ex <- pharmaversesdtm::ex
pc <- pharmaversesdtm::pc
vs <- pharmaversesdtm::vs
lb <- pharmaversesdtm::lb

adsl <- pharmaverseadam::adsl

ex <- convert_blanks_to_na(ex)
pc <- convert_blanks_to_na(pc)
vs <- convert_blanks_to_na(vs)
lb <- convert_blanks_to_na(lb)

## ----r------------------------------------------------------------------------
# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTSDTM, TRT01P, TRT01A)

pc_dates <- pc %>%
  # Join ADSL with PC (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Derive analysis date/time
  # Impute missing time to 00:00:00
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = PCDTC,
    time_imputation = "00:00:00"
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST,
    NFRLT = if_else(PCTPTNUM < 0, 0, PCTPTNUM), .after = USUBJID
  )

## ----r------------------------------------------------------------------------
# ---- Get dosing information ----

ex_dates <- ex %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Keep records with nonzero dose
  filter(EXDOSE > 0) %>%
  # Add time and set missing end date to start date
  # Impute missing time to 00:00:00
  # Note all times are missing for dosing records in this example data
  # Derive Analysis Start and End Dates
  derive_vars_dtm(
    new_vars_prefix = "AST",
    dtc = EXSTDTC,
    time_imputation = "00:00:00"
  ) %>%
  derive_vars_dtm(
    new_vars_prefix = "AEN",
    dtc = EXENDTC,
    time_imputation = "00:00:00"
  ) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 1,
    NFRLT = case_when(
      VISITDY == 1 ~ 0,
      TRUE ~ 24 * VISITDY
    )
  ) %>%
  # Set missing end dates to start date
  mutate(AENDTM = case_when(
    is.na(AENDTM) ~ ASTDTM,
    TRUE ~ AENDTM
  )) %>%
  # Derive dates from date/times
  derive_vars_dtm_to_dt(exprs(ASTDTM)) %>%
  derive_vars_dtm_to_dt(exprs(AENDTM))

## ----r------------------------------------------------------------------------
ex_exp <- ex_dates %>%
  create_single_dose_dataset(
    dose_freq = EXDOSFRQ,
    start_date = ASTDT,
    start_datetime = ASTDTM,
    end_date = AENDT,
    end_datetime = AENDTM,
    nominal_time = NFRLT,
    lookup_table = dose_freq_lookup,
    lookup_column = CDISC_VALUE,
    keep_source_vars = exprs(
      STUDYID, USUBJID, EVID, EXDOSFRQ, EXDOSFRM,
      NFRLT, EXDOSE, EXDOSU, EXTRT, ASTDT, ASTDTM, AENDT, AENDTM,
      VISIT, VISITNUM, VISITDY,
      TRT01A, TRT01P, DOMAIN, EXSEQ, !!!adsl_vars
    )
  ) %>%
  # Derive AVISIT based on nominal relative time
  # Derive AVISITN to nominal time in whole days using integer division
  # Define AVISIT based on nominal day
  mutate(
    AVISITN = NFRLT %/% 24 + 1,
    AVISIT = paste("Day", AVISITN),
    ADTM = ASTDTM,
    DRUG = EXTRT
  ) %>%
  # Derive dates and times from datetimes
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ASTDTM)) %>%
  derive_vars_dtm_to_tm(exprs(AENDTM))

## ----r------------------------------------------------------------------------
# ---- Find first dose per treatment per subject ----
# ---- Join with ADPPK data and keep only subjects with dosing ----

adppk_first_dose <- pc_dates %>%
  derive_vars_merged(
    dataset_add = ex_exp,
    filter_add = (!is.na(ADTM)),
    new_vars = exprs(FANLDTM = ADTM, EXDOSE_first = EXDOSE),
    order = exprs(ADTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID, DRUG)
  ) %>%
  filter(!is.na(FANLDTM)) %>%
  # Derive AVISIT based on nominal relative time
  # Derive AVISITN to nominal time in whole days using integer division
  # Define AVISIT based on nominal day
  mutate(
    AVISITN = NFRLT %/% 24 + 1,
    AVISIT = paste("Day", AVISITN),
  )

## ----r------------------------------------------------------------------------
# ---- Find previous dose  ----

adppk_prev <- adppk_first_dose %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(ADTM),
    new_vars = exprs(
      ADTM_prev = ADTM, EXDOSE_prev = EXDOSE, AVISIT_prev = AVISIT,
      AENDTM_prev = AENDTM
    ),
    join_vars = exprs(ADTM),
    join_type = "all",
    filter_add = NULL,
    filter_join = ADTM > ADTM.join,
    mode = "last",
    check_type = "none"
  )

## ----r------------------------------------------------------------------------
adppk_nom_prev <- adppk_prev %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(NFRLT),
    new_vars = exprs(NFRLT_prev = NFRLT),
    join_vars = exprs(NFRLT),
    join_type = "all",
    filter_add = NULL,
    filter_join = NFRLT > NFRLT.join,
    mode = "last",
    check_type = "none"
  )

## ----r------------------------------------------------------------------------
adppk_aprlt <- bind_rows(adppk_nom_prev, ex_exp) %>%
  group_by(USUBJID, DRUG) %>%
  mutate(
    FANLDTM = min(FANLDTM, na.rm = TRUE),
    min_NFRLT = min(NFRLT, na.rm = TRUE),
    maxdate = max(ADT[EVID == 0], na.rm = TRUE), .after = USUBJID
  ) %>%
  arrange(USUBJID, ADTM) %>%
  ungroup() %>%
  filter(ADT <= maxdate) %>%
  # Derive Actual Relative Time from First Dose (AFRLT)
  derive_vars_duration(
    new_var = AFRLT,
    start_date = FANLDTM,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Reference Dose (APRLT)
  derive_vars_duration(
    new_var = APRLT,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive APRLT
  mutate(
    APRLT = case_when(
      EVID == 1 ~ 0,
      is.na(APRLT) ~ AFRLT,
      TRUE ~ APRLT
    ),
    NPRLT = case_when(
      EVID == 1 ~ 0,
      is.na(NFRLT_prev) ~ NFRLT - min_NFRLT,
      TRUE ~ NFRLT - NFRLT_prev
    )
  )

## ----r------------------------------------------------------------------------
# ---- Derive Analysis Variables ----
# Derive actual dose DOSEA and planned dose DOSEP,
# Derive AVAL and DV

adppk_aval <- adppk_aprlt %>%
  mutate(
    # Derive Actual Dose
    DOSEA = case_when(
      EVID == 1 ~ EXDOSE,
      is.na(EXDOSE_prev) ~ EXDOSE_first,
      TRUE ~ EXDOSE_prev
    ),
    # Derive Planned Dose
    DOSEP = case_when(
      TRT01P == "Xanomeline High Dose" ~ 81,
      TRT01P == "Xanomeline Low Dose" ~ 54,
      TRT01P == "Placebo" ~ 0
    ),
    # Derive PARAMCD
    PARAMCD = case_when(
      EVID == 1 ~ "DOSE",
      TRUE ~ PCTESTCD
    ),
    ALLOQ = PCLLOQ,
    # Derive CMT
    CMT = case_when(
      EVID == 1 ~ 1,
      PCSPEC == "PLASMA" ~ 2,
      TRUE ~ 3
    ),
    # Derive BLQFL/BLQFN
    BLQFL = case_when(
      PCSTRESC == "<BLQ" ~ "Y",
      TRUE ~ "N"
    ),
    BLQFN = case_when(
      PCSTRESC == "<BLQ" ~ 1,
      TRUE ~ 0
    ),
    AMT = case_when(
      EVID == 1 ~ EXDOSE,
      TRUE ~ NA_real_
    ),
    # Derive DV and AVAL
    DV = PCSTRESN,
    DVID = PCTESTCD,
    AVAL = DV,
    DVL = case_when(
      DV != 0 ~ log(DV),
      TRUE ~ NA_real_
    ),
    # Derive MDV
    MDV = case_when(
      EVID == 1 ~ 1,
      is.na(DV) ~ 1,
      TRUE ~ 0
    ),
    AVALU = case_when(
      EVID == 1 ~ NA_character_,
      TRUE ~ PCSTRESU
    ),
    RLTU = "h",
    USTRESC = PCSTRESC,
    UDTC = format_ISO8601(ADTM),
    II = if_else(EVID == 1, 1, 0),
    SS = if_else(EVID == 1, 1, 0),
    ADDL = 0,
    OCC = 1,
  )

## ----r------------------------------------------------------------------------
# ---- Add ASEQ ----

adppk_aseq <- adppk_aval %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(AFRLT, EVID, CMT),
    check_type = "error"
  ) %>%
  mutate(
    PROJID = DRUG,
    PROJIDN = 1,
    PART = 1,
  )

## ----r------------------------------------------------------------------------
#---- Derive Covariates ----
# Include numeric values for STUDYIDN, USUBJIDN, SEXN, RACEN etc.

covar <- adsl %>%
  create_var_from_codelist(metacore, input_var = STUDYID, out_var = STUDYIDN) %>%
  create_var_from_codelist(metacore, input_var = SEX, out_var = SEXN) %>%
  create_var_from_codelist(metacore, input_var = RACE, out_var = RACEN) %>%
  create_var_from_codelist(metacore, input_var = ETHNIC, out_var = AETHNIC) %>%
  create_var_from_codelist(metacore, input_var = AETHNIC, out_var = AETHNICN) %>%
  create_var_from_codelist(metacore, input_var = ARMCD, out_var = COHORT) %>%
  create_var_from_codelist(metacore, input_var = ARMCD, out_var = COHORTC) %>%
  create_var_from_codelist(metacore, input_var = COUNTRY, out_var = COUNTRYN) %>%
  create_var_from_codelist(metacore, input_var = COUNTRY, out_var = COUNTRYL) %>%
  mutate(
    STUDYIDN = as.numeric(word(USUBJID, 1, sep = fixed("-"))),
    SITEIDN = as.numeric(word(USUBJID, 2, sep = fixed("-"))),
    USUBJIDN = as.numeric(word(USUBJID, 3, sep = fixed("-"))),
    SUBJIDN = as.numeric(SUBJID),
    ROUTE = unique(ex$EXROUTE),
    FORM = unique(ex$EXDOSFRM),
    REGION1 = COUNTRY,
    REGION1N = COUNTRYN,
    SUBJTYPC = "Volunteer",
  ) %>%
  create_var_from_codelist(metacore, input_var = FORM, out_var = FORMN) %>%
  create_var_from_codelist(metacore, input_var = ROUTE, out_var = ROUTEN) %>%
  create_var_from_codelist(metacore, input_var = SUBJTYPC, out_var = SUBJTYP)

## ----r------------------------------------------------------------------------
labsbl <- lb %>%
  filter(LBBLFL == "Y" & LBTESTCD %in% c("CREAT", "ALT", "AST", "BILI")) %>%
  mutate(LBTESTCDB = paste0(LBTESTCD, "BL")) %>%
  select(STUDYID, USUBJID, LBTESTCDB, LBSTRESN)

covar_vslb <- covar %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "HEIGHT",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(HTBL = VSSTRESN)
  ) %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WTBL = VSSTRESN)
  ) %>%
  derive_vars_transposed(
    dataset_merge = labsbl,
    by_vars = exprs(STUDYID, USUBJID),
    key_var = LBTESTCDB,
    value_var = LBSTRESN
  ) %>%
  mutate(
    BMIBL = compute_bmi(height = HTBL, weight = WTBL),
    BSABL = compute_bsa(
      height = HTBL,
      weight = HTBL,
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

adppk_prefinal <- adppk_aseq %>%
  derive_vars_merged(
    dataset_add = select(covar_vslb, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  arrange(STUDYIDN, USUBJIDN, AFRLT, EVID) %>%
  # Add RECSEQ
  # Exclude records if needed
  mutate(
    RECSEQ = row_number(),
    EXCLFCOM = "None"
  ) %>%
  create_var_from_codelist(metacore, input_var = DVID, out_var = DVIDN) %>%
  create_var_from_codelist(metacore, input_var = EXCLFCOM, out_var = EXCLF)

## ----r------------------------------------------------------------------------
adppk <- adppk_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

adppk_xpt <- adppk %>%
  xportr_type(metacore, domain = "ADPPK") %>% # Coerce variable type to match spec
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adppk.xpt")) # Write xpt v5 transport file

