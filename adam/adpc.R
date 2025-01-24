## ----r echo=TRUE, message=FALSE-----------------------------------------------
# Load Packages
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)
library(metacore)
library(metatools)
library(xportr)
library(pharmaversesdtm)
library(pharmaverseadam)
library(reactable)

## ----r echo=TRUE--------------------------------------------------------------
# ---- Load Specs for Metacore ----

metacore <- spec_to_metacore("./metadata/pk_spec.xlsx") %>%
  select_dataset("ADPC")

## ----r------------------------------------------------------------------------
# ---- Load source datasets ----
# Load PC, EX, VS, LB and ADSL
data("pc")
data("ex")
data("vs")

data("adsl")

ex <- convert_blanks_to_na(ex)
pc <- convert_blanks_to_na(pc)
vs <- convert_blanks_to_na(vs)

## ----r------------------------------------------------------------------------
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
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST,
    NFRLT = if_else(PCTPTNUM < 0, 0, PCTPTNUM), .after = USUBJID
  )

reactable(pc_dates)

## ----r------------------------------------------------------------------------
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

reactable(ex_dates)

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
  derive_vars_dtm_to_tm(exprs(AENDTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

reactable(ex_exp)

## ----r------------------------------------------------------------------------
adpc_first_dose <- pc_dates %>%
  derive_vars_merged(
    dataset_add = ex_exp,
    filter_add = (EXDOSE > 0 & !is.na(ADTM)),
    new_vars = exprs(FANLDTM = ADTM),
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

reactable(adpc_first_dose)

## ----r------------------------------------------------------------------------
adpc_prev <- adpc_first_dose %>%
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

adpc_next <- adpc_prev %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(ADTM),
    new_vars = exprs(
      ADTM_next = ADTM, EXDOSE_next = EXDOSE, AVISIT_next = AVISIT,
      AENDTM_next = AENDTM
    ),
    join_vars = exprs(ADTM),
    join_type = "all",
    filter_add = NULL,
    filter_join = ADTM <= ADTM.join,
    mode = "first",
    check_type = "none"
  )

reactable(adpc_prev)

## ----r------------------------------------------------------------------------
adpc_nom_prev <- adpc_next %>%
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

adpc_nom_next <- adpc_nom_prev %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(NFRLT),
    new_vars = exprs(NFRLT_next = NFRLT),
    join_vars = exprs(NFRLT),
    join_type = "all",
    filter_add = NULL,
    filter_join = NFRLT <= NFRLT.join,
    mode = "first",
    check_type = "none"
  )

reactable(adpc_nom_prev)

## ----r------------------------------------------------------------------------
adpc_arrlt <- bind_rows(adpc_nom_next, ex_exp) %>%
  group_by(USUBJID, DRUG) %>%
  mutate(
    FANLDTM = min(FANLDTM, na.rm = TRUE),
    min_NFRLT = min(NFRLT_prev, na.rm = TRUE),
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
  # Derive Actual Relative Time from Reference Dose (ARRLT)
  derive_vars_duration(
    new_var = ARRLT,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Next Dose (AXRLT not kept)
  derive_vars_duration(
    new_var = AXRLT,
    start_date = ADTM_next,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  mutate(
    ARRLT = case_when(
      EVID == 1 ~ 0,
      is.na(ARRLT) ~ AXRLT,
      TRUE ~ ARRLT
    ),
    # Derive Reference Dose Date
    PCRFTDTM = case_when(
      EVID == 1 ~ ADTM,
      is.na(ADTM_prev) ~ ADTM_next,
      TRUE ~ ADTM_prev
    )
  ) %>%
  # Derive dates and times from datetimes
  derive_vars_dtm_to_dt(exprs(FANLDTM)) %>%
  derive_vars_dtm_to_tm(exprs(FANLDTM)) %>%
  derive_vars_dtm_to_dt(exprs(PCRFTDTM)) %>%
  derive_vars_dtm_to_tm(exprs(PCRFTDTM))

reactable(adpc_arrlt)

## ----r------------------------------------------------------------------------
# Derive Nominal Relative Time from Reference Dose (NRRLT)

adpc_nrrlt <- adpc_arrlt %>%
  mutate(
    NRRLT = case_when(
      EVID == 1 ~ 0,
      is.na(NFRLT_prev) ~ NFRLT - min_NFRLT,
      TRUE ~ NFRLT - NFRLT_prev
    ),
    NXRLT = case_when(
      EVID == 1 ~ 0,
      TRUE ~ NFRLT - NFRLT_next
    )
  )

## ----r------------------------------------------------------------------------
adpc_aval <- adpc_nrrlt %>%
  mutate(
    PARCAT1 = PCSPEC,
    ATPTN = case_when(
      EVID == 1 ~ 0,
      TRUE ~ PCTPTNUM
    ),
    ATPT = case_when(
      EVID == 1 ~ "Dose",
      TRUE ~ PCTPT
    ),
    ATPTREF = case_when(
      EVID == 1 ~ AVISIT,
      is.na(AVISIT_prev) ~ AVISIT_next,
      TRUE ~ AVISIT_prev
    ),
    # Derive baseline flag for pre-dose records
    ABLFL = case_when(
      ATPT == "Pre-dose" ~ "Y",
      TRUE ~ NA_character_
    ),
    # Derive BASETYPE
    BASETYPE = paste(ATPTREF, "Baseline"),

    # Derive Actual Dose
    DOSEA = case_when(
      EVID == 1 ~ EXDOSE,
      is.na(EXDOSE_prev) ~ EXDOSE_next,
      TRUE ~ EXDOSE_prev
    ),
    # Derive Planned Dose
    DOSEP = case_when(
      TRT01P == "Xanomeline High Dose" ~ 81,
      TRT01P == "Xanomeline Low Dose" ~ 54
    ),
    DOSEU = "mg",
  ) %>%
  # Derive relative time units
  mutate(
    FRLTU = "h",
    RRLTU = "h",
    # Derive PARAMCD
    PARAMCD = coalesce(PCTESTCD, "DOSE"),
    ALLOQ = PCLLOQ,
    # Derive AVAL
    AVAL = case_when(
      EVID == 1 ~ EXDOSE,
      PCSTRESC == "<BLQ" & NFRLT == 0 ~ 0,
      PCSTRESC == "<BLQ" & NFRLT > 0 ~ 0.5 * ALLOQ,
      TRUE ~ PCSTRESN
    ),
    AVALU = case_when(
      EVID == 1 ~ EXDOSU,
      TRUE ~ PCSTRESU
    ),
    AVALCAT1 = if_else(PCSTRESC == "<BLQ", PCSTRESC, prettyNum(signif(AVAL, digits = 3))),
  ) %>%
  # Add SRCSEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = coalesce(PCSEQ, EXSEQ)
  )

reactable(adpc_aval)

## ----r------------------------------------------------------------------------
dtype <- adpc_aval %>%
  filter(NFRLT > 0 & NXRLT == 0 & EVID == 0 & !is.na(AVISIT_next)) %>%
  select(-PCRFTDT, -PCRFTTM) %>%
  # Re-derive variables in for DTYPE copy records
  mutate(
    ABLFL = NA_character_,
    ATPTREF = AVISIT_next,
    ARRLT = AXRLT,
    NRRLT = NXRLT,
    PCRFTDTM = ADTM_next,
    DOSEA = EXDOSE_next,
    BASETYPE = paste(AVISIT_next, "Baseline"),
    ATPT = "Pre-dose",
    ATPTN = NFRLT,
    ABLFL = "Y",
    DTYPE = "COPY"
  ) %>%
  derive_vars_dtm_to_dt(exprs(PCRFTDTM)) %>%
  derive_vars_dtm_to_tm(exprs(PCRFTDTM))

reactable(dtype)

## ----r------------------------------------------------------------------------
adpc_dtype <- bind_rows(adpc_aval, dtype) %>%
  arrange(STUDYID, USUBJID, BASETYPE, ADTM, NFRLT) %>%
  mutate(
    # Derive MRRLT, ANL01FL and ANL02FL
    MRRLT = if_else(ARRLT < 0, 0, ARRLT),
    ANL01FL = "Y",
    ANL02FL = if_else(is.na(DTYPE), "Y", NA_character_),
  )

reactable(adpc_dtype)

## ----r------------------------------------------------------------------------
# ---- Derive BASE and Calculate Change from Baseline ----

adpc_base <- adpc_dtype %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, PARCAT1, BASETYPE),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  )

adpc_chg <- derive_var_chg(adpc_base)

## ----r------------------------------------------------------------------------
# ---- Add ASEQ ----

adpc_aseq <- adpc_chg %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(ADTM, BASETYPE, EVID, AVISITN, ATPTN, PARCAT1, DTYPE),
    check_type = "error"
  ) %>%
  # Derive PARAM and PARAMN using metatools
  create_var_from_codelist(metacore, input_var = PARAMCD, out_var = PARAM) %>%
  create_var_from_codelist(metacore, input_var = PARAMCD, out_var = PARAMN)

## ----r------------------------------------------------------------------------
#---- Derive additional baselines from VS ----

adpc_baselines <- adpc_aseq %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "HEIGHT",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(HTBL = VSSTRESN, HTBLU = VSSTRESU)
  ) %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WTBL = VSSTRESN, WTBLU = VSSTRESU)
  ) %>%
  mutate(
    BMIBL = compute_bmi(height = HTBL, weight = WTBL),
    BMIBLU = "kg/m^2"
  )

reactable(adpc_baselines)

## ----r------------------------------------------------------------------------
# ---- Add all ADSL variables ----

# Add all ADSL variables
adpc_prefinal <- adpc_baselines %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r------------------------------------------------------------------------
# Apply metadata and perform associated checks ----
adpc <- adpc_prefinal %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore) %>% # Check all variables specified are present and no more
  check_ct_data(metacore) %>% # Checks all variables with CT only contain values within the CT
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

reactable(adpc)

## ----r------------------------------------------------------------------------
dir <- tempdir() # Change to whichever directory you want to save the dataset in

adpc_xpt <- adpc %>%
  xportr_type(metacore, domain = "ADPC") %>% # Coerce variable type to match spec
  xportr_length(metacore) %>% # Assigns SAS length from a variable level metadata
  xportr_label(metacore) %>% # Assigns variable label from metacore specifications
  xportr_format(metacore) %>% # Assigns variable format from metacore specifications
  xportr_df_label(metacore) %>% # Assigns dataset label from metacore specifications
  xportr_write(file.path(dir, "adpc.xpt")) # Write xpt v5 transport file

