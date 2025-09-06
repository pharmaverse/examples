## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(sdtm.oak)
library(pharmaverseraw)
library(dplyr)

# Read in input data
vs_raw <- pharmaverseraw::vs_raw
dm <- pharmaversesdtm::dm

## ----r------------------------------------------------------------------------
vs_raw <- vs_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "vitals"
  )

## ----r, echo = TRUE-----------------------------------------------------------
# TODO: update file path when the file is moved to the correct folder
study_ct <- read.csv("metadata/sdtm_ct.csv")

## ----r------------------------------------------------------------------------
# Map topic variable SYSBP
vs_sysbp <-
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "SYS_BP",
    tgt_var = "VSTESTCD",
    tgt_val = "SYSBP",
    ct_spec = study_ct,
    ct_clst = "C66741"
  ) %>%
  # Filter for records where VSTESTCD is not empty.
  # Only these records need qualifier mappings.
  dplyr::filter(!is.na(.data$VSTESTCD))

## ----r------------------------------------------------------------------------
## Map Rest of the Variables
vs_sysbp <- vs_sysbp %>%
  # Map VSTEST using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "SYS_BP",
    tgt_var = "VSTEST",
    tgt_val = "Systolic Blood Pressure",
    ct_spec = study_ct,
    ct_clst = "C67153",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRES using assign_no_ct algorithm
  assign_no_ct(
    raw_dat = vs_raw,
    raw_var = "SYS_BP",
    tgt_var = "VSORRES",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "SYS_BP",
    tgt_var = "VSORRESU",
    tgt_val = "mmHg",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSPOS using assign_ct algorithm
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "SUBPOS",
    tgt_var = "VSPOS",
    ct_spec = study_ct,
    ct_clst = "C71148",
    id_vars = oak_id_vars()
  )

## ----r------------------------------------------------------------------------
# Repeat Map Topic and Map Rest
vs_diabp <-
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "DIA_BP",
    tgt_var = "VSTESTCD",
    tgt_val = "DIABP",
    ct_spec = study_ct,
    ct_clst = "C66741"
  ) %>%
  dplyr::filter(!is.na(.data$VSTESTCD)) %>%
  # Map VSTEST using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "DIA_BP",
    tgt_var = "VSTEST",
    tgt_val = "Diastolic Blood Pressure",
    ct_spec = study_ct,
    ct_clst = "C67153",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRES using assign_no_ct algorithm
  assign_no_ct(
    raw_dat = vs_raw,
    raw_var = "DIA_BP",
    tgt_var = "VSORRES",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "DIA_BP",
    tgt_var = "VSORRESU",
    tgt_val = "mmHg",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSPOS using assign_ct algorithm
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "SUBPOS",
    tgt_var = "VSPOS",
    ct_spec = study_ct,
    ct_clst = "C71148",
    id_vars = oak_id_vars()
  )

# Map topic variable PULSE and its qualifiers.
vs_pulse <-
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "PULSE",
    tgt_var = "VSTESTCD",
    tgt_val = "PULSE",
    ct_spec = study_ct,
    ct_clst = "C66741"
  ) %>%
  dplyr::filter(!is.na(.data$VSTESTCD)) %>%
  # Map VSTEST using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "PULSE",
    tgt_var = "VSTEST",
    tgt_val = "Pulse Rate",
    ct_spec = study_ct,
    ct_clst = "C67153",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRES using assign_no_ct algorithm
  assign_no_ct(
    raw_dat = vs_raw,
    raw_var = "PULSE",
    tgt_var = "VSORRES",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "PULSE",
    tgt_var = "VSORRESU",
    tgt_val = "beats/min",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSPOS using assign_ct algorithm
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "SUBPOS",
    tgt_var = "VSPOS",
    ct_spec = study_ct,
    ct_clst = "C71148",
    id_vars = oak_id_vars()
  )

# Map topic variable TEMP from raw variable TEMP and its qualifiers.
vs_temp <-
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.TEMP",
    tgt_var = "VSTESTCD",
    tgt_val = "TEMP",
    ct_spec = study_ct,
    ct_clst = "C66741"
  ) %>%
  dplyr::filter(!is.na(.data$VSTESTCD)) %>%
  # Map VSTEST using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.TEMP",
    tgt_var = "VSTEST",
    tgt_val = "Temperature",
    ct_spec = study_ct,
    ct_clst = "C67153",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRES using assign_no_ct algorithm
  assign_no_ct(
    raw_dat = vs_raw,
    raw_var = "IT.TEMP",
    tgt_var = "VSORRES",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.TEMP",
    tgt_var = "VSORRESU",
    tgt_val = "F",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSLOC from TEMPLOC using assign_ct
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "IT.TEMP_LOC",
    tgt_var = "VSLOC",
    ct_spec = study_ct,
    ct_clst = "C74456",
    id_vars = oak_id_vars()
  ) %>%
  # Create VSSTRESC by converting VSORRES from F to C
  mutate(VSSTRESC = as.character(sprintf("%.2f", (as.numeric(VSORRES) - 32) * 5/9))) %>%
  # Map VSSTRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.TEMP",
    tgt_var = "VSSTRESU",
    tgt_val = "C",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  )

# Map topic variable HEIGHT from raw variable IT.HEIGHT_VSORRRES and its qualifiers.
vs_height <-
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.HEIGHT_VSORRES",
    tgt_var = "VSTESTCD",
    tgt_val = "HEIGHT",
    ct_spec = study_ct,
    ct_clst = "C66741"
  ) %>%
  dplyr::filter(!is.na(.data$VSTESTCD)) %>%
  # Map VSTEST using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.HEIGHT_VSORRES",
    tgt_var = "VSTEST",
    tgt_val = "Height",
    ct_spec = study_ct,
    ct_clst = "C67153",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRES using assign_no_ct algorithm
  assign_no_ct(
    raw_dat = vs_raw,
    raw_var = "IT.HEIGHT_VSORRES",
    tgt_var = "VSORRES",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.HEIGHT_VSORRES",
    tgt_var = "VSORRESU",
    tgt_val = "in",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  ) %>%
  # Create VSSRESC by converting VSORRES from in to cm
  mutate(VSSTRESC = as.character(sprintf("%.2f", as.numeric(VSORRES) * 2.54))) %>%
  # Map VSSTRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.HEIGHT_VSORRES",
    tgt_var = "VSSTRESU",
    tgt_val = "cm",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  )

# Map topic variable WEIGHT from raw variable IT.WEIGHT and its qualifiers.
vs_weight <-
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.WEIGHT",
    tgt_var = "VSTESTCD",
    tgt_val = "WEIGHT",
    ct_spec = study_ct,
    ct_clst = "C66741"
  ) %>%
  dplyr::filter(!is.na(.data$VSTESTCD)) %>%
  # Map VSTEST using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.WEIGHT",
    tgt_var = "VSTEST",
    tgt_val = "Weight",
    ct_spec = study_ct,
    ct_clst = "C67153",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRES using assign_no_ct algorithm
  assign_no_ct(
    raw_dat = vs_raw,
    raw_var = "IT.WEIGHT",
    tgt_var = "VSORRES",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSORRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.WEIGHT",
    tgt_var = "VSORRESU",
    tgt_val = "LB",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  ) %>%
  # Create VSSTRESC by converting VSORRES from LB to KG
  mutate(VSSTRESC = as.character(sprintf("%.2f", as.numeric(VSORRES) / 2.20462))) %>%
  # Map VSSTRESU using hardcode_ct algorithm
  hardcode_ct(
    raw_dat = vs_raw,
    raw_var = "IT.WEIGHT",
    tgt_var = "VSSTRESU",
    tgt_val = "kg",
    ct_spec = study_ct,
    ct_clst = "C66770",
    id_vars = oak_id_vars()
  )

## ----r------------------------------------------------------------------------
# Combine all the topic variables into a single data frame.
vs_combined <- dplyr::bind_rows(
  vs_diabp, vs_height, vs_pulse,
  vs_sysbp, vs_temp, vs_weight
)

## ----r------------------------------------------------------------------------
# Map qualifiers common to all topic variables
vs <- vs_combined %>%
  # Map VSDTC using assign_ct algorithm
  assign_datetime(
    raw_dat = vs_raw,
    raw_var = c("VTLD"),
    tgt_var = "VSDTC",
    raw_fmt = c(list(c("d-m-y", "dd-mmm-yyyy")))
  ) %>%
  # Map VSTPT from TMPTC using assign_ct
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "TMPTC",
    tgt_var = "VSTPT",
    ct_spec = study_ct,
    ct_clst = "TPT",
    id_vars = oak_id_vars()
  ) %>%
  # Map VSTPTNUM from TMPTC using assign_ct
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "TMPTC",
    tgt_var = "VSTPTNUM",
    ct_spec = study_ct,
    ct_clst = "TPTNUM",
    id_vars = oak_id_vars()
  ) %>%
  # Map VISIT from INSTANCE using assign_ct
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISIT",
    ct_spec = study_ct,
    ct_clst = "VISIT",
    id_vars = oak_id_vars()
  ) %>%
  # Map VISITNUM from INSTANCE using assign_ct
  assign_ct(
    raw_dat = vs_raw,
    raw_var = "INSTANCE",
    tgt_var = "VISITNUM",
    ct_spec = study_ct,
    ct_clst = "VISITNUM",
    id_vars = oak_id_vars()
  )

## ----r------------------------------------------------------------------------
vs <- vs %>%
  dplyr::mutate(
    STUDYID = "CDISCPILOT01",
    DOMAIN = "VS",
    VSCAT = "VITAL SIGNS",
    USUBJID = paste0("01", "-", .data$patient_number),
    VSSTRESC = ifelse(is.na(VSSTRESC), VSORRES, VSSTRESC),
    VSSTRESN = as.numeric(VSSTRESC),
    VSSTRESU = ifelse(is.na(VSSTRESU), VSORRESU, VSSTRESU),
    VSELTM = ifelse(is.na(VSTPT), NA, paste0("PT", readr::parse_number(VSTPT), "M")),
    VSTPTREF = ifelse(is.na(VSPOS), NA, paste("PATIENT", VSPOS))
  ) %>%
  arrange(USUBJID, VSTESTCD, as.numeric(VISITNUM), as.numeric(VSTPTNUM)) %>%
  derive_seq(tgt_var = "VSSEQ",
             rec_vars = c("USUBJID", "VSTESTCD")) %>%
  derive_study_day(
    sdtm_in = .,
    dm_domain = dm,
    tgdt = "VSDTC",
    refdt = "RFXSTDTC",
    study_day_var = "VSDY"
  ) %>%
  dplyr::select("STUDYID", "DOMAIN", "USUBJID", "VSSEQ", "VSTESTCD", "VSTEST", "VSPOS", "VSORRES", "VSORRESU", "VSSTRESC", "VSSTRESN", "VSSTRESU", "VSLOC", "VISITNUM", "VISIT", "VSDTC", "VSDY", "VSTPT", "VSTPTNUM", "VSELTM", "VSTPTREF")

