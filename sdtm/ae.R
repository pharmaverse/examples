## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(sdtm.oak)
library(pharmaverseraw)
library(dplyr)

ae_raw <- pharmaverseraw::ae_raw

## ----r eval=TRUE, echo=FALSE--------------------------------------------------
dm <- pharmaversesdtm::dm

## ----r------------------------------------------------------------------------
ae_raw <- ae_raw %>%
  generate_oak_id_vars(
    pat_var = "PATNUM",
    raw_src = "ae_raw"
  )

## ----r, echo = TRUE-----------------------------------------------------------
# TODO: update file path when the file is moved to the correct folder
study_ct <- read.csv("metadata/sdtm_ct.csv")

## ----r------------------------------------------------------------------------
ae <-
  # Derive topic variable
  # Map AETERM using assign_no_ct, raw_var=IT.AETERM, tgt_var=AETERM
  assign_no_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AETERM",
    tgt_var = "AETERM",
    id_vars = oak_id_vars()
  )

## ----r eval=TRUE--------------------------------------------------------------
ae <- ae %>%
  # Map AEOUT using assign_ct, raw_var=AEOUTCOME, tgt_var=AEOUT
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "AEOUTCOME",
    tgt_var = "AEOUT",
    ct_spec = study_ct,
    ct_clst = "C66768",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESEV using assign_no_ct, raw_var=IT.AESEV, tgt_var=AESEV
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AESEV",
    tgt_var = "AESEV",
    ct_spec = study_ct,
    ct_clst = "C66769",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESER using assign_no_ct, raw_var=IT.AESER, tgt_var=AESER
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AESER",
    tgt_var = "AESER",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AEACN using assign_no_ct, raw_var=IT.AEACN, tgt_var=AEACN
  assign_no_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AEACN",
    tgt_var = "AEACN",
    id_vars = oak_id_vars()
  ) %>%
  # Map AEREL using assign_ct, raw_var=IT.AEREL, tgt_var=AEREL
  # User-added codelist is in the ct,
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AEREL",
    tgt_var = "AEREL",
    ct_spec = study_ct,
    ct_clst = "AEREL",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESCAN using assign_ct, raw_var=AESCAN, tgt_var=AESCAN
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "AESCAN",
    tgt_var = "AESCAN",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESCNO using assign_ct, raw_var=AESCNO, tgt_var=AESCNO
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "AESCNO",
    tgt_var = "AESCONG",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AEDIS using assign_ct, raw_var=AEDIS, tgt_var=AEDIS
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "AEDIS",
    tgt_var = "AESDISAB",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESDTH using assign_ct, raw_var=IT.AESDTH, tgt_var=AESDTH
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AESDTH",
    tgt_var = "AESDTH",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESHOSP using assign_ct, raw_var=IT.AESHOSP, tgt_var=AESHOSP
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AESHOSP",
    tgt_var = "AESHOSP",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESLIFE using assign_ct, raw_var=IT.AESLIFE, tgt_var=AESLIFE
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "IT.AESLIFE",
    tgt_var = "AESLIFE",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AESOD using assign_ct, raw_var=AESOD, tgt_var=AESOD
  assign_ct(
    raw_dat = ae_raw,
    raw_var = "AESOD",
    tgt_var = "AESOD",
    ct_spec = study_ct,
    ct_clst = "C66742",
    id_vars = oak_id_vars()
  ) %>%
  # Map AEDTC using assign_datetime, raw_var=AEDTCOL
  assign_datetime(
    raw_dat = ae_raw,
    raw_var = "AEDTCOL",
    tgt_var = "AEDTC",
    raw_fmt = c("m/d/y")
  ) %>%
  # Map AESTDTC using assign_datetime, raw_var=IT.AESTDAT
  assign_datetime(
    raw_dat = ae_raw,
    raw_var = "IT.AESTDAT",
    tgt_var = "AESTDTC",
    raw_fmt = c("m/d/y"),
    id_vars = oak_id_vars()
  ) %>%
  # Map AEENDTC using assign_datetime, raw_var=IT.AEENDAT
  assign_datetime(
    raw_dat = ae_raw,
    raw_var = "IT.AEENDAT",
    tgt_var = "AEENDTC",
    raw_fmt = c("m/d/y"),
    id_vars = oak_id_vars()
  )

## ----r------------------------------------------------------------------------
ae <- ae %>%
  dplyr::mutate(
    STUDYID = ae_raw$STUDY,
    DOMAIN = "AE",
    USUBJID = paste0("01-", ae_raw$PATNUM),
    AELLT = ae_raw$AELLT,
    AELLTCD = ae_raw$AELLTCD,
    AEDECOD = ae_raw$AEDECOD,
    AEPTCD = ae_raw$AEPTCD,
    AEHLT = ae_raw$AEHLT,
    AEHLTCD = ae_raw$AEHLTCD,
    AEHLGT = ae_raw$AEHLGT,
    AEHLGTCD = ae_raw$AEHLGTCD,
    AEBODSYS = ae_raw$AEBODSYS,
    AEBDSYCD = ae_raw$AEBDSYCD,
    AESOC = ae_raw$AESOC,
    AESOCCD = ae_raw$AESOCCD,
    AETERM = toupper(AETERM)
  ) %>%
  derive_seq(
    tgt_var = "AESEQ",
    rec_vars = c("USUBJID", "AETERM")
  ) %>%
  derive_study_day(
    sdtm_in = .,
    dm_domain = dm,
    tgdt = "AESTDTC",
    refdt = "RFXSTDTC",
    study_day_var = "AESTDY"
  ) %>%
  derive_study_day(
    sdtm_in = .,
    dm_domain = dm,
    tgdt = "AEENDTC",
    refdt = "RFXENDTC",
    study_day_var = "AEENDY"
  ) %>%
  select(
    "STUDYID", "DOMAIN", "USUBJID", "AESEQ", "AETERM", "AELLT", "AELLTCD", "AEDECOD", "AEPTCD", "AEHLT", "AEHLTCD", "AEHLGT",
    "AEHLGTCD", "AEBODSYS", "AEBDSYCD", "AESOC", "AESOCCD", "AESEV", "AESER", "AEACN", "AEREL", "AEOUT", "AESCAN", "AESCONG",
    "AESDISAB", "AESDTH", "AESHOSP", "AESLIFE", "AESOD", "AEDTC", "AESTDTC", "AEENDTC", "AESTDY", "AEENDY"
  )
