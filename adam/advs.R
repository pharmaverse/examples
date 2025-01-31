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

# Read in input data
adsl <- pharmaverseadam::adsl
vs <- pharmaversesdtm::vs

vs <- convert_blanks_to_na(vs)

## ----r echo=TRUE--------------------------------------------------------------
# ---- Load Specs for Metacore ----
metacore <- spec_to_metacore(
  path = "./metadata/safety_specs.xlsx",
  where_sep_sheet = FALSE,
  quiet = TRUE
) %>%
  select_dataset("ADVS")

## ----r------------------------------------------------------------------------
# Select required ADSL variables
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

# Join ADSL variables with VS
advs <- vs %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r------------------------------------------------------------------------
# Calculate ADT, ADY
advs <- advs %>%
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC,
    # Below arguments are default values and not necessary to add in our case
    highest_imputation = "n", # means no imputation is performed on partial/missing dates
    flag_imputation = "auto" # To automatically create ADTF variable when highest_imputation is "Y", "M" or "D"
  ) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ADT)
  )

## ----r eval=TRUE, include=FALSE-----------------------------------------------
param_lookup <- tibble::tribble(
  ~VSTESTCD, ~PARAMCD,                            ~PARAM, ~PARAMN,
  "SYSBP",    "SYSBP", " Systolic Blood Pressure (mmHg)",        1,
  "DIABP",    "DIABP", "Diastolic Blood Pressure (mmHg)",        2,
  "PULSE",    "PULSE",          "Pulse Rate (beats/min)",        3,
  "WEIGHT",  "WEIGHT",                     "Weight (kg)",        4,
  "HEIGHT", "HEIGHT",                      "Height (cm)",        5,
  "TEMP",     "TEMP",                  "Temperature (C)",        6,
  "MAP",       "MAP",    "Mean Arterial Pressure (mmHg)",        7,
  "BMI",       "BMI",          "Body Mass Index(kg/m^2)",        8,
  "BSA",       "BSA",           "Body Surface Area(m^2)",        9
)
attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"

## ----r------------------------------------------------------------------------
advs <- advs %>%
  # Add PARAMCD only - add PARAM etc later
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = exprs(PARAMCD),
    by_vars = exprs(VSTESTCD),
    # Below arguments are default values and not necessary to add in our case
    print_not_mapped = TRUE # Printing whether some parameters are not mapped
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>% 
  mutate(
    AVAL = VSSTRESN,
    AVALU = VSSTRESU
  ) 

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>% 
  derive_param_map(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM, AVALU), # Other variables than the defined ones here won't be populated
    set_values_to = exprs(PARAMCD = "MAP"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT),
    # Below arguments are default values and not necessary to add in our case
    sysbp_code = "SYSBP",
    diabp_code = "DIABP",
    hr_code = NULL
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>% 
  derive_param_computed(
    by_vars = exprs(STUDYID, USUBJID, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    parameters = "WEIGHT",
    set_values_to = exprs(
      AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      PARAMCD = "BMI",
      AVALU = "kg/m^2"
    ),
    constant_parameters = c("HEIGHT"),
    constant_by_vars = exprs(USUBJID)
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>% 
  derive_param_bsa(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    method = "Mosteller",
    set_values_to = exprs(
      PARAMCD = "BSA",
      AVALU = "m^2"
    ),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT),
    constant_by_vars = exprs(USUBJID),
    # Below arguments are default values and not necessary to add in our case
    height_code = "HEIGHT",
    weight_code = "WEIGHT"
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>%
  mutate(
    ATPTN = VSTPTNUM,
    ATPT = VSTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(case_when(
      VISIT == "BASELINE" ~ "0",
      str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", "")),
      TRUE ~ NA_character_
    ))
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- derive_summary_records(
  dataset = advs,
  dataset_add = advs, # Observations from the specified dataset are going to be used to calculate and added as new records to the input dataset.
  by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT, ADY, AVALU),
  filter_add = !is.na(AVAL),
  set_values_to = exprs(
    AVAL = mean(AVAL),
    DTYPE = "AVERAGE"
  )
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- derive_var_ontrtfl(
  advs,
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  filter_pre_timepoint = toupper(AVISIT) == "BASELINE" # Observations as not on-treatment
)

## ----r include=FALSE----------------------------------------------------------
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "SYSBP",      90,    130,    70,   140,
  "DIABP",      60,     80,    40,    90,
  "PULSE",      60,    100,    40,   110,
  "TEMP",     36.5,   37.5,    35,    38
)

advs <- derive_vars_merged(
  advs,
  dataset_add = range_lookup,
  by_vars = exprs(PARAMCD)
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- derive_var_anrind(
  advs,
  # Below arguments are default values and not necessary to add in our case
  signif_dig = get_admiral_option("signif_digits"),
  use_a1hia1lo = FALSE
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- derive_basetype_records(
  dataset = advs,
  basetypes = exprs(
    "LAST: AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
    "LAST: AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
    "LAST: AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
    "LAST" = is.na(ATPTN)
  )
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- restrict_derivation(
  advs,
  derivation = derive_var_extreme_flag,
  args = params(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = exprs(ADT, VISITNUM, VSSEQ),
    new_var = ABLFL,
    mode = "last", # Determines of the first or last observation is flagged
    # Below arguments are default values and not necessary to add in our case
    true_value = "Y"
  ),
  filter = (!is.na(AVAL) &
              ADT <= TRTSDT & !is.na(BASETYPE) & is.na(DTYPE)
  )
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- derive_var_base(
  advs,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
  source_var = AVAL,
  new_var = BASE,
  # Below arguments are default values and not necessary to add in our case
  filter = ABLFL == "Y"
) 

advs <- derive_var_base(
  advs,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
  source_var = ANRIND,
  new_var = BNRIND
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- restrict_derivation(
  advs,
  derivation = derive_var_chg,
  filter = AVISITN > 0
)

advs <- restrict_derivation(
  advs,
  derivation = derive_var_pchg,
  filter = AVISITN > 0
)

## ----r eval=TRUE--------------------------------------------------------------
advs <-   restrict_derivation(
  advs,
  derivation = derive_var_extreme_flag,
  args = params(
    new_var = ANL01FL,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
    order = exprs(ADT, AVAL),
    mode = "last", # Determines of the first or last observation is flagged - As seen while deriving ABLFL
    # Below arguments are default values and not necessary to add in our case
    true_value = "Y"
  ),
  filter = !is.na(AVISITN) & ONTRTFL == "Y"
)

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>% 
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- derive_var_obs_number(
  advs,
  new_var = ASEQ,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
  check_type = "error"
)

## ----r eval=TRUE--------------------------------------------------------------
avalcat_lookup <- exprs(
  ~PARAMCD,  ~condition,   ~AVALCAT1, ~AVALCA1N,
  "HEIGHT",  AVAL > 140,   ">140 cm",         1,
  "HEIGHT", AVAL <= 140, "<= 140 cm",         2
)

advs <- advs %>%
  derive_vars_cat(
    definition = avalcat_lookup,
    by_vars = exprs(PARAMCD)
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>% 
  create_var_from_codelist(
    metacore,
    input_var = PARAMCD,
    out_var = PARAM,
    decode_to_code = FALSE # input_var is the code column of the codelist
  ) %>%
  create_var_from_codelist(
    metacore,
    input_var = PARAMCD,
    out_var = PARAMN
  )

## ----r eval=TRUE--------------------------------------------------------------
advs <- advs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

## ----r, message=FALSE, warning=FALSE------------------------------------------
dir <- tempdir() # Specify the directory for saving the XPT file

# Apply metadata and perform checks
advs_prefinal <- advs %>%
  drop_unspec_vars(metacore) %>% # Drop unspecified variables from specs
  check_variables(metacore, dataset_name = "ADVS") %>% # Check all variables specified are present and no more
  order_cols(metacore) %>% # Orders the columns according to the spec
  sort_by_key(metacore) # Sorts the rows by the sort keys

# Apply apply labels, formats, and export the dataset to an XPT file.
advs_final <- advs_prefinal %>% 
  xportr_type(metacore) %>%
  xportr_length(metacore) %>%
  xportr_label(metacore) %>%
  xportr_format(metacore, domain = "ADVS") %>%
  xportr_df_label(metacore, domain = "ADVS") %>% 
  xportr_write(file.path(dir, "advs.xpt"), metadata = metacore, domain = "ADVS")
