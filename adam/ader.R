## ----r echo=TRUE, message=FALSE-----------------------------------------------
# Load Packages
library(admiral)
library(admiralonco)
# pharmaverseadam contains example datasets generated from the CDISC pilot
# project SDTM ran through admiral templates
library(pharmaverseadam)
library(dplyr)
library(lubridate)
library(stringr)
library(purrr)
library(rlang)
library(metacore)
library(metatools)
library(xportr)

## ----r echo=TRUE, message=FALSE-----------------------------------------------
# ---- Load Specs for Metacore ----
metacore <- spec_to_metacore("./metadata/pk_spec.xlsx", quiet = TRUE) %>%
  select_dataset("ADER")

