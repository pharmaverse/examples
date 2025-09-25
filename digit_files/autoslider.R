## ----r setup------------------------------------------------------------------
library(autoslider.core)
library(dplyr)

# 1. Load ALL necessary packages
library(rtables) # For append_topleft()
library(assertthat) # For assert_that() you had issues with before
library(tern)

# define path to the yml files
spec_file <- file.path("metadata/autoslideR_spec.yml")

filters_file <- file.path("metadata/autoslideR_filters.yml")
# load all filters
filters::load_filters(filters_file, overwrite = TRUE)

## ----r, table-----------------------------------------------------------------
# read data
data <- list(
  "adsl" = pharmaverseadam::adsl %>%
    mutate(
      FASFL = SAFFL, # add FASFL for illustrative purpose for t_pop_slide
      # DISTRTFL DCSREAS is needed for t_ds_slide but is missing in example data
      DISTRTFL = sample(c("Y", "N"), size = length(TRT01A), replace = TRUE, prob = c(.1, .9)),
      DCSREAS = sample(
        c(
          "DEATH", "LACK OF EFFICACY", "PROTOCOL VIOLATION",
          "WITHDRAWAL BY PARENT/GUARDIAN", "ADVERSE EVENT"
        ),
        size = length(TRT01A), replace = TRUE, prob = c(.1, .2, .3, .3, .1)
      )
    ),
  "adae" = pharmaverseadam::adae,
  "adtte" = pharmaverseadam::adtte_onco,
  "adrs" = pharmaverseadam::adrs_onco,
  "adlb" = pharmaverseadam::adlb
)

# create outputs based on the specs and the functions
outputs <- spec_file %>%
  read_spec() %>%
  # we can also filter for specific programs:
  filter_spec(., program %in% c("t_dm_slide")) %>%
  # these filtered specs are now piped into the generate_outputs function.
  # this function also requires the data
  generate_outputs(datasets = data) %>%
  # now we decorate based on the specs, i.e. add footnotes and titles
  decorate_outputs(
    version_label = NULL
  )

## ----r------------------------------------------------------------------------
outputs$t_dm_slide_SE

## ----r------------------------------------------------------------------------
outputs %>%
  generate_slides(
    outfile = "presentation.pptx",
    template = file.path(system.file(package = "autoslider.core"), "/theme/basic.pptx"),
    table_format = autoslider_format
  )
