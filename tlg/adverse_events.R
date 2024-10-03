## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(pharmaverseadam)
library(tern)
library(dplyr)

adsl <- adsl %>%
  df_explicit_na()

adae <- adae %>%
  df_explicit_na()

## ----r preproc----------------------------------------------------------------
adae <- adae %>%
  var_relabel(
    AEBODSYS = "MedDRA System Organ Class",
    AEDECOD = "MedDRA Preferred Term"
  ) %>%
  filter(SAFFL == "Y")

# Define the split function
split_fun <- drop_split_levels

## ----r table------------------------------------------------------------------
lyt <- basic_table(show_colcounts = TRUE) %>%
  split_cols_by(var = "ACTARM") %>%
  add_overall_col(label = "All Patients") %>%
  analyze_num_patients(
    vars = "USUBJID",
    .stats = c("unique", "nonunique"),
    .labels = c(
      unique = "Total number of patients with at least one adverse event",
      nonunique = "Overall total number of events"
    )
  ) %>%
  split_rows_by(
    "AEBODSYS",
    child_labels = "visible",
    nested = FALSE,
    split_fun = split_fun,
    label_pos = "topleft",
    split_label = obj_label(adae$AEBODSYS)
  ) %>%
  summarize_num_patients(
    var = "USUBJID",
    .stats = c("unique", "nonunique"),
    .labels = c(
      unique = "Total number of patients with at least one adverse event",
      nonunique = "Total number of events"
    )
  ) %>%
  count_occurrences(
    vars = "AEDECOD",
    .indent_mods = -1L
  ) %>%
  append_varlabels(adae, "AEDECOD", indent = 1L)

result <- build_table(lyt, df = adae, alt_counts_df = adsl)

result

