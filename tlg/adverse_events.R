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

## ----r message=FALSE, warning=FALSE-------------------------------------------
library(pharmaverseadam) # for clinical trial data
library(dplyr) # for data manipulation
library(cards) # for creating analysis result displays
library(tfrmt) # for formatting tables in R

## ----r------------------------------------------------------------------------
# Filter to include only subjects marked as part of the safety population
adsl <- pharmaverseadam::adsl |>
  filter(SAFFL == "Y")
# Load adverse event data
adae <- pharmaverseadam::adae |>
  filter(SAFFL == "Y" & TRTEMFL == "Y")

## ----r------------------------------------------------------------------------
# Create an ARD that stacks hierarchical data of adverse events
# Grouping by treatment, system organ class, and preferred term
ae_ard <- ard_stack_hierarchical(
  data = adae,
  by = TRT01A, # Note: by variables must be present in the denominator dataset
  variables = c(AEBODSYS, AETERM),
  statistic = ~ c("n", "p"), # Calculate count and percentage
  denominator = adsl,
  id = USUBJID,
  over_variables = TRUE,
  overall = TRUE
)

# Filter adae and adsl with trt01a set to "Total" and create a new ARD for the total column
adae2 <- adae |>
  mutate(TRT01A = "Total")
adsl2 <- adsl |>
  mutate(TRT01A = "Total")

ae2_ard <- ard_stack_hierarchical(
  data = adae2,
  by = TRT01A, # Note: by variables must be present in the denominator dataset
  variables = c(AEBODSYS, AETERM),
  denominator = adsl2,
  statistic = ~ c("n", "p"),
  id = USUBJID,
  over_variables = TRUE,
  overall = TRUE
) |>
  filter(group2 == "TRT01A" | variable == "TRT01A") # filter to stats we need

## ----r------------------------------------------------------------------------
ae3_ard <- bind_ard(ae_ard, ae2_ard) |>
  # reshape the data
  shuffle_card(fill_hierarchical_overall = "ANY EVENT") |>
  # transform group-level freqs/pcts into a singular "bigN" row
  prep_big_n(vars = "TRT01A") |>
  # for nested variables, fill any missing values with "ANY EVENT"
  prep_hierarchical_fill(vars = c("AEBODSYS", "AETERM"), fill = "ANY EVENT") |>
  mutate(TRT01A = ifelse(TRT01A == "Overall TRT01A", "Total", TRT01A))

# create ordering columns, sort by AEBODSYS
ordering_aebodsys <- ae3_ard |>
  filter(TRT01A == "Total", stat_name == "n", AETERM == "ANY EVENT") |>
  arrange(desc(stat)) |>
  mutate(ord1 = row_number()) |>
  select(AEBODSYS, ord1)

# sort by AETERM after AEBODSYS order
ordering_aeterm <- ae3_ard |>
  filter(TRT01A == "Total", stat_name == "n") |>
  group_by(AEBODSYS) |>
  arrange(desc(stat)) |>
  mutate(ord2 = row_number()) |>
  select(AEBODSYS, AETERM, ord2)

# join on our ordering columns and keep required columns
ae4_ard <- ae3_ard |>
  full_join(ordering_aebodsys, by = "AEBODSYS") |>
  full_join(ordering_aeterm, by = c("AEBODSYS", "AETERM")) |>
  select(AEBODSYS, AETERM, ord1, ord2, stat, stat_name, TRT01A)

## ----r------------------------------------------------------------------------
AE_T01 <- tfrmt_n_pct(
  n = "n", pct = "p",
  pct_frmt_when = frmt_when(
    "==1" ~ frmt("(100%)"),
    ">=0.995" ~ frmt("(>99%)"),
    "==0" ~ frmt(""),
    "<=0.01" ~ frmt("(<1%)"),
    "TRUE" ~ frmt("(xx.x%)", transform = ~ . * 100)
  )
) |>
  tfrmt(
    group = AEBODSYS,
    label = AETERM,
    param = stat_name,
    value = stat,
    column = TRT01A,
    sorting_cols = c(ord1, ord2),
    col_plan = col_plan(
      "System Organ Class
            Preferred Term" = AEBODSYS, Placebo, `Xanomeline High Dose`, `Xanomeline Low Dose`,
      -ord1, -ord2
    ),
    row_grp_plan = row_grp_plan(row_grp_structure(
      group_val = ".default", element_block(post_space = " ")
    )),
    big_n = big_n_structure(param_val = "bigN", n_frmt = frmt(" (N=xx)"))
  ) |>
  print_to_gt(ae4_ard)

AE_T01
