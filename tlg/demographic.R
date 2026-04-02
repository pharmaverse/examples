## ----r preproc----------------------------------------------------------------
library(dplyr)

# Create categorical variables, remove screen failures, and assign column labels
adsl <- pharmaverseadam::adsl |>
  filter(!ACTARM %in% "Screen Failure") |>
  mutate(
    SEX = case_match(SEX, "M" ~ "MALE", "F" ~ "FEMALE"),
    AGEGR1 =
      case_when(
        between(AGE, 18, 40) ~ "18-40",
        between(AGE, 41, 64) ~ "41-64",
        AGE > 64 ~ ">=65"
      ) |>
        factor(levels = c("18-40", "41-64", ">=65"))
  ) |>
  labelled::set_variable_labels(
    AGE = "Age (yr)",
    AGEGR1 = "Age group",
    SEX = "Sex",
    RACE = "Race"
  )

## ----r gtsummary-table--------------------------------------------------------
library(cards)
library(gtsummary)
theme_gtsummary_compact() # reduce default padding and font size for a gt table

# build the ARD with the needed summary statistics using {cards}
ard <-
  ard_stack(
    adsl,
    ard_continuous(variables = AGE),
    ard_categorical(variables = c(AGEGR1, SEX, RACE)),
    .by = ACTARM, # split results by treatment arm
    .attributes = TRUE # optionally include column labels in the ARD
  )

# use the ARD to create a demographics table using {gtsummary}
tbl_ard_summary(
  cards = ard,
  by = ACTARM,
  include = c(AGE, AGEGR1, SEX, RACE),
  type = AGE ~ "continuous2",
  statistic = AGE ~ c("{N}", "{mean} ({sd})", "{median} ({p25}, {p75})", "{min}, {max}")
) |>
  bold_labels() |>
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n}") |> # add Ns to header
  modify_footnote(everything() ~ NA) # remove default footnote

## ----r gtsummary-ard----------------------------------------------------------
# build demographics table directly from a data frame
tbl <- adsl |> tbl_summary(by = ACTARM, include = c(AGE, AGEGR1, SEX, RACE))

# extract ARD from table object
gather_ard(tbl)[[1]] |> select(-gts_column) # removing column so ARD fits on page

## ----r rtables-setup, message=FALSE, warning=FALSE, results='hold'------------
library(tern)

adsl2 <- adsl |>
  df_explicit_na()

## ----r rtables-table----------------------------------------------------------
vars <- c("AGE", "AGEGR1", "SEX", "RACE")
var_labels <- c(
  "Age (yr)",
  "Age group",
  "Sex",
  "Race"
)

lyt <- basic_table(show_colcounts = TRUE) |>
  split_cols_by(var = "ACTARM") |>
  add_overall_col("All Patients") |>
  analyze_vars(
    vars = vars,
    var_labels = var_labels
  )

result <- build_table(lyt, adsl2)

result

## ----r tfrmt-table------------------------------------------------------------
library(cards)
library(forcats)
library(tfrmt)

# build the ARD with the needed summary statistics using {cards}
ard <-
  ard_stack(
    adsl,
    ard_continuous(
      variables = AGE,
      statistic = ~ continuous_summary_fns(c("N", "mean", "sd", "min", "max"))
    ),
    ard_categorical(variables = c(AGEGR1, SEX, RACE)),
    .by = ACTARM, # split results by treatment arm
    .overall = TRUE,
    .total_n = TRUE
  )

# tidy the ARD for use in {tfrmt}
ard_tbl <-
  ard |>
  # reshape the data
  shuffle_card(fill_overall = "Total") |>
  # transform group-level freqs/pcts into a singular "bigN" row
  prep_big_n(vars = "ACTARM") |>
  # consolidate vars into a single variable column
  prep_combine_vars(vars = c("AGE", "AGEGR1", "SEX", "RACE")) |>
  # coalesce categorical levels + continuous stats into a "label"
  prep_label() |>
  group_by(ACTARM, stat_variable) |>
  mutate(across(c(variable_level, label), ~ ifelse(stat_name == "N", "n", .x))) |>
  ungroup() |>
  unique() |>
  # sorting
  mutate(
    ord1 = fct_inorder(stat_variable) |> fct_relevel("SEX", after = 0) |> as.numeric(),
    ord2 = ifelse(label == "n", 1, 2)
  ) |>
  # relabel the variables
  mutate(stat_variable = case_when(
    stat_variable == "AGE" ~ "Age (YEARS) at First Dose",
    stat_variable == "AGEGR1" ~ "Age Group (YEARS) at First Dose",
    stat_variable == "SEX" ~ "Sex",
    stat_variable == "RACE" ~ "High Level Race",
    .default = stat_variable
  )) |>
  # drop variables not needed
  select(ACTARM, stat_variable, label, stat_name, stat, ord1, ord2) |>
  # remove duplicates (extra denominators per variable level)
  unique()

# create a demographics table using {tfrmt}
DM_T01 <- tfrmt(
  group = stat_variable,
  label = label,
  param = stat_name,
  value = stat,
  column = ACTARM,
  sorting_cols = c(ord1, ord2),
  body_plan = body_plan(
    frmt_structure(group_val = ".default", label_val = ".default", frmt("xxx")),
    frmt_structure(
      group_val = ".default", label_val = ".default",
      frmt_combine("{n} ({p}%)",
        n = frmt("xxx"),
        p = frmt("xx", transform = ~ . * 100)
      )
    )
  ),
  big_n = big_n_structure(param_val = "bigN", n_frmt = frmt(" (N=xx)")),
  col_plan = col_plan(
    -starts_with("ord")
  ),
  col_style_plan = col_style_plan(
    col_style_structure(col = c("Placebo", "Xanomeline High Dose", "Xanomeline Low Dose", "Total"), align = "left")
  ),
  row_grp_plan = row_grp_plan(
    row_grp_structure(group_val = ".default", element_block(post_space = " "))
  )
) |>
  print_to_gt(ard_tbl)

DM_T01
