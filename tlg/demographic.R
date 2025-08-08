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

