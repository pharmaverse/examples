## ----r preproc----------------------------------------------------------------
library(dplyr)

# Create categorical variables, remove scren failures, and assign column labels
adsl <- pharmaverseadam::adsl |>
  filter(!ACTARM %in% "Screen Failure") |>
  mutate(
    SEX = case_when(
      SEX == "M" ~ "Male",
      SEX == "F" ~ "Female"
    ) |>
      factor(),
    AGEGR1 =
      case_when(
        between(AGE, 18, 40) ~ "18-40",
        between(AGE, 41, 64) ~ "41-64",
        AGE > 64 ~ ">=65"
      ) |>
      factor(levels = c("18-40", "41-64", ">=65")
      )
  ) |> 
  labelled::set_variable_labels(AGE = "Age (yr)",
                                AGEGR1 = "Age group",
                                SEX = "Sex",
                                RACE = "Race")

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
    .by = ACTARM,      # split results by treatment arm
    .missing = TRUE,   # add information about missingness rates
  )

# use the ARD to create a demographics table
tbl_ard_summary(
  cards = ard, 
  by = ACTARM, 
  include = c(AGE, AGEGR1, SEX, RACE),
  type = AGE ~ "continuous2",
  statistic = AGE ~ c("{mean} ({sd})", 
                      "{median} ({p25}, {p75})", 
                      "{min}, {max}"),
  missing = "always"
) |> 
  bold_labels() |> 
  modify_header(all_stat_cols() ~ "**{level}**  \nN = {n}") |> # add Ns to header
  modify_footnote(everything() ~ NA) # remove default footnote

## ----r gtsummary-ard----------------------------------------------------------
tbl <- adsl |> 
  tbl_summary(
    by = ACTARM, 
    include = c(AGE, AGEGR1, SEX, RACE),
    # display summary stats for AGE on multiple rows
    type = AGE ~ "continuous2",
    statistic = AGE ~ c("{mean} ({sd})", 
                        "{median} ({p25}, {p75})", 
                        "{min}, {max}"),
    missing = "always"
  ) |> 
  bold_labels() |> 
  modify_footnote(everything() ~ NA) # remove default footnote

gather_ard(tbl)

## ----r rtables-setup, message=FALSE, warning=FALSE, results='hold'------------
library(tern)

adsl2 <- adsl %>%
  df_explicit_na()

## ----r rtables-table----------------------------------------------------------
vars <- c("AGE", "AGEGR1", "SEX", "RACE")
var_labels <- c(
  "Age (yr)",
  "Age group",
  "Sex",
  "Race"
)

lyt <- basic_table(show_colcounts = TRUE) %>%
  split_cols_by(var = "ACTARM") %>%
  add_overall_col("All Patients") %>%
  analyze_vars(
    vars = vars,
    var_labels = var_labels
  )

result <- build_table(lyt, adsl2)

result

