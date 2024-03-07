## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(pharmaverseadam)
library(tern)
library(dplyr)

adsl <- adsl %>%
  df_explicit_na()

## ----r preproc----------------------------------------------------------------
# Create categorical variables
adsl <- adsl %>%
  mutate(
    SEX = factor(case_when(
      SEX == "M" ~ "Male",
      SEX == "F" ~ "Female",
      SEX == "U" ~ "Unknown",
      SEX == "UNDIFFERENTIATED" ~ "Undifferentiated"
    )),
    AGEGR1 = factor(
      case_when(
        between(AGE, 18, 40) ~ "18-40",
        between(AGE, 41, 64) ~ "41-64",
        AGE > 64 ~ ">=65"
      ),
      levels = c("18-40", "41-64", ">=65")
    )
  )

## ----r table------------------------------------------------------------------
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

result <- build_table(lyt, adsl)

result

