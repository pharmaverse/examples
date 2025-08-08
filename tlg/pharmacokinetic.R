## ----r setup, message=FALSE, warning=FALSE, results='hold'--------------------
library(pharmaverseadam)
library(tern)
library(dplyr)
library(ggplot2)
library(nestcolor)
library(rlistings)

# Read data from pharmaverseadam
adpc <- pharmaverseadam::adpc
adsl <- pharmaverseadam::adsl

# Use tern::df_explicit_na() to end encode missing values as categorical
adsl <- adsl %>%
  df_explicit_na()

adpc <- adpc %>%
  df_explicit_na()

# For ADPC keep only concentration records and treated subjects
# Keep only plasma records for this example
# Remove DTYPE = COPY records with ANL02FL == "Y"
adpc <- adpc %>%
  filter(PARAMCD != "DOSE" & TRT01A != "Placebo" & PARCAT1 == "PLASMA" & ANL02FL == "Y")

## ----r table------------------------------------------------------------------
# Setting up the data for table
adpc_t <- adpc %>%
  mutate(
    NFRLT = as.factor(NFRLT),
    AVALCAT1 = as.factor(AVALCAT1),
    NOMTPT = as.factor(paste(NFRLT, "/", PCTPT))
  ) %>%
  select(NOMTPT, ACTARM, VISIT, AVAL, PARAM, AVALCAT1)

adpc_t$NOMTPT <- factor(
  adpc_t$NOMTPT,
  levels = levels(adpc_t$NOMTPT)[order(as.numeric(gsub(".*?([0-9\\.]+).*", "\\1", levels(adpc_t$NOMTPT))))]
)

# Row structure
lyt_rows <- basic_table() %>%
  split_rows_by(
    var = "ACTARM",
    split_fun = drop_split_levels,
    split_label = "Treatment Group",
    label_pos = "topleft"
  ) %>%
  add_rowcounts(alt_counts = TRUE) %>%
  split_rows_by(
    var = "VISIT",
    split_fun = drop_split_levels,
    split_label = "Visit",
    label_pos = "topleft"
  ) %>%
  split_rows_by(
    var = "NOMTPT",
    split_fun = drop_split_levels,
    split_label = "Nominal Time (hr) / Timepoint",
    label_pos = "topleft",
    child_labels = "hidden"
  )

lyt <- lyt_rows %>%
  analyze_vars_in_cols(
    vars = c("AVAL", "AVALCAT1", rep("AVAL", 8)),
    .stats = c("n", "n_blq", "mean", "sd", "cv", "geom_mean", "geom_cv", "median", "min", "max"),
    .formats = c(
      n = "xx.", n_blq = "xx.", mean = format_sigfig(3), sd = format_sigfig(3), cv = "xx.x", median = format_sigfig(3),
      geom_mean = format_sigfig(3), geom_cv = "xx.x", min = format_sigfig(3), max = format_sigfig(3)
    ),
    .labels = c(
      n = "n", n_blq = "Number\nof\nLTRs/BLQs", mean = "Mean", sd = "SD", cv = "CV (%) Mean",
      geom_mean = "Geometric Mean", geom_cv = "CV % Geometric Mean", median = "Median", min = "Minimum", max = "Maximum"
    ),
    na_str = "NE",
    .aligns = "decimal"
  )

result <- build_table(lyt, df = adpc_t, alt_counts_df = adsl) %>% prune_table()

# Decorating
main_title(result) <- "Summary of PK Concentrations by Nominal Time and Treatment: PK Evaluable"
subtitles(result) <- c(
  "Protocol: xxxxx",
  paste("Analyte: ", unique(adpc_t$PARAM)),
  paste("Treatment:", unique(adpc_t$ACTARM))
)
main_footer(result) <- "NE: Not Estimable"

result

## ----r graph------------------------------------------------------------------
# Keep only treated subjects for graph
adsl_f <- adsl %>%
  filter(SAFFL == "Y" & TRT01A != "Placebo")

# Set titles and footnotes
use_title <- "Plot of Mean (+/- SD) Plasma Concentrations Over Time by Treatment, \nPK Evaluable Patients"
use_subtitle <- "Analyte:"
use_footnote <- "Program: \nOutput:"

result <- g_lineplot(
  df = adpc,
  variables = control_lineplot_vars(
    x = "NFRLT",
    y = "AVAL",
    group_var = "ARM",
    paramcd = "PARAM",
    y_unit = "AVALU",
    subject_var = "USUBJID"
  ),
  alt_counts_df = adsl_f,
  position = ggplot2::position_dodge2(width = 0.5),
  y_lab = "Concentration",
  y_lab_add_paramcd = FALSE,
  y_lab_add_unit = TRUE,
  interval = "mean_sdi",
  whiskers = c("mean_sdi_lwr", "mean_sdi_upr"),
  title = use_title,
  subtitle = use_subtitle,
  caption = use_footnote,
  ggtheme = theme_nest()
)

plot <- result + theme(plot.caption = element_text(hjust = 0))
plot

## ----r listing----------------------------------------------------------------
# Get value of Analyte
analyte <- unique(adpc$PARAM)

# Select columns for listing
out <- adpc %>%
  select(ARM, USUBJID, VISIT, NFRLT, AFRLT, AVALCAT1)

# Add descriptive labels
var_labels(out) <- c(
  ARM = "Treatment Group",
  USUBJID = "Subject ID",
  VISIT = "Visit",
  NFRLT = paste0("Nominal\nSampling\nTime (", adpc$RRLTU[1], ")"),
  AFRLT = paste0("Actual Time\nFrom First\nDose (", adpc$RRLTU[1], ")"),
  AVALCAT1 = paste0("Concentration\n(", adpc$AVALU[1], ")")
)

# Create listing
lsting <- as_listing(
  out,
  key_cols = c("ARM", "USUBJID", "VISIT"),
  disp_cols = names(out),
  default_formatting = list(
    all = fmt_config(align = "left"),
    numeric = fmt_config(
      format = "xx.xx",
      na_str = " ",
      align = "right"
    )
  ),
  main_title = paste(
    "Listing of",
    analyte,
    "Concentration by Treatment Group, Subject and Nominal Time, PK Population\nProtocol: xxnnnnn"
  ),
  subtitles = paste("Analyte:", analyte)
)

head(lsting, 28)

