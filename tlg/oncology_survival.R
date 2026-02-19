## ----r setup------------------------------------------------------------------
library(pharmaverseadam)
library(ggsurvfit)
library(ggplot2)
library(gtsummary)
library(dplyr)

# ── Read data ──────────────────────────────────────────────────────────────────
adtte_onco <- pharmaverseadam::adtte_onco

# Overview of available endpoints and their event rates
adtte_onco |>
  group_by(PARAMCD, PARAM) |>
  summarise(
    N           = n(),
    N_events    = sum(CNSR == 0),
    Pct_events  = round(100 * mean(CNSR == 0), 1),
    Median_AVAL = round(median(AVAL), 2),
    .groups     = "drop"
  )

## ----r filter-pfs-------------------------------------------------------------
# ── Filter to PFS endpoint ─────────────────────────────────────────────────────
adtte_pfs <- adtte_onco |>
  filter(PARAMCD == "PFS")

# Preview key variables
adtte_pfs |>
  select(USUBJID, PARAM, PARAMCD, AVAL, CNSR, EVNTDESC, CNSDTDSC) |>
  head(5)

## ----r km-plot----------------------------------------------------------------
# ── Fit the survival model ─────────────────────────────────────────────────────
# survfit2() is preferred over survfit() as it tracks the calling environment,
# enabling clean legend labels and p-value computation downstream.
km_fit <- survfit2(Surv_CNSR(AVAL, CNSR) ~ 1, data = adtte_pfs)

# ── Build the plot ─────────────────────────────────────────────────────────────
km_fit |>
  ggsurvfit(linewidth = 1) +
  # ── Explicit color scales prevent the site theme from stripping colour to B&W ─
  scale_color_manual(values = c("#2c7bb6")) +
  scale_fill_manual(values = c("#2c7bb6")) +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = c("n.risk", "cum.event", "cum.censor"),
    stats_label = list(
      n.risk     = "At Risk",
      cum.event  = "Events (cum.)",
      cum.censor = "Censored (cum.)"
    ),
    # ── Slightly smaller label text keeps row labels from squeezing number columns
    theme = theme_risktable_default(axis.text.y.size = 9, plot.title.size = 9)
  ) +
  add_quantile(
    y_value   = 0.5, # median survival guideline
    color     = "gray40",
    linewidth = 0.75,
    linetype  = "dashed"
  ) +
  add_censor_mark(shape = 3, size = 2) +
  scale_ggsurvfit() +
  labs(
    title = paste0(unique(adtte_pfs$PARAM), "\nKaplan-Meier Estimate"),
    x = paste0("Time (", unique(adtte_pfs$AVALU), ")"),
    y = "Progression-Free Survival Probability",
    caption = paste0(
      "Analysis dataset: ADTTE_ONCO  |  PARAMCD: ", unique(adtte_pfs$PARAMCD),
      "\nCensored observations marked with '+'"
    )
  ) +
  theme_ggsurvfit_default() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

## ----r median-table-classic---------------------------------------------------
# ── Median survival with 95% CI ────────────────────────────────────────────────
tbl_survfit(
  km_fit,
  probs        = 0.5, # 50th percentile = median
  label_header = "**Median (95% CI)**"
) |>
  modify_caption(
    paste0(
      "**Table 1. Median Progression-Free Survival**",
      "\nADTTE_ONCO  |  PARAMCD: ", unique(adtte_pfs$PARAMCD)
    )
  ) |>
  bold_labels()

## ----r prob-table-classic-----------------------------------------------------
# ── Survival probability at selected time points ───────────────────────────────
# AVAL in adtte_onco is in years; express months as fractions of a year
tbl_survfit(
  km_fit,
  times        = c(0.25, 0.5, 0.75, 1), # 3, 6, 9, 12 months in years
  label_header = "**PFS Probability (95% CI)**"
) |>
  modify_header(
    label = "**Time Point**"
  ) |>
  modify_table_body(
    ~ .x |>
      mutate(label = recode(label,
        "0.25" = "3 months",
        "0.5"  = "6 months",
        "0.75" = "9 months",
        "1"    = "12 months"
      ))
  ) |>
  modify_caption(
    paste0(
      "**Table 2. Progression-Free Survival Probability at Selected Time Points**",
      "\nADTTE_ONCO  |  PARAMCD: ", unique(adtte_pfs$PARAMCD)
    )
  ) |>
  bold_labels()

## ----r cnsr-note--------------------------------------------------------------
# # ✗  Error-prone: requires manual recoding of CNSR
# survival::Surv(adtte_pfs$AVAL, 1 - adtte_pfs$CNSR)
#
# # ✓  Correct CDISC-aware approach
# ggsurvfit::Surv_CNSR(adtte_pfs$AVAL, adtte_pfs$CNSR)

## ----r endpoint-note----------------------------------------------------------
# # Overall Survival — expect few events; median may not be estimable
# adtte_os  <- adtte_onco |> filter(PARAMCD == "OS")
#
# # Duration of Response — responders only; smaller N than OS/PFS
# # Note: admiralonco uses PARAMCD = "RSD", not "DOR"
# adtte_rsd <- adtte_onco |> filter(PARAMCD == "RSD")

## ----r strata-note------------------------------------------------------------
# survfit2(Surv_CNSR(AVAL, CNSR) ~ ARM, data = adtte_pfs) |>
#   ggsurvfit(linewidth = 1) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_fill_brewer(palette  = "Dark2") +
#   add_confidence_interval() +
#   add_risktable(
#     theme = theme_risktable_default(axis.text.y.size = 9, plot.title.size = 9)
#   ) +
#   add_pvalue(location = "annotation") +
#   scale_ggsurvfit()

## ----r km-plot-adtte----------------------------------------------------------
# ── ggsurvfit::adtte — four-arm breast cancer PFS trial ────────────────────────
arm_labels <- c(
  "V"    = "Vinorelbine (V)",
  "T"    = "Docetaxel (T)",
  "T->V" = "Docetaxel → Vinorelbine",
  "T+V"  = "Docetaxel + Vinorelbine"
)

survfit2(Surv_CNSR(AVAL, CNSR) ~ STR01, data = ggsurvfit::adtte) |>
  ggsurvfit(linewidth = 1) +
  scale_color_brewer(palette = "Dark2", labels = arm_labels) +
  scale_fill_brewer(palette = "Dark2", labels = arm_labels) +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = "n.risk",
    stats_label = list(n.risk = "At Risk"),
    theme = theme_risktable_default(axis.text.y.size = 9, plot.title.size = 9)
  ) +
  add_quantile(
    y_value   = 0.5,
    color     = "gray40",
    linewidth = 0.75,
    linetype  = "dashed"
  ) +
  add_censor_mark(shape = 3, size = 1.5) +
  add_pvalue(location = "annotation", x = 4.5) +
  scale_ggsurvfit() +
  labs(
    title = "Progression-Free Survival by Treatment Arm",
    x = "Time (Years)",
    y = "Progression-Free Survival Probability",
    caption = paste0(
      "Dataset: ggsurvfit::adtte  |  HER2+ breast cancer Phase III trial",
      "\nCensored observations marked with '+'"
    )
  ) +
  theme_ggsurvfit_default() +
  theme(
    plot.caption  = element_text(hjust = 0, size = 8),
    legend.title  = element_blank()
  )
