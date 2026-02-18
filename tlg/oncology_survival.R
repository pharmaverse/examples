## ----r setup, message=FALSE, warning=FALSE------------------------------------
library(pharmaverseadam)
library(ggsurvfit)
library(gtsummary)
library(dplyr)

# ── Read data ──────────────────────────────────────────────────────────────────
adtte <- pharmaverseadam::adtte_onco

# Keep a single endpoint for the basic example: Overall Survival (OS).
# Adjust PARAMCD to the endpoint of interest for other analyses.
adtte_os <- adtte |>
  filter(PARAMCD == "OS")

# Preview key variables
adtte_os |>
  select(USUBJID, PARAM, PARAMCD, AVAL, CNSR, EVNTDESC, CNSDTDSC) |>
  head(5)

## ----r km-plot, message=FALSE, warning=FALSE, fig.width=9, fig.height=7-------
# ── Fit the survival model ─────────────────────────────────────────────────────
# survfit2() is preferred over survfit() as it tracks the calling environment,
# enabling clean legend labels and p-value computation downstream.
km_fit <- survfit2(Surv_CNSR(AVAL, CNSR) ~ 1, data = adtte_os)

# ── Build the plot ─────────────────────────────────────────────────────────────
km_fit |>
  ggsurvfit(linewidth = 1) +
  add_confidence_interval() +
  add_risktable(
    risktable_stats = c("n.risk", "cum.event", "cum.censor"),
    stats_label = list(
      n.risk = "At Risk",
      cum.event = "Events (cum.)",
      cum.censor = "Censored (cum.)"
    )
  ) +
  add_quantile(
    y_value   = 0.5, # median survival guideline
    color     = "gray40",
    linewidth = 0.75,
    linetype  = "dashed"
  ) +
  add_censor_mark(shape = 3, size = 2) +
  scale_ggsurvfit(x_scales = list(breaks = 0:6)) +
  labs(
    title = paste0(unique(adtte_os$PARAM), "\nKaplan-Meier Estimate"),
    x = "Time (Years)",
    y = "Overall Survival Probability",
    caption = paste0(
      "Analysis dataset: ADTTE  |  PARAMCD: ", unique(adtte_os$PARAMCD),
      "\nCensored observations marked with '+'"
    )
  ) +
  theme_ggsurvfit_default() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

## ----r median-table, message=FALSE, warning=FALSE-----------------------------
# ── Median survival with 95% CI ────────────────────────────────────────────────
tbl_survfit(
  km_fit,
  probs = 0.5, # 50th percentile = median
  label_header = "**Median (95% CI)**"
) |>
  modify_caption(
    paste0(
      "**Table 1. Median Overall Survival**",
      "\nADTTE  |  PARAMCD: ", unique(adtte_os$PARAMCD)
    )
  ) |>
  bold_labels()

## ----r prob-table, message=FALSE, warning=FALSE-------------------------------
# ── Survival probability at selected time points ───────────────────────────────
tbl_survfit(
  km_fit,
  times        = c(0.5, 1, 1.5, 2), # time in years (AVAL units in adtte)
  label_header = "**Survival Probability (95% CI)**"
) |>
  modify_header(
    label = "**Time Point (Years)**"
  ) |>
  modify_caption(
    paste0(
      "**Table 2. Overall Survival Probability at Selected Time Points**",
      "\nADTTE  |  PARAMCD: ", unique(adtte_os$PARAMCD)
    )
  ) |>
  bold_labels()

## ----r cnsr-note, eval=FALSE--------------------------------------------------
# # ✗  Error-prone: requires manual recoding of CNSR
# survival::Surv(adtte_os$AVAL, 1 - adtte_os$CNSR)
#
# # ✓  Correct CDISC-aware approach
# ggsurvfit::Surv_CNSR(adtte_os$AVAL, adtte_os$CNSR)

## ----r strata-note, eval=FALSE------------------------------------------------
# # Merge treatment from ADSL
# adsl <- pharmaverseadam::adsl
# adtte_os_arm <- adtte_os |>
#   left_join(select(adsl, USUBJID, TRT01A), by = "USUBJID")
#
# # Stratified KM
# survfit2(Surv_CNSR(AVAL, CNSR) ~ TRT01A, data = adtte_os_arm) |>
#   ggsurvfit(linewidth = 1) +
#   add_confidence_interval() +
#   add_risktable() +
#   scale_ggsurvfit()
