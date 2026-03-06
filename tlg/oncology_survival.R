## ----r setup------------------------------------------------------------------
library(pharmaverseadam)
library(ggsurvfit)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(forcats)
library(broom)
library(survival)

# ── Read data ──────────────────────────────────────────────────────────────────
adtte_onco <- pharmaverseadam::adtte_onco

# ── ADRS_ONCO: tumor response data ───────────────────────────────────────────
adrs_onco <- pharmaverseadam::adrs_onco |>
  filter(ARMCD != "Scrnfail")

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
# ── PFS endpoint ────────────────────────────────────
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
    x = "Time (Days)",
    y = "Progression-Free Survival Probability",
    caption = paste0(
      "Analysis dataset: ADTTE_ONCO  |  PARAMCD: ", unique(adtte_pfs$PARAMCD),
      "\nCensored observations marked with '+'"
    )
  ) +
  theme_ggsurvfit_default() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

## ----r median-table-----------------------------------------------------------
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

## ----r prob-table-------------------------------------------------------------
# ── Survival probability at selected time points ───────────────────────────────
# AVAL in adtte_onco is in days
tbl_survfit(
  km_fit,
  times        = c(30, 60, 90, 180), # 1, 2, 3, 6 months in days
  label_header = "**PFS Probability (95% CI)**"
) |>
  modify_header(
    label = "**Time Point**"
  ) |>
  modify_table_body(
    ~ .x |>
      mutate(label = recode(label,
        "30"  = "1 month",
        "60"  = "2 months",
        "90"  = "3 months",
        "180" = "6 months"
      ))
  ) |>
  modify_caption(
    paste0(
      "**Table 2. Progression-Free Survival Probability at Selected Time Points**",
      "\nADTTE_ONCO  |  PARAMCD: ", unique(adtte_pfs$PARAMCD)
    )
  ) |>
  bold_labels()

## ----r bor-setup--------------------------------------------------------------
# ── Filter to CBOR parameter ──────────────────────────────────
# ANL01FL == "Y" restricts to the primary analysis flag records.
adrs_bor <- adrs_onco |>
  filter(PARAMCD == "CBOR" & ANL01FL == "Y") |>
  mutate(
    # Order AVALC from best to worst response for table display
    AVALC = fct_relevel(
      AVALC,
      "CR", "PR", "SD", "NON-CR/NON-PD", "PD", "NE", "MISSING"
    )
  )

## ----r bor-table--------------------------------------------------------------
# ── Best Overall Response table by treatment arm ───────────────────────────────
adrs_bor |>
  tbl_summary(
    by = ARM,
    include = AVALC,
    label = list(AVALC = "Best Overall Response"),
    statistic = list(AVALC = "{n} ({p}%)"),
    digits = list(AVALC = list(0, 1))
  ) |>
  add_overall(last = TRUE) |>
  add_n() |>
  bold_labels() |>
  modify_header(label = "**Response**") |>
  modify_caption(
    paste0(
      "**Table 3. Confirmed Best Overall Response (RECIST 1.1)**",
      "\nADRS_ONCO  |  PARAMCD: CBOR  |  ANL01FL = Y"
    )
  )

## ----r orr-inline-------------------------------------------------------------
# ── ORR: proportion with CR or PR, by arm ─────────────────────────────────────
adrs_bor |>
  summarise(
    .by    = ARM,
    n_resp = sum(AVALC %in% c("CR", "PR"), na.rm = TRUE),
    n_tot  = n(),
    orr    = round(100 * n_resp / n_tot, 1)
  ) |>
  mutate(label = paste0(n_resp, "/", n_tot, " (", orr, "%)")) |>
  select(ARM, ORR = label) |>
  knitr::kable(caption = "Overall Response Rate (CR + PR) by Arm")

## ----r cnsr-note--------------------------------------------------------------
# # ✗  Error-prone: requires manual recoding of CNSR
# survival::Surv(adtte_pfs$AVAL, 1 - adtte_pfs$CNSR)
# 
# # ✓  Correct CDISC-aware approach
# ggsurvfit::Surv_CNSR(adtte_pfs$AVAL, adtte_pfs$CNSR)

## ----r endpoint-note----------------------------------------------------------
# # Overall Survival — expect few events; median may not be estimable
# adtte_os <- adtte_onco |> filter(PARAMCD == "OS")
# 
# # Duration of Response — responders only; smaller N than OS/PFS
# # Note: admiralonco uses PARAMCD = "RSD", not "DOR"
# adtte_rsd <- adtte_onco |> filter(PARAMCD == "RSD")

## ----r strata-note------------------------------------------------------------
# # With two arms (ARM), add_pvalue() computes and annotates a log-rank test p-value.
# # Not applicable for single-arm fits (~ 1) — only add when comparing groups.
# survfit2(Surv_CNSR(AVAL, CNSR) ~ ARM, data = adtte_pfs) |>
#   ggsurvfit(linewidth = 1) +
#   scale_color_brewer(palette = "Dark2") +
#   scale_fill_brewer(palette = "Dark2") +
#   add_confidence_interval() +
#   add_risktable(
#     theme = theme_risktable_default(axis.text.y.size = 9, plot.title.size = 9)
#   ) +
#   add_pvalue(location = "annotation") +
#   scale_ggsurvfit()

## ----r km-plot-adtte----------------------------------------------------------
# ── ggsurvfit::adtte — four-arm breast cancer PFS trial ────────────────────────
# TRT01P = planned treatment; STR01 = hormone receptor status (not treatment)
survfit2(Surv_CNSR(AVAL, CNSR) ~ TRT01P, data = ggsurvfit::adtte) |>
  ggsurvfit(linewidth = 1) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
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
    x = "Time (Years)", # ggsurvfit::adtte AVAL is in years
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

## ----r forest-setup-----------------------------------------------------------
# ── Restrict to two arms ───────────────────────────────────────────────────────
adtte_2arm <- ggsurvfit::adtte |>
  filter(TRT01PN %in% c(1, 2)) |>
  mutate(TRT01P = factor(TRT01P,
    levels = c("tablemab x 52 weeks", "vismab x 52 weeks")
  ))

# ── Helper: fit Cox and return a one-row HR summary ──────────────────────────
cox_hr <- function(data, subgroup, level) {
  fit <- coxph(Surv_CNSR(AVAL, CNSR) ~ TRT01P, data = data)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    mutate(
      subgroup = subgroup,
      level    = level,
      n        = nrow(data),
      n_events = sum(data$CNSR == 0)
    )
}

# ── Build rows: Overall + by STR01L + by STR02L ───────────────────────────────
str01_rows <- lapply(unique(adtte_2arm$STR01L), function(lv) {
  cox_hr(filter(adtte_2arm, STR01L == lv), "Hormone Receptor Status", lv)
})

str02_rows <- lapply(unique(adtte_2arm$STR02L), function(lv) {
  cox_hr(filter(adtte_2arm, STR02L == lv), "Prior Radiotherapy", lv)
})

forest_data <- bind_rows(
  cox_hr(adtte_2arm, "Overall", "Overall"),
  bind_rows(str01_rows),
  bind_rows(str02_rows)
) |>
  mutate(
    label = ifelse(level == "Overall", "Overall", paste0("  ", level)),
    label = factor(label, levels = rev(unique(label))),
    hr_label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, conf.low, conf.high)
  )

## ----r forest-plot------------------------------------------------------------
ggplot(forest_data, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_pointrange(
    aes(xmin = conf.low, xmax = conf.high),
    color = "#2c7bb6",
    linewidth = 0.75,
    fatten = 4
  ) +
  geom_text(
    aes(x = 3.5, label = hr_label),
    hjust = 1, size = 3
  ) +
  scale_x_log10(
    limits = c(0.3, 4),
    breaks = c(0.5, 1, 2),
    labels = c("0.5", "1", "2")
  ) +
  facet_grid(subgroup ~ ., scales = "free_y", space = "free") +
  labs(
    title    = "Subgroup Forest Plot: PFS Hazard Ratio",
    subtitle = "vismab x 52 weeks vs. tablemab x 52 weeks (reference)",
    x        = "Hazard Ratio (log scale)  |  \u2190 Favours tablemab   Favours vismab \u2192",
    y        = NULL,
    caption  = "Dataset: ggsurvfit::adtte  |  Cox proportional hazards model  |  95% CI"
  ) +
  theme_bw() +
  theme(
    strip.background    = element_rect(fill = "gray90"),
    strip.text          = element_text(face = "bold"),
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_blank(),
    plot.caption        = element_text(hjust = 0, size = 8)
  )

