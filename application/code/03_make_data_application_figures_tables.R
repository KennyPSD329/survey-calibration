# ==============================================================================
# 03_make_data_application_figures_tables.R
#
# Reads processed data and model results, generates all primary manuscript
# figures and supplementary tables, and saves them to figures/primary/ and
# results/.
#
# Run from the application/ root directory after Scripts 01 and 02.
#
# Primary outputs (figures/primary/):
#   realdata_ci_plot_all.{pdf,png}           - Main CI plot (paper Figure)
#   DXA_bodyfat_by_NHANES_cycle_boxplot.{pdf,png} - Supplement: DXA by cycle
#   BMI_distribution_NHANES_NHIS.{pdf,png}  - Supplement: BMI distributions
#
# Primary table outputs (results/):
#   eda_diagnostics/sample_summary_table.csv
#   eda_diagnostics/tab_eda_sample_summary.tex
#   CI_estimates_table.csv
# ==============================================================================

rm(list = ls())
library(haven)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(survey)

if (!file.exists("data/processed/nhanes_analysis.rds")) {
  stop("Run 01_prepare_data_application_data.R first.")
}
if (!file.exists("results/CI_table_tidy.rds")) {
  stop("Run 02_run_data_application_models.R first.")
}
if (!file.exists("data/processed/nhanes_4cycles.rds")) {
  stop("Run 01_prepare_data_application_data.R first (nhanes_4cycles.rds missing).")
}

dir.create("figures/primary", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/extra",   showWarnings = FALSE, recursive = TRUE)
dir.create("results/eda_diagnostics", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Figure 1 (Primary): Confidence interval plot
# ==============================================================================
cat("Generating CI plot...\n")

CI_table_tidy <- readRDS("results/CI_table_tidy.rds")

method_order <- c(
  "Benchmark: NHANES 01-06",
  "Internal: NHANES 03-04",
  "NHANES pooled calibration",
  "NHANES external calibration",
  "NHIS pooled calibration",
  "NHIS external calibration"
)

variable_order <- c(
  "Age", "Female", "NH Black", "Hispanic", "NH Other",
  "Former smoker", "Current smoker", "Some-day drinker",
  "Heavy drinker", "Physical activity", "Total fat"
)

variable_labels <- c(
  "Age"              = "Age\n(x 10^{-2})",
  "Female"           = "Female\n(x 10^{-1})",
  "NH Black"         = "NH Black\n(x 10^{-1})",
  "Hispanic"         = "Hispanic\n(x 10^{-1})",
  "NH Other"         = "NH Other\n(x 10^{-1})",
  "Former smoker"    = "Former smoker\n(x 10^{-1})",
  "Current smoker"   = "Current smoker\n(x 10^{-1})",
  "Some-day drinker" = "Some-day drinker\n(x 10^{-1})",
  "Heavy drinker"    = "Heavy drinker\n(x 10^{-2})",
  "Physical activity"= "Physical activity\n(x 10^{-2})",
  "Total fat"        = "Total fat\n(x 10^{-6})"
)

vline_df <- data.frame(
  variable   = factor(setdiff(variable_order, c("Age", "Current smoker")),
                      levels = variable_order),
  xintercept = 0
)

ci_plot_data <- CI_table_tidy %>%
  filter(method %in% method_order, variable %in% variable_order) %>%
  mutate(
    method   = factor(method,   levels = rev(method_order)),
    variable = factor(variable, levels = variable_order)
  )

method_colors <- c(
  "Benchmark: NHANES 01-06"    = "black",
  "Internal: NHANES 03-04"     = "grey45",
  "NHANES pooled calibration"  = "#0072B2",
  "NHANES external calibration"= "#56B4E9",
  "NHIS pooled calibration"    = "#009E73",
  "NHIS external calibration"  = "#7CAE00"
)

p_ci <- ggplot(ci_plot_data, aes(y = method, x = estimate, color = method)) +
  geom_vline(
    data = vline_df,
    aes(xintercept = xintercept),
    linetype = "dashed", linewidth = 0.45, color = "grey65",
    inherit.aes = FALSE
  ) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.22, linewidth = 0.85) +
  geom_point(size = 3.1) +
  facet_wrap2(
    ~variable, scales = "free_x", ncol = 2,
    labeller = labeller(variable = variable_labels)
  ) +
  facetted_pos_scales(
    x = list(
      variable == "Age" ~ scale_x_continuous(
        limits = c(6.5, 15), breaks = c(7.5, 10, 12.5, 15),
        expand = expansion(mult = c(0.02, 0.02))
      ),
      variable == "Current smoker" ~ scale_x_continuous(
        limits = c(7, 21), breaks = c(8, 12, 16, 20),
        expand = expansion(mult = c(0.02, 0.02))
      )
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  scale_color_manual(values = method_colors, breaks = method_order) +
  labs(
    x     = "Point estimate and 95% confidence interval",
    y     = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 12.5, face = "bold"),
    legend.key.width   = unit(1.4, "lines"),
    legend.key.height  = unit(1.0, "lines"),
    strip.background   = element_rect(fill = "grey92", color = "grey65"),
    strip.text         = element_text(face = "bold", size = 15),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.40),
    axis.text.y        = element_text(size = 11.5, face = "bold"),
    axis.text.x        = element_text(size = 11.5, face = "bold"),
    axis.title.x       = element_text(size = 13.5, face = "bold"),
    panel.spacing      = unit(1.15, "lines"),
    plot.margin        = margin(t = 5, r = 8, b = 5, l = 8)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figures/primary/realdata_ci_plot_all.pdf", p_ci,
       width = 8.8, height = 10.2, units = "in")
ggsave("figures/primary/realdata_ci_plot_all.png", p_ci,
       width = 8.8, height = 10.2, units = "in", dpi = 300)
cat("  Saved figures/primary/realdata_ci_plot_all.{pdf,png}\n")

# ==============================================================================
# Figure 2 (Primary): DXA total body fat by NHANES cycle and mortality status
#
# Shows distributions across all 4 cycles (1999-2000 through 2005-2006),
# grouped by 10-year all-cause mortality.  Demonstrates the distributional
# difference in NHANES 1999-2000 that motivates its exclusion from the main
# calibration analysis.  Corresponds to Fat.png in the original analysis.
# ==============================================================================
cat("Generating DXA body fat by NHANES cycle boxplot...\n")

nhanes_4cyc <- readRDS("data/processed/nhanes_4cycles.rds")

nhanes_4cyc$cycle_label <- factor(
  nhanes_4cyc$cycle,
  levels = c(1, 2, 3, 4),
  labels = c("1999-2000", "2001-2002", "2003-2004", "2005-2006")
)
nhanes_4cyc$Mortality <- factor(nhanes_4cyc$mortality,
                                levels = c(0, 1), labels = c("Alive", "Dead"))

p_fat <- ggplot(nhanes_4cyc, aes(x = cycle_label, y = DXDTOFAT, fill = Mortality)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  coord_cartesian(
    ylim = quantile(nhanes_4cyc$DXDTOFAT, c(0.05, 0.95), na.rm = TRUE)
  ) +
  scale_fill_manual(values = c("Alive" = "#F8766D", "Dead" = "#00BFC4")) +
  labs(
    x    = "NHANES Cycles",
    y    = "Total Body Fat (g)",
    fill = "Mortality"
  ) +
  theme_minimal() +
  theme(
    axis.title  = element_text(size = 18),
    axis.text   = element_text(size = 16),
    legend.title= element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave("figures/primary/DXA_bodyfat_by_NHANES_cycle_boxplot.pdf", p_fat,
       width = 7, height = 5.5)
ggsave("figures/primary/DXA_bodyfat_by_NHANES_cycle_boxplot.png", p_fat,
       width = 7, height = 5.5, dpi = 300)
cat("  Saved figures/primary/DXA_bodyfat_by_NHANES_cycle_boxplot.{pdf,png}\n")

# ==============================================================================
# Figure 3 (Primary): BMI distribution comparison (NHANES measured, self-reported, NHIS)
# ==============================================================================
cat("Generating BMI distribution plot...\n")

safe_xpt <- function(path, vars) {
  tryCatch({
    df <- read_xpt(path)
    df[, intersect(vars, names(df)), drop = FALSE]
  }, error = function(e) {
    warning(sprintf("Cannot read %s: %s", path, conditionMessage(e)))
    NULL
  })
}

load_nhanes_bmi_cycle <- function(demo_path, bmx_path, whq_path, label) {
  demo <- safe_xpt(demo_path, c("SEQN", "SDMVPSU", "SDMVSTRA", "WTMEC2YR", "RIDAGEYR"))
  bmx  <- safe_xpt(bmx_path,  c("SEQN", "BMXBMI"))
  whq  <- safe_xpt(whq_path,  c("SEQN", "WHD010", "WHD020"))
  if (is.null(demo) || is.null(bmx)) return(NULL)
  df <- left_join(demo, bmx, by = "SEQN")
  if (!is.null(whq)) df <- left_join(df, whq, by = "SEQN")
  df$cycle <- label
  df
}

bmi_01 <- load_nhanes_bmi_cycle(
  "data/raw/nhanes/01/DEMO_B.XPT", "data/raw/nhanes/01/BMX_B.XPT",
  "data/raw/nhanes/01/WHQ_B.XPT",  "2001-2002")
bmi_03 <- load_nhanes_bmi_cycle(
  "data/raw/nhanes/03/DEMO_C.XPT", "data/raw/nhanes/03/BMX_C.XPT",
  "data/raw/nhanes/03/WHQ_C.XPT",  "2003-2004")
bmi_05 <- load_nhanes_bmi_cycle(
  "data/raw/nhanes/05/DEMO_D.XPT", "data/raw/nhanes/05/BMX_D.XPT",
  "data/raw/nhanes/05/WHQ_D.XPT",  "2005-2006")

# Fallback to lowercase filenames
if (is.null(bmi_01)) bmi_01 <- load_nhanes_bmi_cycle(
  "data/raw/nhanes/01/DEMO_B.xpt", "data/raw/nhanes/01/BMX_B.xpt",
  "data/raw/nhanes/01/WHQ_B.xpt",  "2001-2002")
if (is.null(bmi_03)) bmi_03 <- load_nhanes_bmi_cycle(
  "data/raw/nhanes/03/DEMO_C.xpt", "data/raw/nhanes/03/BMX_C.xpt",
  "data/raw/nhanes/03/WHQ_C.xpt",  "2003-2004")
if (is.null(bmi_05)) bmi_05 <- load_nhanes_bmi_cycle(
  "data/raw/nhanes/05/DEMO_D.xpt", "data/raw/nhanes/05/BMX_D.xpt",
  "data/raw/nhanes/05/WHQ_D.xpt",  "2005-2006")

nhanes_bmi <- bind_rows(Filter(Negate(is.null), list(bmi_01, bmi_03, bmi_05)))
n_bmi_cycles <- length(unique(nhanes_bmi$cycle))
nhanes_bmi$wt_mec <- nhanes_bmi$WTMEC2YR / n_bmi_cycles
nhanes_bmi$WHD020[!is.na(nhanes_bmi$WHD020) & nhanes_bmi$WHD020 > 1000] <- NA
nhanes_bmi$WHD010[!is.na(nhanes_bmi$WHD010) & nhanes_bmi$WHD010 > 1000] <- NA
nhanes_bmi$BMI_self_raw <- nhanes_bmi$WHD020 * 0.45359237 / (nhanes_bmi$WHD010 * 0.0254)^2

nhanes_meas <- nhanes_bmi %>%
  filter(!is.na(RIDAGEYR), RIDAGEYR >= 40, RIDAGEYR <= 69,
         !is.na(wt_mec), wt_mec > 0,
         !is.na(BMXBMI), BMXBMI >= 10, BMXBMI <= 80)

nhanes_self_bmi <- nhanes_bmi %>%
  filter(!is.na(RIDAGEYR), RIDAGEYR >= 40, RIDAGEYR <= 69,
         !is.na(wt_mec), wt_mec > 0,
         !is.na(BMI_self_raw), BMI_self_raw >= 10, BMI_self_raw <= 80)

nhis_raw <- tryCatch(
  read.csv("data/raw/nhis/nhis_2003_2004.csv", stringsAsFactors = FALSE),
  error = function(e) { warning("Cannot read nhis_2003_2004.csv"); NULL }
)

plot_df <- NULL
if (!is.null(nhis_raw)) {
  nhis_raw$WEIGHT[nhis_raw$WEIGHT %in% c(996,997,998,999)] <- NA
  nhis_raw$WEIGHT[!is.na(nhis_raw$WEIGHT) & (nhis_raw$WEIGHT < 50 | nhis_raw$WEIGHT > 299)] <- NA
  nhis_raw$HEIGHT[nhis_raw$HEIGHT %in% c(95,96,97,98,99)] <- NA
  nhis_raw$HEIGHT[!is.na(nhis_raw$HEIGHT) & (nhis_raw$HEIGHT < 36 | nhis_raw$HEIGHT > 90)] <- NA
  nhis_raw$BMI_self_raw <- (nhis_raw$WEIGHT * 0.45359237) / ((nhis_raw$HEIGHT * 0.0254)^2)
  nhis_raw$wt <- nhis_raw$SAMPWEIGHT / 2
  nhis_bmi <- nhis_raw %>%
    filter(!is.na(AGE), AGE >= 40, AGE <= 69,
           !is.na(wt), wt > 0,
           !is.na(BMI_self_raw), BMI_self_raw >= 10, BMI_self_raw <= 80)
  if ("PREGNANTNOW" %in% names(nhis_bmi) && "SEX" %in% names(nhis_bmi))
    nhis_bmi <- nhis_bmi %>% filter(!(SEX == 2 & PREGNANTNOW %in% c(2,7,8,9)))

  nhanes_meas$wt_norm     <- nhanes_meas$wt_mec     / sum(nhanes_meas$wt_mec)
  nhanes_self_bmi$wt_norm <- nhanes_self_bmi$wt_mec / sum(nhanes_self_bmi$wt_mec)
  nhis_bmi$wt_norm        <- nhis_bmi$wt            / sum(nhis_bmi$wt)

  plot_df <- bind_rows(
    data.frame(BMI = nhanes_meas$BMXBMI,           wt = nhanes_meas$wt_norm,
               source = "NHANES measured",      stringsAsFactors = FALSE),
    data.frame(BMI = nhanes_self_bmi$BMI_self_raw,  wt = nhanes_self_bmi$wt_norm,
               source = "NHANES self-reported", stringsAsFactors = FALSE),
    data.frame(BMI = nhis_bmi$BMI_self_raw,         wt = nhis_bmi$wt_norm,
               source = "NHIS self-reported",   stringsAsFactors = FALSE)
  )
}

if (!is.null(plot_df)) {
  src_levels <- c("NHANES measured", "NHANES self-reported", "NHIS self-reported")
  plot_df$source <- factor(plot_df$source, levels = src_levels)
  bmi_colors <- c(
    "NHANES measured"       = "#2166AC",
    "NHANES self-reported"  = "#1A9850",
    "NHIS self-reported"    = "#D01C8B"
  )

  p_bmi <- ggplot(plot_df, aes(x = BMI, weight = wt, color = source, fill = source)) +
    geom_density(alpha = 0.15, linewidth = 0.85) +
    scale_color_manual(values = bmi_colors, name = NULL) +
    scale_fill_manual( values = bmi_colors, name = NULL) +
    coord_cartesian(xlim = c(15, 65)) +
    labs(
      x = "BMI (kg/m²)",
      y = "Weighted density",
      title = "Weighted BMI Distributions: NHANES vs NHIS",
      subtitle = "NHANES cycles 2001-2006 · NHIS 2003-2004 · Adults age 40-69"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.key.width = unit(1.5, "cm"),
      plot.title       = element_text(face = "bold")
    )

  ggsave("figures/primary/BMI_distribution_NHANES_NHIS.pdf", p_bmi, width = 8, height = 5.5)
  ggsave("figures/primary/BMI_distribution_NHANES_NHIS.png", p_bmi, width = 8, height = 5.5, dpi = 300)
  cat("  Saved figures/primary/BMI_distribution_NHANES_NHIS.{pdf,png}\n")
} else {
  cat("  Skipped BMI distribution: NHIS file not available.\n")
}

# ==============================================================================
# Table 1: EDA sample summary (weighted, survey-design-aware)
# Covers all 4 NHANES cycles + NHIS 2003-2004; used in Supplementary Material.
# Logic adapted from eda_supplement_diagnostics.R.
# ==============================================================================
cat("Generating EDA sample summary table...\n")

nhis_summ   <- readRDS("data/processed/nhis_analysis.rds")
nhanes_4cyc <- readRDS("data/processed/nhanes_4cycles.rds")
nhis_summ$female <- as.integer(as.character(nhis_summ$RIAGENDR) == "Female")
nhanes_4cyc$female <- as.integer(as.character(nhanes_4cyc$RIAGENDR) == "Female")

make_summary <- function(df, src_label) {
  ds <- tryCatch(
    svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weight = ~WTINT2YR, nest = TRUE, data = df),
    error = function(e) svydesign(id = ~1, weight = ~WTINT2YR, data = df)
  )

  n_total  <- nrow(df)
  n_deaths <- sum(df$mortality, na.rm = TRUE)
  pct_mort   <- round(100 * as.numeric(svymean(~mortality, ds, na.rm = TRUE)), 1)
  pct_female <- round(100 * as.numeric(svymean(~female,    ds, na.rm = TRUE)), 1)

  r_tbl <- svymean(~race4, ds, na.rm = TRUE)
  rn    <- names(r_tbl)
  get_r <- function(pat) {
    idx <- grep(pat, rn)
    if (length(idx)) round(100 * as.numeric(r_tbl[idx]), 1) else NA
  }
  pct_nhblack  <- get_r("NH-Black")
  pct_hispanic <- get_r("Hispanic")
  pct_nhother  <- get_r("NH-Other")

  n_hispanic    <- sum(df$race4 == "Hispanic", na.rm = TRUE)
  n_nhother     <- sum(df$race4 == "NH-Other",  na.rm = TRUE)
  mort_hispanic <- sum(df$mortality[df$race4 == "Hispanic"], na.rm = TRUE)
  mort_nhother  <- sum(df$mortality[df$race4 == "NH-Other"],  na.rm = TRUE)

  s_tbl <- svymean(~as.factor(smk3), ds, na.rm = TRUE)
  sn    <- names(s_tbl)
  get_s <- function(suf) {
    idx <- grep(paste0(suf, "$"), sn)
    if (length(idx)) round(100 * as.numeric(s_tbl[idx]), 1) else NA
  }
  pct_former  <- get_s("1")
  pct_current <- get_s("2")
  n_former    <- sum(df$smk3 == 1, na.rm = TRUE)
  n_current   <- sum(df$smk3 == 2, na.rm = TRUE)
  mort_former  <- sum(df$mortality[df$smk3 == 1], na.rm = TRUE)
  mort_current <- sum(df$mortality[df$smk3 == 2], na.rm = TRUE)

  d_tbl <- svymean(~as.factor(hvy.drk), ds, na.rm = TRUE)
  dn    <- names(d_tbl)
  get_d <- function(suf) {
    idx <- grep(paste0(suf, "$"), dn)
    if (length(idx)) round(100 * as.numeric(d_tbl[idx]), 1) else NA
  }
  pct_social <- get_d("1")
  pct_heavy  <- get_d("2")

  data.frame(
    Source        = src_label,
    n             = n_total,
    Deaths_n      = n_deaths,
    Deaths_pct    = pct_mort,
    Female_pct    = pct_female,
    NHBlack_pct   = pct_nhblack,
    Hispanic_n    = n_hispanic,
    Hispanic_pct  = pct_hispanic,
    Hispanic_mort = mort_hispanic,
    NHOther_n     = n_nhother,
    NHOther_pct   = pct_nhother,
    NHOther_mort  = mort_nhother,
    Former_n      = n_former,
    Former_pct    = pct_former,
    Former_mort   = mort_former,
    Current_n     = n_current,
    Current_pct   = pct_current,
    Current_mort  = mort_current,
    Social_pct    = pct_social,
    Heavy_pct     = pct_heavy,
    stringsAsFactors = FALSE
  )
}

s99   <- make_summary(nhanes_4cyc %>% filter(cycle == 1), "NHANES 1999-2000")
s01   <- make_summary(nhanes_4cyc %>% filter(cycle == 2), "NHANES 2001-2002")
s03   <- make_summary(nhanes_4cyc %>% filter(cycle == 3), "NHANES 2003-2004 (internal)")
s05   <- make_summary(nhanes_4cyc %>% filter(cycle == 4), "NHANES 2005-2006")
snhis <- make_summary(nhis_summ,                           "NHIS 2003-2004")

eda_tab <- rbind(s99, s01, s03, s05, snhis)
write.csv(eda_tab, "results/eda_diagnostics/sample_summary_table.csv", row.names = FALSE)
cat("  Saved results/eda_diagnostics/sample_summary_table.csv\n")

# --- LaTeX version (two-panel table for Supplementary Material) ---
fmt_pct <- function(x) paste0(x, "\\%")

panA_header <- paste(
  "\\multicolumn{9}{l}{\\textit{Panel A: Sample size, mortality, and race/ethnicity}} \\\\",
  "\\midrule",
  "Data source & $n$ & \\makecell{Deaths\\\\$n$ (\\%)} &",
  "\\makecell{Female\\\\(\\%)} & \\makecell{NH-Black\\\\(\\%)} &",
  "\\makecell{Hispanic\\\\$n$ (\\%)} & \\makecell{Hispanic\\\\deaths} &",
  "\\makecell{NH-Other\\\\$n$ (\\%)} & \\makecell{NH-Other\\\\deaths} \\\\",
  sep = "\n"
)

panA_rows <- sapply(seq_len(nrow(eda_tab)), function(i) {
  r   <- eda_tab[i, ]
  src <- if (i == 1) paste0("\\textit{", r$Source, "}$^{\\dagger}$") else r$Source
  paste(src, "&", r$n, "&",
        paste0(r$Deaths_n, " (", r$Deaths_pct, "\\%)"), "&",
        fmt_pct(r$Female_pct), "&", fmt_pct(r$NHBlack_pct), "&",
        paste0(r$Hispanic_n, " (", r$Hispanic_pct, "\\%)"), "&", r$Hispanic_mort, "&",
        paste0(r$NHOther_n, " (", r$NHOther_pct, "\\%)"), "&", r$NHOther_mort, "\\\\")
})

panB_header <- paste(
  "\\\\[4pt]",
  "\\multicolumn{9}{l}{\\textit{Panel B: Smoking and alcohol use}} \\\\",
  "\\midrule",
  "Data source & \\makecell{Former smk\\\\$n$ (\\%)} & \\makecell{Former\\\\deaths} &",
  "\\makecell{Current smk\\\\$n$ (\\%)} & \\makecell{Current\\\\deaths} &",
  "\\makecell{Social\\\\drinker (\\%)} & \\makecell{Heavy\\\\drinker (\\%)} & & \\\\",
  sep = "\n"
)

panB_rows <- sapply(seq_len(nrow(eda_tab)), function(i) {
  r   <- eda_tab[i, ]
  src <- if (i == 1) paste0("\\textit{", r$Source, "}$^{\\dagger}$") else r$Source
  paste(src, "&",
        paste0(r$Former_n, " (", r$Former_pct, "\\%)"), "&", r$Former_mort, "&",
        paste0(r$Current_n, " (", r$Current_pct, "\\%)"), "&", r$Current_mort, "&",
        fmt_pct(r$Social_pct), "&", fmt_pct(r$Heavy_pct), "& & \\\\")
})

ltx_lines <- c(
  "% -----------------------------------------------------------------------",
  "% tab_eda_sample_summary.tex",
  "% Insert via \\input{tab_eda_sample_summary} in Supplementary Materials",
  "% -----------------------------------------------------------------------",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\caption{Unweighted sample sizes and weighted covariate distributions",
  "  for each data source after applying the inclusion criteria of the main analysis",
  "  (age $\\geq 40$, valid DXA measurement, valid self-reported BMI $\\geq 18.5$,",
  "  height 120--220\\,cm).",
  "  Percentages are weighted using the respective survey sampling weights.",
  "  For Hispanic and NH-Other subgroups, and for former and current smokers,",
  "  unweighted counts ($n$) and the number of deaths within each category are",
  "  also shown to facilitate interpretation of the confidence-interval plots.",
  "  $\\dagger$\\,NHANES 1999--2000 is excluded from the main analysis owing to",
  "  distributional differences in DXA-measured total body fat",
  "  (see Figure~\\ref{fig:bodyft-boxplot}).}",
  "\\label{tab:eda-sample-summary}",
  "\\begin{tabular}{lrrrrrrrr}",
  "\\toprule",
  panA_header,
  panA_rows,
  "\\addlinespace",
  panB_header,
  panB_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(ltx_lines, "results/eda_diagnostics/tab_eda_sample_summary.tex")
cat("  Saved results/eda_diagnostics/tab_eda_sample_summary.tex\n")

# ==============================================================================
# Table 2: CI estimates table (formatted for paper)
# ==============================================================================
cat("Generating CI estimates table...\n")

ci_table_wide <- readRDS("results/CI_table_tidy.rds") %>%
  filter(method %in% c("Benchmark: NHANES 01-06",
                        "Internal: NHANES 03-04",
                        "NHANES pooled calibration",
                        "NHANES external calibration",
                        "NHIS pooled calibration",
                        "NHIS external calibration")) %>%
  mutate(ci_str = sprintf("%.2f (%.2f, %.2f)", estimate, lower, upper)) %>%
  select(variable, method, ci_str) %>%
  tidyr::pivot_wider(names_from = method, values_from = ci_str)

write.csv(ci_table_wide, "results/CI_estimates_table.csv", row.names = FALSE)
cat("  Saved results/CI_estimates_table.csv\n")

cat("\nScript 03 complete. All primary figures and tables generated.\n")
cat("  Primary figures: figures/primary/\n")
cat("  Extra/archived:  figures/extra/  (populate manually if needed)\n")
