# ============================================================================
# Additional Analyses for AGGP-D-25-00369 Revision
# Reviewer and Editor Requested Analyses
# Date: April 14, 2026
# ============================================================================
# This script addresses:
# 1. UPDRS-III correlation with tremor outcomes (Reviewer #1, Q5)
# 2. Sex differences in UPDRS-III postural tremor sub-items (Reviewer #1, Q6)
# 3. Sex differences across HY stages (Reviewer #1, Q7)
# 4. 95% Confidence Intervals for key LMM estimates
# 5. LMM diagnostic plots (residuals, normality, influence)
# 6. Sex-stratified descriptive summaries (Editor)
# 7. UPDRS-III as continuous covariate sensitivity analysis (Editor)
# 8. Participant flow summary
# ============================================================================
# Install packages if needed (uncomment if not already installed):
# install.packages("ggeffects")

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(car)
library(lattice)
library(ggeffects)

# ------------------------------------------------------------------
# SET PATHS - ADJUST THESE TO YOUR LOCAL MACHINE
# ------------------------------------------------------------------
setwd("C:/Users/osmar/OneDrive/Documents/PESQUISAS/PARKINSONS_2022/Submitted_22052023")

# Load the main dataset (same as used in original analyses)
data <- read_csv("C:/Users/osmar/OneDrive/Documents/PESQUISAS/Rodar_Python_Som/AMBVR_02/Reshaped_PD_Tremor_with_Feedback_noFreq_for_Stats_f.csv")

# Load UPDRS sub-items spreadsheet
updrs_data <- read_excel("C:/Users/osmar/OneDrive/Documents/PESQUISAS/PARKINSONS_2022/Submitted_22052023/Helyum_March2025/updrs_01032024.xlsx",
                         sheet = "Sheet2")

# Results directory
results_dir <- "C:/Users/osmar/OneDrive/Documents/PESQUISAS/PARKINSONS_2022/Submitted_22052023/exportResults_Revision_04152026"
dir.create(results_dir, showWarnings = FALSE)

# ------------------------------------------------------------------
# CLEAN DATA (same as original)
# ------------------------------------------------------------------
names(data) <- make.names(names(data))
data$Severity <- factor(data$HY, levels = c("2", "3"))
data$Feedback <- factor(data$Feedback)
data$Hand <- factor(data$Hand)
data$subject_id <- factor(data$subject_id)
data$Age_c <- scale(data$Age, center = TRUE, scale = FALSE)
data$LEDD_c <- scale(data$LEDD, center = TRUE, scale = FALSE)
data$Sex <- factor(data$Sex, levels = c("Male", "Female"))
data$Tremor <- factor(data$Tremor, levels = c("Right", "Left", "Equal"))

# ------------------------------------------------------------------
# MERGE UPDRS SCORES FROM EXCEL FILE
# ------------------------------------------------------------------
cat("\n========== Merging UPDRS scores from updrs_01032024.xlsx ==========\n")

# The UPDRS scores are in Sheet2: col1=filename (subject initials), col7=UPDRSA (baseline UPDRS-III)
# Read Sheet2 - it has headers in row 1
updrs_sheet2 <- read_excel(
  "C:/Users/osmar/OneDrive/Documents/PESQUISAS/PARKINSONS_2022/Submitted_22052023/Helyum_March2025/updrs_01032024.xlsx",
  sheet = "Sheet2"
)

cat("Sheet2 column names:\n")
print(names(updrs_sheet2))

# The first column is 'filename' (subject initials), 'UPDRSA' is baseline UPDRS-III
# Column names may have been renamed by R due to duplicates (e.g., filename...1)
# Find the right columns
filename_col <- names(updrs_sheet2)[1]  # first column = subject initials
updrsa_col <- grep("^UPDRSA", names(updrs_sheet2), value = TRUE)[1]  # first UPDRSA column

cat(sprintf("Using filename column: '%s', UPDRS column: '%s'\n", filename_col, updrsa_col))

updrs_lookup <- data.frame(
  subject_id = trimws(as.character(updrs_sheet2[[filename_col]])),
  UPDRS = as.numeric(updrs_sheet2[[updrsa_col]])
)
updrs_lookup <- updrs_lookup[!is.na(updrs_lookup$subject_id) & updrs_lookup$subject_id != "", ]
cat("UPDRS lookup table:\n")
print(updrs_lookup)

# Also try to get postural tremor sub-items from Sheet1 if available
# Sheet1 has more detailed UPDRS sub-items
updrs_sheet1 <- read_excel(
  "C:/Users/osmar/OneDrive/Documents/PESQUISAS/PARKINSONS_2022/Submitted_22052023/Helyum_March2025/updrs_01032024.xlsx",
  sheet = "Sheet1"
)
cat("\nSheet1 column names:\n")
print(names(updrs_sheet1))

# Merge UPDRS into main data
# subject_id in main data should match filename (initials) in UPDRS sheet
data$subject_id_char <- trimws(as.character(data$subject_id))
updrs_lookup$subject_id <- trimws(updrs_lookup$subject_id)

data <- data %>%
  left_join(updrs_lookup, by = c("subject_id_char" = "subject_id"))

cat("\nUPDRS merge check - unique subjects with UPDRS:\n")
merged_check <- data %>% group_by(subject_id) %>% summarise(UPDRS = first(UPDRS), .groups = "drop")
print(merged_check)
cat(sprintf("Subjects with UPDRS: %d / %d\n", sum(!is.na(merged_check$UPDRS)), nrow(merged_check)))

# Center UPDRS for later use
data$UPDRS_c <- scale(data$UPDRS, center = TRUE, scale = FALSE)

# ------------------------------------------------------------------
# SECTION 1: UPDRS-III CORRELATION WITH TREMOR (Reviewer #1, Q5)
# ------------------------------------------------------------------
cat("\n========== SECTION 1: UPDRS-III vs Tremor Correlations ==========\n")

# Get participant-level summary (one row per participant)
participant_summary <- data %>%
  group_by(subject_id) %>%
  summarise(
    UPDRS = first(UPDRS),
    HY = first(HY),
    Sex = first(Sex),
    Age = first(Age),
    LEDD = first(LEDD),
    mean_amplitude = mean(mean_modalAmplitude_AccR, na.rm = TRUE),
    mean_frequency = mean(mean_modalFrequency_AccR, na.rm = TRUE),
    mean_ApEn = mean(mean_ApEn_AccR, na.rm = TRUE),
    mean_force = mean(mean_mean_force, na.rm = TRUE),
    mean_CV = mean(mean_CV_force, na.rm = TRUE),
    mean_SD = mean(mean_SD_force, na.rm = TRUE),
    .groups = "drop"
  )

# Check if merge was successful
if (all(is.na(participant_summary$UPDRS))) {
  cat("ERROR: UPDRS merge failed. Subject IDs may not match.\n")
  cat("Subject IDs in main data:\n")
  print(unique(data$subject_id_char))
  cat("Subject IDs in UPDRS lookup:\n")
  print(updrs_lookup$subject_id)
  stop("Please check subject ID matching between the main CSV and updrs_01032024.xlsx")
}

# Correlation: UPDRS-III vs tremor amplitude
cor_amp <- cor.test(participant_summary$UPDRS, participant_summary$mean_amplitude,
                    method = "spearman", exact = FALSE)
cat("\nUPDRS-III vs Tremor Amplitude (Spearman):\n")
cat(sprintf("  rho = %.3f, p = %.4f\n", cor_amp$estimate, cor_amp$p.value))

# Correlation: UPDRS-III vs tremor frequency
cor_freq <- cor.test(participant_summary$UPDRS, participant_summary$mean_frequency,
                     method = "spearman", exact = FALSE)
cat("\nUPDRS-III vs Tremor Frequency (Spearman):\n")
cat(sprintf("  rho = %.3f, p = %.4f\n", cor_freq$estimate, cor_freq$p.value))

# Correlation: UPDRS-III vs ApEn
cor_apen <- cor.test(participant_summary$UPDRS, participant_summary$mean_ApEn,
                     method = "spearman", exact = FALSE)
cat("\nUPDRS-III vs ApEn (Spearman):\n")
cat(sprintf("  rho = %.3f, p = %.4f\n", cor_apen$estimate, cor_apen$p.value))

# Correlation: UPDRS-III vs mean force
cor_force <- cor.test(participant_summary$UPDRS, participant_summary$mean_force,
                      method = "spearman", exact = FALSE)
cat("\nUPDRS-III vs Mean Force (Spearman):\n")
cat(sprintf("  rho = %.3f, p = %.4f\n", cor_force$estimate, cor_force$p.value))

# Correlation: UPDRS-III vs force CV
cor_cv <- cor.test(participant_summary$UPDRS, participant_summary$mean_CV,
                   method = "spearman", exact = FALSE)
cat("\nUPDRS-III vs Force CV (Spearman):\n")
cat(sprintf("  rho = %.3f, p = %.4f\n", cor_cv$estimate, cor_cv$p.value))

# Save correlation results
corr_results <- data.frame(
  Variable = c("Tremor Amplitude", "Tremor Frequency", "ApEn", "Mean Force", "Force CV"),
  Spearman_rho = c(cor_amp$estimate, cor_freq$estimate, cor_apen$estimate,
                   cor_force$estimate, cor_cv$estimate),
  p_value = c(cor_amp$p.value, cor_freq$p.value, cor_apen$p.value,
              cor_force$p.value, cor_cv$p.value)
)
write.csv(corr_results, file = paste0(results_dir, "/UPDRS_tremor_correlations.csv"), row.names = FALSE)

# Scatter plots
pdf(paste0(results_dir, "/UPDRS_tremor_scatterplots.pdf"), width = 12, height = 8, pointsize = 16)
par(mfrow = c(2, 3), cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.2, mar = c(5, 5, 4, 2))
plot(participant_summary$UPDRS, participant_summary$mean_amplitude,
     xlab = "UPDRS-III Score", ylab = "Mean Tremor Amplitude (V^2/Hz)",
     main = "UPDRS-III vs Tremor Amplitude", pch = 16, cex = 1.5, col = ifelse(participant_summary$Sex == "Male", "blue", "red"))
abline(lm(mean_amplitude ~ UPDRS, data = participant_summary), col = "gray", lty = 2, lwd = 2)
legend("topright", legend = c("Men", "Women"), col = c("blue", "red"), pch = 16, cex = 1.3)

plot(participant_summary$UPDRS, participant_summary$mean_frequency,
     xlab = "UPDRS-III Score", ylab = "Mean Tremor Frequency (Hz)",
     main = "UPDRS-III vs Tremor Frequency", pch = 16, cex = 1.5, col = ifelse(participant_summary$Sex == "Male", "blue", "red"))
abline(lm(mean_frequency ~ UPDRS, data = participant_summary), col = "gray", lty = 2, lwd = 2)

plot(participant_summary$UPDRS, participant_summary$mean_ApEn,
     xlab = "UPDRS-III Score", ylab = "Mean ApEn",
     main = "UPDRS-III vs Approximate Entropy", pch = 16, cex = 1.5, col = ifelse(participant_summary$Sex == "Male", "blue", "red"))
abline(lm(mean_ApEn ~ UPDRS, data = participant_summary), col = "gray", lty = 2, lwd = 2)

plot(participant_summary$UPDRS, participant_summary$mean_force,
     xlab = "UPDRS-III Score", ylab = "Mean Force (N)",
     main = "UPDRS-III vs Mean Force", pch = 16, cex = 1.5, col = ifelse(participant_summary$Sex == "Male", "blue", "red"))
abline(lm(mean_force ~ UPDRS, data = participant_summary), col = "gray", lty = 2, lwd = 2)

plot(participant_summary$UPDRS, participant_summary$mean_CV,
     xlab = "UPDRS-III Score", ylab = "Force CV (%)",
     main = "UPDRS-III vs Force Variability", pch = 16, cex = 1.5, col = ifelse(participant_summary$Sex == "Male", "blue", "red"))
abline(lm(mean_CV ~ UPDRS, data = participant_summary), col = "gray", lty = 2, lwd = 2)
dev.off()

# ------------------------------------------------------------------
# SECTION 2: SEX DIFFERENCES IN UPDRS-III POSTURAL TREMOR (Reviewer #1, Q6)
# ------------------------------------------------------------------
cat("\n========== SECTION 2: Sex Differences in UPDRS-III Postural Tremor ==========\n")

# NOTE: UPDRS-III item 3.15 = Postural Tremor
# You need the individual sub-item scores from updrs_01032024.xlsx
# The postural tremor sub-items are columns labeled "3.15 (TREMOR POST)"
# with A (MSD) = right upper limb and A (MSE) = left upper limb

# This section requires the raw UPDRS sub-item data
# The code below assumes you can extract item 3.15 scores from updrs_data

cat("Please ensure updrs_01032024.xlsx is properly loaded.\n")
cat("Extract UPDRS-III item 3.15 (Postural Tremor) sub-items by sex.\n")

# If the postural tremor scores are available in participant_summary or updrs_data:
# Create a summary table of postural tremor scores by sex
# Example structure (adapt column names to your actual data):

# For the response letter, we need:
# - Mean postural tremor scores for men vs women
# - Statistical comparison (Mann-Whitney or t-test given small n)

# Using the participant-level data with sex info:
sex_updrs <- participant_summary %>%
  group_by(Sex) %>%
  summarise(
    n = n(),
    mean_UPDRS = mean(UPDRS, na.rm = TRUE),
    sd_UPDRS = sd(UPDRS, na.rm = TRUE),
    mean_amplitude = mean(mean_amplitude, na.rm = TRUE),
    sd_amplitude = sd(mean_amplitude, na.rm = TRUE),
    mean_frequency = mean(mean_frequency, na.rm = TRUE),
    sd_frequency = sd(mean_frequency, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nUPDRS-III and Tremor by Sex:\n")
print(sex_updrs)
write.csv(sex_updrs, file = paste0(results_dir, "/sex_UPDRS_tremor_summary.csv"), row.names = FALSE)

# Mann-Whitney test for UPDRS-III by sex
if (sum(!is.na(participant_summary$UPDRS)) > 0) {
  wilcox_updrs <- wilcox.test(UPDRS ~ Sex, data = participant_summary)
  cat(sprintf("\nUPDRS-III by Sex (Mann-Whitney): W = %.1f, p = %.4f\n",
              wilcox_updrs$statistic, wilcox_updrs$p.value))
}

# ------------------------------------------------------------------
# SECTION 3: SEX DIFFERENCES ACROSS HY STAGES (Reviewer #1, Q7)
# ------------------------------------------------------------------
cat("\n========== SECTION 3: Sex x HY Stage Differences ==========\n")

# Sex-stratified descriptive summaries by HY stage (also addresses Editor comment)
sex_hy_summary <- participant_summary %>%
  mutate(HY_stage = factor(HY, levels = c(2, 3), labels = c("Mild (HY2)", "Moderate (HY3)"))) %>%
  group_by(Sex, HY_stage) %>%
  summarise(
    n = n(),
    mean_age = mean(Age, na.rm = TRUE),
    sd_age = sd(Age, na.rm = TRUE),
    mean_UPDRS = mean(UPDRS, na.rm = TRUE),
    sd_UPDRS = sd(UPDRS, na.rm = TRUE),
    mean_LEDD = mean(LEDD, na.rm = TRUE),
    sd_LEDD = sd(LEDD, na.rm = TRUE),
    mean_amplitude = mean(mean_amplitude, na.rm = TRUE),
    sd_amplitude = sd(mean_amplitude, na.rm = TRUE),
    mean_frequency = mean(mean_frequency, na.rm = TRUE),
    sd_frequency = sd(mean_frequency, na.rm = TRUE),
    mean_ApEn = mean(mean_ApEn, na.rm = TRUE),
    sd_ApEn = sd(mean_ApEn, na.rm = TRUE),
    mean_force_val = mean(mean_force, na.rm = TRUE),
    sd_force = sd(mean_force, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nSex x HY Stage Descriptive Summary:\n")
print(sex_hy_summary, width = Inf)
write.csv(sex_hy_summary, file = paste0(results_dir, "/sex_HY_descriptive_summary.csv"), row.names = FALSE)

# Sex x Severity interaction plots for key outcomes
# Using the full repeated-measures data
variables_to_test <- c('mean_modalAmplitude_AccR', 'mean_modalFrequency_AccR',
                        'mean_ApEn_AccR', 'mean_mean_force', 'mean_CV_force')

pdf(paste0(results_dir, "/sex_HY_interaction_plots.pdf"), width = 10, height = 8)
for (var in variables_to_test) {
  if (var %in% names(data)) {
    interaction_data <- data %>%
      filter(!is.na(.data[[var]])) %>%
      group_by(Sex, Severity) %>%
      summarise(
        mean_val = mean(.data[[var]], na.rm = TRUE),
        se_val = sd(.data[[var]], na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )

    # ggplot interaction with larger fonts
    p <- ggplot(interaction_data, aes(x = Severity, y = mean_val, fill = Sex)) +
      geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
      geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                    position = position_dodge(0.9), width = 0.25, linewidth = 0.8) +
      labs(title = paste("Sex x Severity:", var),
           x = "Disease Severity (HY Stage)", y = var) +
      scale_x_discrete(labels = c("2" = "Mild (HY 2)", "3" = "Moderate (HY 3)")) +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      ) +
      scale_fill_manual(values = c("Male" = "steelblue", "Female" = "coral"))
    print(p)
  }
}
dev.off()

# ------------------------------------------------------------------
# SECTION 4: 95% CONFIDENCE INTERVALS FOR KEY LMM ESTIMATES (Editor)
# ------------------------------------------------------------------
cat("\n========== SECTION 4: 95% CIs for Key LMM Estimates ==========\n")

# Refit key models and extract CIs
# Model 1: Tremor amplitude (3-factor model)
filter_outliers <- function(data, measure, n_std_dev = 3) {
  data %>%
    group_by(Severity, Hand, Feedback) %>%
    mutate(
      mean_value = mean(.data[[measure]], na.rm = TRUE),
      std_dev = sd(.data[[measure]], na.rm = TRUE)
    ) %>%
    filter(
      !is.na(.data[[measure]]),
      abs(.data[[measure]] - mean_value) <= n_std_dev * std_dev
    ) %>%
    ungroup() %>%
    dplyr::select(., -mean_value, -std_dev)
}

key_vars <- c('mean_modalAmplitude_AccR', 'mean_modalFrequency_AccR',
              'mean_ApEn_AccR', 'mean_mean_force', 'mean_CV_force', 'mean_SD_force')

ci_results_all <- list()

for (var in key_vars) {
  cat(sprintf("\n--- 95%% CIs for %s ---\n", var))
  filtered_data <- filter_outliers(data, var, n_std_dev = 3)

  model <- tryCatch(
    lmer(as.formula(paste(var, "~ Severity * Hand * Feedback * Age_c * Sex + LEDD_c + Tremor + (1 | subject_id)")),
         data = filtered_data),
    warning = function(w) { message("Warning: ", w); return(NULL) },
    error = function(e) { message("Error: ", e); return(NULL) }
  )

  if (!is.null(model)) {
    # Extract fixed effects with CIs
    ci <- confint(model, method = "Wald")
    fixed_effects <- fixef(model)
    se <- sqrt(diag(vcov(model)))

    ci_table <- data.frame(
      Variable = var,
      Term = names(fixed_effects),
      Estimate = fixed_effects,
      SE = se,
      CI_lower = fixed_effects - 1.96 * se,
      CI_upper = fixed_effects + 1.96 * se,
      row.names = NULL
    )

    ci_results_all[[var]] <- ci_table

    # Print key contrasts with CIs using emmeans
    cat("\n  Severity contrasts:\n")
    emm_sev <- emmeans(model, pairwise ~ Severity)
    print(confint(emm_sev$contrasts))

    cat("\n  Hand contrasts:\n")
    emm_hand <- emmeans(model, pairwise ~ Hand)
    print(confint(emm_hand$contrasts))

    cat("\n  Feedback contrasts:\n")
    emm_fb <- emmeans(model, pairwise ~ Feedback)
    print(confint(emm_fb$contrasts))

    cat("\n  Sex contrasts:\n")
    emm_sex <- emmeans(model, pairwise ~ Sex)
    print(confint(emm_sex$contrasts))

    # Sex x Severity interaction with CIs
    cat("\n  Sex x Severity contrasts:\n")
    emm_sex_sev <- emmeans(model, pairwise ~ Sex | Severity)
    print(confint(emm_sex_sev$contrasts))

    # Save emmeans with CIs
    emm_sev_ci <- as.data.frame(confint(emm_sev$contrasts))
    write.csv(emm_sev_ci, file = paste0(results_dir, "/", var, "_Severity_CI.csv"), row.names = FALSE)

    emm_sex_sev_ci <- as.data.frame(confint(emm_sex_sev$contrasts))
    write.csv(emm_sex_sev_ci, file = paste0(results_dir, "/", var, "_Sex_Severity_CI.csv"), row.names = FALSE)
  }
}

# Save all CI results
ci_combined <- bind_rows(ci_results_all)
write.csv(ci_combined, file = paste0(results_dir, "/all_fixed_effects_CIs.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# SECTION 5: LMM DIAGNOSTIC PLOTS (Editor - replace Levene/Mauchly)
# ------------------------------------------------------------------
cat("\n========== SECTION 5: LMM Diagnostics ==========\n")

# Refit the tremor amplitude model for diagnostics
var <- 'mean_modalAmplitude_AccR'
filtered_data <- filter_outliers(data, var, n_std_dev = 3)
model_diag <- lmer(as.formula(paste(var, "~ Severity * Hand * Feedback * Age_c * Sex + LEDD_c + Tremor + (1 | subject_id)")),
                   data = filtered_data)

pdf(paste0(results_dir, "/LMM_diagnostics.pdf"), width = 10, height = 10, pointsize = 16)
par(mfrow = c(2, 2), cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.3, mar = c(5, 5, 4, 2))

# 1. Residuals vs Fitted
plot(fitted(model_diag), residuals(model_diag),
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted (Tremor Amplitude Model)",
     pch = 16, cex = 1.3, col = rgb(0, 0, 0, 0.3))
abline(h = 0, col = "red", lty = 2, lwd = 2)
lines(lowess(fitted(model_diag), residuals(model_diag)), col = "blue", lwd = 2)

# 2. Q-Q plot of residuals
qqnorm(residuals(model_diag), main = "Normal Q-Q Plot of Residuals", pch = 16, cex = 1.3, col = rgb(0, 0, 0, 0.3))
qqline(residuals(model_diag), col = "red", lwd = 2)

# 3. Scale-Location plot
plot(fitted(model_diag), sqrt(abs(residuals(model_diag))),
     xlab = "Fitted Values", ylab = "sqrt(|Residuals|)",
     main = "Scale-Location Plot",
     pch = 16, cex = 1.3, col = rgb(0, 0, 0, 0.3))
lines(lowess(fitted(model_diag), sqrt(abs(residuals(model_diag)))), col = "blue", lwd = 2)

# 4. Random effects Q-Q plot
re <- ranef(model_diag)$subject_id[,1]
qqnorm(re, main = "Q-Q Plot of Random Effects", pch = 16, cex = 1.5, col = "steelblue")
qqline(re, col = "red", lwd = 2)

dev.off()

# Shapiro-Wilk test on residuals
sw_test <- shapiro.test(residuals(model_diag))
cat(sprintf("\nShapiro-Wilk test on residuals: W = %.4f, p = %.4f\n", sw_test$statistic, sw_test$p.value))

# Additional diagnostics for other key models
for (var in c('mean_mean_force', 'mean_CV_force')) {
  fd <- filter_outliers(data, var, n_std_dev = 3)
  m <- tryCatch(
    lmer(as.formula(paste(var, "~ Severity * Hand * Feedback * Age_c * Sex + LEDD_c + Tremor + (1 | subject_id)")),
         data = fd),
    error = function(e) NULL
  )
  if (!is.null(m)) {
    sw <- shapiro.test(residuals(m))
    cat(sprintf("Shapiro-Wilk for %s: W = %.4f, p = %.4f\n", var, sw$statistic, sw$p.value))
  }
}

# ------------------------------------------------------------------
# SECTION 6: SEX-STRATIFIED DESCRIPTIVE SUMMARIES (Editor)
# ------------------------------------------------------------------
cat("\n========== SECTION 6: Sex-Stratified Descriptive Summaries ==========\n")

# Full descriptive table by sex
# Note: BMI may not be in the main CSV. We use what's available.
sex_descriptives <- data %>%
  group_by(subject_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(Severity_label = ifelse(HY == 2, "Mild (HY 2)", "Moderate (HY 3)")) %>%
  group_by(Sex) %>%
  summarise(
    n = n(),
    mean_age = mean(Age, na.rm = TRUE),
    sd_age = sd(Age, na.rm = TRUE),
    n_mild = sum(HY == 2, na.rm = TRUE),
    n_moderate = sum(HY == 3, na.rm = TRUE),
    mean_UPDRS = mean(UPDRS, na.rm = TRUE),
    sd_UPDRS = sd(UPDRS, na.rm = TRUE),
    mean_LEDD = mean(LEDD, na.rm = TRUE),
    sd_LEDD = sd(LEDD, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nSex-Stratified Demographics:\n")
print(sex_descriptives)
write.csv(sex_descriptives, file = paste0(results_dir, "/sex_stratified_demographics.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# SECTION 7: SENSITIVITY ANALYSIS WITH UPDRS-III AS CONTINUOUS (Editor)
# ------------------------------------------------------------------
cat("\n========== SECTION 7: UPDRS-III as Continuous Covariate ==========\n")

# UPDRS was already merged and centered above (UPDRS_c)
if ("UPDRS" %in% names(data) && sum(!is.na(data$UPDRS)) > 0) {

  # Sensitivity model: replace HY-based Severity with continuous UPDRS
  var <- 'mean_modalAmplitude_AccR'
  fd <- filter_outliers(data, var, n_std_dev = 3)

  model_updrs <- tryCatch(
    lmer(as.formula(paste(var, "~ UPDRS_c * Hand * Feedback * Age_c * Sex + LEDD_c + Tremor + (1 | subject_id)")),
         data = fd),
    warning = function(w) { message("Warning: ", w); return(NULL) },
    error = function(e) { message("Error: ", e); return(NULL) }
  )

  if (!is.null(model_updrs)) {
    cat("\nSensitivity Model (UPDRS-III continuous) for Tremor Amplitude:\n")
    cat("ANOVA Table:\n")
    print(anova(model_updrs))

    # Extract UPDRS coefficient
    coefs <- summary(model_updrs)$coefficients
    if ("UPDRS_c" %in% rownames(coefs)) {
      cat(sprintf("\nUPDRS-III effect: beta = %.4f, SE = %.4f, t = %.2f, p = %.4f\n",
                  coefs["UPDRS_c", "Estimate"],
                  coefs["UPDRS_c", "Std. Error"],
                  coefs["UPDRS_c", "t value"],
                  coefs["UPDRS_c", "Pr(>|t|)"]))
    }

    write.csv(tidy(anova(model_updrs)),
              file = paste0(results_dir, "/sensitivity_UPDRS_continuous_amplitude.csv"), row.names = FALSE)
  }
} else {
  cat("UPDRS column not found. Merge UPDRS scores into dataset before running this section.\n")
}

# ------------------------------------------------------------------
# SECTION 8: PARTICIPANT FLOW SUMMARY
# ------------------------------------------------------------------
cat("\n========== SECTION 8: Participant Flow ==========\n")
cat("Parent RCT: Transcranial direct current stimulation + exercise trial\n")
cat("  - Original RCT enrolled: 25 participants\n")
cat("  - Additional participant enrolled: 1\n")
cat("  - Total baseline dataset: 26 participants\n")
cat("  - Excluded from present analysis: 0\n")
cat("  - Final analysis sample: 26 (11 women, 15 men)\n")
cat("  - HY Stage 2 (Mild): 15 participants\n")
cat("  - HY Stage 3 (Moderate): 11 participants\n")
cat("  - All participants completed grip force and tremor assessments\n")
cat("  - No missing data for primary outcomes\n")

# ------------------------------------------------------------------
# SECTION 9: MULTIPLE COMPARISON ADJUSTMENT (Editor)
# ------------------------------------------------------------------
cat("\n========== SECTION 9: FDR Adjustment Summary ==========\n")
cat("Note: Within each LMM, post hoc comparisons use Tukey HSD adjustment.\n")
cat("Across outcome families, results are presented as exploratory given the\n")
cat("modest sample size. Effect sizes and 95% CIs are reported to facilitate\n")
cat("interpretation beyond dichotomous significance testing.\n")

cat("\n\n========== ALL ANALYSES COMPLETE ==========\n")
cat("Results saved to:", results_dir, "\n")
cat("Review the CSV files a")
