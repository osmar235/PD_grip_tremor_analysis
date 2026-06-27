# Load necessary libraries
library(arm)
library(readxl)
library(dplyr)
library(ggplot2)
library(lme4)
library(broom.mixed)
library(tidyverse)
library(ggeffects)
library(sjPlot)
library(boot)
library(emmeans)
library(lmerTest)

# Function to filter outliers more than n standard deviations away from the mean
filter_outliers <- function(data, measure, n_std_dev) {
  data %>%
    group_by(Severity, Hand, Feedback, Freq) %>%
    mutate(
      mean_value = mean(.data[[measure]], na.rm = TRUE),  # Calculate mean
      std_dev = sd(.data[[measure]], na.rm = TRUE)  # Calculate standard deviation
    ) %>%
    filter(
      !is.na(.data[[measure]]),  # Ensure non-missing data points
      abs(.data[[measure]] - mean_value) <= n_std_dev * std_dev  # Filter within n std devs
    ) %>%
    ungroup() %>%
    #select(-mean_value, -std_dev)  # Remove temporary columns
    dplyr::select(., -mean_value, -std_dev)
}

# Set working directory and load data
data_path <- file.path("data", "Reshaped_PD_Tremor_with_Feedback_Freq_for_Stats_f.csv")

if (!file.exists(data_path)) {
  stop(
    "Data file not found. The participant-level dataset is not publicly available ",
    "due to consent restrictions. Please request the dataset from the corresponding author ",
    "and place it in the /data folder before running this script."
  )
}

data <- readr::read_csv(data_path)

# Clean column names and convert relevant columns to factors
names(data) <- make.names(names(data))
data$Severity <- factor(data$Severity)
data$Severity <- factor(data$HY)
#data$Feedback <- factor(data$Feedback)
data$Hand <- factor(data$Hand)
data$Freq <- factor(data$Freq)  # Add Freq as a factor
data$subject_id <- factor(data$subject_id)
# Center continuous variables
data$Age_c <- scale(data$Age, center = TRUE, scale = FALSE)
#data$HY_c <- scale(data$HY, center = TRUE, scale = FALSE)
data$LEDD_c <- scale(data$LEDD, center = TRUE, scale = FALSE)

# Ensure categorical variables are factors
data$Sex <- factor(data$Sex, levels = c("Male", "Female"))
data$Tremor <- factor(data$Tremor, levels = c("Right", "Left", "Equal"))
#data$HY <- factor(data$HY, levels = c("2", "3"))

# Define the list of variables to test
#variables <- c('mean_media_power_AccRxFORCE','mean_media_Rsq_AccRxFORCE','mean_media_power_AccR',
#        'mean_media_power_force','mean_media_Npower_AccR','mean_media_Npower_force','mean_band_AccR','mean_band_force','mean_Nband_AccR','mean_Nband_force')

variables <- c('mean_media_Rsq_AccRxFORCE')


# Filter data to remove rows with NA values in key variables
filtered_data <- data %>% filter(!is.na(mean_media_Npower_force))

# Directory to save the plots and results
results_dir <- file.path("results", "main_4factor_models")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Initialize lists to store results
anova_results <- list()
mixed_model_results <- list()
emmeans_results <- list()

# Initialize a list to store summary statistics
summary_stats <- list()

# Loop through each variable and generate plots, perform stats
for (var in variables) {
  measure <- var
  print(paste("Analyzing", measure))
  
  # Filter outliers for each variable
  filtered_data <- filter_outliers(data, measure, n_std_dev = 1)

  # Build Mixed Effects Model including Freq
  model <- tryCatch(
    lmer(as.formula(paste(measure, "~ Severity * Hand * Feedback * Freq + Age_c + Sex + LEDD_c + Tremor + (1 | subject_id)")), data = filtered_data),
    warning = function(w) { message("Warning: ", w); return(NULL) },
    error = function(e) { message("Error: ", e); return(NULL) }
  )
  
  # Only proceed if model was successfully fit
  if (!is.null(model)) {
    # Perform ANOVA on the model and store results
    anova_result <- tidy(anova(model))
    # New code (uses Kenward-Roger)
    #anova_result <- tidy(anova(model, ddf = "Kenward-Roger"))
    anova_results[[measure]] <- anova_result
    
    # Save detailed mixed model results
    mixed_model_results[[measure]] <- summary(model)
    
    # New code
    #emmeans_Severity <- emmeans(model, pairwise ~ Severity, lmer.df = "kenward-roger")
    #emmeans_Severity_Hand <- emmeans(model, pairwise ~ Severity * Hand, lmer.df = "kenward-roger")
    # ... and so on for all other emmeans calls

    # Perform pairwise comparisons using emmeans
    emmeans_Severity <- emmeans(model, pairwise ~ Severity)
    emmeans_Hand <- emmeans(model, pairwise ~ Hand)
    emmeans_Feedback <- emmeans(model, pairwise ~ Feedback)
    emmeans_Freq <- emmeans(model, pairwise ~ Freq)
    
    # Two-way interactions
    emmeans_Severity_Hand <- emmeans(model, pairwise ~ Severity * Hand)
    emmeans_Severity_Feedback <- emmeans(model, pairwise ~ Severity * Feedback)
    emmeans_Hand_Feedback <- emmeans(model, pairwise ~ Hand * Feedback)
    emmeans_Hand_Freq <- emmeans(model, pairwise ~ Hand | Freq)
    emmeans_Feedback_Freq <- emmeans(model, pairwise ~ Feedback | Freq)
    emmeans_Severity_Freq <- emmeans(model, pairwise ~ Severity * Freq)
    
    # Three-way interactions
    emmeans_Severity_Hand_Feedback <- emmeans(model, pairwise ~ Severity * Hand * Feedback)
    emmeans_Severity_Hand_Freq <- emmeans(model, pairwise ~ Severity * Hand * Freq)
    emmeans_Hand_Feedback_Freq <- emmeans(model, pairwise ~ Hand * Feedback * Freq)
    emmeans_Severity_Feedback_Freq <- emmeans(model, pairwise ~  Severity * Feedback | Freq )
    emmeans_Severity_Hand_Freq <- emmeans(model, pairwise ~  Severity * Hand | Freq )
    emmeans_Severity_Sex_Freq <- emmeans(model, pairwise ~  Severity * Sex | Freq )
    
    # Four-way interaction
    emmeans_Severity_Hand_Feedback_Freq <- emmeans(model, pairwise ~ Severity * Hand * Feedback * Freq)

    # Store the results
    emmeans_results[[measure]] <- list(
      Severity = as.data.frame(emmeans_Severity$contrasts),
      Hand = as.data.frame(emmeans_Hand$contrasts),
      Feedback = as.data.frame(emmeans_Feedback$contrasts),
      Freq = as.data.frame(emmeans_Freq$contrasts),
      Severity_Hand = as.data.frame(emmeans_Severity_Hand$contrasts),
      Severity_Feedback = as.data.frame(emmeans_Severity_Feedback$contrasts),
      Hand_Feedback = as.data.frame(emmeans_Hand_Feedback$contrasts),
      Hand_Freq = as.data.frame(emmeans_Hand_Freq$contrasts),
      Feedback_Freq = as.data.frame(emmeans_Feedback_Freq$contrasts),
      Severity_Freq = as.data.frame(emmeans_Severity_Freq$contrasts),
      Severity_Hand_Feedback = as.data.frame(emmeans_Severity_Hand_Feedback$contrasts),
      Severity_Hand_Freq = as.data.frame(emmeans_Severity_Hand_Freq$contrasts),
      Hand_Feedback_Freq = as.data.frame(emmeans_Hand_Feedback_Freq$contrasts),
      Severity_Feedback_Freq = as.data.frame(emmeans_Severity_Feedback_Freq$contrasts),
      Severity_Sex_Freq = as.data.frame(emmeans_Severity_Sex_Freq$contrasts),
      Severity_Hand_Feedback_Freq = as.data.frame(emmeans_Severity_Hand_Feedback_Freq$contrasts)
    )

    # Save the pairwise comparisons for each variable
    for (effect in names(emmeans_results[[measure]])) {
      write.csv(emmeans_results[[measure]][[effect]], file = paste0(results_dir, "/", measure, "_", effect, "_pairwise_comparisons.csv"), row.names = FALSE)
    }

    # Save mean, standard deviation, and standard error for each main effect and interaction
    summary_stats_main <- list(
      Severity = filtered_data %>%
        group_by(Severity) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Hand = filtered_data %>%
        group_by(Hand) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Feedback = filtered_data %>%
        group_by(Feedback) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Freq = filtered_data %>%
        group_by(Freq) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      HY = filtered_data %>%
        group_by(HY) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      

        )
    )
    
    summary_stats_two_way <- list(
      Severity_Hand = filtered_data %>%
        group_by(Severity, Hand) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Severity_Feedback = filtered_data %>%
        group_by(Severity, Feedback) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Hand_Feedback = filtered_data %>%
        group_by(Hand, Feedback) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Hand_Freq = filtered_data %>%
        group_by(Hand, Freq) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Feedback_Freq = filtered_data %>%
        group_by(Feedback, Freq) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Severity_Freq = filtered_data %>%
        group_by(Severity, Freq) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        )
    )

    
    
    summary_stats_three_way_Severity_Hand_Freq <- filtered_data %>%
      group_by(Severity, Hand, Freq) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )  

    summary_stats_three_way_SHF <- filtered_data %>%
      group_by(Severity, Hand, Feedback) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )

    summary_stats_three_way_SFF <- filtered_data %>%
      group_by(Severity, Feedback, Freq) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )  

    summary_stats_three_way_SSF <- filtered_data %>%
      group_by(Severity, Sex, Freq) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )  
    
    summary_stats_four_way <- filtered_data %>%
      group_by(Severity, Hand, Feedback, Freq) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )
    
    # Save main effect, two-way, three-way, and four-way summary stats to CSV
    write.csv(summary_stats_main$Severity, file = paste0(results_dir, "/", measure, "_Severity_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Hand, file = paste0(results_dir, "/", measure, "_Hand_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Feedback, file = paste0(results_dir, "/", measure, "_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Freq, file = paste0(results_dir, "/", measure, "_Freq_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$HY, file = paste0(results_dir, "/", measure, "_HY_summary_stats.csv"), row.names = FALSE)


    write.csv(summary_stats_two_way$Severity_Hand, file = paste0(results_dir, "/", measure, "_Severity_Hand_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Severity_Feedback, file = paste0(results_dir, "/", measure, "_Severity_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Hand_Feedback, file = paste0(results_dir, "/", measure, "_Hand_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Severity_Freq, file = paste0(results_dir, "/", measure, "_Severity_Freq_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Hand_Freq, file = paste0(results_dir, "/", measure, "_Hand_Freq_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Feedback_Freq, file = paste0(results_dir, "/", measure, "_Feedback_Freq_summary_stats.csv"), row.names = FALSE)

    write.csv(summary_stats_three_way_SHF, file = paste0(results_dir, "/", measure, "_Severity_Hand_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_three_way_SFF, file = paste0(results_dir, "/", measure, "_Severity_Feedback_Freq_summary_stats.csv"), row.names = FALSE)
    #write.csv(summary_stats_three_way_SSF, file = paste0(results_dir, "/", measure, "_Severity_Sex_Freq_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_four_way, file = paste0(results_dir, "/", measure, "_Severity_Hand_Feedback_Freq_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_three_way_Severity_Hand_Freq, file = paste0(results_dir, "/", measure, "_Severity_Hand_Freq_summary_stats.csv"), row.names = FALSE)


    
    # # Save mixed model results
    mixed_model_result <- summary(model)
    mixed_model_results[[measure]] <- mixed_model_result

      # Manual input: Specify the level of the Severity and the Hand where the asterisk should go
  asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_hand_level <- "D"  # Set "L" or "D" for Left or Dominant Hand
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually
  
  # Specify custom names for the levels
  custom_severity_names <- c("Mild", "Moderate")  # Example names
  custom_hand_names <- c("Right", "Left")  # Adjust for Hand factor

  # Manually specify axis labels
  x_axis_label <- "Severity Levels"
  #y_axis_label <- paste("Mean", measure)
  y_axis_label <- "Wavelet Acc Power  (v^2/Hz)"

  # Box Plot (Severity x Hand interaction) with manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Severity, y = .data[[measure]], fill = Hand)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of Severity and Hand"),
         x = x_axis_label, y = y_axis_label) +
         scale_y_continuous(limits = c(min(filtered_data[[measure]]), max(filtered_data[[measure]]) * 0.1)) + # Adjust y-axis limits
    theme_minimal() +
    scale_x_discrete(labels = custom_severity_names) +  # Custom x-axis labels
    scale_fill_discrete(labels = custom_hand_names) +  # Custom fill labels for Hand
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    ) +
    # Manually place asterisks above specific Severity and Hand levels
    geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & Hand == asterisk_hand_level),
              aes(x = Severity, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              position = position_dodge(width = 0.9), size = 8, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asteriskHAND.png"), plot = box_plot)


    #PLOT com FREQ

   # Manual input: Specify the levels of the Severity and Frequencies where the asterisks should go
  asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_freq_levels <- c("1")#c("2", "4")  # You can specify multiple frequency levels here
  asterisk_y_adjust <- 2  # Adjust the height of the asterisk manually

  # Specify custom names for the levels
  custom_severity_names <- c("Mild", "Moderate")  # Example names
  custom_freq_names <- c("0-1", "1-4", "4-8", "8-16")  # Adjust for Frequency levels

  # Manually specify axis labels
  x_axis_label <- "Frequency (Hz)"
  y_axis_label <- "X-Wavelet Power (V^2/Hz)"
  #y_axis_label <- "Acc Wavelet Power (V^2/Hz)"
  #y_axis_label <- "Acc Normalized Wavelet Power (%)"

  # Box Plot (Severity x Freq interaction) with manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Severity)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of Severity and Freq"),
        x = x_axis_label, y = y_axis_label) +
    theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +  # Custom x-axis labels
    scale_fill_discrete(labels = custom_severity_names) +  # Custom fill labels for Severity
    #scale_y_continuous(limits = c(min(filtered_data[[measure]]), max(filtered_data[[measure]]) * 0.1)) + # Adjust y-axis limits
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    )

  # Loop through each frequency level where an asterisk should be placed
  for (freq_level in asterisk_freq_levels) {
    # Add an asterisk at each specified frequency level
    box_plot <- box_plot +
      geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & Freq == freq_level),
                #aes(x = Freq, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
                aes(x = Freq, y = 0.15, label = "*"), #manual adjustment
                #aes(x = Freq, y = 20, label = "*"), #manual adjustment
                position = position_dodge(width = 0.9), size = 10, color = "black")
  }

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_testebox_plot_with_manual_asterisk_SeverityxFreq.png"), plot = box_plot)




  
    }


  #PLOT HAND FEEDBACK

  # Manual input: Specify the level of the Severity and the Hand where the asterisk should go
  asterisk_feedback_level <- "NFB"  # Asterisk only for 'Without Feedback'
  asterisk_hand_level <- "Dominant"  # Set "L" or "D" for Left or Dominant Hand
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually
  
  # Specify custom names for the levels
  custom_feedback_names <- c("With Feedback", "Without Feedback")  # Adjust for Feedback levels
  custom_hand_names <- c("Dominant", "Non-Dominant")  # Adjust for Hand factor

  # Manually specify axis labels
  x_axis_label <- "Feedback Levels"
  #y_axis_label <- paste("Mean", measure)
  y_axis_label <- "Wavelet Acc Power  (v^2/Hz)"

  # Box Plot (Feedback x Hand interaction) with manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Feedback, y = .data[[measure]], fill = Hand)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of Feedback and Hand"),
         x = x_axis_label, y = y_axis_label) +
         #scale_y_continuous(limits = c(min(filtered_data[[measure]]), max(filtered_data[[measure]]) * 0.1)) + # Adjust y-axis limits
    theme_minimal() +
    scale_x_discrete(labels = custom_feedback_names) +  # Custom x-axis labels
    scale_fill_discrete(labels = custom_hand_names) +  # Custom fill labels for Hand
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    ) +
    # Manually place asterisks above specific Severity and Hand levels
    geom_text(data = filtered_data %>% filter(Feedback == asterisk_feedback_level & Hand == asterisk_hand_level),
              aes(x = Feedback, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              position = position_dodge(width = 0.9), size = 8, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asteriskFeedbackxHAND.png"), plot = box_plot)





      # Manual input: Specify the levels of the Severity, Feedback, and Frequencies where the asterisk should go
  asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_freq_level <- "1"  # The frequency level where the asterisk should be placed
  asterisk_feedback_level <- "NFB"  # Asterisk only for 'Without Feedback'
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually

  # Specify custom names for the levels
  custom_severity_names <- c("Mild", "Moderate")  # Example names
  custom_freq_names <- c("0-1", "1-4", "4-8", "8-16")  # Adjust for Frequency levels
  custom_feedback_names <- c("With Feedback", "Without Feedback")  # Adjust for Feedback levels

  # Manually specify axis labels
  x_axis_label <- "Frequency (Hz)"
  #y_axis_label <- "Wavelet Coherence (0-1)" mean_media_Npower_force
  y_axis_label <- "Force Wavelet Normalized Power (%)"
  #y_axis_label <- "Wavelet Coherence (0-1)"
 
  # Box Plot (Severity x Freq interaction) with Feedback facet and manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Severity)) +
    #geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    geom_boxplot(
      outlier.color = "red",
      position = position_dodge(width = 0.8),
      color = "black",
      alpha = 1
    ) +
    geom_jitter(
      aes(color = Severity),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.6,
      size = 2,
      show.legend = FALSE
    ) +
    labs(title = paste("Interaciton between Severity, Feedback, and Freq on FWNP"),
        x = x_axis_label, y = y_axis_label) +
    facet_wrap(~ Feedback, labeller = labeller(Feedback = custom_feedback_names)) +  # Facet by Feedback
    theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +  # Custom x-axis labels for Frequencies
    scale_fill_discrete(labels = custom_severity_names) +  # Custom fill labels for Severity
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")  # Adjust facet label size
    )

  # Add an asterisk only for the 'Without Feedback' condition and specified frequency
  box_plot <- box_plot +
    geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & 
                                              Freq == asterisk_freq_level & 
                                              Feedback == asterisk_feedback_level),
              #aes(x = Freq, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              aes(x = Freq, y = 90, label = "*"),
              position = position_dodge(width = 0.9), size = 12, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asterisk_Severity_Feedback_Freq_01122024.png"), plot = box_plot)


  
# Define custom colors for Severity levels
  severity_colors <- c("Mild" = "#1f77b4",    # Blue color
                      "Moderate" = "#ff7f0e")  # Orange color

  # Box Plot with custom colors
  box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Severity)) +
    geom_boxplot(
      outlier.color = "red",
      position = position_dodge(width = 0.8),
      color = "black",
      alpha = 1
    ) +
    geom_jitter(
      aes(color = Severity),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.6,
      size = 2,
      show.legend = FALSE
    ) +
    labs(
      #title = "Interaction between Severity, Feedback, and Frequency on FWNP",
      title = "Interaction between Severity, Feedback, and Frequency on FWNP",
      x = x_axis_label,
      y = y_axis_label
    ) +
    facet_wrap(~ Feedback, labeller = labeller(Feedback = custom_feedback_names)) +
    theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +
    scale_fill_manual(values = severity_colors, labels = custom_severity_names) +
    scale_color_manual(values = severity_colors, labels = custom_severity_names) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")
    )

  # Add an asterisk only for the 'Without Feedback' condition and specified frequency
  box_plot <- box_plot +
    geom_text(
      data = filtered_data %>% filter(
        Severity == asterisk_severity_level &
          Freq == asterisk_freq_level &
          Feedback == asterisk_feedback_level
      ),
      aes(x = Freq, y = .5, label = "*"),
      position = position_dodge(width = 0.9),
      size = 12,
      color = "black"
    )

  # Save the box plot with asterisks
  ggsave(
    filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asterisk_Severity_Feedback_Freq_colors.png"),
    plot = box_plot
  )



#PLOT HAND FEEDBACK

  # Manual input: Specify the level of the Severity and the Hand where the asterisk should go
  asterisk_feedback_level <- "NFB"  # Asterisk only for 'Without Feedback'
  asterisk_hand_level <- "Dominant"  # Set "L" or "D" for Left or Dominant Hand
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually
  
  # Specify custom names for the levels
  custom_feedback_names <- c("With Feedback", "Without Feedback")  # Adjust for Feedback levels
  custom_hand_names <- c("Dominant", "Non-Dominant")  # Adjust for Hand factor

  # Manually specify axis labels
  x_axis_label <- "Feedback Levels"
  #y_axis_label <- paste("Mean", measure)
  y_axis_label <- "Wavelet Acc Power  (v^2/Hz)"

  # Box Plot (Feedback x Hand interaction) with manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Feedback, y = .data[[measure]], fill = Hand)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of Feedback and Hand"),
         x = x_axis_label, y = y_axis_label) +
         #scale_y_continuous(limits = c(min(filtered_data[[measure]]), max(filtered_data[[measure]]) * 0.1)) + # Adjust y-axis limits
    theme_minimal() +
    scale_x_discrete(labels = custom_feedback_names) +  # Custom x-axis labels
    scale_fill_discrete(labels = custom_hand_names) +  # Custom fill labels for Hand
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    ) +
    # Manually place asterisks above specific Severity and Hand levels
    geom_text(data = filtered_data %>% filter(Feedback == asterisk_feedback_level & Hand == asterisk_hand_level),
              aes(x = Feedback, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              position = position_dodge(width = 0.9), size = 8, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asteriskFeedbackxHAND.png"), plot = box_plot)





      # Manual input: Specify the levels of the Severity, Feedback, and Frequencies where the asterisk should go
  asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_freq_level <- c("1","2")  # The frequency level where the asterisk should be placed
  asterisk_feedback_level <- "NFB"  # Asterisk only for 'Without Feedback'
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually

  # Specify custom names for the levels
  custom_severity_names <- c("Mild", "Moderate")  # Example names
  custom_freq_names <- c("0-1", "1-4", "4-8", "8-16")  # Adjust for Frequency levels
  custom_feedback_names <- c("With Feedback", "Without Feedback")  # Adjust for Feedback levels

  # Manually specify axis labels
  x_axis_label <- "Frequency (Hz)"
  #y_axis_label <- "Wavelet Coherence (0-1)" mean_media_Npower_force
  y_axis_label <- "Force Wavelet Normalized Power (%)"
  #y_axis_label <- "Wavelet Coherence (0-1)"
 
  # Box Plot (Severity x Freq interaction) with Feedback facet and manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Severity)) +
    #geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    geom_boxplot(
      outlier.color = "red",
      position = position_dodge(width = 0.8),
      color = "black",
      alpha = 1
    ) +
    geom_jitter(
      aes(color = Severity),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.6,
      size = 2,
      show.legend = FALSE
    ) +
    labs(title = paste("Interaciton between Severity, Feedback, and Freq on FWNP"),
        x = x_axis_label, y = y_axis_label) +
    facet_wrap(~ Feedback, labeller = labeller(Feedback = custom_feedback_names)) +  # Facet by Feedback
    theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +  # Custom x-axis labels for Frequencies
    scale_fill_discrete(labels = custom_severity_names) +  # Custom fill labels for Severity
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")  # Adjust facet label size
    )

  # Add an asterisk only for the 'Without Feedback' condition and specified frequency
  box_plot <- box_plot +
    geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & 
                                              Freq == asterisk_freq_level & 
                                              Feedback == asterisk_feedback_level),
              #aes(x = Freq, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              aes(x = Freq, y = 90, label = "*"),
              position = position_dodge(width = 0.9), size = 12, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asterisk_Severity_Feedback_Freq_01122024.png"), plot = box_plot)


  
# Define custom colors for Severity levels
  severity_colors <- c("Mild" = "#1f77b4",    # Blue color
                      "Moderate" = "#ff7f0e")  # Orange color

  # Box Plot with custom colors
  box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Severity)) +
    geom_boxplot(
      outlier.color = "red",
      position = position_dodge(width = 0.8),
      color = "black",
      alpha = 1
    ) +
    geom_jitter(
      aes(color = Severity),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.6,
      size = 2,
      show.legend = FALSE
    ) +
    labs(
      #title = "Interaction between Severity, Feedback, and Frequency on FWNP",
      title = "Interaction between Severity, Feedback, and Frequency on COH",
      x = x_axis_label,
      y = "Wavelet Coherence"
    ) +
    facet_wrap(~ Feedback, labeller = labeller(Feedback = custom_feedback_names)) +
    theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +
    scale_fill_manual(values = severity_colors, labels = custom_severity_names) +
    scale_color_manual(values = severity_colors, labels = custom_severity_names) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")
    )

  # Add an asterisk only for the 'Without Feedback' condition and specified frequency
  box_plot <- box_plot +
    geom_text(
      data = filtered_data %>% filter(
        Severity == asterisk_severity_level &
          Freq == asterisk_freq_level &
          Feedback == asterisk_feedback_level
      ),
      aes(x = Freq, y = .5, label = "*"),
      position = position_dodge(width = 0.9),
      size = 12,
      color = "black"
    )

  # Save the box plot with asterisks
  ggsave(
    filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asterisk_Severity_Feedback_Freq_colors_COH.png"),
    plot = box_plot
  )






severity_colors <- c("Mild" = "#1f77b4",    # Blue color
                      "Moderate" = "#ff7f0e")  # Orange color

  # Manual input: Specify the levels of the Severity, Hand, and Frequencies where the asterisk should go
asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate"
asterisk_hand_level <- "Dominant"  # Set to "Dominant" or "Non-dominant" (match your factor levels)
asterisk_freq_level <- "3"  # The frequency level where the asterisk should be placed
asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually

# Specify custom names for the levels
custom_severity_names <- c("Mild", "Moderate")  # Example names
custom_freq_names <- c("0-1", "1-4", "4-8", "8-16")  # Adjust for Frequency levels
custom_hand_names <- c("Dominant", "Non-dominant")  # Adjust for Hand factor levels

# Manually specify axis labels
x_axis_label <- "Frequency (Hz)"
y_axis_label <- "Force Wavelet Normalized Power (%)"

# Box Plot (Severity x Hand x Freq interaction) with manual asterisk placement
box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Severity)) +
  geom_boxplot(
      outlier.color = "red",
      position = position_dodge(width = 0.8),
      color = "black",
      alpha = 1
    ) +
    geom_jitter(
      aes(color = Severity),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.6,
      size = 2,
      show.legend = FALSE
    ) +
  labs(title = paste("Interaction between Severity, Hand, and Frequency on FWNP"),
       x = x_axis_label, y = y_axis_label) +
  facet_wrap(~ Hand, labeller = labeller(Hand = custom_hand_names)) +  # Facet by Hand
  theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +
    scale_fill_manual(values = severity_colors, labels = custom_severity_names) +
    scale_color_manual(values = severity_colors, labels = custom_severity_names) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")  # Adjust facet label size
  )

# Add an asterisk only for the specified Hand, Severity, and Frequency levels
box_plot <- box_plot +
  geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & 
                                            Freq == asterisk_freq_level & 
                                            Hand == asterisk_hand_level),
            #aes(x = Freq, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
            aes(x = Freq, y = 89, label = "*"),
            position = position_dodge(width = 0.9), size = 12, color = "black")

# Save the box plot with asterisks
ggsave(filename = paste0(results_dir, "/", measure, "_Zbox_plot_with_manual_asterisk_Severity_Hand_Freq_12012024.png"), plot = box_plot)


 # Manual input: Specify the levels of the Severity, Hand, and Feedback where the asterisk should go
asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate"
asterisk_hand_level <- "Dominant"  # Set to "Dominant" or "Non-dominant" (match your factor levels)
asterisk_Feedback_level <- "NFB"  # The frequency level where the asterisk should be placed
asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually

# Specify custom names for the levels
custom_severity_names <- c("Mild", "Moderate")  # Example names
custom_feedback_names <- c("With Feedback", "Without Feedback")  # Adjust for Feedback levels
custom_hand_names <- c("Dominant", "Non-dominant")  # Adjust for Hand factor levels

# Manually specify axis labels
x_axis_label <- "Feedback"
y_axis_label <- "Force Wavelet Normalized Power (%)"

# Box Plot (Severity x Hand x feedback interaction) with manual asterisk placement
box_plot <- ggplot(filtered_data, aes(x = Feedback, y = .data[[measure]], fill = Severity)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
  labs(title = paste("Box Plot of Severity, Hand, and Feedback"),
       x = x_axis_label, y = y_axis_label) +
  facet_wrap(~ Hand, labeller = labeller(Hand = custom_hand_names)) +  # Facet by Hand
  theme_minimal() +
  scale_x_discrete(labels = custom_feedback_names) +  # Custom x-axis labels for feedback
  scale_fill_discrete(labels = custom_severity_names) +  # Custom fill labels for Severity
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 18, face = "bold")  # Adjust facet label size
  )

# Add an asterisk only for the specified Hand, Severity, and feedback levels
box_plot <- box_plot +
  geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & 
                                            Feedback == asterisk_feedback_level & 
                                            Hand == asterisk_hand_level),
            #aes(x = Feedback, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
            aes(x = Feedback, y = 89, label = "*"),
            position = position_dodge(width = 0.9), size = 8, color = "black")

# Save the box plot with asterisks
ggsave(filename = paste0(results_dir, "/", measure, "_Zbox_plot_with_manual_asterisk_Severity_Hand_Feedback.png"), plot = box_plot)


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure your data frame is named 'filtered_data' and contains the necessary variables:
# 'Severity', 'Freq', and 'mean_media_Npower_AccR'

# Convert variables to factors with specified levels
filtered_data$Severity <- factor(filtered_data$Severity, levels = c("Mild", "Moderate"))
filtered_data$Freq <- factor(filtered_data$Freq, levels = c("1", "2", "3", "4"))

# Specify custom labels for frequency bands
custom_freq_names <- c("0–1", "1–4", "4–8", "8–16")
names(custom_freq_names) <- c("1", "2", "3", "4")

# Specify colors for severity levels
severity_colors <- c("Mild" = "#1f77b4", "Moderate" = "#ff7f0e")

# Create the box plot
interaction_plot <- ggplot(filtered_data, aes(x = Freq, y = mean_media_Npower_AccR, fill = Severity)) +
  geom_boxplot(
    outlier.color = "red",
    position = position_dodge(width = 0.8),
    color = "black",
    alpha = 1
  ) +
  geom_jitter(
    aes(color = Severity),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 0.6,
    size = 1.5,
    show.legend = FALSE
  ) +
  labs(
    title = "Interaction between Severity and Frequency on AWNP",
    x = "Frequency Band (Hz)",
    y = "Acceleration Wavelet Normalized Power (%)",
    fill = "Disease Severity"
  ) +
  theme_minimal() +
  scale_x_discrete(labels = custom_freq_names) +
  scale_fill_manual(values = severity_colors) +
  scale_color_manual(values = severity_colors) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# Calculate the maximum y-value for the Moderate group at Freq 2 to position the asterisk
max_y <- filtered_data %>%
  filter(Severity == "Moderate", Freq == "2") %>%
  summarise(max_y = max(mean_media_Npower_AccR, na.rm = TRUE)) %>%
  pull(max_y)

# Adjust y position for the asterisk
asterisk_y <- 20  # Adjust this value as needed

# Create annotation data frame
annot_data <- data.frame(
  Freq = factor("2", levels = c("1", "2", "3", "4")),
  Severity = factor("Moderate", levels = c("Mild", "Moderate")),
  y = asterisk_y,
  label = "*"
)

# Add the asterisk to the plot
interaction_plot <- interaction_plot +
  geom_text(
    data = annot_data,
    aes(x = Freq, y = y, label = label, group = Severity),
    position = position_dodge(width = 0.8),
    color = "black",
    size = 12,
    vjust = -0.5
  )

# Save the plot
ggsave(
  filename = paste0(results_dir, "/mean_media_Npower_AccR_severity_freq_boxplot.png"),
  plot = interaction_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Display the plot
print(interaction_plot)


# Manual input: Specify the levels of the Severity, Feedback, and Frequencies where the asterisk should go
  asterisk_hand_level <- "Dominant"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_freq_level <- "1"  # The frequency level where the asterisk should be placed
  asterisk_feedback_level <- "NFB"  # Asterisk only for 'Without Feedback'
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually

  # Specify custom names for the levels
  custom_hand_names <- c("Dominant", "Non-dominant")  # Adjust for Hand factor levels
  custom_freq_names <- c("0-1", "1-4", "4-8", "8-16")  # Adjust for Frequency levels
  custom_feedback_names <- c("With Feedback", "Without Feedback")  # Adjust for Feedback levels

  # Manually specify axis labels
  x_axis_label <- "Frequency (Hz)"
  #y_axis_label <- "Wavelet Coherence (0-1)" mean_media_Npower_force
  y_axis_label <- "Force Wavelet Normalized Power (%)"
  # Box Plot (hand x Freq interaction) with Feedback facet and manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Freq, y = .data[[measure]], fill = Hand)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of hand, Feedback, and Freq"),
        x = x_axis_label, y = y_axis_label) +
    facet_wrap(~ Feedback, labeller = labeller(Feedback = custom_feedback_names)) +  # Facet by Feedback
    theme_minimal() +
    scale_x_discrete(labels = custom_freq_names) +  # Custom x-axis labels for Frequencies
    scale_fill_discrete(labels = custom_hand_names) +  # Custom fill labels for Severity
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18, face = "bold")  # Adjust facet label size
    )

  # Add an asterisk only for the 'Without Feedback' condition and specified frequency
  box_plot <- box_plot +
    geom_text(data = filtered_data %>% filter(Hand == asterisk_hand_level & 
                                              Freq == asterisk_freq_level & 
                                              Feedback == asterisk_feedback_level),
              aes(x = Freq, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              position = position_dodge(width = 0.9), size = 8, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asterisk_Hand_Feedback_Freq.png"), plot = box_plot)




}

# Save all ANOVA and mixed model results
if (length(anova_results) > 0) {
  anova_results_df <- do.call(rbind, lapply(names(anova_results), function(name) {
    anova_results[[name]] %>% mutate(measure = name)
  }))
  write.csv(anova_results_df, file = paste0(results_dir, "/anova_results.csv"), row.names = FALSE)
}

if (length(mixed_model_results) > 0) {
  mixed_model_results_df <- do.call(rbind, lapply(names(mixed_model_results), function(name) {
    tibble(
      measure = name,
      term = rownames(mixed_model_results[[name]]$coefficients),
      estimate = mixed_model_results[[name]]$coefficients[, "Estimate"],
      std.error = mixed_model_results[[name]]$coefficients[, "Std. Error"],
      df = if ("df" %in% colnames(mixed_model_results[[name]]$coefficients)) {
        mixed_model_results[[name]]$coefficients[, "df"]
      } else {
        NA
      },
      t.value = mixed_model_results[[name]]$coefficients[, "t value"],
      p.value = if ("Pr(>|t|)" %in% colnames(mixed_model_results[[name]]$coefficients)) {
        mixed_model_results[[name]]$coefficients[, "Pr(>|t|)"]
      } else {
        NA
      }
    )
  }))
  write.csv(mixed_model_results_df, file = paste0(results_dir, "/mixed_model_results.csv"), row.names = FALSE)
}
