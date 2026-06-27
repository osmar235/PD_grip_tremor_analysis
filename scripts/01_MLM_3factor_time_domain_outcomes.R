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
filter_outliers <- function(data, measure, n_std_dev = 5) {
  data %>%
    group_by(Severity, Hand, Feedback) %>%
    mutate(
      mean_value = mean(.data[[measure]], na.rm = TRUE),  # Calculate mean
      std_dev = sd(.data[[measure]], na.rm = TRUE)  # Calculate standard deviation
    ) %>%
    filter(
      !is.na(.data[[measure]]),  # Ensure non-missing data points
      abs(.data[[measure]] - mean_value) <= n_std_dev * std_dev  # Filter within n std devs
    ) %>%
    ungroup() %>%
    dplyr::select(., -mean_value, -std_dev)  # Remove temporary columns
}

# Set working directory and load data
data_path <- file.path("data", "Reshaped_PD_Tremor_with_Feedback_noFreq_for_Stats_f.csv")

if (!file.exists(data_path)) {
  stop(
    "Data file not found. The participant-level dataset is not publicly available ",
    "due to consent restrictions. Please request the dataset from the corresponding author ",
    "and place it in the /data folder before running this script."
  )
}

data <- readr::read_csv(data_path)

# Clean column names and convert relevant columns to factors
#Age 	Sex	HY	LEDD	Tremor	Dominance
# Load necessary libraries
library(lme4)
library(car)


names(data) <- make.names(names(data))
#data$Severity <- factor(data$Severity)
data$Severity <- factor(data$HY, levels = c("2", "3"))
data$Feedback <- factor(data$Feedback)
data$Hand <- factor(data$Hand)
data$subject_id <- factor(data$subject_id)

# Center continuous variables
data$Age_c <- scale(data$Age, center = TRUE, scale = FALSE)
#data$HY_c <- scale(data$HY, center = TRUE, scale = FALSE)
data$LEDD_c <- scale(data$LEDD, center = TRUE, scale = FALSE)

# Ensure categorical variables are factors
data$Sex <- factor(data$Sex, levels = c("Male", "Female"))
data$Tremor <- factor(data$Tremor, levels = c("Right", "Left", "Equal"))



# Define the list of variables to test (remove AccGroup and Freq)
#variables <- c('mean_modalAmplitude_AccR', 'mean_mean_force', 'mean_rms_force')  # Add other variables of interest here
#variables <- c('mean_mean_force','norm_mean_force','mean_SD_force','mean_CV_force','mean_rms_force',
#                    'mean_modalFrequency_AccR','mean_modalAmplitude_AccR','mean_ApEn_AccR','mean_SampEnt_AccR')# List of variables to analyze
variables <- c('mean_ApEn_AccR')

# Filter data to remove rows with NA values in key variables
#filtered_data <- data %>% filter(!is.na(mean_modalAmplitude_AccR))

# Directory to save the plots and results
results_dir <- file.path("results", "main_3factor_models")
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
  filtered_data <- filter_outliers(data, measure, n_std_dev = 3)

  # Build Mixed Effects Model
  model <- tryCatch(
    lmer(as.formula(paste(measure, "~ Severity * Hand * Feedback * Age_c * Sex + LEDD_c + Tremor + (1 | subject_id)")), data = filtered_data),
    warning = function(w) { message("Warning: ", w); return(NULL) },
    error = function(e) { message("Error: ", e); return(NULL) }
  )
  
  #simple_slopes <- sim_slopes(model, pred = "Age_c", modx = "Hand", mod2 = "Severity")
    #simple_slopes <- sim_slopes(model, pred = "Age_c", modx = "Sex")

    # Print the default output which shows estimates, SE, t-values, degrees of freedom and p-values.
    #print(simple_slopes)  

  #model <- lmer(
  #  measure ~ Severity * Hand * Feedback + Age_c + Sex + HY_c + LEDD_c + Tremor + (1 | subject_id),
  #  data = filtered_data
  #)

  # Summarize the model
  summary(model)

  # Check for multicollinearity
  #vif_model <- lm(measure ~ Severity * Hand * Feedback + Age_c + Sex + HY_c + LEDD_c + Tremor, data = filtered_data)
  #vif(vif_model)

  # Only proceed if model was successfully fit
  if (!is.null(model)) {
    # Perform ANOVA on the model and store results
    anova_result <- tidy(anova(model))
    anova_results[[measure]] <- anova_result
    
    # Save detailed mixed model results
    mixed_model_results[[measure]] <- summary(model)
    
    # Perform pairwise comparisons using emmeans
    emmeans_Severity <- emmeans(model, pairwise ~ Severity)
    emmeans_Hand <- emmeans(model, pairwise ~ Hand)
    emmeans_Feedback <- emmeans(model, pairwise ~ Feedback)
    emmeans_Tremor <- emmeans(model, pairwise ~ Tremor)
    
    # Two-way interactions
    emmeans_Severity_Hand <- emmeans(model,  pairwise ~ Hand | Severity)
    emmeans_Severity_Feedback <- emmeans(model, pairwise ~ Severity * Feedback)
    emmeans_Hand_Feedback <- emmeans(model, pairwise ~ Hand * Feedback)
    emmeans_Severity_Sex <- emmeans(model, pairwise ~ Severity * Sex)
    emmeans_Hand_Sex <- emmeans(model, pairwise ~ Hand * Sex)
    # Three-way interaction
    emmeans_Severity_Hand_Sex <- emmeans(model, pairwise ~ Severity * Hand | Sex)
    emmeans_Severity_Sex_Feedback <- emmeans(model, pairwise ~ Severity * Sex * Feedback)

    # Store the results
    emmeans_results[[measure]] <- list(
      Severity = as.data.frame(emmeans_Severity$contrasts),
      Hand = as.data.frame(emmeans_Hand$contrasts),
      Feedback = as.data.frame(emmeans_Feedback$contrasts),
      Tremor = as.data.frame(emmeans_Tremor$contrasts),
      Severity_Hand = as.data.frame(emmeans_Severity_Hand$contrasts),
      Severity_Feedback = as.data.frame(emmeans_Severity_Feedback$contrasts),
      Hand_Feedback = as.data.frame(emmeans_Hand_Feedback$contrasts),
      Severity_Sex = as.data.frame(emmeans_Severity_Sex$contrasts),
      Hand_Sex = as.data.frame(emmeans_Hand_Sex$contrasts),
      Severity_Hand_Sex = as.data.frame(emmeans_Severity_Hand_Sex$contrasts),
      Severity_Sex_Feedback = as.data.frame(emmeans_Severity_Sex_Feedback$contrasts)
    )

    # Save the pairwise comparisons for each variable
    for (effect in names(emmeans_results[[measure]])) {
      write.csv(emmeans_results[[measure]][[effect]], file = paste0(results_dir, "/", measure, "_", effect, "_pairwise_comparisons.csv"), row.names = FALSE)
    }

    # Save mean, standard deviation, and standard error for each main effect and interaction
    summary_stats_main <- list(
      Age_c = filtered_data %>%
        group_by(Age_c) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
        Sex = filtered_data %>%
        group_by(Sex) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
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
      Tremor = filtered_data %>%
        group_by(Tremor) %>%
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
      Sex_Hand = filtered_data %>%
        group_by(Sex, Hand) %>%
        summarise(
          mean = mean(.data[[measure]], na.rm = TRUE),
          sd = sd(.data[[measure]], na.rm = TRUE),
          se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
        ),
      Sex_Severity = filtered_data %>%
        group_by(Sex, Severity) %>%
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
        )
    )
    
    summary_stats_three_way1 <- filtered_data %>%
      group_by(Severity, Hand, Feedback) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )
    
    summary_stats_three_way2 <- filtered_data %>%
      group_by(Severity, Hand, Sex) %>%
      summarise(
        mean = mean(.data[[measure]], na.rm = TRUE),
        sd = sd(.data[[measure]], na.rm = TRUE),
        se = sd(.data[[measure]], na.rm = TRUE) / sqrt(n())
      )

    # --- NEW CODE BLOCK TO EXTRACT AND SAVE AGE RESULTS ---

    # 1. Get the summary of the model's fixed effects coefficients
    model_summary_coef <- summary(model)$coefficients

    # 2. Check if 'Age_c' is in the model and extract its results
    if ("Age_c" %in% rownames(model_summary_coef)) {
      
      # Extract the row for Age_c
      age_results <- model_summary_coef["Age_c", ]
      
      # 3. Format the results into your desired string format
      #    Using sprintf() for clean formatting of numbers
      age_formatted_string <- sprintf(
        "age (β = %.3f per year, SE = %.3f, t(%.1f) = %.2f, p = %.4f)",
        age_results["Estimate"],
        age_results["Std. Error"],
        age_results["df"],
        age_results["t value"],
        age_results["Pr(>|t|)"]
      )
      
      # Print the string to the console to see it during the run
      print(age_formatted_string)
      
      # 4. Save the formatted string to a text file
      write(age_formatted_string, file = paste0(results_dir, "/", measure, "_age_results.txt"))
      
    } else {
      print("Age_c not found in the model summary.")
    }

    # --- END OF NEW CODE BLOCK ---
    # Save main effect, two-way, and three-way summary stats to CSV
    write.csv(summary_stats_main$Age_c, file = paste0(results_dir, "/", measure, "_Age_c_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Tremor, file = paste0(results_dir, "/", measure, "_Tremor_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Sex, file = paste0(results_dir, "/", measure, "_Sex_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Severity, file = paste0(results_dir, "/", measure, "_Severity_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Hand, file = paste0(results_dir, "/", measure, "_Hand_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_main$Feedback, file = paste0(results_dir, "/", measure, "_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Sex_Hand, file = paste0(results_dir, "/", measure, "_Sex_Hand_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Sex_Severity, file = paste0(results_dir, "/", measure, "_Sex_Severity_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Severity_Hand, file = paste0(results_dir, "/", measure, "_Severity_Hand_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Severity_Feedback, file = paste0(results_dir, "/", measure, "_Severity_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_two_way$Hand_Feedback, file = paste0(results_dir, "/", measure, "_Hand_Feedback_summary_stats.csv"), row.names = FALSE)

    write.csv(summary_stats_three_way1, file = paste0(results_dir, "/", measure, "_Severity_Hand_Feedback_summary_stats.csv"), row.names = FALSE)
    write.csv(summary_stats_three_way2, file = paste0(results_dir, "/", measure, "_Severity_Hand_Sex_summary_stats.csv"), row.names = FALSE)
    
    # # Save mixed model results
    mixed_model_result <- summary(model)
    mixed_model_results[[measure]] <- mixed_model_result

    # Manual input: Specify the level of the Severity and the Hand where the asterisk should go
  asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_hand_level <- "Dominant"  # Set "L" or "D" for Left or Dominant Hand
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually
  
  # Specify custom names for the levels
  custom_severity_names <- c("Mild", "Moderate")  # Example names
  custom_hand_names <- c("Dominant", "Non-Dominant")  # Adjust for Hand factor

  # Manually specify axis labels
  x_axis_label <- "Severity Levels"
  #y_axis_label <- paste("Mean", measure)
  y_axis_label <- "Mean Force (N)"

  # Box Plot (Severity x Hand interaction) with manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Severity, y = .data[[measure]], fill = Hand)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of Severity and Hand"),
         x = x_axis_label, y = y_axis_label) +
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
  ggsave(filename = paste0(results_dir, "/", measure, "_box_plot_with_manual_asterisk.png"), plot = box_plot)



 # Manual input: Specify the level of the Severity and the Hand where the asterisk should go
  asterisk_severity_level <- "Moderate"  # Can set to "Mild" or "Moderate" (match your factor levels)
  asterisk_feedback_level <- "FB"  # Set "L" or "D" for Left or Dominant Hand
  asterisk_y_adjust <- 0.1  # Adjust the height of the asterisk manually
  
  # Specify custom names for the levels
  custom_severity_names <- c("Mild", "Moderate")  # Example names
  custom_feedback_names <- c("With", "Without")  # Adjust for Hand factor

  # Manually specify axis labels
  x_axis_label <- "Severity Levels"
  #y_axis_label <- paste("Mean", measure)
  y_axis_label <- "Mean Modal Frequency (Hz)"
  #y_axis_label <- "CV of Force (%)"

  # Box Plot (Severity x Hand interaction) with manual asterisk placement
  box_plot <- ggplot(filtered_data, aes(x = Severity, y = .data[[measure]], fill = Feedback)) +
    geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
    labs(title = paste("Box Plot of Severity and Feedback"),
         x = x_axis_label, y = y_axis_label) +
    theme_minimal() +
    scale_x_discrete(labels = custom_severity_names) +  # Custom x-axis labels
    scale_fill_discrete(labels = custom_feedback_names) +  # Custom fill labels for Hand
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16)
    ) +
    # Manually place asterisks above specific Severity and Hand levels
    geom_text(data = filtered_data %>% filter(Severity == asterisk_severity_level & Feedback == asterisk_feedback_level),
              aes(x = Severity, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
              position = position_dodge(width = 0.9), size = 8, color = "black")

  # Save the box plot with asterisks
  ggsave(filename = paste0(results_dir, "/", measure, "_FeedbackxSeverity_box_plot_with_manual_asterisk.png"), plot = box_plot)
}

    # Specify axis labels and custom names for levels
    x_axis_label <- "Feedback Levels"
    y_axis_label <- "Std of Force (N)"
    custom_feedback_names <- c("With", "Without")  # Example custom labels
    # Boxplot for visualization
    box_plot <- ggplot(data, aes(x = Feedback, y = .data[[measure]], fill = Feedback)) +
        geom_boxplot(outlier.color = "red", outlier.shape = 16, position = position_dodge(width = 0.9)) +
        labs(title = "Box Plot of Std of Force by Feedback",
            x = x_axis_label, y = y_axis_label) +
        scale_x_discrete(labels = custom_feedback_names) +  # Custom x-axis labels
        #scale_y_continuous(limits = c(min(data[[measure]]), max(data[[measure]]) * 0.001)) + # Adjust y-axis limits
        theme_minimal() +
        theme(
            plot.title = element_text(size = 20, face = "bold"),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 16),
            legend.position = "none"
    )
    # Optional: Manually place asterisks if needed based on significance
    asterisk_feedback_level <- "FB"  # Can set to "Mild" or "Moderate"
    asterisk_y_adjust <- 5  # Adjust the height of the asterisk manually

    # If the test result is significant, place the asterisk
    #if (test_result$p_value < 0.05) {
    #box_plot_with_asterisk <- box_plot + 
    #    geom_text(data = data %>% filter(Feedback == asterisk_feedback_level),
    #            aes(x = Feedback, y = max(.data[[measure]]) + asterisk_y_adjust, label = "*"),
    #            size = 8, color = "black")
    
    # Show where the plot with asterisks will be saved
    #asterisk_plot_file_path <- paste0(results_dir, "/zMainEffect_box_plot_with_asterisk.png")
    #print(paste0("Saving box plot with asterisk to: ", asterisk_plot_file_path))
    
    # Save the box plot with the asterisk
    #ggsave(filename = asterisk_plot_file_path, plot = box_plot_with_asterisk)
    #}

  # Load necessary libraries
  library(ggplot2)
  library(dplyr)

  # Ensure that your data frame is named 'filtered_data' and contains the necessary variables
  # 'Severity', 'Sex', and 'mean_modalAmplitude_AccR'

  # Convert variables to factors if they are not already
  filtered_data$Severity <- factor(filtered_data$Severity, levels = c("Mild", "Moderate"))
  filtered_data$Sex <- factor(filtered_data$Sex, levels = c("Male", "Female"))

  # Create summary statistics for each group (optional, for adding mean points)
  summary_stats <- filtered_data %>%
    group_by(Severity, Sex) %>%
    summarise(
      mean_value = mean(mean_modalAmplitude_AccR, na.rm = TRUE),
      se = sd(mean_modalAmplitude_AccR, na.rm = TRUE) / sqrt(n())
    )

  # Create the box plot
  interaction_plot <- ggplot(filtered_data, aes(x = Severity, y = mean_modalAmplitude_AccR, fill = Sex)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    geom_jitter(
      aes(color = Sex),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      alpha = 0.6,
      size = 2,
      show.legend = FALSE
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 18,
      size = 3,
      position = position_dodge(width = 0.8),
      color = "black"
    ) +
    labs(
      title = "Interaction between Severity and Sex on Acc MMFA",
      x = "Disease Severity",
      y = expression(paste("Acc Mean Modal Frequency Amplitude (", V^2, "/Hz)")),
      fill = "Sex"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size  = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
    scale_color_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e"))

  # Save the plot
  ggsave(
    filename = paste0(results_dir, "/mean_modalAmplitude_AccR_severity_sex_boxplot.png"),
    plot = interaction_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  # Display the plot
  print(interaction_plot)


}

# Load necessary libraries
library(ggplot2)
library(lme4)
library(ggeffects)
library(dplyr)

# Assuming 'filtered_data' is your dataset after processing
# Fit the linear mixed-effects model
#model <- lmer(
#  mean_mean_force ~ Severity * Hand * Age_c + Sex + LEDD_c + Tremor + (1 | subject_id),
#  data = filtered_data
#)

# Generate predicted values for plotting
# Use the ggpredict function to obtain predictions with confidence intervals
library(ggeffects)
pred_data <- ggpredict(model, terms = c("Age_c [all]", "Hand", "Severity"))

# Convert centered Age back to actual Age for plotting
mean_age <- mean(data$Age, na.rm = TRUE)
pred_data$x <- pred_data$x + mean_age

# Plot the interaction
interaction_plot <- ggplot(pred_data, aes(x = x, y = predicted, color = group, linetype = facet)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    title = "Predicted Mean Force by Age, Hand Dominance, and Disease Severity",
    x = "Age (years)",
    y = "Predicted Mean Force (N)",
    color = "Hand Dominance",
    fill = "Hand Dominance",
    linetype = "Disease Severity"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "dashed"))

# Save the interaction plot
ggsave(
  filename = paste0(results_dir, "/mean_mean_force_interaction_plot.png"),
  plot = interaction_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Optionally, display the plot
print(interaction_plot)

# Run ggpredict with Sex included
pred_data <- ggpredict(model, terms = c("Age_c [all]", "Hand", "Severity", "Sex"))

# Check the structure of pred_data
print("a partir daqui")
str(pred_data)
head(pred_data)

# Load necessary libraries
library(ggplot2)
library(lme4)
library(ggeffects)
library(dplyr)

# Assuming 'filtered_data' is your dataset after processing
# Fit the linear mixed-effects model with the four-way interaction
#model <- lmer(
#  mean_mean_force ~ Severity * Hand * Age_c * Sex + LEDD_c + Tremor + (1 | subject_id),
#  data = filtered_data
#)

# Generate predicted values for plotting
# Use the ggpredict function to obtain predictions with confidence intervals
pred_data <- ggpredict(
  model,
  terms = c("Age_c [all]", "Hand", "Severity", "Sex"),
  type = "fe"  # Fixed effects only
)

# Convert centered Age back to actual Age for plotting
mean_age <- mean(data$Age, na.rm = TRUE)
pred_data$x <- pred_data$x + mean_age

# Adjust variable names based on pred_data structure
pred_data$Hand <- factor(pred_data$group, levels = c("Dominant", "Non-Dominant"))
pred_data$Severity <- factor(pred_data$facet, levels = c("Mild", "Moderate"))
pred_data$Sex <- factor(pred_data$panel, levels = c("Male", "Female"))

# Now, proceed to plot
interaction_plot <- ggplot(pred_data, aes(x = x, y = predicted, color = Sex, linetype = Hand)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Sex), alpha = 0.1, color = NA) +
  labs(
    title = "Predicted Mean Force by Age, Hand Dominance, Sex, and Disease Severity",
    x = "Age (years)",
    y = "Predicted Mean Force (N)",
    color = "Sex",
    fill = "Sex",
    linetype = "Hand Dominance"
  ) +
  facet_wrap(~ Severity) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  scale_linetype_manual(values = c("solid", "dashed"))

# Save the interaction plot
ggsave(
  filename = paste0(results_dir, "/mean_force_four_way_interaction_plot.png"),
  plot = interaction_plot,
  width = 12,
  height = 8,
  dpi = 300
)

# Optionally, display the plot
print(interaction_plot)



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

### --- ADDITIONAL GRAPHS --- ###

# (1) Plot for the main effect of Age_c
# Generate predicted values along the range of Age_c
pred_age <- ggpredict(model, terms = "Age_c [all]")
# Convert the centered Age back to the original scale by adding the mean Age
mean_age <- mean(data$Age, na.rm = TRUE)
pred_age$x <- pred_age$x + mean_age

age_plot <- ggplot(pred_age, aes(x = x, y = predicted)) +
  geom_line(size = 1.2, color = "darkgreen") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "darkgreen", alpha = 0.2) +
  labs(
    title = "Predicted Outcome by Age",
    x = "Age (years)",
    y = "Predicted Outcome"
  ) +
  theme_minimal()

ggsave(filename = paste0(results_dir, "/", measure, "_predicted_by_age.png"),
       plot = age_plot, width = 8, height = 6, dpi = 300)
print(age_plot)

# (2) Plot for the Hand x Sex interaction
# Here we request predictions only for Sex and Hand.
# With two terms, the first term ("Sex") becomes the x-axis variable and the second ("Hand") is used to distinguish groups.
pred_hand_sex <- ggpredict(model, terms = c("Sex", "Hand"))

hand_sex_plot <- ggplot(pred_hand_sex, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(
    title = "Predicted Outcome by Sex and Hand",
    x = "Sex",
    y = "Predicted Outcome",
    color = "Hand"
  ) +
  theme_minimal()

ggsave(filename = paste0(results_dir, "/", measure, "_predicted_by_hand_sex.png"),
       plot = hand_sex_plot, width = 8, height = 6, dpi = 300)
print(hand_sex_plot)

# (3) Plot for the Severity x Hand x Sex interaction
# Request predictions for three terms; with three terms, ggpredict returns:
#   - x: the first term ("Severity")
#   - group: the second term ("Hand")
#   - facet: the third term ("Sex")
pred_sev_hand_sex <- ggpredict(model, terms = c("Severity", "Hand", "Sex"))

# Create the plot: here we map x to Severity, color to Hand, and use facet_wrap on the facet variable (Sex).
sev_hand_sex_plot <- ggplot(pred_sev_hand_sex, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(
    title = "Predicted Outcome by Severity, Hand, and Sex",
    x = "Severity",
    y = "Predicted Outcome",
    color = "Hand"
  ) +
  facet_wrap(~ facet) +
  theme_minimal()

ggsave(filename = paste0(results_dir, "/", measure, "_predicted_by_severity_hand_sex.png"),
       plot = sev_hand_sex_plot, width = 8, height = 6, dpi = 300)
print(sev_hand_sex_plot)


# Convert Severity so that 2 becomes "Mild" and 3 becomes "Moderate"
filtered_data$Severity <- factor(filtered_data$Severity, levels = c("2", "3"), labels = c("Mild", "Moderate"))

# --- Three-Way Interaction Plot using filtered_data (Severity × Hand × Sex) --- #

# Manual settings and aesthetics:
severity_colors   <- c("2" = "#1f77b4", "3" = "#ff7f0e")
custom_hand_names <- c("Dominant", "Non-Dominant")
asterisk_severity_level <- "3"   # Level where an asterisk should be added.
asterisk_hand_level     <- "Dominant"    # Hand condition for the asterisk.
asterisk_sex            <- "Female"      # Only add the asterisk for the Female subgroup.
asterisk_y_adjust       <- 0.1           # Vertical offset for the asterisk.

# Create the three-way interaction boxplot:
three_way_plot <- ggplot(filtered_data, aes(x = Severity, y = .data[[measure]], fill = Hand)) +
  geom_boxplot(outlier.color = "red",
               position = position_dodge(width = 0.8),
               color = "black", alpha = 1) +
  geom_jitter(aes(color = Hand),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 2, show.legend = FALSE) +
  labs(title = "Interaction between Severity, Hand, and Sex on Measure",
       x = "Severity", 
       y = paste("Mean", measure)) +
  facet_wrap(~ Sex, labeller = labeller(Sex = c("Male" = "Male", "Female" = "Female"))) +
  theme_minimal() +
  scale_fill_manual(values = severity_colors, labels = custom_hand_names) +
  scale_color_manual(values = severity_colors, labels = custom_hand_names) +
  theme(plot.title   = element_text(size = 20, face = "bold"),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 16),
        strip.text   = element_text(size = 18, face = "bold"))

# Add an asterisk for the specified condition (Female, Moderate severity, Dominant hand):
three_way_plot <- three_way_plot +
  geom_text(
    data = filtered_data %>% 
      filter(Severity == asterisk_severity_level,
             Hand == asterisk_hand_level,
             Sex == asterisk_sex),
    aes(x = Severity, y = max(.data[[measure]], na.rm = TRUE) + asterisk_y_adjust, label = "*"),
    position = position_dodge(width = 0.8),
    size = 12, color = "black"
  )

# Save and display the plot:
ggsave(filename = paste0(results_dir, "/", measure, "_three_way_interaction_filtered_data_with_asterisk.png"),
       plot = three_way_plot, width = 10, height = 8, dpi = 300)
print(three_way_plot)

pred_sev_hand_sex <- ggpredict(model, terms = c("Severity", "Hand", "Sex"))

# Recode Severity: "2" becomes "Mild" and "3" becomes "Moderate"
pred_sev_hand_sex$x <- factor(pred_sev_hand_sex$x, levels = c("2", "3"), labels = c("Mild", "Moderate"))


library(ggplot2)
library(dplyr)

# Manual settings:
severity_colors   <- c("Mild" = "#1f77b4", "Moderate" = "#ff7f0e")
custom_hand_names <- c("Dominant", "Non-Dominant")
custom_severity_names2 <- c("Mild", "Moderate")
asterisk_severity_level <- "Moderate"
asterisk_hand_level     <- "Dominant"
asterisk_sex            <- "Female"
asterisk_y_adjust       <- 0.1

# Create the boxplot with jitter:
three_way_plot <- ggplot(filtered_data, aes(x = Severity, y = .data[[measure]], fill = Hand)) +
  geom_boxplot(outlier.color = "red",
               position = position_dodge(width = 0.8),
               color = "black", alpha = 1) +
  geom_jitter(aes(color = Hand),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 2, show.legend = FALSE) +
  labs(title = "Interaction between Severity, Hand, and Sex on Measure",
       x = "Severity", y = paste("Mean", measure)) +
  facet_wrap(~ Sex, labeller = labeller(Sex = c("Male" = "Male", "Female" = "Female"))) +
  theme_minimal() +
  scale_x_discrete(labels = custom_severity_names2) +  # Custom x-axis labels
  scale_fill_manual(values = severity_colors, labels = custom_hand_names) +
  scale_color_manual(values = severity_colors, labels = custom_hand_names) +
  theme(plot.title   = element_text(size = 20, face = "bold"),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 16),
        strip.text   = element_text(size = 18, face = "bold"))

# Add an asterisk for the condition: (Female, Moderate severity, Dominant hand)
three_way_plot <- three_way_plot +
  geom_text(
    data = filtered_data %>% 
      filter(Severity == asterisk_severity_level,
             Hand == asterisk_hand_level,
             Sex == asterisk_sex),
    aes(x = Severity, y = max(.data[[measure]], na.rm = TRUE) + asterisk_y_adjust, label = "*"),
    position = position_dodge(width = 0.8),
    size = 12, color = "black"
  )

# Save and display:
ggsave(filename = paste0(results_dir, "/", measure, "_three_way_interaction_filtered_data_with_asterisk_02.png"),
       plot = three_way_plot, width = 10, height = 8, dpi = 300)
print(three_way_plot)
