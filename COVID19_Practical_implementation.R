#Title: A Comprehensive Guide to Selecting the Right Modeling Strategy for Explanatory and Predictive Data Analysis
#Analysis: COVID-19 Immunologic Profile Analysis
# Comprehensive script including EDA, preprocessing, and modeling

# Load required packages
library(tidyverse)    # Data manipulation and visualization
library(caret)       # Machine learning framework
library(glmnet)      # Regularized regression
library(randomForest) # Random forest implementation
library(GGally)      # Advanced plotting
library(ggcorrplot)  # Correlation visualization
library(pROC)        # ROC curve analysis
library(factoextra)  # PCA visualization
library(mice)        # Missing data imputation
library(VIM)         # Missing data visualization
library(gridExtra)   # For arranging multiple plots
library(tidyr)      #Data manipulation
library(glmmLasso)  # For mixed-effects lasso
library(lme4)       # For mixed model formulas
library(kableExtra) # For the analysis report

# Custom white theme function
theme_white <- function() {
  theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid = element_line(color = "gray90"),
      panel.border = element_blank()
    )
}
#Step 1: Defining the Study Goal: 
# Developing a prediction model to identify biomarkers for Covid19-disease progress.
#------------------------------------------------------------------------------------
# Step 2: Understanding the Data 
#Data Loading and Initial Exploration 
data <- readxl::read_excel("Immunologic profiles in patients with COVID-19.xlsx")
# 1. Data  Preprocessing --------------------------------------------------- 
#1) Basic Data Overview
# View the structure of the dataset
str(data) #some cytokines include "OOR", showing as character:
# 2) Data cleaning:
# Dealing with "OOR":
data[data == "OOR <"] <- NA        # replace 'OOR <' with NA
cytokine_cols <- 6:30  # Columns from 6 to 30 (cytokines)   # Convert the columns that should be numeric (cytokines) from character to numeric
data[, cytokine_cols] <- sapply(data[, cytokine_cols], as.numeric)
min_detectable <- min(data[, cytokine_cols], na.rm = TRUE)  # Find the minimum detectable value (smallest value excluding NA)
imputation_value <- min_detectable - 0.01  # Subtract a small amount, e.g., 0.01
data[, cytokine_cols] <- lapply(data[, cytokine_cols], function(x) {
  ifelse(is.na(x), imputation_value, x)
})  # Impute 'OOR <' values (now NA) with a small amount less than the minimum detectable value

str(data) # The data now has the 'OOR <' values replaced with the imputation value
# Convert varaiables to factor
data$Patient  <- factor(data$Patient )
data$SEX <- factor(data$SEX)
data$week <- factor(data$week)
data<- na.omit(data) 

# 2. Exploratory Data Analysis (EDA) --------------------------------------
# Summary statistics
summary(data)
# Check for missing values
colSums(is.na(data))
# Clean column names
colnames(data) <- make.names(colnames(data))
# Convert Death to factor (our outcome variable)
data$Death <- as.factor(data$Death)
# 1) Missing data analysis
missing_plot <- aggr(data, numbers = TRUE, sortVars = TRUE)
print(missing_plot) # no missing values

# 2) Distribution of key variables
# Bar plots for categorical variables (e.g., week, sex)
sex_plot <- ggplot(data, aes(x = factor(SEX))) + 
  geom_bar(fill = "lightgreen", alpha = 0.7) + 
  labs(title = "Distribution of Sex", x = "SEX", y = "Count") +
  theme_minimal()
ggsave("Distribution_Sex.png", sex_plot, width = 8, height = 6, dpi = 300, bg = "white")

# Bar plots for categorical variables (e.g., week, sex)
death_plot <- ggplot(data, aes(x = factor(Death))) + 
  geom_bar(fill = "lightblue", alpha = 0.7) + 
  labs(title = "Distribution of Death", x = "Death", y = "Count") +
  theme_minimal()

ggsave("Distribution_death.png", death_plot, width = 8, height = 6, dpi = 300, bg = "white")

# Plot distributions of all immunologic markers using Colorblind-friendly palette
dist_plots <-data %>% 
  select(-Patient, -AGE, -SEX, -Death, -week) %>% 
  gather(key = "Marker", value = "Value") %>% 
  mutate(Value = as.numeric(Value)) %>%  # make sure all cytokines are numeric
  ggplot(aes(x = Value, fill = Marker)) +
  geom_histogram(color = "white", alpha = 0.8, bins = 30) +
  facet_wrap(~Marker, scales = "free") +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Immunologic Markers",
       x = "Measurement Value",
       y = "Frequency") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(dist_plots)

ggsave("Distribution_all_cytokines.png", dist_plots, width = 8, height = 6, dpi = 300, bg = "white")

# 3) Normalize/scale numeric features (except outcome and IDs)
preprocess_params <- preProcess(data %>% select(-Death, -Patient, -SEX), 
                                method = c("center", "scale"))
processed_data <- predict(preprocess_params, data)

library(DataExplorer)
create_report(data, output_file = "report.pdf", output_format = "pdf_document")
create_report(data, output_file = "data_report.html", output_dir = getwd(), report_title = "Covid Data Report")

#------------------------------------------------------------------------------------
#Step 3: Variable Selection 
# Feature Filtering ------------------------------------------------------------------------------------
# 1) Correlation analysis
#Calculate correlation matrix (excluding non-numeric columns)
cor_matrix <- cor(data %>% select(-Patient, -Death, -week,-Death,-SEX), 
                  use = "pairwise.complete.obs")

cor_plot <- ggcorrplot(cor_matrix, hc.order = TRUE, 
                       outline.color = "white",
                       ggtheme = ggplot2::theme_minimal,
                       colors = c("#6D9EC1", "white", "red"),
                       title = "Correlation Matrix of Immunologic Markers")

print(cor_plot) # 6 correlated factors

#Remove highly correlated features (threshold > 0.85)

# Find and remove correlated variables
highly_correlated <- findCorrelation(cor_matrix, 
                                     cutoff = 0.85,  # Adjust threshold as needed
                                     names = TRUE,  # Return column names
                                     exact = TRUE)  # More precise calculation

# Verify removed variables
cat("Removing", length(highly_correlated), "highly correlated variables:\n")
print(highly_correlated) 
# as "TNF.alpha" is an important antinflammatory that we dont want to lose it.
#control which variable to remove:
# Get full correlation pairs above threshold
high_cor_pairs <- which(abs(cor_matrix) > 0.85 & upper.tri(cor_matrix), arr.ind = TRUE)

# Create a tibble of correlated pairs
cor_tibble <- tibble(
  var1 = rownames(cor_matrix)[high_cor_pairs[,1]],
  var2 = colnames(cor_matrix)[high_cor_pairs[,2]],
  cor_value = cor_matrix[high_cor_pairs]
) %>% arrange(desc(abs(cor_value)))

# Manually review and select variables to remove
print(cor_tibble)

# I will then manually specify which to remove
vars_to_remove <- c("BTLA","CCL4.MIP.1.beta","IL.10")  # Example - use your actual variables
filtered_data <- data %>% select(-all_of(vars_to_remove))

# Verify new dimensions
cat("\nOriginal variables:", ncol(data), 
    "\nFiltered variables:", ncol(filtered_data), 
    "\nRemoved:", ncol(data) - ncol(filtered_data), "variables\n")

# Visualize remaining correlations
filtered_cor <- cor(filtered_data %>% select_if(is.numeric), 
                    use = "pairwise.complete.obs")
ggcorrplot(filtered_cor, 
           hc.order = TRUE, 
           type = "lower",
           colors = c("#6D9EC1", "white", "red"),
           title = "Filtered Correlation Matrix")


# 2) Feature Extraction --------------------------------------------------
# Principal Component Analysis (PCA) 
# Normalize/scale numeric features (except outcome and IDs)
preprocess_params <- preProcess(filtered_data %>% select(-Death, -Patient, -SEX), 
                                method = c("center", "scale"))
processed_data <- predict(preprocess_params, filtered_data)

# Perform PCA on immunologic markers
# Convert Death to factor with descriptive labels
processed_data$Death <- factor(processed_data$Death, 
                               levels = c(0, 1), 
                               labels = c("Survived", "Died"))

# Perform PCA on immunologic markers (excluding metadata)
pca_data <- processed_data %>% 
  select(-Death, -Patient, -SEX, -week, -AGE) %>% 
  select_if(is.numeric)  # Ensure only numeric columns are included

pca_result <- prcomp(pca_data, scale. = TRUE)

# Enhanced PCA visualization
pca_plot <- fviz_pca_ind(pca_result,
                         col.ind = processed_data$Death,  # Use the factor with new labels
                         palette = c("darkgreen", "#FC4E07"),
                         addEllipses = TRUE,
                         #ellipse.type = "confidence",
                         #ellipse.level = 0.95,  # 95% confidence ellipses
                         legend.title = "Mortality Status",  # More descriptive
                         title = "PCA of Immunologic Markers by Mortality Status",
                         ggtheme = theme_white()) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank())

# Display plots
print(pca_plot)
# Scree plot with percentage variance explained
scree_plot <- fviz_eig(pca_result, 
                       addlabels = TRUE, 
                       barfill = "#4E84C4",
                       barcolor = "#4E84C4",
                       linecolor = "#293352",
                       main = "Scree Plot: Variance Explained by PCA Components",
                       ggtheme = theme_white()) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

# Display plots
print(scree_plot)

# Save high-quality versions
#ggsave("pca_by_mortality.png", pca_plot, width = 8, height = 6, dpi = 300, bg = "white")
#ggsave("pca_scree_plot.png", scree_plot, width = 8, height = 6, dpi = 300, bg = "white")

#---------------------------------------------------------------------------------------------------------
# Step 4: Model Assumption 
# perform 3 models: Regularized Logistic Regression, RF, and GLMMLasso ------------------------------------------------------------
# 1) Data preperation:
#Data splitting (75% training, 25% testing)
# Convert Death to factor with explicit level names
processed_data <- predict(preprocess_params, filtered_data)
processed_data$Death <- factor(processed_data$Death, levels = c(0, 1), labels = c("Survived", "Died"))
set.seed(129)
trainIndex <- createDataPartition(processed_data$Death, p = .75, 
                                  list = FALSE, 
                                  times = 1)

train_data <- processed_data[trainIndex, ]
test_data <- processed_data[-trainIndex, ]

# 2) Regularized Logistic Regression (Elastic Net)
# Set up cross-validation
ctrl <- trainControl(method = "cv",
                     number = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")

# Train model
logit_model <- train(Death ~ . - Patient,
                     data = train_data,
                     method = "glmnet",
                     family = "binomial",
                     trControl = ctrl,
                     tuneLength = 10,
                     metric = "ROC")

# Print model results
print(logit_model)
plot(logit_model, main = "Regularized Logistic Regression Tuning")
summary(logit_model)
logit_model$bestTune
logit_model$results


# Variable Importance
# Logistic regression coefficients
logit_coef <- varImp(logit_model, scale = FALSE)
plot(logit_coef, top = 15, main = "Top 15 Important Variables (Logistic Regression)")
# CD40, CTLA-4,TIM-3, CCL2-JE/MCP-1, AGE


# 3) Random Forest
rf_model <- train(Death ~ . - Patient,
                  data = train_data,
                  method = "rf",
                  trControl = ctrl,
                  tuneLength = 5,
                  importance = TRUE,
                  metric = "ROC")

print(rf_model)
plot(rf_model, main = "Random Forest Tuning")


# Random Forest importance
rf_importance <- varImp(rf_model, scale = FALSE)
plot(rf_importance, top = 15, main = "Top 15 Important Variables (Random Forest)")
#CD40, CTLA-4, AGE, IL-6, TIM-3, CCL2-JE/MCP-1, LAG-3, CD-27, GITR, TNF-alpha

# 4) GLMMLasso 
# here, we will use the original data as glmmlasso fuction perform scale and feature selection via regularization
data$Death <- as.numeric(as.character(data$Death)) # Convert back to 0/1
# prediction model:
set.seed(129)  # Set seed for reproducibility
# Split the data into training (80%) and test (20%) sets
trainIndex <- createDataPartition(data$Death, p = .8, 
                                  list = FALSE, 
                                  times = 1)
trainData <- data[trainIndex, ]
testData  <- data[-trainIndex, ]
trainData <- as.data.frame(trainData)
#  1) choose the best lambda:
# Manual k-fold cross-validation was performed using the glmmLasso function to find the optimal lambda value, taking into account fixed and random effects.
# since glmmLasso package lacks built-in cross-validation capabilities
# conventional cross-validation methods in glmnet caret packages unsuitable (doesnt account for random effects).
# Set up the k-fold cross-validation
k <- 10  # Number of folds
set.seed(123)  # For reproducibility

# Define lambda values to test (you can adjust this range)
lambda_values <- seq(0.001, 1, by = 0.05)

# Create k-folds
folds <- createFolds(trainData$Death, k = k, list = TRUE)

# Initialize variable to store results
cv_errors <- numeric(length(lambda_values))

# Perform k-fold cross-validation
for (i in 1:length(lambda_values)) {
  lambda_value <- lambda_values[i]
  fold_errors <- numeric(k)
  
  for (j in 1:k) {
    # Create training and validation sets for the current fold
    train_indices <- unlist(folds[-j])
    valid_indices <- unlist(folds[j])
    
    train_fold <- trainData[train_indices, ]
    valid_fold <- trainData[valid_indices, ]
    train_fold <- na.omit(train_fold)
    # Fit the glmmLasso model on the training fold
    model <- glmmLasso(
      fix = Death ~ AGE + SEX + week + BTLA + CD27 + CD28 + `TIM.3` + HVEM + CD40 +
        GITR + `LAG.3` + `TLR.2` + GITRL + `PD.1` + `CTLA.4` + `CD80.B7.1` +
        `CD86.B7.2` + `PD.L1` + `PD.L2` + ICOS + `CCL2.JE.MCP.1` +
        CCL4.MIP.1.beta +`CXCL10.IP.10.CRG.2` + `GM.CSF` + `IL.6` +
        `IL.8.CXCL8` + `IL.10` + `TNF.alpha`,
      rnd = list(Patient = ~1),  # Random effect for 'Patient'
      data = train_fold,
      lambda = lambda_value,
      family = binomial()
    )
    
    # Predict on the validation fold
    predictions <- predict(model, newdata = valid_fold, type = "response")
    predicted_classes <- ifelse(predictions > 0.5, 1, 0)
    
    # Calculate the error (misclassification rate or another metric)
    fold_errors[j] <- mean(predicted_classes != valid_fold$Death)
  }
  
  # Average error across all folds for the current lambda
  cv_errors[i] <- mean(fold_errors)
}

# Find the lambda value that minimizes cross-validation error
best_lambda <- lambda_values[which.min(cv_errors)]
cat("Best lambda value:", best_lambda)
# Also using lambda value within 1 SD of the mean
# Calculate the mean and standard deviation of cv_errors
mean_cv_error <- mean(cv_errors)
sd_cv_error <- sd(cv_errors)

# Create a data frame
cv_results <- data.frame(lambda = lambda_values, AUC = auc_scores)
# Simulated example

best_lambda <- lambda_values[which.min(cv_errors)]
lambda_1se <- max(lambda_values[cv_errors <= min(cv_errors) + sd_cv_error])

library(ggplot2)

# Plotting
ggplot(data = data.frame(lambda = lambda_values, error = cv_errors),
       aes(x = lambda, y = error)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = best_lambda, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = lambda_1se, linetype = "dotted", color = "darkred") +
  labs(title = "GLMMLasso Tuning Curve",
       subtitle = "Cross-validated error vs. Lambda",
       x = expression(lambda),
       y = "Cross-Validation Error") +
  theme_minimal()



# Define the threshold for acceptable lambda values (within 1 SD of the mean)
acceptable_lambda_indices <- which(cv_errors <= (mean_cv_error + sd_cv_error))

# Find the best lambda value within the acceptable range
best_lambda_within_sd <- lambda_values[acceptable_lambda_indices][which.min(cv_errors[acceptable_lambda_indices])]

cat("Best lambda value:", best_lambda_within_sd)

# Step 2: Fit the GLMM with Lasso (glmmLasso) on the training data
lambda_value <- best_lambda  #0.101 # # using the estimated lambda value 
glmmLasso_model <- glmmLasso(
  fix = Death ~ AGE + SEX + week + BTLA + CD27 + CD28 + `TIM.3` + HVEM + CD40 + 
    GITR + `LAG.3` + `TLR.2` + GITRL + `PD.1` + `CTLA.4` + `CD80.B7.1` + 
    `CD86.B7.2` + `PD.L1` + `PD.L2` + ICOS + `CCL2.JE.MCP.1` + 
    `CCL4.MIP.1.beta` + `CXCL10.IP.10.CRG.2` + `GM.CSF` + `IL.6` + 
    `IL.8.CXCL8` + `IL.10` + `TNF.alpha`,
  rnd = list(Patient = ~1),  # Random effect for 'Patient'
  data = trainData,
  lambda = lambda_value,  # Regularization parameter
  family = binomial()  # Binary outcome
)

summary(glmmLasso_model)
# identifying the top biomarkers:
non_zero_coefs <- glmmLasso_model$coefficients[glmmLasso_model$coefficients != 0]

# Display the non-zero coefficients and their corresponding predictors
print(non_zero_coefs)
# Rank the predictors by absolute value of the coefficients (importance)
ranked_predictors <- sort(abs(non_zero_coefs), decreasing = TRUE)

# Display the ranked predictors with their coefficients
print(ranked_predictors)
# Step 3: Make predictions on the test set
testData <- as.data.frame(testData)
testPredictions <- predict(glmmLasso_model, newdata = testData, type = "response")

# Convert probabilities to binary outcome (threshold 0.5 for classification)
testPredBinary <- ifelse(testPredictions > 0.5, 1, 0)

# Step 5: Evaluate the model 1) Confusion Matrix
confMatrix <- confusionMatrix(as.factor(testPredBinary), as.factor(testData$Death))
print(confMatrix)

# 2) Calculate AUC for better evaluation
roc_curve <- roc(testData$Death, testPredictions)
auc(roc_curve)

# 3) Create ROC curve object
roc_curve <- roc(testData$Death, testPredictions)

# Plot smoothed ROC curve
plot(smooth(roc_curve),  # This will smooth the curve
     col = "darkred",        # Color of the curve
     main = "ROC Curve for GLMM with Lasso",
     xlab = "1 - Specificity",
     ylab = "Sensitivity",
     lwd = 2)             # Line width

# Add AUC to the plot
auc_value <- auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), col = "darkred", lwd = 2)

# Print AUC value for reference
print(auc_value)

# GLMMLasso Variable Importance 
# Extract non-zero coefficients
glmm_model <- glmmLasso_model
glmm_coef <- coef(glmm_model)
glmm_important <- glmm_coef[glmm_coef != 0 & !names(glmm_coef) %in% c("(Intercept)", "Patient.ID")]

# Create importance plot
glmm_imp_plot <- data.frame(
  Variable = names(glmm_important),
  Importance = abs(glmm_important)
) %>%
  arrange(desc(Importance)) %>%
  ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "#1E88E5") +
  coord_flip() +
  labs(title = "GLMM Lasso Variable Importance",
       x = "Variable",
       y = "Coefficient Magnitude") +
  theme_minimal()

print(glmm_imp_plot)

#---------------------------------------------------------------------------------------------------------
# Step 5: Model Evaluation ----------------------------------------------------
# 1) Logistic Regression Evaluation
logit_pred <- predict(logit_model, newdata = test_data, type = "prob")
logit_roc <- roc(test_data$Death, logit_pred[,2]) #0.714

# Predict on test set
pred_logit <- predict(logit_model, newdata = test_data)

# Create confusion matrix
confusionMatrix(pred_logit, test_data$Death, positive = "Died")


# 2) Random Forest Evaluation
rf_pred <- predict(rf_model, newdata = test_data, type = "prob")
rf_roc <- roc(test_data$Death, rf_pred[,2]) # 0.818

# Confusion Matrix:
# Predict on test set
pred_rf <- predict(rf_model, newdata = test_data)

# Create confusion matrix
confusionMatrix(pred_rf, test_data$Death, positive = "Died")



# 3)GLMMLasso Evaluation
roc_curve <- roc(testData$Death, testPredictions)

#  ROC comparison plot
roc_comparison <- ggroc(list(Logistic = logit_roc, 
                             RandomForest = rf_roc,
                             GLMMLasso = roc_curve)) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
  labs(title = "ROC Curve Comparison for GLMMLasso, Logistic, RF",
       color = "Model Type") +
  theme_minimal()

print(roc_comparison)

ggsave("roc_3models.png", roc_comparison, width = 8, height = 6, dpi = 300, bg = "white")


# Compute AUC for each model
auc_logit <- auc(logit_roc)
auc_rf <- auc(rf_roc)
auc_glmm <- auc(roc_curve)

# Create named ROC list
roc_list <- list(
  Logistic = logit_roc,
  RandomForest = rf_roc,
  GLMMLasso = roc_curve
)

# Generate ROC plot with ggroc
roc_comparison <- ggroc(roc_list, legacy.axes = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve Comparison: Regularized Logistic, Random Forest, and GLMMLasso",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Model"
  ) +
  annotate("text", x = 0.65, y = 0.25,
           label = paste0("AUC (Logistic) = ", round(auc_logit, 3), "\n",
                          "AUC (RF) = ", round(auc_rf, 3), "\n",
                          "AUC (GLMMLasso) = ", round(auc_glmm, 3)),
           hjust = 0, size = 4) +
  theme_minimal(base_size = 14)

# Print plot
print(roc_comparison)

ggsave("auc_3models.png", roc_comparison, width = 8, height = 6, dpi = 300, bg = "white")

# 7. Generate Final Report -----------------------------------------------
# Create a list of all results
results <- list(
  missing_data_plot = missing_plot,
  distribution_plots = dist_plots,
  correlation_plot = cor_plot,
  pca_plot = pca_plot,
  scree_plot = scree_plot,
  logistic_model = logit_model,
  random_forest_model = rf_model,
  glmmLasso_model = glmmLasso_model,
  roc_comparison = roc_comparison,
  highly_correlated_vars = names(highly_correlated),
  preprocess_params = preprocess_params
)

# Save results for reporting
saveRDS(results, "covid_immunologic_analysis_results.rds")

# 8. Create Publication-Quality Plots ------------------------------------
# Actual vs Predicted plot
pred_actual_data <- data.frame(
  Actual = test_data$Death,
  Logistic = logit_pred[,2],
  RandomForest = rf_pred[,2],
  GLMMLasso  = glmm_pred
)
# Create Combined Actual vs Predicted Plot ---------------------------
combined_pred_plot <- pred_actual_data %>%
  gather(key = "Model", value = "Prediction", -Actual) %>%
  ggplot(aes(x = as.numeric(Actual)-1, y = Prediction, color = Model)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("#E46726", "#6D9EC1", "#1E88E5")) +
  labs(title = "Actual vs Predicted Mortality by GLMMLasoo, Logistic, RandomForest",
       x = "Actual Mortality (0 = Survived, 1 = Died)",
       y = "Predicted Probability of Death") +
  theme_minimal()

print(combined_pred_plot)


# Save all plots to files
ggsave("correlation_matrix.png", cor_plot, width = 12, height = 10, bg="white")
ggsave("pca_plot.png", pca_plot, width = 10, height = 8, bg="white")
ggsave("roc_comparison.png", roc_comparison, width = 10, height = 8, bg="white")
ggsave("actual_vs_predicted.png", combined_pred_plot, width = 12, height = 8, bg="white")

models_summary %>%
  kable("html", caption = "Model Comparison") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  add_header_above(c(" " = 1, "Logistic" = 1, "RF" = 1, "GLMMLasso" = 1)) %>%
  footnote(general = "Test set performance evaluated on 20% holdout data.")

