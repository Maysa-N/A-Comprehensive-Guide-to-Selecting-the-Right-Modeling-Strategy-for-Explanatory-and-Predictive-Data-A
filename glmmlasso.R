# Oct 27, Maysa Niazy
# Practical implementation for Guidelines paper
#using COVID-19 data found in:https://doi.org/10.7910/DVN/TAS2PH 

library(glmmLasso)
library(caret)
library(tidyverse)  # For data manipulation and visualization
library(dplyr)
library(pROC)
library(car)
library(corrplot)   # For correlation matrix plot
library(ggplot2)    # For data visualization
library(DataExplorer) # To generate automated EDA reports
setwd("/Users/mniazy/Downloads")
# Load the data:
data <- readxl::read_excel("Immunologic profiles in patients with COVID-19.xlsx")
data <- as.data.frame(data)

#
#Step 1: Goal: Prediction model
# Step 2: Data preperation 
# 1) Basic Data Overview
# View the structure of the dataset
str(data)
#Convert varaiables to factor
data$Patient  <- factor(data$Patient )
data$SEX <- factor(data$SEX)
data$week <- factor(data$week)
data<- na.omit(data)
# 2- data cleaning:
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
# 2) Exploratory Data Analysis
# Summary statistics
summary(data)
# Check for missing values
colSums(is.na(data))
# 2.1 Univariate Analysis
# Histograms for continuous variables (e.g., age, cytokine levels)
data %>% 
  select_if(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) + 
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) + 
  facet_wrap(~key, scales = "free_x") + 
  theme_minimal() +
  labs(title = "Histograms of Continuous Variables", x = "Value", y = "Frequency")

# Bar plots for categorical variables (e.g., week, sex)
ggplot(data, aes(x = factor(SEX))) + 
  geom_bar(fill = "lightgreen", alpha = 0.7) + 
  labs(title = "Distribution of Sex", x = "SEX", y = "Count") +
  theme_minimal()

# 3) Bivariate Analysis
# Boxplots to compare continuous variables (e.g., cytokines) across the binary outcome 'death'
ggplot(data, aes(x = factor(Death), y = `CCL2/JE/MCP-1`)) + 
  geom_boxplot(fill = "orange", alpha = 0.7) + 
  labs(title = "MCP-1 Levels by Death Outcome", x = "Death (0=No, 1=Yes)", y = "MCP-1")
# Boxplots to compare continuous variables (e.g., cytokines) across the binary outcome 'death'
ggplot(data, aes(x = factor(Death), y = BTLA)) + 
  geom_boxplot(fill = "orange", alpha = 0.7) + 
  labs(title = "BTLA Levels by Death Outcome", x = "Death (0=No, 1=Yes)", y = "BTLA")

# 2.3 Correlation Matrix for Numeric Variables
# Compute the correlation matrix
numeric_vars <- data %>% select_if(is.numeric)
corr_matrix <- cor(numeric_vars, use = "complete.obs")

# Visualize the correlation matrix
corrplot(corr_matrix, method = "color", type = "upper", tl.cex = 0.7, number.cex = 0.7, title = "Correlation Matrix")

# 3) Check Multicollinearity
# Variance Inflation Factor (VIF)
vif_model <- lm(BTLA ~ ., data = data)  # Example linear model for VIF calculation
vif(vif_model)
# Standardization for columns 6 to 30
data <- data %>%
  mutate(across(6:30, ~ ( . - mean(.)) / sd(.)))

# Final check before model fitting
summary(covid_data_scaled)
# Automated EDA Report
#create_report(data)
# Step 3: Variable Selection:
# Filtering selection based on correlations
# Check Multicollinearity
# Variance Inflation Factor (VIF)
vif_model <- lm(BTLA ~ ., data = data)  # Example linear model for VIF calculation
vif(vif_model)

# Step 4: Model assumptions
###################################
# prediction model:
set.seed(129)  # Set seed for reproducibility

# Split the data into training (80%) and test (20%) sets
trainIndex <- createDataPartition(data$Death, p = .8, 
                                  list = FALSE, 
                                  times = 1)
trainData <- data[trainIndex, ]
testData  <- data[-trainIndex, ]

#  1) choose the best lambda:
# Manual k-fold cross-validation using the glmmLasso function is used to find the optimal lambda value, taking into account fixed and random effects.
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
    
    # Fit the glmmLasso model on the training fold
    model <- glmmLasso(
      fix = Death ~ AGE + SEX + week + BTLA + CD27 + CD28 + `TIM-3` + HVEM + CD40 +
        GITR + `LAG-3` + `TLR-2` + GITRL + `PD-1` + `CTLA-4` + `CD80/B7-1` +
        `CD86/B7-2` + `PD-L1` + `PD-L2` + ICOS + `CCL2/JE/MCP-1` +
        `CCL4/MIP-1 beta` + `CXCL10/IP-10/CRG-2` + `GM-CSF` + `IL-6` +
        `IL-8/CXCL8` + `IL-10` + `TNF-alpha`,
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
# or using lambda value within 1 SD of the mean
# Calculate the mean and standard deviation of cv_errors
mean_cv_error <- mean(cv_errors)
sd_cv_error <- sd(cv_errors)

# Define the threshold for acceptable lambda values (within 1 SD of the mean)
acceptable_lambda_indices <- which(cv_errors <= (mean_cv_error + sd_cv_error))

# Find the best lambda value within the acceptable range
best_lambda_within_sd <- lambda_values[acceptable_lambda_indices][which.min(cv_errors[acceptable_lambda_indices])]

cat("Best lambda value:", best_lambda_within_sd)


# Step 2: Fit the GLMM with Lasso (glmmLasso) on the training data
lambda_value <- best_lambda  #0.101 # # using the estimated lambda value 
glmmLasso_model <- glmmLasso(
  fix = Death ~ AGE + SEX + week + BTLA + CD27 + CD28 + `TIM-3` + HVEM + CD40 + 
    GITR + `LAG-3` + `TLR-2` + GITRL + `PD-1` + `CTLA-4` + `CD80/B7-1` + 
    `CD86/B7-2` + `PD-L1` + `PD-L2` + ICOS + `CCL2/JE/MCP-1` + 
    `CCL4/MIP-1 beta` + `CXCL10/IP-10/CRG-2` + `GM-CSF` + `IL-6` + 
    `IL-8/CXCL8` + `IL-10` + `TNF-alpha`,
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

# Plot the smoothed ROC curve
plot(smooth(roc_curve),  # This will smooth the curve
     col = "blue",        # Color of the curve
     main = "ROC Curve for GLMM with Lasso",
     xlab = "1 - Specificity",
     ylab = "Sensitivity",
     lwd = 2)             # Line width

# Add AUC to the plot
auc_value <- auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), col = "blue", lwd = 2)

# Print AUC value for reference
print(auc_value)


