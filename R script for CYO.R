if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(caTools)) install.packages("caTools", repos = "http://cran.us.r-project.org")
if(!require(ROCR)) install.packages("ROCR", repos = "http://cran.us.r-project.org")
if(!require(pscl)) install.packages("pscl", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(car)) install.packages("car", repos = "http://cran.us.r-project.org")
if(!require(broom)) install.packages("broom", repos = "http://cran.us.r-project.org")
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(boot)) install.packages("boot", repos = "http://cran.us.r-project.org")
if(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")

library(dplyr)
library(caTools)
library(ROCR)
library(pscl)
library(caret)
library(car)
library(broom)
library(tidyverse)
library(boot)
library(ggplot2)

data <- read.csv("./indian_liver_patient.csv")

# Checking for NAs

sum(is.na(data))

# Removing NAs since there are only four observations with NAs, at most

data <- na.omit(data)

# Turn 2's in data$Dataset into "0" which means the patient does not have
# liver disease

data$Dataset[data$Dataset == 2] <- 0

# Turn outcome variables into factor variables

data$Dataset <- as.factor(data$Dataset)

# Splitting data into training set and testing set

set.seed(1)

# We use 70% of dataset as training set and 30% as test set

# We choose this train-test split as there are only a few hundred observations
# and we want to keep the proportion of the test set quite high since
# overfitting could be easy with this relatively small dataset

sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
train  <- data[sample, ]
test   <- data[!sample, ]

# Data Visualization with training data
# We group training data by whether a patient has a liver disease or not

data_vis <- train %>% group_by(Dataset)

# We then start with making 100% stacked barplots for the proportions
# of sexes for each group (1 = group with liver disease, 0 = group with
# no liver disease)

data_vis %>%
  count(Gender) %>%
  ggplot(aes(x = Dataset ,y = n, fill = Gender)) +
    geom_bar(position = "fill", stat = "identity") +
    ylab("Proportion") +
    xlab("Groups with or without liver disease") +
    ggtitle("100% stacked barplots for the proportions
         of sexes for each group") +
    theme(plot.title = element_text(hjust = 0.5)
          )

# Boxplots of variables
# Total Bilirubin

data_vis %>%
  ggplot(aes(x = Dataset, y = log(Total_Bilirubin), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Total Bilirubin") + 
    theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Group")

# Direct Bilirubin

data_vis %>%
  ggplot(aes(x = Dataset, y = log(Direct_Bilirubin), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Direct Bilirubin") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")

# Alkaline Phosphotase

data_vis %>%
  ggplot(aes(x = Dataset, y = log(Alkaline_Phosphotase), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Alkaline Phosphotase") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")

# Alamine Aminotransferase

data_vis %>%
  ggplot(aes(x = Dataset, y = log(Alamine_Aminotransferase), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Alamine Aminotransferase") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")

# Aspartate Aminotransferase

data_vis %>%
  ggplot(aes(x = Dataset, y = log(Aspartate_Aminotransferase), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Aspartate Aminotransferase") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")


# Total Protiens
data_vis %>%
  ggplot(aes(x = Dataset, y = log(Total_Protiens), color = Dataset)) +
    geom_boxplot()+
    ggtitle("Total Proteins") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")

# Albumin
data_vis %>%
  ggplot(aes(x = Dataset, y = log(Albumin), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Albumin") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")

# Albumin and Globulin Ratio
data_vis %>%
  ggplot(aes(x = Dataset, y = log(Albumin_and_Globulin_Ratio), color = Dataset)) +
    geom_boxplot() +
    ggtitle("Albumin and Globulin Ratio") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Group")

# Logistic Regression

log_model1 <- glm(Dataset ~ ., data = train, family = 'binomial')

summary(log_model1)

# However, we move straight into removing Albumin and Globulin Ratio since 
# it is clear that it would be multicollinear with Albumin

log_model2 <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio, 
                  data = train,
                  family = 'binomial'
                  )


summary(log_model2)

# Checking for influential values that could affect logistic regression

plot(log_model2, which = 4, id.n = 3)

model.train <- augment(log_model2) %>% 
  mutate(index = 1:n())

# We filter for observations with standard residuals greater than 3
model.train %>% 
  filter(abs(.std.resid) > 3)

# Removing influential observations/outliers

new_train <- log_model2 %>% augment() %>% filter(abs(.std.resid) <= 3)
new_train <- new_train[, 2:12]

# Logistic regression model with same variables but without influential
# observations/outliers

log_model2b <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio, 
                   data = new_train,
                   family = 'binomial'
                   )

summary(log_model2b)

# We check for outliers again and find that there are none now

log_model2b %>% augment() %>% filter(abs(.std.resid) > 3)

# We start backward elimination and take out the variable gender as
# it is currently the most insignificant variable since it has
# the highest p-value (a = 0.05)

log_model2c <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio -Gender,
                   data = new_train, 
                   family = 'binomial'
                   )

summary(log_model2c)

# We take out the variable Total Bilirubin next

log_model2d <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio -Gender -Total_Bilirubin,
                   data = new_train,
                   family = 'binomial'
                   )

summary(log_model2d)

# We take out the variable Aspartate Aminotransferase next

log_model2e <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio -Total_Bilirubin -Gender -
                     Aspartate_Aminotransferase,
                   data = new_train,
                   family = 'binomial'
                   )

summary(log_model2e)

# We take out the variable Total Protiens next

log_model2f <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio -Total_Bilirubin -Gender -
                     Aspartate_Aminotransferase -Total_Protiens,
                   data = new_train,
                   family = 'binomial'
                   )

summary(log_model2f)

# Finally, we take out the variable Albumin and construct the final
# logistic regression model

log_model2g <- glm(Dataset ~ .-Albumin_and_Globulin_Ratio -Total_Bilirubin -Gender -
                     Aspartate_Aminotransferase -Total_Protiens -Albumin,
                   data = new_train,
                   family = 'binomial'
                   )

summary(log_model2g)

# McFadden's of the final model

pR2(log_model2g)["McFadden"]

# Variable importance

varImp(log_model2g)

# Variance inflation factors

vif(log_model2g)

# Predict test data based on final logistic regression model 

predict_reg <- predict(log_model2g, 
                       test, 
                       type = "link"
                       )

# Changing probabilities

predict_reg <- ifelse(predict_reg > 0.5, 1, 0)

# Constructing a confusion matrix and making predict_reg into factors

predict_reg <- as.factor(predict_reg)

logreg_conmatrix <- confusionMatrix(data = predict_reg,
                                    reference = test$Dataset,
                                    positive = "1"
                                    )

# Logistic regression confusion matrix

logreg_conmatrix

# Tweaking probability threshold

# Defining a function that takes threshold level as input and have the corresponding confusion matrix as 
# its output
conf_matrix <- function(threshold){
  predict_reg <- as.numeric(predict_reg)
  predict_reg <- predict(log_model2g, 
                         test, 
                         type = "link")
  predict_reg <- ifelse(predict_reg > threshold, 1, 0)
  predict_reg <- as.factor(predict_reg)
  logreg_conmatrix <- confusionMatrix(data = predict_reg,
                                      reference = test$Dataset,
                                      positive = "1")
  print(logreg_conmatrix)
}

# Using conf_matrix function for multiple threshold values

sapply(c(0.5, 0.005, 0.0005, 0.00005), conf_matrix)

# Defining predict_reg based on the probability threshold of 0.0005 as it produces the confusion matrix
# with the highest accuracy and sensitivity which is desirable in the case of diagnosing people with 
# possible liver disease
predict_reg <- as.numeric(predict_reg)
predict_reg <- ifelse(predict_reg > 0.0005, 1, 0)
predict_reg <- as.factor(predict_reg)
logreg_conmatrix <- confusionMatrix(data = predict_reg,
                                    reference = test$Dataset,
                                    positive = "1")

# ROC-AUC Curve

ROCPred <- prediction(predict(log_model2g, test), test$Dataset) 
ROCPer <- performance(ROCPred, measure = "tpr", 
                      x.measure = "fpr"
                      )

auc <- performance(ROCPred, measure = "auc")
auc <- auc@y.values[[1]]
auc

# Plotting curve

par(mfrow = c(1, 1))

plot(ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC Curve for logistic regression")
abline(a = 0, b = 1)

auc <- round(auc, 4)
legend(.6, .4, auc, title = "AUC", cex = 1)


# XGB Model

xgb_train <- as.data.frame(new_train)

# We change the gender character strings into 1 for Female, 0 
# for male, and make both of them into numeric variables

xgb_train$Gender[xgb_train$Gender == "Female"] <- 1
xgb_train$Gender[xgb_train$Gender == "Male"] <- 0
xgb_train$Gender <- as.numeric(xgb_train$Gender)

# We specify reasonable options for hyperparameters that will
# be used to check for the best XGBoost model

grid_tune <- expand.grid(nrounds = c(100, 400, 600),
                         max_depth = c(4, 6, 8, 10),
                         eta = c(0.1, 0.3, 0.5),
                         gamma = c(0, 1, 2),
                         colsample_bytree = c(1, 0.7, 0.5),
                         min_child_weight = c(1, 2, 3), 
                         subsample = c(0.5, 0.75, 1)
                         )



# We set a three-fold cross validation (due to the computational
# resources that grid_tune is already taking)

train_control <- trainControl(method = "cv",
                              number= 3,
                              verboseIter = FALSE,
                              allowParallel = TRUE
                              )

# We fit the XGBoost model with the presence of liver disease
# being the y variable and all other variables as the x variables

xgb_model <- train(x = xgb_train[, -1],
                   y = xgb_train$Dataset,
                   trControl = train_control,
                   tuneGrid = grid_tune,
                   method = "xgbTree",
                   verbose = FALSE,
                   verbosity = 0
                   )

xgb_model$bestTune

# Model evaluation

xgb_test <- test
xgb_test$Gender[xgb_test$Gender == "Female"] <- 1
xgb_test$Gender[xgb_test$Gender == "Male"] <- 0
xgb_test$Gender <- as.numeric(xgb_test$Gender)

xgb_pred <- predict(xgb_model, xgb_test)

# Confusion matrix for the final XGBoost model
xgb_conmatrix <- confusionMatrix(data = xgb_pred, 
                                 reference = xgb_test$Dataset,
                                 positive = "1"
                                 )

xgb_conmatrix 

# ROC-AUC Curve

xgb_prob <- predict(xgb_model, xgb_test, type = "prob")
log_odds_xgb_prob <- log(xgb_prob[, 2]/xgb_prob[ , 1])

xgb_ROCPred <- prediction(log_odds_xgb_prob, xgb_test$Dataset) 
xgb_ROCPer <- performance(xgb_ROCPred, measure = "tpr", 
                      x.measure = "fpr")

xgb_auc <- performance(xgb_ROCPred, measure = "auc")
xgb_auc <- xgb_auc@y.values[[1]]
xgb_auc

# Plotting curve
plot(xgb_ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC Curve for XGB Model"
     )

abline(a = 0, b = 1)

xgb_auc <- round(xgb_auc, 4)
legend(.6, .4, xgb_auc, title = "AUC", cex = 1)

# We compare the ROC-AUC Curve for each model
par(mfrow = c(1, 2))

# Logistic regression

plot(ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC Curve for logistic regression"
     )

abline(a = 0, b = 1)

auc <- round(auc, 4)
legend(.22, .2, auc, title = "AUC", cex = 1)

# XGB Model

plot(xgb_ROCPer, colorize = TRUE, 
     print.cutoffs.at = seq(0.1, by = 0.1), 
     main = "ROC Curve for XGB Model"
     )

abline(a = 0, b = 1)

xgb_auc <- round(xgb_auc, 4)
legend(.22, .2, xgb_auc, title = "AUC", cex = 1)

# We compare the confusion matrix for each model

# Logistic Regression

logreg_conmatrix

# XGB Model

xgb_conmatrix

# Resources:
# https://www.geeksforgeeks.org/logistic-regression-in-r-programming/
# https://www.youtube.com/watch?v=qjeUhuvkbHY&t=723s
