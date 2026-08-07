# SKRYPT 04: Modeling 

library(caret)
library(randomForest)
library(kernlab)
library(pheatmap)
library(dplyr)

set.seed(777)

load("results/models/Boruta_Results.RData")

train_index <- createDataPartition(ml_data_final$Target, p = 0.75, list = FALSE)
train_set <- ml_data_final[train_index, ]
test_set  <- ml_data_final[-train_index, ]

# trening
fit_control <- trainControl(
  method = "repeatedcv",
  number = 5,           # 5-fold 
  repeats = 3,          
  classProbs = TRUE,
  savePredictions = "final"
)

message("\n> Trenowanie modeli.")

# Random Forest 
message("Trenowanie Random Forest.")
set.seed(123)
model_rf <- train(
  Target ~ ., data = train_set, 
  method = "rf", 
  trControl = fit_control, 
  metric = "Accuracy",
  ntree = 500,
  tuneLength = 5
)
message("   -> Random Forest: done!")

# k-NN 
message("Trenowanie k-NN.")
set.seed(123)
model_knn <- train(
  Target ~ ., data = train_set, 
  method = "knn", 
  trControl = fit_control, 
  metric = "Accuracy",
  preProcess = c("center", "scale"),
  tuneLength = 10
)
message("   -> k-NN: done!")

# SVM
message("Trenowanie SVM (Linear).")
set.seed(123)
model_svm_linear <- tryCatch({
  train(
    Target ~ ., data = train_set, 
    method = "svmLinear", 
    trControl = fit_control, 
    metric = "Accuracy",
    preProcess = c("center", "scale"),
    tuneLength = 5
  )
}, error = function(e) { 
  return(NULL) 
})

if(!is.null(model_svm_linear)) message("   -> SVM: done!")

message("koniec treningu")

# zestawienie wynikow
models_list <- list(Random_Forest = model_rf, kNN = model_knn)
if(!is.null(model_svm_linear)) models_list$SVM <- model_svm_linear

results <- resamples(models_list)

message("\n wyniki treningowe (Cross-Validation) ---")
print(summary(results))

message("\n weryfikacja na zbiorze testowym")

evaluate_model <- function(model, name) {
  if(is.null(model)) return(NULL)
  pred <- predict(model, test_set)
  cm <- confusionMatrix(pred, test_set$Target)
  acc <- round(cm$overall['Accuracy'] * 100, 2)
  kappa <- round(cm$overall['Kappa'], 2)
  message(sprintf("%-15s : Accuracy = %5.2f%%, Kappa = %4.2f", name, acc, kappa))
  return(cm)
}

cm_rf <- evaluate_model(model_rf, "Random Forest")
cm_knn <- evaluate_model(model_knn, "k-NN")
if(!is.null(model_svm_linear)) cm_svm <- evaluate_model(model_svm_linear, "SVM")

# zapis wynikow
if(!dir.exists("results/models")) dir.create("results/models", recursive = TRUE)
save(models_list, results, cm_rf, file = "results/models/Final_Models.RData")

message("\n>>> Analiza zakończona. <<<")
