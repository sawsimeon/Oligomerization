library(protr)
OliFP_100 <- protr::readFASTA("OliFP_100.Fasta")
OliFP_99 <- protr::readFASTA("OliFP_99.Fasta")
OliFP_95 <- protr::readFASTA("OliFP_95.Fasta")
OliFP_90 <- protr::readFASTA("OliFP_90.Fasta")

OliFP_100 <- OliFP_100[(sapply(OliFP_100, protcheck))]
OliFP_99 <- OliFP_99[(sapply(OliFP_99, protcheck))]
OliFP_95 <- OliFP_95[(sapply(OliFP_95, protcheck))]
OliFP_90 <- OliFP_90[(sapply(OliFP_90, protcheck))]

acc_models <- function(x, model){
  library(protr)
  library(caret)
  library(RWeka)
  library(stringr)
  data <- protr::readFASTA(x)
  data <- data[(sapply(data, protr::protcheck))]
  acc  = t(sapply(data, protr::extractAAC))
  Oligomerization <- row.names(acc)
  Dimer <- stringr::str_extract(Oligomerization, "Dimer")
  Monomer <- stringr::str_extract(Oligomerization, "Monomer")
  Tetramer <- stringr::str_extract(Oligomerization, "Tetramer")
  Oligomerization <- paste(Dimer, Monomer, Tetramer)
  Oligomerization <- sapply(Oligomerization, gsub, pattern = "NA", replacement = "")
  Oligomerization <- stringr::str_trim(Oligomerization, side = "both")
  Oligomerization <- ifelse(Oligomerization == "Monomer", "Monomer", "Oligomer")
  Oligomerization <- as.factor(Oligomerization)
  acc <- suppressWarnings(data.frame(acc))
  row.names(acc) <- NULL
  data <- cbind(Oligomerization, acc)
  results <- vector("list", 100)
  if (model == "train") {
  for (i in 1:100) {
  trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
  train <- data[trainIndex,]
  test <- data[-trainIndex,]
  model <- RWeka::J48(Oligomerization~., data = train)
  summary <- summary(model)
  confusionmatrix <- summary$confusionMatrix
  results[[i]] <- as.numeric(confusionmatrix)
  }
  } else if (model == "cv") {
    for (i in 1:100) {
      trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
      train <- data[trainIndex,]
      test <- data[-trainIndex,]
      model_train <- RWeka::J48(Oligomerization ~., data = train)
      eval_j48 <- RWeka::evaluate_Weka_classifier(model_train, numFolds = 10, 
                                                  complexity = FALSE, seed = 1,
                                                  class = TRUE)
      confusionmatrix <- eval_j48$confusionMatrix
      results[[i]] <- as.numeric(confusionmatrix)
    }
  } else if (model == "test") {
    for (i in 1:100) {
    trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
    train <- data[trainIndex,]
    test <- data[-trainIndex,]
    model_train <- RWeka::J48(Oligomerization ~., data = train)
    eval_j48 <- RWeka::evaluate_Weka_classifier(model_train, newdata = test, 
                                                complexity = FALSE, seed = 1,
                                                class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
    }
  }
  return(results)
}

dpc_models <- function(x, model){
  library(protr)
  library(caret)
  library(RWeka)
  library(stringr)
  data <- protr::readFASTA(x)
  data <- data[(sapply(data, protr::protcheck))]
  dpc  = t(sapply(data, protr::extractDC))
  Oligomerization <- row.names(dpc)
  Dimer <- stringr::str_extract(Oligomerization, "Dimer")
  Monomer <- stringr::str_extract(Oligomerization, "Monomer")
  Tetramer <- stringr::str_extract(Oligomerization, "Tetramer")
  Oligomerization <- paste(Dimer, Monomer, Tetramer)
  Oligomerization <- sapply(Oligomerization, gsub, pattern = "NA", replacement = "")
  Oligomerization <- stringr::str_trim(Oligomerization, side = "both")
  Oligomerization <- ifelse(Oligomerization == "Monomer", "Monomer", "Oligomer")
  Oligomerization <- as.factor(Oligomerization)
  dpc <- suppressWarnings(data.frame(dpc))
  row.names(dpc) <- NULL
  data <- cbind(Oligomerization, dpc)
  results <- vector("list", 100)
  if (model == "train") {
    for (i in 1:100) {
      trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
      train <- data[trainIndex,]
      test <- data[-trainIndex,]
      model <- RWeka::J48(Oligomerization~., data = train)
      summary <- summary(model)
      confusionmatrix <- summary$confusionMatrix
      results[[i]] <- as.numeric(confusionmatrix)
    }
  } else if (model == "cv") {
    for (i in 1:100) {
      trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
      train <- data[trainIndex,]
      test <- data[-trainIndex,]
      model_train <- RWeka::J48(Oligomerization ~., data = train)
      eval_j48 <- RWeka::evaluate_Weka_classifier(model_train, numFolds = 10, 
                                                  complexity = FALSE, seed = 1,
                                                  class = TRUE)
      confusionmatrix <- eval_j48$confusionMatrix
      results[[i]] <- as.numeric(confusionmatrix)
    }
  } else if (model == "test") {
    for (i in 1:100) {
      trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
      train <- data[trainIndex,]
      test <- data[-trainIndex,]
      model_train <- RWeka::J48(Oligomerization ~., data = train)
      eval_j48 <- RWeka::evaluate_Weka_classifier(model_train, newdata = test, 
                                                  complexity = FALSE, seed = 1,
                                                  class = TRUE)
      confusionmatrix <- eval_j48$confusionMatrix
      results[[i]] <- as.numeric(confusionmatrix)
    }
  }
  return(results)
}




mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}




acc_results <- function(x, model, descriptors, ...) {
  if (descriptors == "ACC") {
  ok <- acc_models(x, model) }
  else if (descriptors == "DPC") {
    ok <- dpc_models(x, model) }
  results <- data.frame(ok)
  rm(ok)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  rm(m)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}



