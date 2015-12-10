library(protr)
OliFP_100 <- protr::readFASTA("OliFP_100.Fasta")
OliFP_99 <- protr::readFASTA("OliFP_99.Fasta")
OliFP_95 <- protr::readFASTA("OliFP_95.Fasta")
OliFP_90 <- protr::readFASTA("OliFP_90.Fasta")

OliFP_100 <- OliFP_100[(sapply(OliFP_100, protcheck))]
OliFP_99 <- OliFP_99[(sapply(OliFP_99, protcheck))]
OliFP_95 <- OliFP_95[(sapply(OliFP_95, protcheck))]
OliFP_90 <- OliFP_90[(sapply(OliFP_90, protcheck))]

acc_models <- function(x){
  library(protr)
  library(caret)
  library(RWeka)
  library(stringr)
  data <- protr::readFASTA(x)
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
  trainIndex <- caret::createDataPartition(data$Oligomerization, p = 0.8, list = FALSE)
  train <- data[trainIndex,]
  test <- data[-trainIndex,]
  model <- Rweka::J48(Oligomerization~., data = train)
  results <- summary(model)
  return(results)
}



Dimer <- str_extract(Oligomerization, "Dimer")
Monomer <- str_extract(Oligomerization, "Monomer")
Tetramer <- str_extract(Oligomerization, "Tetramer")
Oligomerization <- paste(Dimer, Monomer, Tetramer)
Oligomerization <- sapply(Oligomerization, gsub, pattern = "NA", replacement = "")
Oligomerization <- str_trim(Oligomerization, side = "both")
Oligomerization <- ifelse(Oligomerization == "Monomer", "Monomer", "Oligomer")

Oligomerization <- as.factor(Oligomerization)
acc  = t(sapply(OliFP_100, extractAAC))
acc <- suppressWarnings(data.frame(acc))
row.names(acc) <- NULL
data <- cbind(Oligomerization, acc)

data <- ifelse(Oligomerization == "Monomer", "Monomer", "Oligomer")




model <- J48(Oligomerization~., data = data)

#Oligomerization <- data.frame(Oligomerization = row.names(acc))
Oligomerization <- row.names(acc)
for (i in Oligomerization){
  data <- string
}





library(stringr)
word.list = str_split(Oligomerization, '\\s+')
new_list <- unlist(word.list)
Dimer <- new_list[grep("Dimer", new_list)]
Dimer <- data.frame(Dimer)
Monomer <- new_list[grep("Monomer", new_list)]
Tetramer <- new_list[grep("Tetramer", new_list)]
Dimer$Oligomerization <- Dimer
Monomer$Oligomerization <- Monomer
Tetramer$Oligomerization <- Tetramer