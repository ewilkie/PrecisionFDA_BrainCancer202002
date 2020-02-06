############################################
## Phase 1 - Sub Challenge 2 - Submission ##
############################################

pacman::p_load(ggplot2, reshape, dplyr,caret,DMwR,purrr,pROC,doParallel,data.table,statmod, gam)

fits <- readRDS("/Users/ewilkie/Documents/PrecisionFDA/Remote/Fits_GE_permutation_v2.rds")
summary(fits)

## chosen model
fit <- fits[[1]]$fits[[4]]

## need to get variables 
selv <- fit$finalModel$xNames
r <- gsub("`", "", selv)

Phase1_GE <- fread("../Data/sc1_Phase1_GE_FeatureMatrix.txt", sep="\t", header=T)
Phase1_GE <- as.data.frame(Phase1_GE)
rownames(Phase1_GE) <- Phase1_GE[,1]
Phase1_GE <- Phase1_GE[,-1]
gene.names <- gsub("\\.", "-",colnames(Phase1_GE))
colnames(Phase1_GE) <- gene.names

data <- Phase1_GE[,which(colnames(Phase1_GE) %in% r)]


y <- read.table("../Data/sc1_Phase1_GE_Outcome.txt", sep="\t", header=T, row.names=1)
## change outcome variable coding
prm <- y$SURVIVAL_STATUS
prm[which(prm == 0)] <- "Alive"
prm[which(prm == 1)] <- "Dead"

indata <- cbind(prm,data)
colnames(indata)[1] <- "Class"
dim(indata)

## predict data
pred <- predict(fit, newdata=indata, type = "prob")

## find which ones do and don't match
comb <- cbind(pred, indata$Class)

check <- list()
for (i in 1:nrow(comb)){
  ind <- which.max(comb[i,1:2])
  check[[i]] <- names(ind)
}

check_vec <- unlist(check)
merg <- cbind(comb, check_vec)

## classes don't match
wrong_class <- merg[which(as.vector(merg[,3]) != as.vector(merg[,4])),]

##b) Short listed features selected by the model
v <- varImp(fit)$importance
ver <- v[which(v[,1] !=0 & v[,2] != 0 ),]
paste(rownames(ver) , collapse=",")

## d) Confusion matrix
Confusion.matrix <- confusionMatrix(as.factor(check_vec), as.factor(indata$Class))
print(Confusion.matrix$table)

## e) Overall accuracy
error_rate <- nrow(wrong_class) / length(y$SURVIVAL_STATUS)
OA <- 1- error_rate
print(OA)

## f) Specificity 
spec <- Confusion.matrix$byClass[2]
print(spec)

#g) Sensitivity 
sens <- Confusion.matrix$byClass[1]
print(sens)
## Sensitivity 
## 0.6046512

## h) Area under the curve (AUC)
## plot ROC with AUC label
plot(roc(indata$Class,predict(fit, indata, type = "prob")[,"Dead"]), print.auc=TRUE)
