############################################
## Phase 1 - Sub Challenge 2 - Submission ##
############################################

fits <- readRDS("/Users/ewilkie/Documents/PrecisionFDA/Remote/Fits_CN_permutation_v2.rds")
summary(fits)

## chosen model
fit <- fits[[1]]$fits[[1]]

## need to get variables 
selv <- summary(fit, plot=F)$var
r <- gsub("`", "", selv)


## Cleaning data
Phase1_CN_out <- read.table("../Data/sc2_Phase1_CN_FeatureMatrix.txt", sep="\t", header=T, row.names=1)

x <- Phase1_CN_out
for(col in 1:ncol(x)){
  minv <- min(x[-which(x[,col] == 0),col], na.rm=T)
  x[which(x[,col] == 0),col] <- minv
}

x[is.na(x)] <- min(x)
x[x == Inf] <- min(x)

notin <- r[-which(r %in% colnames(x))]
length(notin)
colnames(x)[which(colnames(x) %in% notin)]

data <- x[,which(colnames(x) %in% r)]

y <- read.table("../Data/sc2_Phase1_CN_Outcome.txt", sep="\t", header=T, row.names=1)
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
b <- which(summary(fit, plot=F)$rel.inf != 0)
selv <- summary(fit, plot=F)$var[b]
r <- gsub("`", "", selv)
paste(r, collapse=",")

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

## h) Area under the curve (AUC)
## plot ROC with AUC label
plot(roc(indata$Class,predict(fit, indata, type = "prob")[,"Dead"]), print.auc=TRUE)
