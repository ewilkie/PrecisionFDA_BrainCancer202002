############################################
## Phase 1 - Sub Challenge 3 - Submission ##
############################################

fits <- readRDS("/Users/ewilkie/Documents/PrecisionFDA/Remote/Fits_GE_CN_permutation_v2.rds")
summary(fits)

## chosen model
fit <- fits[[1]]$fits[[4]]

## need to get variables 
selv <- fit$finalModel$xNames
r <- gsub("`", "", selv)

## read in patient data
Phase1_CN_GE_out <- read.table("../Data/sc3_Phase1_CN_GE_FeatureMatrix.txt", sep="\t", header=T, row.names=1)
dim(Phase1_CN_GE_out) 

cyto.ind <- grep("^X[0-9]", colnames(Phase1_CN_GE_out))
cyto <- Phase1_CN_GE_out[,cyto.ind]
dim(cyto)

x <- cyto
for(col in 1:ncol(x)){
  minv <- min(x[-which(x[,col] == 0),col], na.rm=T)
  x[which(x[,col] == 0),col] <- minv
}

x[x == NA] <- min(x)
x[x == Inf] <- min(x)

## check for SD
sd <- apply(x, 2 , sd)
x <- x[,-which(sd < 0.1)]

## combine data

exp <- Phase1_CN_GE_out[,-cyto.ind]
colnames(exp) <- gsub("\\.", "-",colnames(exp))
d <- cbind(x, exp)

notin <- r[-which(r %in% colnames(d))]
length(notin)
colnames(x)[which(colnames(x) %in% notin)]

## select variables
data <- d[,which(colnames(d) %in% r)]
y <- read.table("../Data/sc3_Phase1_CN_GE_Outcome.txt", sep="\t", header=T, row.names=1)
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

## h) Area under the curve (AUC)
## plot ROC with AUC label
plot(roc(indata$Class,predict(fit, indata, type = "prob")[,"Dead"]), print.auc=TRUE)
