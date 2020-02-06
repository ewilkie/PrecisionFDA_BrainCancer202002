
ela <- Sys.time()

#library(pacmac)
pacman::p_load(ggplot2, reshape, dplyr,caret,DMwR,purrr,pROC,doParallel,data.table,statmod, gam)

## Parallele processin
## 2GB per core * number of cvs 10 = 20GBs total 5 clusters 
cl <- makePSOCKcluster(12)
registerDoParallel(cl)

Phase1_GE <- fread("../Data/sc1_Phase1_GE_FeatureMatrix.txt", sep="\t", header=T)
Phase1_GE <- as.data.frame(Phase1_GE)
rownames(Phase1_GE) <- Phase1_GE[,1]
Phase1_GE <- Phase1_GE[,-1]
gene.names <- gsub("\\.", "-",colnames(Phase1_GE))
colnames(Phase1_GE) <- gene.names

## gene expression
hu <- read.csv("../Data/Sc1_GE_names.csv", sep=",", header=T)

## extract data
x <- Phase1_GE[,which(colnames(Phase1_GE) %in% as.vector(hu$Input))]

zv <- apply(x, 2, function(x) length(unique(x)) == 1)
x <- x[, !zv]
varlrds <- list()
varlrds$zv <- zv
descrCorr <- cor(x)
highCorr <- findCorrelation(descrCorr, 0.3)
indata <- x[,-highCorr]
dim(indata)

varlrds$higcorr <- colnames(x)[highCorr]
saveRDS(varlrds, "../Output/Varlsel_GE.rds")

## read in outcome
y <- read.table("../Data/sc1_Phase1_GE_Outcome.txt", sep="\t", header=T, row.names=1)

set.seed(5627)
options(expressions = 5e5)

## first item in loop is non-permuted sample
rds_out <- list()
for(i in 1:3){
  print(i)
  if(i == 1){
    prm <- y$SURVIVAL_STATUS
  }else{
    prm <- y[sample(nrow(y)),]
  }
  ## Important variable selection
  sbftest <- sbf(x=indata,y=prm, sbfControl = sbfControl(functions = caretSBF,verbose = TRUE, method = "cv", number = 5, saveDetails = TRUE))
  
  selvar <- sbftest$optVariables
  print("number of variables selected")
  print(length(selvar))
  
  rds_out[[i]] <- list(selvar)
  #rds_out[[i]]$selvar <- selvar
  
  seldata <- indata[,which(colnames(indata) %in% selvar)]
  
  ## build classifier control
  ctrl <- trainControl(method = "repeatedcv",number = 10,repeats = 5, summaryFunction = twoClassSummary, classProbs = TRUE)
  
  ## change outcome variable coding
  prm[which(prm == 0)] <- "Alive"
  prm[which(prm == 1)] <- "Dead"
  
  ## format data
  x_sub <- seldata
  data <- cbind(prm, x_sub)
  colnames(data)[1] <- "Class"
  dim(data)
  
  print("Starting models...")
  
  ## train data with gbm, knn and nb
  gbm_fit <- train(Class ~ ., data = data, method = "gbm",verbose = FALSE, metric = "ROC", trControl = ctrl)
  knn_fit <- train(Class ~ ., data = data, method = "knn", metric="ROC", trControl = ctrl)
  nb_fit <- train(Class ~ ., data = data, method = "naive_bayes", metric="ROC", trControl = ctrl)
  glmnet_fit <- train(Class ~ ., data = data, method = "glmnet", metric="ROC", trControl = ctrl)
  
  fits <- list(gbm_fit,glmnet_fit,knn_fit,nb_fit)
  
  rds_out[[i]]$fits <- fits
  
  ## extract and combine relevant data
  roc_data_gbm <- roc(data$Class,predict(gbm_fit, data, type = "prob")[,"Dead"])
  pd_gbm <- cbind(roc_data_gbm$specificities,roc_data_gbm$sensitivities, "GBM")
  colnames(pd_gbm) <- c("specificities","sensitivities", "method")
  
  roc_data_knn <- roc(data$Class,predict(knn_fit, data, type = "prob")[,"Dead"])
  pd_knn <- cbind(roc_data_knn$specificities,roc_data_knn$sensitivities, "KNN")
  colnames(pd_knn) <- c("specificities","sensitivities", "method")
  
  roc_data_nb <- roc(data$Class,predict(nb_fit, data, type = "prob")[,"Dead"])
  pd_nb <- cbind(roc_data_nb$specificities,roc_data_nb$sensitivities, "NB")
  colnames(pd_nb) <- c("specificities","sensitivities", "method")
  
  roc_data_glmnet <- roc(data$Class,predict(glmnet_fit, data, type = "prob")[,"Dead"])
  pd_glmnet <- cbind(roc_data_glmnet$specificities,roc_data_glmnet$sensitivities, "glmnet")
  colnames(pd_glmnet) <- c("specificities","sensitivities", "method")
  
  pdfs <- rbind(pd_gbm,pd_knn,pd_nb,pd_glmnet)
  pdfs <- as.data.frame(pdfs)
  
  ## format
  pdfs$specificities <- as.numeric(as.vector(pdfs$specificities))
  pdfs$sensitivities <- as.numeric(as.vector(pdfs$sensitivities))
  
  ## TO DO annotate with AUC value
  gp <- ggplot(pdfs, aes(x=specificities,y=sensitivities,colour=method)) + geom_point(size=.5) + geom_line(size=1.15) + scale_x_continuous(trans = "reverse") + ggtitle("AUC - Comparisons on Stochastic Gradient Boosting (GBM), K-nearest neighbours (KNN), Naive Bayes (NB) and glmnet")
  outf <- paste("../Plots/SC1_AUC_perm",i, "_v2.pdf", sep="")
  ggsave(outf, gp)
  
  saveRDS(rds_out, file = "../Output/Fits_GE_permutation_v2.rds")
}

