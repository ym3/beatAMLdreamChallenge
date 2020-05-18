library(readr)
library(ggfortify)
library(gplots)
library(RColorBrewer)
library(DESeq2)
library(limma)
library(pls)
library(EnhancedVolcano)
library(pheatmap)
library(NbClust)
library(factoextra)
library(dendextend)
library(corrplot)
library(dplyr)
library(tibble)
library(reshape2)
library(caret)
library(FeatureHashing)
library(xgboost)
library(doParallel)
library(gridExtra)
registerDoParallel(detectCores())
library(BiocParallel)
library(e1071)
library(survival)
library(WGCNA)
library(ape)
registered()
bpparam()

Sys.setenv(OMP_NUM_THREADS=1)
Sys.setenv(OPENBLAS_NUM_THREADS=1)

dnaseq <- read_csv("dnaseq.csv")
rnaseq <- read_csv("rnaseq.csv")
aucs <- read_csv("aucs.csv")
clinical_categorical <- read_csv("clinical_categorical.csv")
clinical_numerical <- read_csv("clinical_numerical.csv")
clinical_categorical_legend <- read_csv("clinical_categorical_legend.csv")
colnames(clinical_categorical) <- make.names(colnames(clinical_categorical))
response <- read_csv("response.csv")

# drop columns with missing values 
clinical_numerical <- clinical_numerical %>% select(-c("%.Blasts.in.PB","WBC.Count"))

# drop non-recurrent genes as they reduce the model fit and cause sparsity
#dnaseq <- dnaseq[dnaseq$Hugo_Symbol %in% names(which(table(dnaseq$Hugo_Symbol) > 2)), ]
#barplot(summary(as.factor(dnaseq$Hugo_Symbol),maxsum=40),las=2)

# prepare rnaseq data
rna_log2counts <- as.matrix(rnaseq[,3:dim(rnaseq)[2]])
rownames(rna_log2counts) <- rnaseq$Gene
rna_counts <- round(2^rna_log2counts)
rownames(rna_counts) <- rownames(rna_log2counts)

# rnaseq cpm check, sum of total cpm must be the same for all samples
colSums(2^rna_log2counts)

# variance stabilizing transformation (VST)
#vstcounts <- vst(rna_counts, blind = FALSE)
#rownames(vstcounts) <- rownames(rna_log2counts)
#rna_log2counts <- vstcounts
#for (i in 1:dim(rna_log2counts)[2]) {
#  par(mfrow=c(1,2))
#  plotMA(rna_log2counts, array = i)
#  abline(h=0,col="grey")
#  plotMA(vstcounts, array = i)
#  abline(h=0,col="grey")
# composition bias is minimal so no need for VST as 
# differential testing also utilises raw counts

  
# caret's training control
fold <- 5
res <- 5
set.seed(42)
seeds <- vector(mode = "list", length = 501)
for(i in 1:500) seeds[[i]] <- sample.int(1000, 600)
seeds[[501]] <- sample.int(1000,1)
train_ctrl <- trainControl(method="repeatedcv",
                           number = fold,
                           repeats = res,
                           savePredictions = TRUE,
                           allowParallel = FALSE,
                           search = "grid",
                           seeds = seeds)

model.fits <- list()
dnaseqGenes <- list()
rnaseqGenes <- list()
rna.preProc <- list()
dev.off()
pdf("res.pdf")

Genes_pvalues <- c()
for (drug in levels(as.factor(aucs$inhibitor))) {
  
  # feature selection (genotypes)
  dnaseq_aucs <- dcast(
    dnaseq, lab_id ~ Hugo_Symbol, fun.aggregate = 
      function(x)(as.numeric(!isEmpty(x))) ) %>% 
      select(-c("NPM1")) %>%
    inner_join(
    as.data.frame(aucs) %>% filter(inhibitor==drug) 
    %>% select(c("lab_id","auc")), by = "lab_id") %>%
    inner_join(
      select(clinical_categorical,c("lab_id","NPM1","FLT3.ITD")), 
      by = "lab_id") 
  
  rownames(dnaseq_aucs)=dnaseq_aucs$lab_id
  dnaseq_aucs <- dnaseq_aucs %>% select(-c("lab_id"))
  colnames(dnaseq_aucs) <- make.names(colnames(dnaseq_aucs))
  
  # keep genes with more than 2 occurances among the samples
  dnaseq_aucs <- dnaseq_aucs[, colSums(dnaseq_aucs) > 1]
  
  # Student's t-test on individual genes
  t.tests <- list()
  symbols <- names(dnaseq_aucs)
  symbols <- symbols[!symbols %in% "auc"]
  
  avg_diff <- apply(dnaseq_aucs, 2, function(x) 
    mean(dnaseq_aucs$auc[which(x==0)])-mean(dnaseq_aucs$auc[which(x==1)]))
  
  for (symbol in symbols){
    t <- t.test(as.formula(paste("auc ~ ", symbol)), 
                data=dnaseq_aucs, var.equal=TRUE, conf.level=0.95)
    t.tests[[symbol]] <- t
    pvalues <- data.frame(p_value = t$p.value,
                          avg_diff = avg_diff[symbol],
                          symbol = symbol)
    Genes_pvalues = rbind(Genes_pvalues,pvalues)
  }
  
  # select statistically significant genes as dnaseq features
  topGenes_dnaseq <- c()
  for (symbol in names(t.tests)){
    if (t.tests[[symbol]]$p.value<0.05){
      topGenes_dnaseq <- c(topGenes_dnaseq, symbol)
    }
  }
  print(drug)
  print(topGenes_dnaseq)
  
  # dimensionality reduction (expression)
  samples_with_auc = as.data.frame(aucs) %>% 
      filter(inhibitor==drug) %>% select(c("lab_id","auc")) 
  
  rnaseq_voom = voom(countdata)$E
  topGenes_rnaseq = data.frame(Gene = rownames(rnaseq_voom[
    order(apply(rnaseq_voom,1,mad), decreasing = T)[1:1000],]))
  
  sigDE_log2counts <- rna_log2counts[which(rownames(rna_log2counts) 
                                           %in% topGenes_rnaseq$Gene),]
  sigDE_log2counts <- sigDE_log2counts %>% subset(select=samples_with_auc$lab_id)

  # one hot encoding of presence/absense mutations
  full_data <- dcast(
    dnaseq, lab_id ~ Hugo_Symbol, fun.aggregate =
      function(x)(as.numeric(!isEmpty(x))) ) %>% select(-c("NPM1"))
  
  # set mutations to 0 for samples with rnaseq data but not in dnaseq.csv
  missingGenos = setdiff(colnames(rna_log2counts),full_data$lab_id)
  mtx <- cbind(missingGenos,matrix(0,length(missingGenos),dim(full_data)[2]-1))
  colnames(mtx) <- c("lab_id",colnames(full_data)[-1])
  full_data = rbind(full_data, mtx)
  
  # append with clinical variables and intersect with aucs
  full_data <- full_data %>%
    inner_join(
      select(clinical_categorical,c("lab_id","NPM1","FLT3.ITD")),
      by = "lab_id") %>%
    select(c("lab_id",topGenes_dnaseq)) %>% 
    #mutate_if(is.integer, as.factor) %>%
    inner_join(
      as.data.frame(aucs) %>% filter(inhibitor==drug) 
      %>% select(c("lab_id","auc")), by = "lab_id") %>%
    inner_join(
      select(clinical_categorical,-c("NPM1","FLT3.ITD")) %>%
       mutate_at(vars(-one_of("lab_id")),as.numeric), by = "lab_id") %>%
    inner_join(clinical_numerical, by = "lab_id")

  # split into training and test sets
  trainIndex <- createDataPartition(y=full_data$auc, times=1, p=0.7, list=F)
  full_data_Train <- full_data #full_data[trainIndex,]
  full_data_Test <- full_data #full_data[-trainIndex,]
  aucs_Train <- full_data[,c("lab_id","auc")] #[trainIndex,c("lab_id","auc")] #[,c("lab_id","auc")] #
  aucs_Test <- full_data[,c("lab_id","auc")] #[-trainIndex,c("lab_id","auc")] #[,c("lab_id","auc")] #
  
  # project top rnaseq genes onto PCA of training set
  sigDE_log2counts_Train <- t(sigDE_log2counts)[
    which(rownames(t(sigDE_log2counts)) %in% full_data_Train$lab_id),]
  sigDE_log2counts_Test <- t(sigDE_log2counts)[
    which(rownames(t(sigDE_log2counts)) %in% full_data_Test$lab_id),]
  
  sigDE_log2counts_Train <- sigDE_log2counts_Train %>% 
    as.data.frame %>% rownames_to_column(var="lab_id")
  sigDE_log2counts_Test <- sigDE_log2counts_Test %>% 
    as.data.frame %>% rownames_to_column(var="lab_id")     
  
  sigDE_log2counts_Train <- sigDE_log2counts_Train %>%
    inner_join( aucs_Train, by = "lab_id")
  sigDE_log2counts_Test <- sigDE_log2counts_Test %>%
    inner_join( aucs_Test, by = "lab_id")
  
  rownames(sigDE_log2counts_Train)=sigDE_log2counts_Train$lab_id
  sigDE_log2counts_Train <- sigDE_log2counts_Train %>% select(-c("lab_id"))
                                                                          
  rownames(sigDE_log2counts_Test)=sigDE_log2counts_Test$lab_id
  sigDE_log2counts_Test <- sigDE_log2counts_Test %>% select(-c("lab_id"))
  
  
  Ncomp <- min(10,dim(sigDE_log2counts_Train)[2]-1)
  sigDE_preProc <- plsr(auc ~ ., data = sigDE_log2counts_Train, 
                        ncomp = Ncomp, scale = TRUE, validation = "CV")
  plot(RMSEP(sigDE_preProc), legendpos = "topright")
  explvar(sigDE_preProc)
  plot(sigDE_preProc, plottype = "scores", comps = 1:min(Ncomp,3))
  Ncomp_randm <- selectNcomp(sigDE_preProc, method = "randomization", plot = TRUE)
  Ncomp_sigma <- selectNcomp(sigDE_preProc, method = "onesigma", plot = TRUE)
  Ncomp<-max(1,Ncomp_randm,Ncomp_sigma)
  print(Ncomp)
  PLS_sigDE_Train <- predict(sigDE_preProc, comps = 1:Ncomp, type = "scores", 
                             newdata = sigDE_log2counts_Train) 
  PLS_sigDE_Test <- predict(sigDE_preProc, comps = 1:Ncomp, type = "scores", 
                            newdata = sigDE_log2counts_Test)
  
  if (Ncomp!=0){
  full_data_Train <- full_data_Train %>%
    inner_join( as.data.frame(PLS_sigDE_Train) %>% select(1:Ncomp) %>%
                  mutate(lab_id=rownames(PLS_sigDE_Train)), by = "lab_id") 
  full_data_Test <- full_data_Test %>%
    inner_join( as.data.frame(PLS_sigDE_Test) %>% select(1:Ncomp) %>%
                  mutate(lab_id=rownames(PLS_sigDE_Test)), by = "lab_id") 
  }
  
  # drop label column and convert to numeric
  full_data_Train <- full_data_Train %>% 
    mutate_at(vars(-one_of("lab_id")),funs(as.numeric(as.character(.))))
  rownames(full_data_Train)=full_data_Train$lab_id
  full_data_Train <- full_data_Train %>% select(-c("lab_id"))
  
  full_data_Test <- full_data_Test %>% 
    mutate_at(vars(-one_of("lab_id")),funs(as.numeric(as.character(.))))
  rownames(full_data_Test)=full_data_Test$lab_id
  full_data_Test <- full_data_Test %>% select(-c("lab_id"))
  
  # mean impute missing values as only a couple NAs
  cM <- colMeans(full_data_Train, na.rm=TRUE)
  indx <- which(is.na(full_data_Train), arr.ind=TRUE)
  if (!isEmpty(indx)){ full_data_Train[indx] <- cM[indx[,2]]}

  cM <- colMeans(full_data_Test, na.rm=TRUE)
  indx <- which(is.na(full_data_Test), arr.ind=TRUE)
  if (!isEmpty(indx)){ full_data_Test[indx] <- cM[indx[,2]]}

  ## drop factors in training set with one level
  fctr <- which(sapply(full_data_Train, is.factor))
  nlevs <- lengths(lapply(full_data_Train[fctr], levels))
  if(length(which(nlevs<2))>0){
    full_data_Train <- full_data_Train[,-fctr[which(nlevs<2)]]
    full_data_Test <- full_data_Test[,-fctr[which(nlevs<2)]]
  }
  
  xgb.tuneGrid = expand.grid( nrounds = c(50,100),
                              max_depth = c(4,6,8),
                              eta = c(0.02,0.1,0.3),
                              gamma = c(0,2),
                              colsample_bytree=c(0.4,0.8),
                              subsample=c(0.5,1.0),
                              min_child_weight = 1)
  model.fit <- train(auc ~ .,
                  data = full_data_Train,
                  method = "xgbTree",
                  #tuneGrid = xgb.tuneGrid,
                  tuneLength = 2,
                  trControl = train_ctrl,
                  importance = TRUE,
                  metric = "RMSE",
                  objective = "reg:squarederror")
  
   print(drug)
   #print(model.fit)

   model.imp <- varImp(model.fit)
   print(model.imp)
   p1 <- theme(axis.text=element_text(size=10,face="bold"),
               axis.title=element_text(size=22,face="bold"),
               plot.title = element_text(size=22))
   plot4 <- ggplot(model.imp, top = 50) +
     labs(title = drug, x='') + p1
  
   cor_cv_pearson <- cor(model.fit$pred$pred,model.fit$pred$obs, method = "pearson")
   cor_cv_spearman <- cor(model.fit$pred$pred,model.fit$pred$obs, method = "spearman")
   test_pred <- predict(model.fit, full_data_Test)
   cor_pearson = cor(full_data_Test$auc, test_pred, method = "pearson")
   cor_spearman = cor(full_data_Test$auc, test_pred, method = "spearman")
   print(c(drug, cor_cv_pearson, cor_cv_spearman))
   print(c(drug, cor_pearson, cor_spearman))
   plot5 <- qplot(full_data_Test$auc, test_pred) +
     xlab("observed auc") +
     ylab("predicted auc") +
     geom_smooth(method = "lm", se = FALSE) +
     annotate("text",label =
                paste("R_ p, s = ", round(cor_pearson,digits=2), ",",
                      round(cor_spearman,digits=2)),
              x=Inf, y =Inf, vjust=1, hjust=1)


  print(plot4)
  print(plot5)
  
  model.fits[[drug]] <- model.fit
  dnaseqGenes[[drug]] <- topGenes_dnaseq
  rnaseqGenes[[drug]] <- topGenes_rnaseq$Gene
  rna.preProc[[drug]] <- sigDE_preProc
}

dev.off()

saveRDS(model.fits, file="model-fits.rds")
saveRDS(dnaseqGenes, file="dnaseqGenes.rds")
saveRDS(rnaseqGenes, file="rnaseqGenes.rds")
saveRDS(rna.preProc, file="rna-preProc.rds")
