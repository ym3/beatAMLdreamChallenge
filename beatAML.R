library(readr)
library(ggfortify)
library(gplots)
library(RColorBrewer)
library(DESeq2)
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

# drop columns with missing values 
clinical_numerical <- clinical_numerical %>% select(-c("%.Blasts.in.PB","WBC.Count"))

# drop non-recurrent genes as they reduce the model fit and cause sparsity
#dnaseq <- dnaseq[dnaseq$Hugo_Symbol %in% names(which(table(dnaseq$Hugo_Symbol) > 2)), ]

# prepare rnaseq data
rna_log2counts <- as.matrix(rnaseq[,3:dim(rnaseq)[2]])
rownames(rna_log2counts) <- rnaseq$Gene
rna_counts <- round(2^rna_log2counts)
rownames(rna_counts) <- rownames(rna_log2counts)

# rnaseq cpm check, sum of total cpm must be the same for all samples
colSums(2^rna_log2counts)

# # variance stabilizing transformation (VST)
# vstcounts <- vst(rna_counts, blind = FALSE)
# rownames(vstcounts) <- rownames(rna_log2counts)
#
# composition bias is minimal so no need for VST as 
# differential testing also utilises raw counts
#for (i in 1:dim(rna_log2counts)[2]) {
#  par(mfrow=c(1,2))
#  plotMA(rna_log2counts, array = i)
#  abline(h=0,col="grey")
#  plotMA(vstcounts, array = i)
#  abline(h=0,col="grey")
#}

# # PCA on log2 counts (rather than vstcounts) better 
# # correspond to PC1=14.4% and PC2=7.8% 
# # from doi:10.1038/ncomms6901
# pcDat <- prcomp(t(rna_log2counts))
# autoplot(pcDat)

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
  for (symbol in symbols){
    t <- t.test(as.formula(paste("auc ~ ", symbol)), 
                data=dnaseq_aucs, var.equal=TRUE, conf.level=0.95)
    t.tests[[symbol]] <- t
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
  
  # standardize auc if supplied to DESeq as continuous variable
  #samples_with_auc$auc <- (samples_with_auc$auc-mean(samples_with_auc$auc))/sd(samples_with_auc$auc)
  #samples_for_DESeq <- intersect(samples_with_auc$lab_id,colnames(rna_counts))
  
  # or make binary class for 20% top (responders) and bottom (non-responders)
  n_top = 20 #round(0.20*dim(samples_with_auc)[1])
  top_aucs <- bind_rows(
      samples_with_auc %>% top_n(n_top,auc) %>% mutate(auc_class="top"),
      samples_with_auc %>% top_n(-n_top,auc) %>% mutate(auc_class="bottom")) %>% 
    mutate_at("auc_class", factor)
  top_aucs <- top_aucs[!duplicated(top_aucs$lab_id),]
  samples_for_DESeq <- intersect(top_aucs$lab_id,colnames(rna_counts))
  countdata <- as.data.frame(rna_counts) %>% subset(select=samples_for_DESeq)
  
  
  if(any(colnames(countdata) != samples_for_DESeq)) {
    stop ("samplesinfo must be in the same order as columns of countData")
  } else {
    ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = top_aucs, #samples_with_auc, ,
                                     design = as.formula(~ auc_class))
  }
  ddsObj <- DESeq(ddsMat, parallel=TRUE)
  res <- results(ddsObj, alpha=0.05)
  res <- as.data.frame(res) %>%
    mutate(Gene=rownames(res)) %>%
    inner_join(
      as.data.frame(rnaseq) %>% select(c("Gene","Symbol")), by = "Gene") %>%
    arrange(padj)
  
  plot1 <- EnhancedVolcano(res,
                lab = res$Symbol, #rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = drug,
                xlim = c(-5, 8))
  print(plot1)

  # extract significant genes with adjusted pvalue<0.05
  if(dim(subset(res,padj<0.05 & log2FoldChange>1))[1]<2) {
    topGenes_rnaseq<-head(res,2) } else {
      topGenes_rnaseq <- res %>% subset(padj<0.05 & log2FoldChange>1) }
    
  # check overlap with Ibrutinib's top genes
  # drug = "Ibrutinib (PCI-32765)"
  # Ibrutinib_topGenes <- #make.names(
  #   c("B3GNT2","RAB43P1","KLHDC3","GYPC","MIR223",
  #     "RP11-333E13.2","CHD5","NAPSB","ERMN","FAM49A","ADAP2",
  #     "TLR2","HNMT","CD86","HSD17B13","RP11-196G18.24","HNRNPA1")#)
  # intersect(Ibrutinib_topGenes,topGenes_rnaseq$Symbol)
  
  sigDE_log2counts <- rna_log2counts[which(rownames(rna_log2counts) 
                                           %in% topGenes_rnaseq$Gene),]
  mat  <- sigDE_log2counts %>% subset(select=samples_for_DESeq)
  mat  <- mat - rowMeans(mat)
  anno <- data.frame(auc_class=top_aucs$auc_class) #auc=samples_with_auc$auc)
  rownames(anno) = samples_for_DESeq
  plot2 <-pheatmap(mat, annotation_col = anno, labels_row = rep("",dim(mat)[1]),
           annotation_legend = TRUE, annotation_names_col = FALSE, 
           main = paste(drug," (n=",dim(mat)[2],", g=",dim(mat)[1],")",sep=""))
  
  sigDE_log2counts <- sigDE_log2counts %>% subset(select=samples_with_auc$lab_id)
  #pcDat <- prcomp(t(sigDE_log2counts))
  
  #plot3 <-autoplot(pcDat, samples_with_auc, colour = 'auc')
  #print(plot3)
  
  # # Sample distances
  # clusts <- NbClust(data = t(sigDE_log2counts), distance = "euclidean", 
  #                   method="complete", index="silhouette", min.nc = 2, max.nc = 10)
  # n.clusters = clusts$Best.nc[["Number_clusters"]] #round(mean(clusts$Best.nc[1,]))
  # print(n.clusters)
  # sampleDists=dist(t(sigDE_log2counts))
  # pheatmap(sampleDists)
  # hclust_avg <- hclust(sampleDists, method="complete")
  # clusters <- cutree(hclust_avg,n.clusters,order_clusters_as_data = FALSE)
  # summary(as.factor(clusters))
  # avg_dend_obj <- as.dendrogram(hclust_avg)
  # avg_col_dend <- color_branches(avg_dend_obj, clusters = clusters)
  # plot(avg_col_dend)

  # one hot encoding of presence/absense of driver genes
  full_data <- dcast(
    dnaseq, lab_id ~ Hugo_Symbol, fun.aggregate =
      function(x)(as.numeric(!isEmpty(x))) ) %>% select(-c("NPM1"))
  
  # set drivers to 0 for samples with rnaseq data but not in dnaseq.csv
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
      %>% select(c("lab_id","auc")), by = "lab_id") #%>%
    # inner_join(
    #   select(clinical_categorical,-c("NPM1","FLT3.ITD")) %>%
    #     mutate_at(vars(-one_of("lab_id")),as.numeric), by = "lab_id") %>%
    # inner_join(clinical_numerical, by = "lab_id")

  # split into training and test sets
  trainIndex <- createDataPartition(y=full_data$auc, times=1, p=0.7, list=F)
  full_data_Train <- full_data[trainIndex,]
  full_data_Test <- full_data[-trainIndex,]
  aucs_Train <- full_data[trainIndex,c("lab_id","auc")] #c("lab_id","auc")] #
  aucs_Test <- full_data[-trainIndex,c("lab_id","auc")] #c("lab_id","auc")] #
  
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
  
  #sigDE_preProc <- preProcess(sigDE_log2counts_Train, method = "pca")
  #PCA_sigDE_Train <- predict(sigDE_preProc, sigDE_log2counts_Train)
  #PCA_sigDE_Test <- predict(sigDE_preProc, sigDE_log2counts_Test)
  
  Ncomp <- min(10,dim(sigDE_log2counts_Train)[2]-1)
  sigDE_preProc <- plsr(auc ~ ., data = sigDE_log2counts_Train, 
                        ncomp = Ncomp, validation = "CV")
  plot(RMSEP(sigDE_preProc), legendpos = "topright")
  explvar(sigDE_preProc)
  plot(sigDE_preProc, plottype = "scores", comps = 1:min(Ncomp,3))
  Ncomp <- Ncomp #selectNcomp(sigDE_preProc, method = "randomization", plot = TRUE)
  PLS_sigDE_Train <- predict(sigDE_preProc, comps = 1:Ncomp, type = "scores", 
                             newdata = sigDE_log2counts_Train) 
  PLS_sigDE_Test <- predict(sigDE_preProc, comps = 1:Ncomp, type = "scores", 
                            newdata = sigDE_log2counts_Test)
  
  # append with PCA of rnaseq differentially expressed genes
  #Ncomp <- min(10,dim(PCA_sigDE_t)[2])
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
  
  # #hash_size = 2^19 #hash.size(full_data)
  # full_data_Train_model <- model.matrix(auc ~ .-1, full_data_Train)
  # full_data_Test_model <- model.matrix(auc ~ .-1, full_data_Test)
  # #full_data_Train_hashed <- hashed.model.matrix(
  # #  auc ~ ., data=full_data_Train, hash.size = hash_size, create.mapping = TRUE)
  # #full_data_Test_hashed <- hashed.model.matrix(
  # #  auc ~ ., data=full_data_Test, hash.size = hash_size, create.mapping = TRUE)
  # 
  # dtrain = xgb.DMatrix(full_data_Train_model, label=full_data_Train$auc)
  # dtest = xgb.DMatrix(full_data_Test_model, label=full_data_Test$auc)
  # watch <- list(train = dtrain, valid = dtest)
  # 
  # xgb.fit <- xgb.train(#booster = "gblinear", nrounds = 10, eta = 0.02,
  #                  booster = "gbtree", max.depth=7, eta=0.02, 
  #                  subsample = 0.7, colsample_bytree = 0.7, nround = 50,
  #                  data = dtrain, objective = "reg:squarederror",
  #                  watchlist = watch)
  # 
  # #feature_names <- vector(mode="character", length=hash_size)
  # #index = as.numeric(hash.mapping(full_data_Train_hashed))
  # #features = names(hash.mapping(full_data_Train_hashed))
  # #feature_names[index] = features
  # 
  # importance_matrix <- xgb.importance(model = xgb.fit)#, feature_names = feature_names)
  # xgb.plot.importance(importance_matrix, top_n = 50, main = drug,
  #                     rel_to_first = TRUE, xlab = "Relative importance")
# 
#   test_pred <- predict(xgb.fit,full_data_Test_model)
#   cor_pearson = cor(full_data_Test$auc, test_pred, method = "pearson")
#   cor_spearman = cor(full_data_Test$auc, test_pred, method = "spearman")
#   print(c(cor_pearson,cor_spearman))
#   plot5 <- qplot(full_data_Test$auc, test_pred) +
#           xlab("observed auc") +
#           ylab("predicted auc") +
#           geom_smooth(method = "lm", se = FALSE) +
#           annotate("text",label =
#                      paste("R_ p, s = ", round(cor_pearson,digits=2), ",",
#                            round(cor_spearman,digits=2)),
#                    x=Inf, y =Inf, vjust=1, hjust=1)
  
  # rr.tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(-3, 3, length = 100))
  # rf.tuneGrid = expand.grid(.mtry=c(1:floor(sqrt(ncol(full_data_Train)))))
  # model.fit <- train(auc ~ .,
  #                    data = full_data_Train,
  #                    method = "parRF", #"glmnet",
  #                    tuneGrid = rf.tuneGrid,
  #                    parms = list(split = "information"),
  #                    trControl = train_ctrl,
  #                    importance = TRUE,
  #                    metric = "RMSE")
  
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
   print(model.fit)

   model.imp <- varImp(model.fit)
   print(model.imp)
   p1 <- theme(axis.text=element_text(size=10,face="bold"),
               axis.title=element_text(size=22,face="bold"),
               plot.title = element_text(size=22))
   plot4 <- ggplot(model.imp, top = 50) +
     labs(title = drug, x='') + p1
  
   test_pred <- predict(model.fit, full_data_Test)
   cor_pearson = cor(full_data_Test$auc, test_pred, method = "pearson")
   cor_spearman = cor(full_data_Test$auc, test_pred, method = "spearman")
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