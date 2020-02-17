library(readr)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(caret)
library(pls)
library(xgboost)

input_dir = "" #"/input/"
output_dir = "" #"/output/"
rds_dir = "" #"/rds/"

out_file=paste0(output_dir,"predictions.csv")
dnaseq_t <- read_csv(paste0(input_dir,"dnaseq.csv"))
rnaseq_t <- read_csv(paste0(input_dir,"rnaseq.csv"))
clinical_categorical_t <- read_csv(paste0(input_dir,"clinical_categorical.csv"))
clinical_numerical_t <- read_csv(paste0(input_dir,"clinical_numerical.csv"))
clinical_categorical_legend_t <- read_csv(paste0(input_dir,"clinical_categorical_legend.csv"))
colnames(clinical_categorical_t) <- make.names(colnames(clinical_categorical_t))

# drop columns with missing values
clinical_numerical_t <- clinical_numerical_t %>% select(-c("%.Blasts.in.PB","WBC.Count"))

# prepare rnaseq data
rna_log2counts_t <- as.matrix(rnaseq_t[,3:dim(rnaseq_t)[2]])
rownames(rna_log2counts_t) <- rnaseq_t$Gene
rna_counts_t <- round(2^rna_log2counts_t)
rownames(rna_counts_t) <- rownames(rna_log2counts_t)

# rnaseq cpm sum of total check
print(colSums(2^rna_log2counts_t))

# load lists of model fits and gene features
model.fits <- readRDS(file = paste0(rds_dir,"model-fits.rds"))
dnaseqGenes <- readRDS(file = paste0(rds_dir,"dnaseqGenes.rds"))
rnaseqGenes <- readRDS(file = paste0(rds_dir,"rnaseqGenes.rds"))
rna.preProc <- readRDS(file = paste0(rds_dir,"rna-preProc.rds"))

output_df <- c()

for (drug_t in names(model.fits)){
  
  my_model <- model.fits[[drug_t]]
  topGenes_dnaseq_t <- dnaseqGenes[[drug_t]]
  topGenes_rnaseq_t <- rnaseqGenes[[drug_t]]
  sigDE_preProc_t <- rna.preProc[[drug_t]]
  
  # PCA of differentially expressed genes
  sigDE_log2counts_t <- rna_log2counts_t[
    which(rownames(rna_log2counts_t) %in% topGenes_rnaseq_t),]
  #pcDat_t <- prcomp(t(sigDE_log2counts_t))
  #PCA_sigDE_t <- predict(sigDE_preProc_t, t(sigDE_log2counts_t))
  PLS_sigDE_t <- predict(sigDE_preProc_t, type = "scores", newdata = t(sigDE_log2counts_t))
  
  full_data_t <- dcast(
    dnaseq_t, lab_id ~ Hugo_Symbol, fun.aggregate =
      function(x)(as.numeric(!isEmpty(x))) ) %>% select(-c("NPM1"))
   
  missingGenos = setdiff(colnames(rna_log2counts_t),full_data_t$lab_id)
  mtx <- cbind(missingGenos,matrix(0,length(missingGenos),dim(full_data_t)[2]-1))
  colnames(mtx) <- c("lab_id",colnames(full_data_t)[-1])
  full_data_t = rbind(full_data_t, mtx)
  
  full_data_t <- full_data_t %>%
    inner_join(
      select(clinical_categorical_t,c("lab_id","NPM1","FLT3.ITD")),
      by = "lab_id") %>%
    select(one_of(c("lab_id",topGenes_dnaseq_t))) %>%
    #mutate_if(is.integer, as.factor) %>%
    #inner_join(
    #  as.data.frame(t(sigDE_log2counts_t)) %>%
    #    mutate(lab_id=rownames(t(sigDE_log2counts_t))), by = "lab_id") %>%
      #as.data.frame(pcDat_t[["x"]]) %>% select(1:min(5,dim(pcDat_t[["x"]])[2])) %>%
      #  mutate(lab_id=rownames(pcDat_t[["x"]])), by = "lab_id") %>%
    inner_join(
      select(clinical_categorical_t,-c("NPM1","FLT3.ITD")) %>%
        mutate_at(vars(-one_of("lab_id")),as.numeric), by = "lab_id") %>%
    inner_join(clinical_numerical_t, by = "lab_id")
  
  full_data_t <- full_data_t %>%
    inner_join( as.data.frame(PLS_sigDE_t) %>% 
                  #select(1:min(10,dim(PCA_sigDE_t)[2])) %>%
                  mutate(lab_id=rownames(PLS_sigDE_t)), by = "lab_id") 
  
  full_data_t[setdiff(topGenes_dnaseq_t,names(full_data_t))] <- 0
  full_data_t[setdiff(rownames(varImp(my_model)$importance),names(full_data_t))] <- 0
  #full_data_t[setdiff(my_model$feature_names,names(full_data_t))] <- 0
  
  full_data_t <- full_data_t %>% 
    mutate_at(vars(-one_of("lab_id")),funs(as.numeric(as.character(.))))
  rownames(full_data_t)=full_data_t$lab_id
  full_data_t <- full_data_t %>% select(-c("lab_id"))
  
  # mean impute missing values as only a couple NAs
  cM <- colMeans(full_data_t, na.rm=TRUE)
  indx <- which(is.na(full_data_t), arr.ind=TRUE)
  if (!isEmpty(indx)){ full_data_t[indx] <- cM[indx[,2]]}
  
  #full_data_t <- full_data_t[, my_model$feature_names]
  #full_data_t_model <- model.matrix(~ .-1, full_data_t)
  
  pred_t <- predict(my_model, newdata = full_data_t) #newdata = full_data_t_model)
  
  aucs_pred <- data.frame(
    lab_id = rownames(full_data_t), #rownames(full_data_t_model),
    inhibitor = rep(drug_t,length(pred_t)),
    auc = pred_t, row.names=NULL)
  output_df = rbind(output_df, aucs_pred)
  
}

write_csv(output_df, out_file)
