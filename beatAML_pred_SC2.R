library(readr)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(limma)
library(survival)

input_dir = "/input/"
output_dir = "/output/"
rds_dir = "/rds/"

out_file=paste0(output_dir,"predictions_response.csv")
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

# variance stabilizing transformation
#vstcounts_t <- vst(rna_counts_t, blind = FALSE)
#rownames(vstcounts_t) <- rownames(rna_log2counts_t)
#rna_log2counts_t <- vstcounts_t

# load the model fit
my_coxph <- readRDS(file = paste0(rds_dir,"coxph-fit.rds"))

# prepare input for predict
countdata_t <- as.data.frame(rna_counts_t)
rnaseq_voom_t = voom(countdata_t)$E
topGenes_rnaseq_t = data.frame(Gene = rownames(rnaseq_voom_t[
  order(apply(rnaseq_voom_t,1,mad), decreasing = T)[1:500],]))
sigDE_log2counts_t <- rna_log2counts_t[which(rownames(rna_log2counts_t) 
                                         %in% topGenes_rnaseq_t$Gene),]
rnaseq_pcDat_t <- prcomp(t(sigDE_log2counts_t))

response_data_t <- clinical_categorical_t %>% 
  inner_join(clinical_numerical_t, by = "lab_id") %>% 
  inner_join(as.data.frame(rnaseq_pcDat_t[["x"]]) %>% 
               select(1:min(5,dim(rnaseq_pcDat_t[["x"]])[2])) %>% 
               mutate(lab_id=rownames(rnaseq_pcDat_t[["x"]])), by = "lab_id")
rownames(response_data_t)=response_data_t$lab_id
response_data_t <- response_data_t %>% select(-c("lab_id"))

#coxph.pred <- survfit(my_coxph, newdata = response_data_t)
coxph.pred <- predict(my_coxph, response_data_t)

output_df <- data.frame(
  lab_id = names(coxph.pred), #names(summary(coxph.pred)$table[,"median"]),
  survival = -1 * coxph.pred, #summary(coxph.pred)$table[,"median"], 
  row.names=NULL)
output_df$survival[which(is.na(output_df$survival))]=0

write_csv(output_df, out_file)
