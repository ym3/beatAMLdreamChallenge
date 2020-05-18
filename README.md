##Tackling AML multiple omics using feature selection and dimensionality reduction

**Dr Yasin Memari, Senior Bioinformatician, University of Cambridge**


**Summary**

Combining mixed numerical and categorical variables of high dimensionality and sparsity in a supervised approach is sensitive to the choice of model and variable selection. Here we show how expression data can be reduced and combined with select genetic markers and clinical information to predict drug sensitivity and response for Acute Myeloid Leukemia.

**Background/Introduction**

For SC1 we use a combination of feature selection of DNAseq mutational data and dimensionality reduction of RNAseq expression data before incorporating the clinical information. Drug sensitivity is modelled using a gradient boosted tree-based approach projecting the combined data onto quantitative ex vivo drug sensitivity (AUC) where parameters of the model are estimated using cross validation. 

For SC2 clinical response is modelled using Cox Proportional Hazard survival analysis.  First, individual covariates are selected based on their statistical significance in a univariate CPH test, then a multivariate CPH is applied to the subset of the significant covariates and used to predict the patients' survival. 


**Methods - SC1**

DNAseq matrix is sparse due to high mutational heterogeneity and exclusivity (Figure 1), therefore we first attempt to overcome the sparsity of the mutation matrix using a univariate feature selection approach.  

${image?fileName=mutations%2Epng&align=Center&scale=50&responsive=true&altText=Figure 1}

Gene symbols where drug response is significantly different between mutant and wild-type samples are selected using a Student’s two-sided t-test. Figure 2 shows p values of individual genes averaged across the drugs versus average AUC difference between the two groups. Upper p-value threshold of 0.05 was applied to select the significant genes for each drug. Key drivers  FLT3-ITD, NPM1, DNMT3A, IDH1 and IDH2 were overrepresented among the significant genes as seen in the figure. 

${image?fileName=Ttest%2Epng&align=Center&scale=55&responsive=true&altText=Figure 2}

For RNAseq, probes with non-zero counts in more that 20% of the samples were retained and mean-variance normalised using voom function in limma package. Several approaches were tested to reduce the dimensionality of the expression matrix, including PCA, PLS and WGCNA of highly variable probes or probes which are differentially expressed between top and bottom 20% drug responders/non-responders (e.g. using DESeq2 in Figure 3).  

${image?fileName=Ibrut%2Epng&align=Center&scale=50&responsive=true&altText=Figure 3}

It was observed that data reduction using supervised Partial Least Squares regression projecting the expression matrix of top variable genes (without DE genes selection) onto AUC values produced the best predictions. Therefore, PLS components of top 1000 variable genes were used in the final submission to reduce the RNAseq data where the one-sigma method of PLS was used to select number of the components.

For model fitting, we used the XGBTree model on the reduced genomic data (RNAseq + DNAseq) combined with clinical categorical and numerical data to predict AUCs. A tree based approach is particularly appropriate as it can handle non-linearity and mixed categorical and numerical features. Fine tuning of the model parameters was performed using 5-fold cross validation in caret package after partitioning the data into 70% training and 30% test sets.  Figure 4 shows model variable importance ranking for an example drug where RNAseq (PLS V1), ageAtDiagnosis, consensus_sex, priorMDS/MPN, FLT3.ITD, NPM1, IDH2 stand out on top in terms of their predictive importance.

${image?fileName=VarImp%2Epng&align=Center&scale=50&responsive=true&altText=Figure 4}

Throughout the analysis missing values of numerical variables in either training or test/validation sets were mean imputed and categorical variables comprising only one level (after partitioning) were dropped from the analysis. After evaluating the model performance using data partitioning, the final submission uses the entire data for training. 

No external data was incorporated in the analysis.


**Methods - SC2**

For predicting drug response (survival rates), we adopt a two-step Cox Proportional Hazard regression approach. Multiple regression is sensitive to the choice of the model covariates, therefore, to increase the statistical power we first select the statistically significant covariates (at 5% p value) using a univariate CPH model, then apply a multivariate CPH to the subset of the significant covariates. In step one, clinical variables age, prior malignancy status as well as RNAseq PCA stand out while DNAseq mutations do not seem to explain much of the survival variability. In step two, the selected covariates are fitted to the data using a multivariate CPH model (Figure 5) which is used to predict the clinical response.

${image?fileName=CoxPH%2Epng&align=Center&scale=40&responsive=true&altText=Figure 5}


**Discussion**

RNAseq, among other data, was found to be the most predictive of drug sensitivity (SC1). In fact, regardless how the expression data was transformed, RNAseq often dominated the variable importance table. XGBTree was preferable to Random Forest which was found to downweigh the categorical variables (DNAseq and clinical categorical) importance. In SC2, age at diagnosis and specimen collection consistently proved to be the most informative as to the drug response.

**References**

Project page: https://www.synapse.org/#!Synapse:syn21561349
Source code: https://github.com/ym3/beatAMLdreamChallenge


**Authors Statement**

Yasin Memari designed and performed all the analysis on behalf of team Signal.
