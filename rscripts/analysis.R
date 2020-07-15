library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(biomaRt)
library(ggplot2)
library(DESeq2)
library(MLSeq)
library(glmnet)
library(pROC)


colors <- c("red","blue")   # set plot colors to be used
SEED <- 1                   # set seed for all randomization

setwd('~/code/earli-study')
source("./rscripts/host_analysis.R")

reports <- read_host_genecounts('./data/earli_paxgene_rna-seq_484/')  # load in host gene count data from IDseq/STAR genecounts
reports_pc <- filter_protein_coding(reports)                          # filter the dataframe to just protein-coding genes

filtered_genecounts <- apply_genecount_qc_filters(reports_pc, 500000, 5000) # require > 500k genecounts / sample across at least 5,000 genes.

#plot_numeric_bar(colSums( filtered_genecounts ))                # plot the sum of protein-coding genes
#plot_numeric_bar(colSums( filtered_genecounts > 0 ))            # plot the number of protein-coding genes with count > 0

# import metadata, filter to include only high-quality samples
metadata <- read.csv("./data/EARLI_metadata_adjudication_IDseq.composite.tsv", sep='\t')
metadata %>% 
  filter(.$HOST_PAXgene_filename %in% colnames(filtered_genecounts)) ->
  metadata

filtered_genecounts %>% 
  dplyr::select(which(colnames(.) %in% metadata$HOST_PAXgene_filename)) -> filtered_genecounts

table(metadata$Group)

# assign training / test / validation groups
set.seed(SEED)
metadata <- assign_classifier_groups(metadata, "Group", c("G1_Sepsis_BldCx_pos"), c("G4_NO_Sepsis"), .6, .2, .2)
table(metadata[,c("Group","classifier")])

### Create ExpressionSet Objects for each of the datasets (train / test / validation) 
eset_train <- get_eset_object(metadata,
                              filtered_genecounts,
                              "HOST_PAXgene_filename",
                              c("Group","classifier"),
                              list(c("G1_Sepsis_BldCx_pos","G4_NO_Sepsis"), c("train")))

eset_test <- get_eset_object(metadata,
                             filtered_genecounts,
                             "HOST_PAXgene_filename",
                             c("Group","classifier"),
                             list(c("G1_Sepsis_BldCx_pos","G4_NO_Sepsis"), c("test")))

eset_valid <- get_eset_object(metadata,
                             filtered_genecounts,
                             "HOST_PAXgene_filename",
                             c("Group","classifier"),
                             list(c("G1_Sepsis_BldCx_pos","G4_NO_Sepsis"), c("valid")))

# eset_valid <- get_eset_object(metadata,
#                               filtered_genecounts,
#                               "HOST_PAXgene_filename",
#                               c("Group","classifier"),
#                               list(c("G2_Sepsis_OtherCx_pos"), c("valid")))

# Run Differential Expression (DEseq2) on TRAINING DATA
dds <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset_train), design= ~ Group)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
res %>% 
  as.data.frame() %>% 
  add_rownames() %>% 
  filter(padj < .001) -> sig_res                 # select only the significant genes
sig_res_sorted <- sig_res[order(sig_res$padj),]  # sort differential expression results by padj
dim(sig_res_sorted)

genes_of_interest <- sig_res_sorted$rowname

# Run DESeq2 normalization on the test and validation datasets (separately from the training dataset)
dds_test <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset_test), design= ~ Group)
dds_test <- estimateSizeFactors(dds_test)
normalized_counts_test <- counts(dds_test, normalized=TRUE)

dds_valid <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset_valid), design= ~ Group)
dds_valid <- estimateSizeFactors(dds_valid)
normalized_counts_valid <- counts(dds_valid, normalized=TRUE)

# IF selecting genes based on differential expression, subset data here.
#training_data <- normalized_counts[sig_res_sorted$rowname,]
#test_data <- normalized_counts_test[sig_res_sorted$rowname,]
#valid_data <- normalized_counts_valid[sig_res_sorted$rowname,]

# otherwise, just provide full normalized count matrices
training_data <- normalized_counts
test_data <- normalized_counts_test
valid_data <- normalized_counts_valid


# separate metadata for each group (training metadata, test metadata, validation metadata)
# TRAINING
metadata %>% 
  filter(classifier == "train") -> 
  metadata.train 
training_data %>%
  as.data.frame() %>% 
  dplyr::select(which(colnames(.) %in% metadata.train[,c("HOST_PAXgene_filename")])) ->
  data.train
metadata %>% 
  filter(classifier == "train") %>% 
  dplyr::select("Group") -> 
  class.train

# TESTING
metadata %>% 
  filter(classifier == "test") -> 
  metadata.test 
test_data %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(.) %in% metadata.test[,c("HOST_PAXgene_filename")])) ->
  data.test
metadata %>% 
  filter(classifier == "test") %>% 
  dplyr::select("Group") -> 
  class.test

# VALIDATION
metadata %>% 
  filter(classifier == "valid") %>% 
  filter(Group %in% c("G1_Sepsis_BldCx_pos","G4_NO_Sepsis")) -> 
  metadata.valid 
valid_data %>%
  as.data.frame() %>%
  dplyr::select(which(colnames(.) %in% metadata.valid[,c("HOST_PAXgene_filename")])) ->
  data.valid
metadata %>% 
  filter(classifier == "valid") %>%
  filter(Group %in% c("G1_Sepsis_BldCx_pos","G4_NO_Sepsis")) %>%
  dplyr::select("Group") -> 
  class.valid


# TRY Random Forest for FEATURE SELECTION (as opposed to differential expression, above)
library(randomForest)

training_data_rf <- as.data.frame(t(data.train))    # set-up the data to work with randomForst pkg
training_data_rf$Group <- droplevels(metadata.train$Group)
names(training_data_rf) <- make.names(names(training_data_rf))

set.seed(SEED) 
rf_classifier = randomForest(Group ~ ., data=training_data_rf, ntree=100, importance=TRUE)  # train RF classifier
rf_classifier
varImpPlot(rf_classifier)
rf_classifier_df <- as.data.frame(rf_classifier$importanceSD)
a <- varImp(rf_classifier)
genes_to_use <- rownames(a[a$G1_Sepsis_BldCx_pos > 0,])   # select only the genes with variable importances > 0

data.train <- data.train[genes_to_use,]   # subset datasets to use the genes identified by RF VarImp
data.test <- data.test[genes_to_use,]
dim(data.train)


# TRAIN CLASSIFIER - BASIC GLMNET
par(mar=c(1,1,1,1))

ALPHA = 1
cvfit = cv.glmnet(t(data.train[metadata.train$HOST_PAXgene_filename,]), as.integer(metadata.train$Group == "G1_Sepsis_BldCx_pos"), alpha = ALPHA)
plot(cvfit)
fit <- glmnet( t(data.train[,as.character(metadata.train$HOST_PAXgene_filename)]), as.integer(metadata.train$Group == "G1_Sepsis_BldCx_pos"), family = c("binomial"),
               alpha = ALPHA,
               lambda = cvfit$lambda.min)

print("# of genes")
print(length(fit$beta[fit$beta > 0]))
print(names(as.matrix(fit$beta)[as.matrix(fit$beta)>0,]))
print(length(fit$beta[fit$beta < 0]))
print(names(as.matrix(fit$beta)[as.matrix(fit$beta)<0,]))

print("TRAIN Performance")
p_train <- predict(fit, newx=t(data.train[,as.character(metadata.train$HOST_PAXgene_filename)]))
plot(p_train, col = metadata.train[metadata.train$HOST_PAXgene_filename %in% colnames(data.train),c("Group")])
pROC::roc(as.numeric(metadata.train[metadata.train$HOST_PAXgene_filename %in% colnames(data.train),c("Group")] == "G1_Sepsis_BldCx_pos"),
          p_train,
          ci = TRUE)
cm = table(class.train$Group == "G1_Sepsis_BldCx_pos", p_train > 0)
cm
par(mar = c(2, 5, 2, 2) + 0.2)
barplot(t(p_train)[1,], col = c("red","blue")[as.integer(metadata.train[metadata.train$HOST_PAXgene_filename %in% colnames(training_data),c("Group")] == "G4_NO_Sepsis") + 1], 
        horiz=TRUE,las=2, cex.names=.8, main = cat("seed: ", SEED, " auc: ", performance$auc ))

print("TEST Performance")
p <- predict(fit, newx=t(data.test[,as.character(metadata.test$HOST_PAXgene_filename)]))
plot(p, col = metadata.test[metadata.test$HOST_PAXgene_filename %in% colnames(data.test),c("Group")])
performance <- pROC::roc(as.numeric(metadata.test[metadata.test$HOST_PAXgene_filename %in% colnames(data.test),c("Group")] == "G1_Sepsis_BldCx_pos"),
          p,
          ci = TRUE)
cm = table(class.test$Group == "G1_Sepsis_BldCx_pos", p > 0)
cm
par(mar = c(2, 5, 2, 2) + 0.2)
barplot(t(p)[1,], col = c("red","blue")[as.integer(metadata.test[metadata.test$HOST_PAXgene_filename %in% colnames(test_data),c("Group")] == "G4_NO_Sepsis") + 1], 
            horiz=TRUE,las=2, cex.names=.8, main = cat("seed: ", SEED, " auc: ", performance$auc ))

# print("VALIDATION Performance")
# p_valid <- predict(fit, newx=t(data.valid[,as.character(metadata.valid$HOST_PAXgene_filename)]))
# plot(p_valid, col = metadata.valid[metadata.valid$HOST_PAXgene_filename %in% colnames(data.valid),c("Group")])
# performance <- pROC::roc(as.numeric(metadata.valid[metadata.valid$HOST_PAXgene_filename %in% colnames(data.valid),c("Group")] == "G1_Sepsis_BldCx_pos"),
#           p_valid,
#           ci = TRUE)
# cm = table(class.valid$Group == "G1_Sepsis_BldCx_pos", p_valid > 0)
# cm
# barplot(t(p_valid)[1,], col = c("red","blue")[as.integer(metadata.valid[metadata.valid$HOST_PAXgene_filename %in% colnames(valid_data),c("Group")] == "G4_NO_Sepsis") + 1], 
#         horiz=TRUE,las=2, cex.names=.8, main = cat("seed: ", SEED, " auc: ", performance$auc ))
# 



####
##### SUPPLEMENTARY ANALYSIS
####


data.trainS4 = DESeqDataSetFromMatrix(countData = data.train, colData = class.train,
                                      design = formula(~Group))
data.testS4 = DESeqDataSetFromMatrix(countData = data.test, colData = class.test,
                                     design = formula(~Group))



# BASIC PCA ANALYSIS
pca_data  <- as.data.frame(t(data.train))
pca_data <- as.data.frame(sapply( pca_data, as.numeric ))
rownames(pca_data) <- rownames(t(data.train))
pca_data <- pca_data[ , which(apply(pca_data, 2, var) != 0)]
pca <- prcomp(pca_data, center=TRUE, scale = TRUE)
pca_result <- pca$x[as.character(metadata.train$HOST_PAXgene_filename),]
plot(pca_result[,1], pca_result[,2], col = colors[as.integer(metadata.train$Group != "G1_Sepsis_BldCx_pos") + 1] )

# project the test data onto the PC components from the training data
s.sc <- scale(as.data.frame(t(data.test)), center= pca$center)
s.pred <- s.sc %*% pca$rotation
points(s.pred[,1], s.pred[,2], col = colors[as.integer(metadata.test$Group != "G1_Sepsis_BldCx_pos") + 1], pch=16)
plot(s.pred[,1], s.pred[,2], col = colors[as.integer(metadata.test$Group != "G1_Sepsis_BldCx_pos") + 1], pch=16)




####
##### MICROBE ANALYSIS
####

source("/Users/katrina.kalantar/Documents/JIRA_TASKS/decontam2/idseqr/R/idseqr.R")
library(reshape2)
library(pheatmap)
library(RColorBrewer)

microbe_reports <- read_microbe_reports('./data/microbe/earliplasma_dna_novaseq_01-2020_785/', tax_level=1, min_nt_alignment_length=70, min_nr_alignment_length=0)
overview <- read.csv("./data/microbe/earlieplasma_dna_novaseq_01-2020_785_overview.csv")

#total mass / ercc mass = total reads / ercc reads
total_input_mass <- (overview$total_reads / overview$total_ercc_reads) * 25
initial_input_mass <- total_input_mass - 25
overview$initial_input_mass_ng <- initial_input_mass/1000
sub_metadata <- metadata[,c("Group","MICROBE_Plasma_DNA.Seq_filename")]
merged <- merge(sub_metadata, overview, by.x=c("MICROBE_Plasma_DNA.Seq_filename"), by.y=c("sample_name"))
# stripchart( list(merged[merged$Group=="G1_Sepsis_BldCx_pos", c("initial_input_mass_ng")], 
#                  merged[merged$Group=="G2_Sepsis_OtherCx_pos", c("initial_input_mass_ng")],
#                  merged[merged$Group=="G3_Sepsis_Cx_neg", c("initial_input_mass_ng")],
#                  merged[merged$Group=="G4_NO_Sepsis", c("initial_input_mass_ng")]), method='jitter', jitter=.1)


water_controls <- as.character(overview[overview$sample_type %in% c("H2Oextraction", "H2Olibcontrol"),c("sample_name")])
microbe_reports_bgc <- filter_background(microbe_reports, water_controls)
sig_microbe_reports_bgc <- microbe_reports_bgc[microbe_reports_bgc$p_val < .01,]

top_taxa <- filter_top_taxa(microbe_reports_bgc, top_tax_per_sample=5)
sig_top_taxa <- top_taxa[top_taxa$p_val < .05,]
#sig_top_taxa <- top_taxa

head(top_taxa[,c("sample_name","name","nt_count","p_val")], 10)
head(sig_top_taxa[,c("sample_name","name","nt_count","p_val")], 10)

mat <- acast(sig_top_taxa, name~sample_name, value.var="nt_rpm")
print(dim(mat))

# SINGLE SAMPLE RBM
par(mar = c(2, 20, 2, 2) + 0.2)
mat_pos <- mat[,as.character(metadata[metadata$Group == "G1_Sepsis_BldCx_pos",c("MICROBE_Plasma_DNA.Seq_filename")])]
id <- 41
print(colnames(mat_pos)[id])
mat_pos <- log10(mat_pos + 1)
values <- sort(mat_pos[,id][!is.na(mat_pos[,id])], decreasing = TRUE)

get_difs <- function(l){
  difs <- list()
  for(i in seq_along(l[1:length(l)-1])){
    difs[[i]] <- l[i] - l[i+1]
  }
  return(unlist(difs))
}
x <- get_difs(values)
top_orgs <- values[1:which.max(x)]
barplot(values, las=2, horiz=TRUE, main = colnames(mat_pos)[id], col= c(rep("red",length(top_orgs)), rep("blue", length(values)-length(top_orgs))))

apply_rbm <- function(df, id){
  log_values <- log10(df + 1)
  values <- sort(log_values[,id][!is.na(log_values[,id])], decreasing = TRUE)
  x <- get_difs(values)
  top_orgs <- values[1:which.max(x)]
  pdf(paste("./outputs/RBM/RBM_",colnames(log_values)[id],".pdf", sep=""), height = 6, width = 6)
  par(mar = c(2, 10, 2, 2) + 0.2)
  barplot(values, las=2, horiz=TRUE, main = colnames(log_values)[id], col= c(rep("red",length(top_orgs)), rep("blue", length(values)-length(top_orgs))), cex.names = .6)
  dev.off()
}

for(i in seq_along(colnames(mat))){
  apply_rbm(mat, i)  
}



# FULL HEATMAP VIEW OF DATA
mat[is.na(mat)] <- 0
logmat <- log10(mat + 1)
pheatmap(logmat)

pos_cases <- logmat[,as.character(metadata[metadata$Group == "G1_Sepsis_BldCx_pos",c("MICROBE_Plasma_DNA.Seq_filename")])]
dim(pos_cases)
pdf("./outputs/top_heatmap.pdf", height = 12, width = 12)
pheatmap(pos_cases[rowSums(pos_cases) > 1,], col = brewer.pal(n = 9, name = "Reds"))
dev.off()




# 
# 
# # methods for classification (using the ML package)
# availableMethods()
# 
# fit.svm <- classify(data = data.trainS4, method = "svmRadial",
#                     preProcessing = "deseq-vst", ref = "G1_Sepsis_BldCx_pos", tuneLength = 10,
#                     control = trainControl(method = "repeatedcv", number = 5,
#                                            repeats = 10, classProbs = TRUE))
# show(fit.svm)
# trained(fit.svm)
# 
# fit <- classify(data = data.trainS4, method = "glmnet", 
#                 preProcessing = "deseq-vst", ref = "G1_Sepsis_BldCx_pos",
#                 control = trainControl(method = "repeatedcv", number = 5,
#                                        repeats = 10, classProbs = TRUE))
# show(fit)
# 
# 
# fit <- classify(data = data.trainS4, method = "voomNSC", classProbs = TRUE,
#                 normalize = "deseq", ref = "G1_Sepsis_BldCx_pos",
#                 control = voomControl(tuneLength = 20))
# trained(fit)
# plot(fit)
# 
# pred <- predict(fit, data.testS4, type = "prob")
# 
# 
# pred <- relevel(pred, ref = "G1_Sepsis_BldCx_pos", type = "prob")
# actual <- relevel(class.test$Group, ref = "G1_Sepsis_BldCx_pos")
# 
# tbl <- table(Predicted = pred, Actual = actual)
# confusionMatrix(tbl, positive = "G1_Sepsis_BldCx_pos")
# 
