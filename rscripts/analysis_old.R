library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(biomaRt)
library(ggplot2)
library(DESeq2)
library(MLSeq)

colors <- c("red","blue")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MLSeq")

setwd('~/code/earli-study')
source("./rscripts/host_analysis.R")

#reports <- read_host_genecounts('./data/test/rna-seq-data/')    # load in host gene count data from IDseq/STAR genecounts
reports <- read_host_genecounts('./data/earli_paxgene_rna-seq_484/')
reports_pc <- filter_protein_coding(reports)                    # filter the dataframe to just protein-coding genes

dim(reports_pc)
barplot(colSums(reports_pc))


# require > 1 million genecounts / sample across at least 10,000 genes.
filtered_genecounts <- apply_genecount_qc_filters(reports_pc, 1000000, 10000)

dim(filtered_genecounts)

plot_numeric_bar(colSums( filtered_genecounts ))                # plot the sum of protein-coding genes
plot_numeric_bar(colSums( filtered_genecounts > 0 ))            # plot the number of protein-coding genes with count > 0

# scatterplot of genecountss x total expressed
plot(colSums( reports_pc ), colSums( reports_pc > 0 ))
points(colSums( filtered_genecounts ), colSums( filtered_genecounts > 0 ), pch=16, col='turquoise')


# import metadata nd filter to include only high-quality samples
metadata <- read.csv("./data/EARLI_metadata_adjudication_IDseq.composite.tsv", sep='\t')
metadata %>% 
  filter(.$HOST_PAXgene_filename %in% colnames(filtered_genecounts)) ->
  metadata

filtered_genecounts %>% 
  dplyr::select(which(colnames(.) %in% metadata$HOST_PAXgene_filename)) -> filtered_genecounts

# assign training / test / validation groups
set.seed(10)
metadata <- assign_classifier_groups(metadata, "Group", c("G1_Sepsis_BldCx_pos"), c("G4_NO_Sepsis"), .5, .25, .25)
table(metadata[,c("Group","classifier")])

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


# run differential expression of TRAINING DATA

dds <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset_train), design= ~ Group)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)#, name="Group_4_NO_Sepsis_vs_1_ssSepsis.BldCx.")
res %>% as.data.frame() %>% add_rownames() %>% filter(padj < .001) -> sig_res
sig_res_sorted <- sig_res[order(sig_res$padj),]

sig_res_sorted_pos <- sig_res_sorted[sig_res_sorted$log2FoldChange > 0,]
sig_res_sorted_neg <- sig_res_sorted[sig_res_sorted$log2FoldChange < 0,]

dim(sig_res_sorted)
head(sig_res_sorted)

plotCounts(dds, gene="LAIR1", intgroup="Group")

normalized_counts[head(sig_res_sorted$rowname, 50),]
library(pheatmap)
vst <- vst(dds)
df <- as.data.frame(colData(vst)[,c("Group")])
rownames(df) <- vst$HOST_PAXgene_filename 
pheatmap(assay(vst)[head(sig_res_sorted$rowname, 250),], cluster_rows=TRUE, annotation_col=df)


#pheatmap(assay(vst)[c(as.character(sig_res_sorted[sig_res_sorted$log2FoldChange > 1,"rowname"]),
#                      as.character(sig_res_sorted[sig_res_sorted$log2FoldChange < -1,"rowname"]))], annotation_col = df)


# run differential expression of TEST DATA
dds_test <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset_test), design= ~ Group)
dds_test <- estimateSizeFactors(dds_test)
normalized_counts_test <- counts(dds_test, normalized=TRUE)
dds_test <- DESeq(dds_test)
resultsNames(dds_test) # lists the coefficients
res_test <- results(dds_test)#, name="Group_4_NO_Sepsis_vs_1_ssSepsis.BldCx.")
res_test %>% as.data.frame() %>% add_rownames() %>% filter(padj < .001) -> sig_res_test
sig_res_test_sorted <- sig_res_test[order(sig_res_test$padj),]


dim(sig_res_test_sorted)
head(sig_res_test_sorted)

plotCounts(dds_test, gene="HLA-DPB1", intgroup="Group")

library(pheatmap)
vst_test <- vst(dds_test)
df_test <- as.data.frame(colData(vst_test)[,c("Group")])
rownames(df_test) <- vst_test$HOST_PAXgene_filename 
pheatmap(assay(vst_test)[head(sig_res_test_sorted$rowname, 50),], cluster_rows=TRUE, annotation_col=df_test)


sig_res_test_sorted$rowname
table(rownames(sig_res_test_sorted) %in% rownames(sig_res_sorted))


plotCounts(dds, gene="ANXA4", intgroup="Group")
plotCounts(dds_test, gene="ANXA4", intgroup="Group")

plotCounts(dds, gene="ARL6IP5", intgroup="Group")
plotCounts(dds_test, gene="ARL6IP5", intgroup="Group")

plotCounts(dds, gene="CTSC", intgroup="Group")
plotCounts(dds_test, gene="CTSC", intgroup="Group")

plotCounts(dds, gene="HRH4", intgroup="Group")
plotCounts(dds_test, gene="HRH4", intgroup="Group")


tr_metric <- colSums(normalized_counts[sig_res_sorted_neg$rowname,])
tr_metric_g1 <-  tr_metric[as.character(eset_train$HOST_PAXgene_filename[eset_train$Group == "G1_Sepsis_BldCx_pos"])]
tr_metric_g4 <-  tr_metric[as.character(eset_train$HOST_PAXgene_filename[eset_train$Group == "G4_NO_Sepsis"])]

te_metric <- colSums(normalized_counts_test[sig_res_sorted_neg$rowname,])
te_metric_g1 <-  te_metric[as.character(eset_test$HOST_PAXgene_filename[eset_test$Group == "G1_Sepsis_BldCx_pos"])]
te_metric_g4 <-  te_metric[as.character(eset_test$HOST_PAXgene_filename[eset_test$Group == "G4_NO_Sepsis"])]

stripchart(list(tr_metric_g1, tr_metric_g4, te_metric_g1, te_metric_g4), method = 'jitter', jitter = .1, vertical = TRUE)
