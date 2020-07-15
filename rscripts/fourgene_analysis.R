library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(biomaRt)
library(ggplot2)
library(DESeq2)
library(MLSeq)
library(scales)

setwd('~/code/earli-study')
source("./rscripts/host_analysis.R")

reports <- read_host_genecounts('./data/earli_paxgene_rna-seq_484/')
reports_pc <- filter_protein_coding(reports)                    # filter the dataframe to just protein-coding genes

dim(reports_pc)
barplot(colSums(reports_pc))

# require > 1 million genecounts / sample across at least 10,000 genes.
filtered_genecounts <- apply_genecount_qc_filters(reports_pc, 1000000, 10000)

dim(filtered_genecounts)


filtered_genecounts[c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),]

# import metadata nd filter to include only high-quality samples
#metadata <- read.csv("./data/EARLI_metadata_adjudication_IDseq.composite.2.tsv", sep='\t')
metadata <- read.csv("./data/metadata_for_included_samples_ALrevised.tsv", sep='\t')
metadata %>% 
  filter(.$HOST_PAXgene_filename %in% colnames(filtered_genecounts)) ->
  metadata

### FILTERING THE RESULTS TO SHOW ONLY THE THREE GROUPS OF INTEREST
metadata %>% 
  filter(.$Group %in% c("1_Sepsis_BldCx_pos","G2_Sepsis_OtherCx_pos","3_Sepsis+Cx-")) ->
  metadata

filtered_genecounts %>% 
  dplyr::select(which(colnames(.) %in% metadata$HOST_PAXgene_filename)) -> filtered_genecounts


# filter out the overlapping samples which were used to derive the four-gene signature
metadata %>% 
  filter(!(.$EARLIStudyId %in% c(1253,1260,1267,1270,1284,1339))) ->
  metadata

metadata$Hospital.death <- as.factor(metadata$Hospital.death)

filtered_genecounts %>% 
  dplyr::select(which(colnames(.) %in% metadata$HOST_PAXgene_filename)) -> filtered_genecounts


sorted_filtered_genecounts <- filtered_genecounts[,as.character(metadata$HOST_PAXgene_filename)]
# sorted_filtered_genecounts[c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),]

rownames(metadata) <- metadata$HOST_PAXgene_filename
eset <- ExpressionSet(assayData = as.matrix(sorted_filtered_genecounts),
                      phenoData = AnnotatedDataFrame(metadata))

dds <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset), design= ~ Group)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
# dds <- DESeq(dds)

#par(mfrow=c(1,4))
plot(metadata$APACHEIII, log10(normalized_counts[c("HLA-DRA"),]), xlab = "APACHE-III", ylab = "log10(HLA-DRA Expression)")
plot(metadata$APACHEIII, log10(normalized_counts[c("HLA-DPA1"),]), xlab = "APACHE-III", ylab = "log10(HLA-DPA1 Expression)")
plot(metadata$APACHEIII, log10(normalized_counts[c("HLA-DPB1"),]), xlab = "APACHE-III", ylab = "log10(HLA-DPB1 Expression)")
plot(metadata$APACHEIII, log10(normalized_counts[c("FUCA1"),]), xlab = "APACHE-III", ylab = "log10(FUCA1 Expression)")

l_dra = lm(as.numeric(normalized_counts[c("HLA-DRA"),]) ~ metadata$APACHEIII)
summary(l_dra)

l_dpa = lm(as.numeric(normalized_counts[c("HLA-DPA1"),]) ~ metadata$APACHEIII)
summary(l_dpa)

l_dpb = lm(as.numeric(normalized_counts[c("HLA-DPB1"),]) ~ metadata$APACHEIII)
summary(l_dpb)

l_fuca = lm(as.numeric(normalized_counts[c("FUCA1"),]) ~ metadata$APACHEIII)
summary(l_fuca)


#plot(metadata$Hospital.death, normalized_counts[c("HLA-DRA"),])
#l = lm(as.numeric(normalized_counts[c("HLA-DRA"),]) ~ as.factor(metadata$Hospital.death))
#summary(l)

table(metadata$Hospital.death == 1)

hd0 <- normalized_counts[,rownames(metadata[metadata$Hospital.death == 0,])]
hd1 <- normalized_counts[,rownames(metadata[metadata$Hospital.death == 1,])]

#hd0 <- assay(vst)[c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),rownames(metadata[metadata$Hospital.death == 0,])]
#hd1 <- assay(vst)[c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),rownames(metadata[metadata$Hospital.death == 1,])]

remove_outliers <- function(data){
  Q0 <- quantile(data, probs=c(.25, .75), na.rm = FALSE)
  print(Q0)
  iqr0 <- IQR(data)
  print(iqr0)
  up0 <-  Q0[2]+1.5*iqr0 # Upper Range  
  print(up0)
  low0 <- Q0[1]-1.5*iqr0 # Lower Range
  print(low0)
  print(length(data[data < low0]))
  print(length(data[data > up0]))
  return(data[data < up0])
}

# # WILCOX HLA-DRA
# boxplot(list(hd0[c("HLA-DRA"),],hd1[c("HLA-DRA"),]))
# wilcox.test(hd0[c("HLA-DRA"),],hd1[c("HLA-DRA"),])
# 
# boxplot(list(remove_outliers(hd0[c("HLA-DRA"),]),remove_outliers(hd1[c("HLA-DRA"),])))
# wilcox.test(remove_outliers(hd0[c("HLA-DRA"),]),remove_outliers(hd1[c("HLA-DRA"),]))
# 
# # WILCOX HLA-DPA1
# boxplot(list(hd0[c("HLA-DPA1"),],hd1[c("HLA-DPA1"),]))
# wilcox.test(hd0[c("HLA-DPA1"),],hd1[c("HLA-DPA1"),])
# 
# boxplot(remove_outliers(hd0[c("HLA-DPA1"),]), remove_outliers(hd1[c("HLA-DPA1"),]))
# wilcox.test(remove_outliers(hd0[c("HLA-DPA1"),]), remove_outliers(hd1[c("HLA-DPA1"),]))
# 
# # WILCOX HLA-DPB1
# 
# boxplot(list(hd0[c("HLA-DPB1"),],hd1[c("HLA-DPB1"),]))
# wilcox.test(hd0[c("HLA-DPB1"),],hd1[c("HLA-DPB1"),])
# 
# boxplot(list(remove_outliers(hd0[c("HLA-DPB1"),]), remove_outliers(hd1[c("HLA-DPB1"),])))
# vioplot(list(remove_outliers(hd0[c("HLA-DPB1"),]), remove_outliers(hd1[c("HLA-DPB1"),])))
# stripchart(list(remove_outliers(hd0[c("HLA-DPB1"),]), remove_outliers(hd1[c("HLA-DPB1"),])), jitter = .1, vertical = TRUE, add=TRUE, method='jitter')
# wilcox.test(remove_outliers(hd0[c("HLA-DPB1"),]), remove_outliers(hd1[c("HLA-DPB1"),]))
# 
# 
# # WILCOX FUCA1
# boxplot(list(hd0[c("FUCA1"),],hd1[c("FUCA1"),]))
# wilcox.test(hd0[c("FUCA1"),],hd1[c("FUCA1"),])
# 
# boxplot(list(remove_outliers(hd0[c("FUCA1"),]), remove_outliers(hd1[c("FUCA1"),])))
# wilcox.test(remove_outliers(hd0[c("FUCA1"),]), remove_outliers(hd1[c("FUCA1"),]))


dds <- DESeqDataSet(makeSummarizedExperimentFromExpressionSet(eset), design= ~ Hospital.death)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
res %>% as.data.frame() %>% add_rownames() %>% filter(padj < .05) -> sig_res

library(pheatmap)
vst <- vst(dds)
df <- as.data.frame(colData(vst)[,c("Hospital.death")])
rownames(df) <- colnames(vst)
pheatmap(assay(vst)[c("HLA-DPA1","HLA-DPB1","FUCA1"),], cluster_rows=TRUE, annotation_col=df)


dim(sig_res_sorted)
sig_res_sorted <- sig_res[order(sig_res$padj),]
sig_res_sorted[sig_res_sorted$rowname %in% c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),]

res[rownames(res) %in% c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),]

dev.off()
par(mar=c(1,1,1,1), mfrow=c(1, 4))
plotCounts(dds, gene="HLA-DRA", intgroup="Hospital.death", pch = 16, col = alpha('black',  .4))
plotCounts(dds, gene="HLA-DPA1", intgroup="Hospital.death", pch = 16, col = alpha('black',  .4))
plotCounts(dds, gene="HLA-DPB1", intgroup="Hospital.death", pch = 16, col = alpha('black',  .4))
plotCounts(dds, gene="FUCA1", intgroup="Hospital.death", pch = 16, col = alpha('black',  .4))

pd <- as.data.frame(cbind(normalized_counts[c("HLA-DRA"),], metadata$Hospital.death == 1))
colnames(pd) <- c("HLA.DRA", "Hospital.death")
p1 <- ggplot(pd, aes(x=factor(Hospital.death), y=HLA.DRA, color=factor(Hospital.death))) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha = .4)+
  scale_color_brewer(palette="Set1") + labs(x = "Hospital Outcome", y = "HLA-DRA") + theme(legend.position="none") + 
  scale_x_discrete(labels=c("0" = "Survival", "1" = "Death"))

pd <- as.data.frame(cbind(normalized_counts[c("HLA-DPA1"),], metadata$Hospital.death == 1))
colnames(pd) <- c("HLA.DPA1", "Hospital.death")
p2 <- ggplot(pd, aes(x=factor(Hospital.death), y=HLA.DPA1, color=factor(Hospital.death))) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha = .4)+
  scale_color_brewer(palette="Set1") + labs(x = "Hospital Outcome", y = "HLA-DPA1") + theme(legend.position="none") + 
  scale_x_discrete(labels=c("0" = "Survival", "1" = "Death"))

pd <- as.data.frame(cbind(normalized_counts[c("HLA-DPB1"),], metadata$Hospital.death == 1))
colnames(pd) <- c("HLA.DPB1", "Hospital.death")
p3 <- ggplot(pd, aes(x=factor(Hospital.death), y=HLA.DPB1, color=factor(Hospital.death))) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha = .4)+
  scale_color_brewer(palette="Set1") + labs(x = "Hospital Outcome", y = "HLA-DPB1") + theme(legend.position="none") + 
  scale_x_discrete(labels=c("0" = "Survival", "1" = "Death"))

pd <- as.data.frame(cbind(normalized_counts[c("FUCA1"),], metadata$Hospital.death == 1))
colnames(pd) <- c("FUCA1", "Hospital.death")
p4 <- ggplot(pd, aes(x=factor(Hospital.death), y=FUCA1, color=factor(Hospital.death))) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha = .4)+
  scale_color_brewer(palette="Set1") + labs(x = "Hospital Outcome", y = "FUCA1") + theme(legend.position="none") + 
  scale_x_discrete(labels=c("0" = "Survival", "1" = "Death"))

library(gridExtra)
g <- grid.arrange(p3, p4, p2, p1, nrow=1)
ggsave(
  "~/Documents/sepsis_study/geneplots.pdf",
  g,
  width = 8,
  height = 3
)

combined_metric <- colSums(normalized_counts[c("HLA-DRA","HLA-DPA1","HLA-DPB1","FUCA1"),])
c_hd0 <- combined_metric[rownames(metadata[metadata$Hospital.death == 0,])]
c_hd1 <- combined_metric[rownames(metadata[metadata$Hospital.death == 1,])]

boxplot(c_hd0, c_hd1, ylab=("sum of four-gene expression (including outliers)"))
stripchart(list(c_hd0, c_hd1), jitter = .1, vertical = TRUE, add=TRUE, method='jitter', col=alpha('blue',.5), pch=16)
wilcox.test(c_hd0, c_hd1)

survival = mean(log2(c_hd0))
death = mean(log2(c_hd1))
foldchange <- death - survival
foldchange


# boxplot(remove_outliers(c_hd0), remove_outliers(c_hd1), ylab=("sum of four-gene expression (without outliers)"))
# stripchart(list(remove_outliers(c_hd0), remove_outliers(c_hd1)), jitter = .1, vertical = TRUE, add=TRUE, method='jitter', col=alpha('blue',.5), pch=16)
# wilcox.test(remove_outliers(c_hd0), remove_outliers(c_hd1))


combined_metric_data <- as.data.frame(cbind(combined_metric, metadata$Hospital.death == 1))
combined_metric_data
colnames(combined_metric_data) <- c("metric", "Hospital.death")

p <- ggplot(combined_metric_data, aes(x=factor(Hospital.death), y=metric, color=factor(Hospital.death))) + 
  geom_boxplot(notch=TRUE, outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2), alpha = .4)+
  scale_color_brewer(palette="Set1")
p + labs(x = "Hospital Outcome", y = "Metric") + theme(legend.position="none") + 
  scale_x_discrete(labels=c("0" = "Survival", "1" = "Death"))
#p # + theme(legend.position="top")
ggsave(
  "~/Documents/sepsis_study/combined_metric_plot.pdf",
  plot = last_plot(),
  width = 3,
  height = 4
)



###
### Output the Metadata
### 
write.table(metadata, "~/Documents/sepsis_study/metadata_for_included_samples.csv", sep=',')



###
### Checking on the outlier sample:
###

outlier <- names(c_hd1[which.max(c_hd1)])
metadata[metadata$HOST_PAXgene_filename == outlier,]

# how many gene counts were there? this is the n-th most gene expression?
summary(colSums(sorted_filtered_genecounts))
which(sort(colSums(sorted_filtered_genecounts)) == sum(sorted_filtered_genecounts[,outlier]))
hist(colSums(sorted_filtered_genecounts), breaks = 20, main = "Distribution of Total Genecounts")
abline(v = sum(sorted_filtered_genecounts[,outlier]), col='purple')

# how many genes were expressed?
summary(colSums(sorted_filtered_genecounts>0))
sum((sorted_filtered_genecounts>0)[,outlier])
which(sort(colSums(sorted_filtered_genecounts > 0)) == sum((sorted_filtered_genecounts>0)[,outlier]))
hist(colSums(sorted_filtered_genecounts > 0), breaks=20, main = "Distribution of Total Expressed Genes")
abline(v=sum((sorted_filtered_genecounts>0)[,outlier]), col = 'purple')

plot(colSums((sorted_filtered_genecounts>0)), colSums((sorted_filtered_genecounts)), 
     xlab = 'total expresesd genes', ylab = 'total gene counts')
points(sum((sorted_filtered_genecounts>0)[,outlier]), sum((sorted_filtered_genecounts)[,outlier]), col = 'red', pch=16)

# this doesn't seem like an outlier due to faulty gene expression