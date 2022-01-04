##### Riddle 2022 differential expression analysis in R using DESeq2
library(DESeq2)
library(dplyr)
library(tidyverse)
library(data.table)

###pairwise analysis  
CountData <- read.csv("/Users/aa9gj/Documents/Riddle_2022_analysis/gene_count_matrix.csv", row.names = 1)
ColData <- read.csv("/Users/aa9gj/Documents/Riddle_2022_analysis/phenotypes.csv", row.names = 1)
# Reorder the columns of the count data
reorder_idx <- match(rownames(ColData),colnames(CountData))
reordered_CountData <- CountData[ , reorder_idx]
all(rownames(ColData) %in% colnames(CountData))


dds <- DESeqDataSetFromMatrix(countData =  reordered_CountData,
                              colData =  ColData,
                              design = ~ Condition)
dds <- estimateSizeFactors(dds)
# Extract the normalized counts
normalized_counts <-counts(dds, normalize = TRUE)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts_mod <- setDT(normalized_counts, keep.rownames = T)[]
normalized_counts_mod <- separate(normalized_counts_mod, rn, into = c("gene_id", "gene_name"), sep = "\\|")
# Transform the normalized counts
vsd <- vst(dds, blind = TRUE)
#Extract the matrix of transformed counts
vsd_mat <- assay(vsd)
#Compute the correlation values between samples
vsd_cor <- cor(vsd_mat, vsd_mat)
#Plot the correlation heatmap
library(pheatmap)
pheatmap(vsd_cor, annotation = select(ColData, Condition))
# Plot the PCA of PC1 and PC2
plotPCA(vsd, intgroup= "Condition")
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
z <- plotPCA(DESeqTransform(se), intgroup= "Condition")
z + geom_label(aes(label = name), size = 2)

# Run the DESeq2 analysis for pairwise comparisons later on
dds <- DESeq(dds, betaPrior = FALSE)
# Plot dispersions
plotDispEsts(dds)

con_vs_treat <- results(dds, 
                        contrast = c("Condition","treatment", "control"), 
                        alpha = 0.05)
con_vs_treat_res <- subset(con_vs_treat, padj < 0.05)

summary(con_vs_treat_res) #upregulated and downregulated in treatment (PTH) vs control
con_vs_treat_res <- as.data.frame(con_vs_treat_res)
con_vs_treat_res <- setDT(con_vs_treat_res, keep.rownames = T)[]
con_vs_treat_res <- separate(con_vs_treat_res, rn, into = c("gene_id", "gene_name"), sep = "\\|")
write.csv(con_vs_treat_res, "/Users/aa9gj/Documents/Riddle_2022_analysis/diff_exp_sig_results", quote = F, row.names = F)

##create heatmaps

fatty_acid <- read.table("/Users/aa9gj/Documents/Riddle_2022_analysis/fatty_acid.txt", header = F)
glycolysis <- read.table("/Users/aa9gj/Documents/Riddle_2022_analysis/glycolysis.txt", header = F)
amino_acid <- read.table("/Users/aa9gj/Documents/Riddle_2022_analysis/amino_acid.txt", header = F)

fattyacid_norm_counts <- inner_join(fatty_acid, normalized_counts_mod, by = c("V1"="gene_name"))
glycolysis_norm_counts <- inner_join(glycolysis, normalized_counts_mod, by = c("V1"="gene_name"))
aminoacid_norm_counts <- inner_join(amino_acid, normalized_counts_mod, by = c("V1"="gene_name"))

# #emailed Ryan about these genes 
# grep("Ascl2", normalized_counts_mod$gene_name) #might have to change genes in the list by ryan to these
# grep("Pfkfb3", normalized_counts_mod$gene_name)# same
# grep("Slc36a4", normalized_counts_mod$gene_name)#same
# x <- grep("Gls", normalized_counts_mod$gene_name)

##make them nicer and send to Ryan
fattyacid_norm_counts <- fattyacid_norm_counts[,-2]
row.names(fattyacid_norm_counts) <- fattyacid_norm_counts$V1
fattyacid_norm_counts[1] <- NULL
pheatmap(fattyacid_norm_counts)
mat <- log2(fattyacid_norm_counts)
is.na(mat)<-sapply(mat, is.infinite)
mat[is.na(mat)]<-0
pheatmap(mat)

glycolysis_norm_counts <- glycolysis_norm_counts[,-2]
row.names(glycolysis_norm_counts) <- glycolysis_norm_counts$V1
glycolysis_norm_counts[1] <- NULL
pheatmap(glycolysis_norm_counts)
mat <- log2(glycolysis_norm_counts)
is.na(mat)<-sapply(mat, is.infinite)
mat[is.na(mat)]<-0
pheatmap(mat)

aminoacid_norm_counts <- aminoacid_norm_counts[,-2]
row.names(aminoacid_norm_counts) <- aminoacid_norm_counts$V1
aminoacid_norm_counts[1] <- NULL
pheatmap(aminoacid_norm_counts)
mat <- log2(aminoacid_norm_counts)
is.na(mat)<-sapply(mat, is.infinite)
mat[is.na(mat)]<-0
pheatmap(mat)
