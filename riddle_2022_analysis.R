##### Riddle 2022 differential expression analysis in R using DESeq2
library(DESeq2)
library(dplyr)
library(tidyverse)
library(data.table)

###pairwise analysis  
CountData <- read.csv("change", row.names = 1)
ColData <- read.csv("", row.names = 1)
ColData$Batch <- as.factor(ColData$Batch)
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
z + geom_label(aes(label = name))

# Run the DESeq2 analysis for pairwise comparisons later on
dds <- DESeq(dds, betaPrior = FALSE)
# Plot dispersions
plotDispEsts(dds)

D0_10_res <- results(dds, 
                     contrast = c("Condition","D0", "D10"), 
                     alpha = 0.05)
D0_10_res <- subset(D0_10_res, padj < 0.05)
D0_10_res_annot <- setDT(as.data.frame(D0_10_res), keep.rownames = TRUE)[]
D0_10_res_annot <- inner_join(D0_10_res_annot, sqanti, by = c("rn" = "isoform"))
D0_D10_res_annot <- as.data.frame(left_join(D0_10_res_annot, gencode_v38_gene, by = c("associated_gene"="gene_id")))
summary(D0_10_res)
write.csv(D0_D10_res_annot, "/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/DE_longread_analysis/DE_D0_10_annot", quote = F)

D0_2_res <- results(dds, 
                    contrast = c("Condition","D0", "D2"), 
                    alpha = 0.05)
D0_2_res <- subset(D0_2_res, padj < 0.05)
D0_2_res_annot <- setDT(as.data.frame(D0_2_res), keep.rownames = TRUE)[]
D0_2_res_annot <- inner_join(D0_2_res_annot, sqanti, by = c("rn" = "isoform"))
D0_D2_res_annot <- as.data.frame(left_join(D0_2_res_annot, gencode_v38_gene, by = c("associated_gene"="gene_id")))
summary(D0_2_res)
write.csv(D0_D2_res_annot, "/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/DE_longread_analysis/DE_D0_2_annot", quote = F)

D0_4_res <- results(dds, 
                    contrast = c("Condition","D0", "D4"), 
                    alpha = 0.05)
D0_4_res <- subset(D0_4_res, padj < 0.05)
D0_4_res_annot <- setDT(as.data.frame(D0_4_res), keep.rownames = TRUE)[]
D0_4_res_annot <- inner_join(D0_4_res_annot, sqanti, by = c("rn" = "isoform"))
D0_D4_res_annot <- as.data.frame(left_join(D0_4_res_annot, gencode_v38_gene, by = c("associated_gene"="gene_id")))
summary(D0_4_res)
write.csv(D0_D4_res_annot, "/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/DE_longread_analysis/DE_D0_4_annot", quote = F)


D2_4_res <- results(dds, 
                    contrast = c("Condition","D2", "D4"), 
                    alpha = 0.05)
D2_4_res <- subset(D2_4_res, padj < 0.05)
D2_4_res_annot <- setDT(as.data.frame(D2_4_res), keep.rownames = TRUE)[]
D2_4_res_annot <- inner_join(D2_4_res_annot, sqanti, by = c("rn" = "isoform"))
D2_D4_res_annot <- as.data.frame(left_join(D2_4_res_annot, gencode_v38_gene, by = c("associated_gene"="gene_id")))
summary(D2_4_res)
write.csv(D2_D4_res_annot, "/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/DE_longread_analysis/DE_D2_4_annot", quote = F)

D2_10_res <- results(dds, 
                     contrast = c("Condition","D2", "D10"), 
                     alpha = 0.05)
D2_10_res <- subset(D2_10_res, padj < 0.05)
D2_10_res_annot <- setDT(as.data.frame(D2_10_res), keep.rownames = TRUE)[]
D2_10_res_annot <- inner_join(D2_10_res_annot, sqanti, by = c("rn" = "isoform"))
D2_D10_res_annot <- as.data.frame(left_join(D2_10_res_annot, gencode_v38_gene, by = c("associated_gene"="gene_id")))
summary(D2_10_res)
write.csv(D2_D10_res_annot, "/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/DE_longread_analysis/DE_D2_10_annot", quote = F)


D4_10_res <- results(dds, 
                     contrast = c("Condition","D4", "D10"), 
                     alpha = 0.05)
D4_10_res <- subset(D4_10_res, padj < 0.05)
D4_10_res_annot <- setDT(as.data.frame(D4_10_res), keep.rownames = TRUE)[]
D4_10_res_annot <- inner_join(D4_10_res_annot, sqanti, by = c("rn" = "isoform"))
D4_D10_res_annot <- as.data.frame(left_join(D4_10_res_annot, gencode_v38_gene, by = c("associated_gene"="gene_id")))
summary(D4_10_res)
write.csv(D4_D10_res_annot, "/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/DE_longread_analysis/DE_D4_10_annot", quote = F)

############################
# Create MA plot
plotMA(D0_2_res)
plotMA(D0_4_res)
plotMA(D0_10_res)
plotMA(D2_4_res)
plotMA(D2_10_res)
plotMA(D4_10_res)

####time series 
##Time series analysis where we find all DE lncRNAs across all timepoints
CountData_time_series <- read.csv("/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/isoseq/cluster_output/all_transcriptome_cupcake.mapped_fl_count.txt", row.names = 1)
ColData_time_series <- read.csv("/project/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/isoseq/cluster_output/phenotypes.csv", row.names = 1)
reorder_idx_time_series <- match(rownames(ColData_time_series),colnames(CountData_time_series))
str(ColData_time_series)
# Reorder the columns of the count data
reordered_CountData_time_series <- CountData_time_series[ , reorder_idx_time_series]
all(rownames(ColData_time_series) %in% colnames(CountData_time_series))


dds_time_series <- DESeqDataSetFromMatrix(countData =  reordered_CountData_time_series,
                                          colData =  ColData_time_series,
                                          design = ~ Condition)
str(dds_time_series$Condition)#for charles when he inquires about why it was changed to factors
dds_time_series <- estimateSizeFactors(dds_time_series)

# Extract the normalized counts
normalized_counts_time_series <- counts(dds_time_series, normalize = TRUE)
# Transform the normalized counts 
vsd_time_series <- vst(dds_time_series, blind = TRUE)
# Extract the matrix of transformed counts
vsd_mat_time_series <- assay(vsd_time_series)
# Compute the correlation values between samples
vsd_cor_time_series <- cor(vsd_mat_time_series, vsd_mat_time_series)
# Plot the PCA of PC1 and PC2
plotPCA(vsd_time_series, intgroup= "Condition")
pheatmap(vsd_cor, annotation = select(ColData, Condition))

## DESeq2 for time series from http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments
dds_time_series <- DESeq(dds_time_series, test="LRT", reduced =~1)
res_time_series <- results(dds_time_series)

res_time_series <- subset(res_time_series, padj < 0.05)
summary(res_time_series)

CountData <- read.csv("/Users/aa9gj/Desktop/days0_10_cupcake.mapped_fl_count.txt", row.names = 1)



