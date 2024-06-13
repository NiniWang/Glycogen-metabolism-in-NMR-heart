# Load necessary library
library("DESeq2")

# Generate a matrix of lengths of each gene in each sample to use for cross-ortholog normalization
# Read the effective lengths of genes from a file
gene_effectiveL <- read.table("merged_{species1}_{species2}_EffectiveLength.output", header = TRUE, sep = "\t", row.names = 1)

# Create a table of reads count and their conditions
# Read the raw count data from a file
raw_count <- read.table("merged_{species1}_{species2}_NumReads.output", header = TRUE, sep = "\t", row.names = 1)

# Define the conditions for each sample
colDate <- data.frame(
  row.names = c("norm_mouse_H_1", "norm_mouse_H_2", "norm_mouse_H_3", "norm_NMR_H_1", "norm_NMR_H_2", "hnorm_NMR_H_3"),
  condition = factor(c("control", "control", "control", "NMR", "NMR", "NMR")),
  levels = c("NMR", "control")
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = raw_count, colData = colDate, design = ~condition)

# Normalize the counts using the effective lengths of genes
dds <- DESeq2::estimateSizeFactors(dds, normMatrix = gene_effectiveL)

# Filter out low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Differential expression analysis
dds <- DESeq(dds)

# Extract results and convert to a data frame
res <- results(dds)
class(res)
res <- as.data.frame(res)

# Add row names as a column and rename columns
res <- cbind(rownames(res), res)
colnames(res) <- c('gene_id', "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# Save the results to a file
write.table(res, "norm_mouseVSnmr_H_all-DESeq2.gene.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, na = '')

# Filter significant results based on adjusted p-value and log2 fold change
resSig <- res[which(res$padj < 0.01 & abs(res$log2FoldChange) > 1.5), ]

# Add a new column to indicate upregulated and downregulated genes
resSig[which(resSig$log2FoldChange > 0), 'up_down'] <- 'up'
resSig[which(resSig$log2FoldChange < 0), 'up_down'] <- 'down'

# Save the significant results to a file
write.table(resSig, "norm_mouseVSnmr_H_DEGs.xls", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, na = '')


