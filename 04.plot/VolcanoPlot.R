library(ggplot2)
raw_count_filt <- read.table("norm_mouseVSnmr_H_all-DESeq2.gene.txt", header = TRUE, sep = "\t")
head(raw_count_filt)

raw_count_filt$Change = as.factor(ifelse(raw_count_filt$padj < 0.01 & abs(raw_count_filt$log2FoldChange) > 1.5,
                           ifelse(raw_count_filt$log2FoldChange > 1.5 ,'UP','DOWN'),'NONE'))

ggplot(data = raw_count_filt, aes(x = log2FoldChange, y = -log10(padj), color = Change)) +
  geom_point(alpha = 0.8, size = 1) + xlim(-20, 20) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  geom_hline(yintercept = -log10(0.01), linetype = 4) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = 4) +
  scale_color_manual(name = "", values = c("coral3", "seagreen4", "grey"), limits = c("UP", "DOWN", "NONE"))

