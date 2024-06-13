raw_count_filt<-read.table("clipboard",header = T,sep="\t")
head(raw_count_filt)
library(ggplot2)
ggplot(data = raw_count_filt, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha=0.8, size = 1)

ggplot(data = raw_count_filt,  aes(x = log2FoldChange, y = -log10(padj), color = up_down, )) +
  geom_point(alpha=0.8, size = 1.5) + xlim(-20,20) +
  theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
  geom_hline(yintercept=-log10(0.01) ,linetype=4) +
  geom_vline(xintercept=c(-1.5,1.5) ,linetype=4 ) +
  scale_color_manual(name = "", values = c("deeppink3", "dodgerblue3", "grey"), limits = c("UP", "DOWN", "NONE")) 
