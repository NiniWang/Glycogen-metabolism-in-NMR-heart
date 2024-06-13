X<-read.table("clipboard",header=T,sep="\t")
y<-X[,2:7]
rownames(y)<-X$gene_name
head(y)
library(pheatmap)

pheatmap(y,scale="row",  cluster_rows =T, cluster_cols =T ,border_color=NA,
         fontsize_row = 10, fontsize_col = 10 ,
         show_rownames = F, show_colnames = T,drop_levels = FALSE)
