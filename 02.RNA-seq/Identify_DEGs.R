getwd()
library("DESeq2")
# Generate a matrix of lengths of each gene in each sample to use for cross-ortholog normalization
gene_effectiveL<-read.table("hypo_L_GeneEffectiveLength.txt",header = T,sep="\t",row.names = 1)
head(gene_effectiveL)

# Create a table of reads count and their conditions
raw_count<-read.table("hypo_L_Readscount.txt",header = T,sep="\t",row.names = 1)
head(raw_count)

colDate<-data.frame(row.names = c("hypo_mouse_L_1","hypo_mouse_L_2","hypo_mouse_L_3","hypo_NMR_L_1","hypo_NMR_L_2","hypo_NMR_L_3"),
                    condition=factor(c("control","control","control","NMR","NMR","NMR")),
                    levels=c("NMR","control"))

dds<-DESeqDataSetFromMatrix(countData = raw_count,colData = colDate,design = ~condition)
dds = DESeq2::estimateSizeFactors(dds, normMatrix = gene_effectiveL)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
######Differential expression analysis####
dds <- DESeq(dds)
res<-results(dds)
class(res)
res<-as.data.frame(res)
res<-cbind(rownames(res),res)
colnames(res)<- c('geng_id',"baseMean" ,"log2FoldChange","lfcSE" ,"stat","pvalue","padj" )
write.table(res,"hypo_mouseVSnmr_L_all-DESeq2.gene.txt",sep = '\t',col.names = T,row.names = F,quote = FALSE,na='')

resSig<-res[which(res$padj<0.01 & abs(res$log2FoldChange)>1.5),]
##新增一列，将log2FoldChange>0标注为up，<0标准为down
resSig[which(resSig$log2FoldChange>0),'up_down']<-'up'
resSig[which(resSig$log2FoldChange<0),'up_down']<-'down'
##保存数据
write.table(resSig,"Norm_L-diff-0.05-FC_1.gene.xls",sep = '\t',
            col.names = T,row.names = F,quote = FALSE,na='')

