
library(pheatmap)

conNum=59                       
treatNum=539                       
expFile="output.txt"        
geneFile="gene.txt"    
setwd("...")


rt=read.table(expFile,header=T,sep="\t",row.names=1,check.names=F)
geneRT=read.table(geneFile,header=F,sep="\t",check.names=F)
hmExp=rt[as.vector(geneRT[,1]),]


hmExp=log2(hmExp+0.1)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(hmExp)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=6,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         scale="row",
         fontsize = 12,
         fontsize_row=10,
         fontsize_col=10)
dev.off()


####