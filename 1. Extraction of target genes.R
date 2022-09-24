#设置工作目录
library(limma)
setwd("...")


rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=data[rowMeans(data)>0,]


gene=read.table("gene.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]

out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="SENexp.txt",sep="\t",quote=F,col.names=F)
