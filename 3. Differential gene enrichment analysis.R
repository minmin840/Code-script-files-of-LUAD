

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot")

setwd("...")          #设置工作目录

library("org.Hs.eg.db")          #引用包
rt=read.table("gene.txt",sep="\t",check.names=F,header=T)    #读取文件
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)    #输出结果

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

rt=read.table("id.txt",sep="\t",header=T,check.names=F)       #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                             #去除基因id为NA的基因
gene=rt$entrezID

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #富集分析
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果

#柱状图
pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#气泡图
pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 30)
dev.off()
ego=read.table("KEGG.txt", header = T,sep="\t",check.names=F)           #读取kegg富集结果文件
go=data.frame(Category = "ALL",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

#读取基因的logFC文件
id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go, genelist)

#绘制KEGG气泡图
pdf(file="KEGGBubble.pdf",width = 10,height = 8)
GOBubble(circ, labels = 3,table.legend =F)
dev.off()

#绘制KEGG圈图
pdf(file="KEGGCircle.pdf",width = 15,height = 9)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10)            #nsub=10中10代表显示KEGG的数据，可修改
dev.off()

#绘制KEGG热图
termNum = 20                                     #限定term数目
geneNum = nrow(genelist)                         #限定基因数目
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="KEGGHeat.pdf",width = 9,height = 5)
GOHeat(chord, nlfc =1, fill.col = c('red', 'white', 'blue'))
dev.off()

#####