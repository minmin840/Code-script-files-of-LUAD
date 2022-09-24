

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")


#install.packages("digest")
#install.packages("GOplot")

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")


#install.packages("digest")
#install.packages("GOplot")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("enrichplot")

setwd("...")          #���ù���Ŀ¼

library("org.Hs.eg.db")          #���ð�
rt=read.table("gene.txt",sep="\t",check.names=F,header=T)    #��ȡ�ļ�
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #�ҳ������Ӧ��id
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)    #������

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

rt=read.table("id.txt",sep="\t",header=T,check.names=F)           #��ȡid.txt�ļ�
rt=rt[is.na(rt[,"entrezID"])==F,]                                 #ȥ������idΪNA�Ļ���
gene=rt$entrezID

#GO��������
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 #���渻�����

#��״ͼ
pdf(file="barplot.pdf",width = 10,height = 8)
barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#����ͼ
pdf(file="bubble.pdf",width = 10,height = 8)
dotplot(kk,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

library(GOplot)


ego=read.table("GO.txt", header = T,sep="\t",check.names=F)           #��ȡkegg��������ļ�
go=data.frame(Category = ego$ONTOLOGY,ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

#��ȡ�����logFC�ļ�
id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go, genelist)

#����GO����ͼ
pdf(file="GOBubble.pdf",width = 10,height = 8)
GOBubble(circ, labels = 3,table.legend =F)
dev.off()

#����GOȦͼ
pdf(file="GOCircle.pdf",width = 11,height = 6)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10)            #nsub=10��10������ʾGO�����ݣ����޸�
dev.off()

#����GO��ͼ
termNum = 20                                     #�޶�term��Ŀ
geneNum = nrow(genelist)                         #�޶�������Ŀ
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="GOHeat.pdf",width = 11,height = 5)
GOHeat(chord, nlfc =1, fill.col = c('red', 'white', 'blue'))
dev.off()
