

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

setwd("C:\\Users\\zm\\Desktop\\�����-LUADϸ��˥��lncRNA\\09KEGG��������")          #���ù���Ŀ¼

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
library(GOplot)

rt=read.table("id.txt",sep="\t",header=T,check.names=F)       #��ȡid.txt�ļ�
rt=rt[is.na(rt[,"entrezID"])==F,]                             #ȥ������idΪNA�Ļ���
gene=rt$entrezID

#kegg��������
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #��������
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                          #���渻�����

#��״ͼ
pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()

#����ͼ
pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 30)
dev.off()
ego=read.table("KEGG.txt", header = T,sep="\t",check.names=F)           #��ȡkegg��������ļ�
go=data.frame(Category = "ALL",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

#��ȡ�����logFC�ļ�
id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go, genelist)

#����KEGG����ͼ
pdf(file="KEGGBubble.pdf",width = 10,height = 8)
GOBubble(circ, labels = 3,table.legend =F)
dev.off()

#����KEGGȦͼ
pdf(file="KEGGCircle.pdf",width = 15,height = 9)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10)            #nsub=10��10������ʾKEGG�����ݣ����޸�
dev.off()

#����KEGG��ͼ
termNum = 20                                     #�޶�term��Ŀ
geneNum = nrow(genelist)                         #�޶�������Ŀ
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="KEGGHeat.pdf",width = 9,height = 5)
GOHeat(chord, nlfc =1, fill.col = c('red', 'white', 'blue'))
dev.off()

#####