

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("scatterplot3d")

library(limma)
library(scatterplot3d)

setwd("...")                        #???ù???Ŀ

myPCA=function(input=null,output=null)
{
		#??ȡ?????????ļ?
		rt=read.table(input,sep="\t",header=T,check.names=F)
		rt=as.matrix(rt)
		rownames(rt)=rt[,1]
		exp=rt[,2:ncol(rt)]
		dimnames=list(rownames(exp),colnames(exp))
		data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
		data=avereps(data)
		data=data[rowMeans(data)>0.5,]
		type=sapply(strsplit(colnames(data),"\\-"),"[",4)
		type=sapply(strsplit(type,""),"[",1)
		type=gsub("2","1",type)
		data=t(data[,type==0])
		rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
		
		#??ȡrisk?????ļ?
		risk=read.table("riskAll.txt",sep="\t",header=T,row.names=1)
		sameSample=intersect(rownames(data),rownames(risk))
		data=data[sameSample,]
		risk=risk[sameSample,]
		group=as.vector(risk[,"risk"])
		
		#PCA????
		data.class <- rownames(data)
		data.pca <- prcomp(data, scale. = TRUE)
		
		#???ӻ?
		color=ifelse(group=="low",3,2)
		pcaPredict=predict(data.pca)
		pdf(file=output,width=5.5,height=5)
		s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
		legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
		dev.off()
}

#ʹ?????л???????��????PCAͼ??????04?ڿ?symbol.txt????ǰĿ¼
myPCA(input="symbol.txt",output="allGene.PCA.pdf")
#ʹ?????????ػ???????��????PCAͼ??????08?ڿ?immuneGeneExp.txt????ǰĿ¼
myPCA(input="diffSENexp.txt",output="diffSENESCENTgene.PCA.pdf")
#ʹ??????????lncRNA????��????PCAͼ??????09?ڿ?immuneLncRNAexp.txt????ǰĿ¼
myPCA(input="lncRNAexp.txt",output="diffSENESCENTgenecoexpressionlncRNA.PCA.pdf")

#ʹ??ģ??????lncRNA????��????PCAͼ??????12?ڿ?risk.txt????ǰĿ¼
risk=read.table("riskAll.txt",sep="\t",header=T,row.names=1)
data=risk[,3:(ncol(risk)-2)]
group=as.vector(risk[,"risk"])
		
#PCA????
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)
		
#???ӻ?
color=ifelse(group=="low",3,2)
pcaPredict=predict(data.pca)
pdf(file="riskGene.PCA.pdf",width=5.5,height=5)
s3d=scatterplot3d(pcaPredict[,1:3], pch = 16, color=color)
legend("top", legend = c("Low risk","High risk"),pch = 16, inset = -0.2, xpd = TRUE, horiz = TRUE,col=c(3,2))
dev.off()


####