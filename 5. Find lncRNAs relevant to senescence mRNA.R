
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
setwd("...")          #设置工作目录

corFilter=0.3                                                        #相关系数过滤标准
pvalueFilter=0.001                                                   #p值过滤标准

#读取lncRNA表达文件,并对数据进行处理
rt = read.table("lncRNA.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
lncRNA=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
lncRNA=avereps(lncRNA)
lncRNA=lncRNA[rowMeans(lncRNA)>0.5,]
group=sapply(strsplit(colnames(lncRNA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
lncRNA=lncRNA[,group==0]

#读取基因表达文件,并对数据进行处理
rt = read.table("diffSENexp.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
ARGgene=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
ARGgene=avereps(ARGgene)
ARGgene=ARGgene[rowMeans(ARGgene)>0.5,]
group=sapply(strsplit(colnames(ARGgene),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
ARGgene=ARGgene[,group==0]

#相关性检验
outTab=data.frame()
for(i in row.names(lncRNA)){
	if(sd(lncRNA[i,])>0.5){
		for(j in row.names(ARGgene)){
			x=as.numeric(lncRNA[i,])
		    y=as.numeric(ARGgene[j,])
			corT=cor.test(x,y)
			cor=corT$estimate
			pvalue=corT$p.value
			if((abs(cor)>corFilter) & (pvalue<pvalueFilter)){
				outTab=rbind(outTab,cbind(gene=j,lncRNA=i,cor,pvalue))
			}
		}
    }
}
write.table(file="corResult.txt",outTab,sep="\t",quote=F,row.names=F)        #输出相关性结果
ARGlncRNA=unique(as.vector(outTab[,"lncRNA"]))
ARGlncRNAexp=lncRNA[ARGlncRNA,]
ARGlncRNAexp=rbind(ID=colnames(ARGlncRNAexp),ARGlncRNAexp)
write.table(ARGlncRNAexp,file="lncRNAexp.txt",sep="\t",quote=F,col.names=F)

