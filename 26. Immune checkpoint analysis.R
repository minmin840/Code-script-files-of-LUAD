
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="symbol.txt"      #表达输入文件
riskFile="riskAll.txt"       #风险输入文件
geneFile="gene.txt"       #免疫检查点的基因文件
setwd("...")     #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
#读取基因文件
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
data=avereps(data)
	
#合并数据
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]

#提取显著差异的基因
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
	if(sd(rt1[,i])<0.001){next}
	wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]

#把数据转换成ggplot2输入文件
rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")
	
#设置比较组
group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
	
#绘制箱线图
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Risk",
				  width=0.8,
				  palette = c("#0066FF", "#FF0000") )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=risk),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
	
#输出图片
pdf(file="checkpoint.diff.pdf", width=8, height=5)
print(boxplot)
dev.off()


##