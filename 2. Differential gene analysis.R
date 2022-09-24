

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

library("limma")

setwd("...")                 #设置工作目录
inputFile="SENexp.txt"                                             #输入文件
fdrFilter=0.05                                                   #fdr临界值
logFCfilter=1                                                     #logFC临界值

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.1,]
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])      
treatNum=length(group[group==0])    
Type=c(rep(1,conNum), rep(2,treatNum))

#读取输入文件
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#差异分析
for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

#绘制热图需要的文件
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffSENexp.txt",sep="\t",col.names=F,quote=F)

#绘制火山图
pdf(file="vol.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="lightsalmon",cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="deepskyblue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图
library(pheatmap)
hmExp=data[as.vector(outDiff[,1]),]
hmExp=log2(hmExp+0.01)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=6,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         #scale = "row",
           color = colorRampPalette(c("deepskyblue", "white","lightsalmon"))(500),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=6,
         fontsize_col=10)
dev.off()

library(ggpubr)


rt=read.table("diffSENexp.txt",sep="\t",header=T,check.names=F,row.names=1)       #读取输入文件
rt=t(rt)
type=c(rep("N",conNum),rep("T",treatNum))

#准备箱线图的输入文件
data=data.frame()
for(i in colnames(rt)){
  data=rbind(data,cbind(expression=log2(rt[,i]+1),gene=i,type))
}
write.table(data,file="diffSENExp1.txt",sep="\t",row.names=F,quote=F)

#绘制箱型图
data=read.table("diffSENExp1.txt",sep="\t",header=T,check.names=F)       #读取箱线图输入文件
p=ggboxplot(data, x="gene", y="expression", color = "type", 
     ylab="Gene expression",
     xlab="",
     palette = c("blue","orange") )
p=p+rotate_x_text(60)
pdf(file="boxplot.pdf",width=11,height=5)                          #输出图片文件
p
dev.off()
#