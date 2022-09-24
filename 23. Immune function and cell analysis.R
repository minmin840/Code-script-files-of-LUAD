

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("reshape2")


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")




#引用包
library(GSVA)
library(limma)
library(GSEABase)
setwd("...")          #设置工作目录

#定义ssGSEA的函数
immuneScore=function(expFile=null, gmtFile=null, socreFile=null){
	#读取表达输入文件，并对输入文件处理
	rt=read.table(expFile, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	mat=avereps(mat)
	mat=mat[rowMeans(mat)>0,]
	
	#读取数据集文件
	geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
	
	#ssgsea分析
	ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
	#定义ssGSEA score矫正函数
	normalize=function(x){
	  return((x-min(x))/(max(x)-min(x)))}
	#对ssGSEA score进行矫正
	ssgseaOut=normalize(ssgseaScore)
	ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
	write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)
}

#ssGSEA分析
immuneScore(expFile="symbol.txt", gmtFile="immune.gmt", socreFile="immScore.TCGA.txt")
#immuneScore(expFile="ICGCsymbol.txt", gmtFile="immune.gmt", socreFile="immScore.ICGC.txt")



#引用包
library(limma)
library(ggpubr)
library(reshape2)

#定义免疫相关性分析函数
scoreCor=function(riskFile=null, scoreFile=null, project=null){
    #读取ssGSEA结果文件
	data=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
	#去除正常样品
	group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
	group=sapply(strsplit(group,""), "[", 1)
	group=gsub("2", "1", group)
	data=t(data[,group==0])
	if(project=="TCGA"){
		rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data)) }
	if(project=="ICGC"){
		rownames(data)=gsub("(.*?)\\-(.*?)\\-.*", "\\2", rownames(data)) }
	data=avereps(data)
	
	#读取风险文件
	risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
	
	#合并数据
	sameSample=intersect(row.names(data),row.names(risk))
	data=data[sameSample,,drop=F]
	risk=risk[sameSample,,drop=F]
	rt=cbind(data,risk[,c("riskScore","risk")])
	rt=rt[,-(ncol(rt)-1)]
	
	#对免疫细胞绘制箱线图
	immCell=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
	          "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
	          "Tfh","Th1_cells","Th2_cells","TIL","Treg")
	rt1=rt[,c(immCell,"risk")]
	data=melt(rt1,id.vars=c("risk"))
	colnames(data)=c("Risk","Type","Score")
	data$Risk=factor(data$Risk, levels=c("low","high"))
	p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
	p=p+rotate_x_text(50)
	p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	#输出图片文件
	pdf(file=paste0(project,".immCell.pdf"), width=7, height=6)
	print(p)
	dev.off()
	
	#对免疫相关功能绘制箱线图
	immFunction=c("APC_co_inhibition","APC_co_stimulation","CCR",
	          "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
	          "MHC_class_I","Parainflammation","T_cell_co-inhibition",
	          "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
	rt1=rt[,c(immFunction,"risk")]
	data=melt(rt1,id.vars=c("risk"))
	colnames(data)=c("Risk","Type","Score")
	data$Risk=factor(data$Risk, levels=c("low","high"))
	p=ggboxplot(data, x="Type", y="Score", color = "Risk",
	     ylab="Score",add = "none",xlab="",palette = c("blue","red") )
	p=p+rotate_x_text(50)
	p=p+stat_compare_means(aes(group=Risk),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	#输出图片文件
	pdf(file=paste0(project,".immFunction.pdf"), width=7, height=6)
	print(p)
	dev.off()
}

#免疫差异分析
scoreCor(riskFile="riskAll.txt", scoreFile="immScore.TCGA.txt", project="TCGA")



####