

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("reshape2")


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")




#���ð�
library(GSVA)
library(limma)
library(GSEABase)
setwd("...")          #���ù���Ŀ¼

#����ssGSEA�ĺ���
immuneScore=function(expFile=null, gmtFile=null, socreFile=null){
	#��ȡ���������ļ������������ļ�����
	rt=read.table(expFile, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	mat=avereps(mat)
	mat=mat[rowMeans(mat)>0,]
	
	#��ȡ���ݼ��ļ�
	geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
	
	#ssgsea����
	ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
	#����ssGSEA score��������
	normalize=function(x){
	  return((x-min(x))/(max(x)-min(x)))}
	#��ssGSEA score���н���
	ssgseaOut=normalize(ssgseaScore)
	ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
	write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)
}

#ssGSEA����
immuneScore(expFile="symbol.txt", gmtFile="immune.gmt", socreFile="immScore.TCGA.txt")
#immuneScore(expFile="ICGCsymbol.txt", gmtFile="immune.gmt", socreFile="immScore.ICGC.txt")



#���ð�
library(limma)
library(ggpubr)
library(reshape2)

#������������Է�������
scoreCor=function(riskFile=null, scoreFile=null, project=null){
    #��ȡssGSEA����ļ�
	data=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
	#ȥ��������Ʒ
	group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
	group=sapply(strsplit(group,""), "[", 1)
	group=gsub("2", "1", group)
	data=t(data[,group==0])
	if(project=="TCGA"){
		rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data)) }
	if(project=="ICGC"){
		rownames(data)=gsub("(.*?)\\-(.*?)\\-.*", "\\2", rownames(data)) }
	data=avereps(data)
	
	#��ȡ�����ļ�
	risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
	
	#�ϲ�����
	sameSample=intersect(row.names(data),row.names(risk))
	data=data[sameSample,,drop=F]
	risk=risk[sameSample,,drop=F]
	rt=cbind(data,risk[,c("riskScore","risk")])
	rt=rt[,-(ncol(rt)-1)]
	
	#������ϸ����������ͼ
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
	#���ͼƬ�ļ�
	pdf(file=paste0(project,".immCell.pdf"), width=7, height=6)
	print(p)
	dev.off()
	
	#��������ع��ܻ�������ͼ
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
	#���ͼƬ�ļ�
	pdf(file=paste0(project,".immFunction.pdf"), width=7, height=6)
	print(p)
	dev.off()
}

#���߲������
scoreCor(riskFile="riskAll.txt", scoreFile="immScore.TCGA.txt", project="TCGA")



####