######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


#引用包
library(limma)
library(ggpubr)
setwd("...")      #设置工作目录

#读取肿瘤突变负荷文件
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#读取风险数据文件
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)
	
#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB=log2(data$TMB+1)
	
#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
#绘制小提琴图
boxplot=ggviolin(data, x="risk", y="TMB", fill="risk",
			      xlab="",
			      ylab="Tumor tmbation burden (log2)",
			      legend.title="",
			      palette = c("#0066FF","#FF0000"),
			      add = "boxplot", add.params = list(fill="white"))+ 
	stat_compare_means(comparisons = my_comparisons)
	
#输出图片
pdf(file="riskTMB.pdf", width=5, height=4.5)
print(boxplot)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

