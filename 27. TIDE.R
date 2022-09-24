######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#引用包
library(limma)
library(ggpubr)
tideFile="TIDE.txt"          #TIDE的打分文件
riskFile="risk.all.txt"      #风险文件
setwd("...")     #设置工作目录

#读取TIDE数据
tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]
row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
	
#合并数据
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)
	
#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制小提琴图
gg1=ggviolin(data, x="risk", y="TIDE", fill = "risk", 
	         xlab="", ylab="TIDE",
	         palette=c("#0066FF","#FF0000"),
	         legend.title="Risk",
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")

#输出图形	
pdf(file="TIDE.pdf", width=6, height=5)
print(gg1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

