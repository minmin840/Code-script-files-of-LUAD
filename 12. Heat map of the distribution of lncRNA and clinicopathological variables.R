#install.packages("pheatmap")


library(pheatmap)           
expFile="riskAll.txt"         
cliFile="clinical.txt"     
outFile="heatmap.pdf"      
var="risk"                  #按照临床性状（此处为风险）对样品排序
setwd("...")                      #设置工作目录
rt=read.table(expFile, sep="\t", header=T, row.names=1, check.names=F)       #读取表达文件
Type=read.table(cliFile, sep="\t", header=T, row.names=1, check.names=F)     #读取临床文件 

#样品取交�?
sameSample=intersect(colnames(rt),row.names(Type))
rt=rt[,sameSample]
Type=Type[sameSample,]
Type=Type[order(Type[,var]),]   #按临床性状排序
rt=rt[,row.names(Type)]

#绘制热图
pdf(outFile,height=5,width=8)
pheatmap(rt, annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,    #是否聚类
         scale="row",   #基因矫正
         show_colnames=F,
         fontsize=7.5,
         fontsize_row=7,
         fontsize_col=5)
dev.off()


