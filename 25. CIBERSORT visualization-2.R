#####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")
#install.packages("pheatmap")
#install.packages("vioplot")


#引用包
library(limma)
library(pheatmap)
library(ggpubr)
library(vioplot)
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
riskFile="riskAll.txt"             #风险文件 
setwd("...")     #设置工作目录


#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(immune[,1:(ncol(immune)-3)])


#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=rownames(risk)[risk[,"risk"]=="low"]
highSample=rownames(risk)[risk[,"risk"]=="high"]

#高低风险组的免疫细胞浸润
lowSameSample=intersect(row.names(data), lowSample)
highSameSample=intersect(row.names(data), highSample)
data=t(data[c(lowSameSample,highSameSample),])
conNum=length(lowSameSample)
treatNum=length(highSameSample)


##########绘制柱状图##########
pdf("barplot.pdf",height=10,width=18)
col=rainbow(nrow(data),s=0.7,v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1], ybottom = -0.01, xright = a1[conNum], ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Low risk",cex=2)
rect(xleft = a1[conNum], ybottom = -0.01, xright =a1[length(a1)] , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"High risk",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


##########绘制小提琴图##########
rt=t(data)
pdf("vioplot.pdf", height=8, width=13.5)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.02),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制小提琴图，低风险用绿色表示，高风险用红色表示
for(i in 1:ncol(rt)){
	  if(sd(rt[1:conNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
	    rt[(conNum+1),i]=0.00001
	  }
	  lowData=rt[1:conNum,i]
	  highData=rt[(conNum+1):(conNum+treatNum),i]
	  vioplot(lowData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(lowData, highData)
	  p=wilcoxTest$p.value
	  mx=max(c(lowData,highData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Low risk", "High risk"),
       lwd=3.5, bty="n", cex=1.2,
       col=c("green","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()


##########绘制热图##########
#风险颜色
ann_colors=list()
bioCol=c("green", "red")
names(bioCol)=c("low", "high")
ann_colors[["Risk"]]=bioCol

#绘制热图
Risk=c(rep("low",conNum), rep("high",treatNum))
names(Risk)=colnames(data)
Risk=as.data.frame(Risk)
pdf(file="heatmap.pdf", width=10, height=6)
data[data>0.4]=0.4
pheatmap(data, 
         annotation=Risk,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("green",1), rep("black",2), rep("red",2)))(50),
         cluster_cols =F,
         show_colnames=F,
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()


###