

#install.packages("pheatmap")


summary(fit)
library(survival)
library(survminer)
library(timeROC)
library(pheatmap)
setwd("...")             #设置工作目录
rt=read.table("riskTrain.txt",sep="\t",header=T,row.names=1,check.names=F)     
rt=rt[order(rt$riskScore),]                                   




riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])


line=rt[,"riskScore"]
line[line>10]=10


pdf(file="riskScore.pdf",width = 10,height = 3.5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("blue",lowLength),
     rep("orange",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("orange","blue"),cex=1.2)
dev.off()

#绘制生存状态图
color=as.vector(rt$fustat)
color[color==1]="orange"
color[color==0]="blue"
pdf(file="survStat.pdf",width = 10,height = 3.5)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("orange","blue"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

#绘制风险热图
rt1=log2(rt[c(3:(ncol(rt)-2))]+0.01)
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 10,height = 3)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         color = colorRampPalette(c("blue", "white", "orange"))(50) )
dev.off()

data=read.table("riskTrain.txt",header=T,sep="\t",check.names=F)
diff=survdiff(Surv(futime, fustat) ~risk,data = data)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = data)
pdf(file="survival.pdf",width=5.5,height=5)
plot(fit,
      lwd=2,
      col=c("red","blue"),
      xlab="Time (year)",
      ylab="Survival rate",
      main=paste("Survival curve (p=", pValue ,")",sep=""),
      mark.time=T)
legend("topright",
      c("High risk", "Low risk"),
     lwd=2,
      col=c("red","blue"))
dev.off()



inputFile="riskTrain.txt"        #输入文件
survFile="survival.pdf"         #生存曲线文件
rocFile="ROC.pdf"               #ROC曲线文件



#读取输入文件
rt=read.table(inputFile,header=T,sep="\t")

#比较高低风险组生存差异，得到显著性p值
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

#绘制生存曲线
surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=5,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
pdf(file=survFile,onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()


###ROC曲线
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file=rocFile,width=5,height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
        c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
          paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
          paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
        col=c("green",'blue','red'),lwd=2,bty = 'n')
dev.off()

