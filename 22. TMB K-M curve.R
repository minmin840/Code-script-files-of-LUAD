######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("survminer")


#引用包
library(survival)
library(survminer)

tmbFile="TMB.txt"            #肿瘤突变负荷文件
riskFile="risk.all.txt"      #风险文件
setwd("...")      #设置工作目录

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #读取肿瘤突变负荷文件

#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

#获取肿瘤突变负荷最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)

#定义生存分析函数
bioSurvival=function(surData=null, outFile=null){
    #比较高低突变负荷组的生存差异,得到显著性的pvalue
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	#print(surv_median(fit))
	
	#绘制生存曲线
	bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 1,
			           palette = bioCol,
			           surv.median.line = "hv",
			           risk.table=F,
			           cumevents=F,
			           risk.table.height=.25)
	#输出图形
	pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
	print(surPlot)
	dev.off()
}

#绘制肿瘤突变负荷的生存曲线
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#绘制肿瘤突变负荷联合病人风险的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

