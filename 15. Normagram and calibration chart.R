#install.packages("rms")

library(rms)
setwd("...")                         #设置工作目录

#列线图绘制
riskFile="indepInput.txt"     #风险输入文件
outFile="Nomogram.pdf"      #输出列线图文件名称
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #读取风险文件
rt=risk[,1:(ncol(risk)-0)]
rt[,"futime"]=rt[,"futime"]/365  
























#数据打包
dd <- datadist(rt)
options(datadist="dd")



#生成函数
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)

#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  
#nomogram可视化
pdf(file=outFile,height=8,width=9)
plot(nom)
dev.off()

#calibration curve
time=1   #预测1年calibration
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="1calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 1-Year OS",ylab="Actual 1-Year OS(proportion)",col="red",sub=F)
dev.off()


#calibration curve
time=2   #预测2年calibration
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="2calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 2-Year OS",ylab="Actual 2-Year OS(proportion)",col="red",sub=F)
dev.off()



#calibration curve
time=3   #预测3年calibration
f <- cph(Surv(futime, fustat) ~age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="3calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 3-Year OS",ylab="Actual 3-Year OS(proportion)",col="red",sub=F)
dev.off()
#calibration curve
time=5   #预测5年calibration
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="5calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 5-Year OS",ylab="Actual 5-Year OS(proportion)",col="red",sub=F)
dev.off()



#
#~!!!