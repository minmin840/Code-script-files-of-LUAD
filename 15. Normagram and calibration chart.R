#install.packages("rms")

library(rms)
setwd("...")                         #���ù���Ŀ¼

#����ͼ����
riskFile="indepInput.txt"     #���������ļ�
outFile="Nomogram.pdf"      #�������ͼ�ļ�����
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        #��ȡ�����ļ�
rt=risk[,1:(ncol(risk)-0)]
rt[,"futime"]=rt[,"futime"]/365  
























#���ݴ��
dd <- datadist(rt)
options(datadist="dd")



#���ɺ���
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)

#����nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  
#nomogram���ӻ�
pdf(file=outFile,height=8,width=9)
plot(nom)
dev.off()

#calibration curve
time=1   #Ԥ��1��calibration
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="1calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 1-Year OS",ylab="Actual 1-Year OS(proportion)",col="red",sub=F)
dev.off()


#calibration curve
time=2   #Ԥ��2��calibration
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="2calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 2-Year OS",ylab="Actual 2-Year OS(proportion)",col="red",sub=F)
dev.off()



#calibration curve
time=3   #Ԥ��3��calibration
f <- cph(Surv(futime, fustat) ~age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="3calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 3-Year OS",ylab="Actual 3-Year OS(proportion)",col="red",sub=F)
dev.off()
#calibration curve
time=5   #Ԥ��5��calibration
f <- cph(Surv(futime, fustat) ~ age	+gender	+stage 	+T	+ N	+riskScore, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=103, B=1000)
pdf(file="5calibration.pdf",height=6,width=8)
plot(cal,xlab="Nomogram-Predicted Probability of 5-Year OS",ylab="Actual 5-Year OS(proportion)",col="red",sub=F)
dev.off()



#
#~!!!