

#install.packages("survival")

library(survival)
pFilter=0.05                                                            

setwd("...")                   

rt=read.table("geneexp.txt",header=T,sep="\t",check.names=F,row.names=1)  
rt$futime=rt$futime/365                                                    

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	   if(sd(rt[,gene])<0.01){
	      next}
	   a=rt[,gene]<=median(rt[,gene])
	   diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	   pValue=1-pchisq(diff$chisq,df=1)
	   fit=survfit(Surv(futime, fustat) ~ a, data = rt)
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	   if((pValue<pFilter) & (coxP<pFilter)){
	         sigGenes=c(sigGenes,gene)
	       	 outTab=rbind(outTab,
	                      cbind(gene=gene,
	                            KM=pValue,
	                            B=coxSummary$coefficients[,"coef"],
	                            SE=coxSummary$coefficients[,"se(coef)"],
	                            HR=coxSummary$conf.int[,"exp(coef)"],
	                            HR.95L=coxSummary$conf.int[,"lower .95"],
	                            HR.95H=coxSummary$conf.int[,"upper .95"],
			                    pvalue=coxP) )
	  }
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)    #Êä³ö»ùÒòºÍpÖµ±í¸ñÎÄ¼þ
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)    #Êä³ö»ùÒòºÍpÖµ±í¸ñÎÄ¼þ
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))


pdf(file="forest.pdf", width = 6,height = 30)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'blue')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#