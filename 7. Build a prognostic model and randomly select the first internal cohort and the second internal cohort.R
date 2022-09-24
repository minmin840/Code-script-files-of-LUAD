
#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("survivalROC")

library(survival)
library(caret)
library(glmnet)
library(survminer)
library(survivalROC)

pFilter=0.05         
setwd("...")                       
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)       
rt[,"futime"]=rt[,"futime"]                                            
rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)


for(i in 1:1000){

		inTrain<-createDataPartition(y=rt[,3],p=0.5,list=F)
		train<-rt[inTrain,]
		test<-rt[-inTrain,]
		trainOut=cbind(id=row.names(train),train)
		testOut=cbind(id=row.names(test),test)
		
		#############µ¥ÒòËØCOX·ÖÎö#############
		outTab=data.frame()
		sigGenes=c("futime","fustat")
		for(i in colnames(train[,3:ncol(train)])){
					 cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
					 coxSummary = summary(cox)
					 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
					 if(coxP<pFilter){
					      a=rt[,i]<=median(rt[,i])
					      diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	                      kmPvalue=1-pchisq(diff$chisq,df=1)
	                      if(kmPvalue<pFilter){
					           sigGenes=c(sigGenes,i)
					           outTab=rbind(outTab,
						              cbind(id=i,
						              HR=coxSummary$conf.int[,"exp(coef)"],
						              HR.95L=coxSummary$conf.int[,"lower .95"],
						              HR.95H=coxSummary$conf.int[,"upper .95"],
						              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
						              )
					        }
					  }
	  }
	  train=train[,sigGenes]
	  test=test[,sigGenes]
	  uniSigExp=train[,sigGenes]
	  uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
	

	  multiCox <- coxph(Surv(futime, fustat) ~ ., data = train)
	  multiCox=step(multiCox,direction = "both")
	  multiCoxSum=summary(multiCox)
		

	  outMultiTab=data.frame()
	  outMultiTab=cbind(
		               coef=multiCoxSum$coefficients[,"coef"],
		               HR=multiCoxSum$conf.int[,"exp(coef)"],
		               HR.95L=multiCoxSum$conf.int[,"lower .95"],
		               HR.95H=multiCoxSum$conf.int[,"upper .95"],
		               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
	
	
		riskScore=predict(multiCox,type="risk",newdata=train)           #
		coxGene=rownames(multiCoxSum$coefficients)
		coxGene=gsub("`","",coxGene)
		outCol=c("futime","fustat",coxGene)
		medianTrainRisk=median(riskScore)
		risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
		trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))

		riskScoreTest=predict(multiCox,type="risk",newdata=test)        #
		riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
		testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest));         
		
		if((length(levels(factor(risk)))==1) | (length(levels(factor(riskTest)))==1)){
				next
		}
		

		diff=survdiff(Surv(futime, fustat) ~risk,data = train)
		pValue=1-pchisq(diff$chisq,df=1)
		roc = survivalROC(Stime=train$futime, status=train$fustat, marker = riskScore, predict.time =5,  method="KM")

		diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
		pValueTest=1-pchisq(diffTest$chisq,df=1)
		rocTest = survivalROC(Stime=test$futime, status=test$fustat, marker = riskScoreTest, predict.time =5,  method="KM")
	
		if((pValue<0.01) & (roc$AUC>0.75) & (pValueTest<0.03) & (rocTest$AUC>0.70)){
	
			   write.table(trainOut,file="train.txt",sep="\t",quote=F,row.names=F)
			   write.table(testOut,file="test.txt",sep="\t",quote=F,row.names=F)
			   write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)
			   write.table(outMultiTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
			   write.table(testRiskOut,file="riskTest.txt",sep="\t",quote=F,row.names=F)
			   write.table(trainRiskOut,file="riskTrain.txt",sep="\t",quote=F,row.names=F)
			   allRiskOut=rbind(trainRiskOut, testRiskOut)
			   write.table(allRiskOut,file="riskAll.txt",sep="\t",quote=F,row.names=F)
			   break
		}
}

#