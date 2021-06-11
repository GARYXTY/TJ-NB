rm(list=ls())
options(stringsAsFactors = F)


library(survival)
library(survminer)
library(rms)

load('./pre/GSE49710/GSE49710.Rdata')
risk<-read.csv('./2.Appraisal/GSE49710_ph.csv')
phenotype_GSE49710$score=risk$status_p
cluster=read.csv('./consensus_cluster/untitled_consensus_cluster.k=2.consensusClass.csv')$cluster
phenotype_GSE49710$cluster=cluster[224:721]

rt=cbind('futime'=as.numeric(phenotype_GSE49710$os_days),
         'fustat'=phenotype_GSE49710$os_bin,
         'age'=phenotype_GSE49710$age.at.diagnosis.ch1,
         'MYCN'=phenotype_GSE49710$mycn.status.ch1,
         'gender'=phenotype_GSE49710$Sex.ch1,
         'stage'=phenotype_GSE49710$inss.stage.ch1,
         'Risk'=phenotype_GSE49710$high.risk.ch1,
         'Probability'=phenotype_GSE49710$score
)
rt=as.data.frame(rt)
rownames(rt)=rownames(phenotype_GSE49710)

rt$futime=as.numeric(rt$futime)
rt$fustat=as.numeric(rt$fustat)
rt$age=as.numeric(rt$age)/365
#rt$age=factor(ifelse(rt$age>1.5,'older','younger'),levels = c('younger','older'))
rt$MYCN=factor(ifelse(rt$MYCN==0,'Non-amplified','Amplified'),levels = c('Non-amplified','Amplified'))
rt$gender=as.factor(rt$gender)
rt$stage=factor(rt$stage,levels = c('1','2','3','4','4S'))
rt$Risk=factor(ifelse(rt$Risk==1,'High','Low'),levels = c('Low','High'))
rt$Probability=factor(ifelse(rt$Probability<0.5,'Low','High'),levels = c('Low','High'))


rt=rt[,-5]

library(regplot)
multiCox <- coxph(Surv(futime, fustat) ~ ., data = rt)
regplot(multiCox,observation=rt[1,], 
        plots = c('density','spikes'),
        points = T,
        title = '',
        failtime = c(1095,1825,3650), 
        prfail = TRUE )
dev.off()





library(rms)
dd<-datadist(rt)
options(datadist="dd")

f <- cph(Surv(futime,fustat)~., data = rt, x=T, y = T, surv = T)
survival <- Survival(f)
survival1 <- function(x)survival(1825,x) 
nom <- nomogram(f, 
                fun = survival1, 
                fun.at = c(0.1,seq(0.1,0.9,by = 0.2), 0.9), 
                funlabel = "5-year-survival")
plot(nom)
dev.off()

calibration <- calibrate(f, cmethod='KM', method="boot", u=1825, m=166, B=250)  

tiff('./5.clinic/result/calibration.tiff',res = 300,width = 1600,height = 1600,compression = 'lzw')
plot(calibration,lwd=2,lty=1,  
     
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),  
     
     xlim=c(0.6,1),ylim=c(0.3,1),  
     
     xlab="Nomogram-Predicted Probability of 5-year-survival",  
     
     ylab="Actual 5-year-survival",  
     
     col=c(rgb(192,98,83,maxColorValue=255)))  
dev.off()

survConcordance(Surv(futime,fustat) ~ predict(f), data = rt) # C-index 0.889
