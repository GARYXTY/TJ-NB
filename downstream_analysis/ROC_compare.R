rm(list = ls())                                    
options(stringsAsFactors = F)

library(survival)
library(timeROC)
#library(survivalROC)

### Zhong et al
######
load('./pre/GSE49711/GSE49711.Rdata')

Zhong=c('ERCC6L','AHCY','STK33','NCAN')
dat_Zhong=dat_GSE49711[Zhong,]
dat_Zhong=as.matrix((dat_Zhong))
dat_Zhong=t(scale(t(dat_Zhong)))

beta_Zhong=c(0.408,0.478,0.345,0.136)

score_Zhong=beta_Zhong%*%dat_Zhong

pd_GSE49711$os_days=pd_GSE49711$os_days/365
pd_GSE49711$Zhong=t(score_Zhong)

roc<-timeROC(T=pd_GSE49711$os_days,
             delta=pd_GSE49711$os_bin,
             marker=pd_GSE49711$Zhong,
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)


x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_Zhong <- data.frame(x,y) 
plotdata_Zhong$group <- "Zhong et al"
######


### Katleen et al
######
load('./pre/GSE49711/GSE49711.Rdata')

Katleen=read.csv('./3.Alternatives/Katleen De Preter.csv')$gene
dat_Katleen=dat_GSE49711[Katleen,]
dat_Katleen=as.matrix((dat_Katleen))
dat_Katleen=t(scale(t(dat_Katleen)))

pd_GSE49711$os_days=pd_GSE49711$os_days/365
df=cbind(futime=pd_GSE49711$os_days,fustat=pd_GSE49711$os_bin,t(dat_Katleen))
df=as.data.frame(df)
df=as.data.frame(lapply(df,as.numeric))

cox <- coxph(Surv(futime, fustat) ~ ., data = df)
coxSummary = summary(cox)
beta_Katleen=coxSummary$coefficients[,"coef"]

score_Katleen=beta_Katleen%*%dat_Katleen


pd_GSE49711$Katleen=t(score_Katleen)

roc<-timeROC(T=pd_GSE49711$os_days,
             delta=pd_GSE49711$os_bin,
             marker=pd_GSE49711$Katleen,
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)


x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_Katleen <- data.frame(x,y) 
plotdata_Katleen$group <- "De Preter et al"
######


### MYCN status
######
load('./pre/GSE49711/GSE49711.Rdata')

pd_GSE49711$os_days=pd_GSE49711$os_days/365
pd_GSE49711=pd_GSE49711[!pd_GSE49711$mycn.status.ch1=="N/A",]
pd_GSE49711$mycn.status.ch1=as.numeric(pd_GSE49711$mycn.status.ch1)
roc<-timeROC(T=pd_GSE49711$os_days,
             delta=pd_GSE49711$os_bin,
             marker=pd_GSE49711$mycn.status.ch1,
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)


x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_MYCN <- data.frame(x,y) 
plotdata_MYCN$group <- "MYCN"
######

### Feng et al
######
rt=read.csv('./pre/GSE49710/GSE49710_phenotype.csv')[,c(13,14)]
rt$p=read.csv('./2.Appraisal/GSE49710_ph.csv')$status_p
rt$os_days=rt$os_days/365


roc<-timeROC(T=rt$os_days,
             delta=rt$os_bin,
             marker=rt$p,
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)

x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_NN <- data.frame(x,y) 
plotdata_NN$group <- "Neural Network"

######


### Plot tROC
######

plotdata <- rbind(plotdata_MYCN,plotdata_Zhong,plotdata_Katleen,plotdata_NN)
plotdata$group=factor(plotdata$group,levels = c('MYCN','Zhong et al','De Preter et al','Neural Network'))

g <- ggplot(plotdata) + 
  geom_path(aes(x = x, y = y,colour = group), size=1) + 
  labs(x="1 - Specificity", y = "Sensitivity") +
  ggpubr::theme_classic2()+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme(axis.text=element_text(size=14,face="bold"),axis.title=element_text(size=16,face="bold"))+
  geom_abline(intercept=0, slope=1,linetype="dashed")+
  scale_colour_hue(name="my legend", labels=c(
    'MYCN amplification: 0.723',
    'Zhong et al: 0.923',
    'De Preter et al: 0.943',
    'DL model: 0.974'
  )
  )+
  theme(legend.title=element_blank())+theme(legend.text = element_text(size = 14,face = "bold"))+
  theme(legend.position=c(.766,.25))
tiff('./3.Alternatives/result/ROC_compare.tiff',res = 300,width = 2400,height = 1500,compression = 'lzw')
g
dev.off()

write.csv(beta_Katleen,'./3.Alternatives/result/Katleen De Preter coef.csv')
######

Score_GSE49710=data.frame('Zhong et al'=t(score_Zhong),
                          'De Preter et al'=t(score_Katleen),
                          'DL model'=rt$p)
write.csv(Score_GSE49710,'./3.Alternatives/result/Score_GSE49710.csv')
