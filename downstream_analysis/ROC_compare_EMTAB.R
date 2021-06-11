rm(list = ls())                                    
options(stringsAsFactors = F)

library(survival)
library(timeROC)


### Zhong et al
######
load('./pre/E-MTAB/E-MTAB.Rdata')

Zhong=c('ERCC6L','AHCY','STK33','NCAN')
dat_Zhong=dat_EMTAB[Zhong,]
dat_Zhong=as.matrix((dat_Zhong))
dat_Zhong=t(scale(t(dat_Zhong)))
beta_Zhong=c(0.408,0.478,0.345,0.136)

score_Zhong=beta_Zhong%*%dat_Zhong


phenotype_EMTAB$Zhong=t(score_Zhong)

roc<-timeROC(T=phenotype_EMTAB$Characteristics.overall.survival.,
             delta=phenotype_EMTAB$Characteristics.overall.survival..1.dead..0.alive.w.o.event.., 
             marker=phenotype_EMTAB$Zhong, 
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)
# AUC:0.43

x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_Zhong <- data.frame(x,y) 
plotdata_Zhong$group <- "Zhong et al"


######


### Katleen et al
######
load('./pre/E-MTAB/E-MTAB.Rdata')

beta_Katleen=read.csv('./3.Alternatives/result/Katleen De Preter coef.csv')

score_Katleen=(beta_Katleen$x)%*%(t(scale(t(dat_EMTAB[beta_Katleen$X,]))))


phenotype_EMTAB$Katleen=t(score_Katleen)

roc<-timeROC(T=phenotype_EMTAB$Characteristics.overall.survival.,
             delta=phenotype_EMTAB$Characteristics.overall.survival..1.dead..0.alive.w.o.event.., 
             marker=phenotype_EMTAB$Katleen, 
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)
# 0.727

x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_Katleen <- data.frame(x,y) 
plotdata_Katleen$group <- "De Preter et al"
######


### MYCN status
######

phenotype_EMTAB$Characteristics.mycn.status.=ifelse(phenotype_EMTAB$Characteristics.mycn.status.=='amplified',1,0)
roc<-timeROC(T=phenotype_EMTAB$Characteristics.overall.survival.,
             delta=phenotype_EMTAB$Characteristics.overall.survival..1.dead..0.alive.w.o.event.., 
             marker=phenotype_EMTAB$Characteristics.mycn.status., 
             cause=1,weighting="marginal",
             times=c(5),
             iid=TRUE)
# 0.658

x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_MYCN <- data.frame(x,y) 
plotdata_MYCN$group <- "MYCN"
######

### Feng et al
######
rt=phenotype_EMTAB[,c(4,6)]

p=read.csv('./2.Appraisal/E_MATAB_ph (1).csv')
rt=rt[p$Unnamed..0,]
rt$p=p$status_p

roc<-timeROC(T=rt$Characteristics.overall.survival.,
             delta=rt$Characteristics.overall.survival..1.dead..0.alive.w.o.event.., 
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
    'MYCN amplification: 0.658',
    'Zhong et al: 0.43',
    'De Preter et al: 0.727',
    'DL model: 0.896'
  )
  )+
  theme(legend.title=element_blank())+theme(legend.text = element_text(size = 14,face = "bold"))+
  theme(legend.position=c(.766,.25))
tiff('./3.Alternatives/result/ROC_compare_EMTAB.tiff',res = 300,width = 2400,height = 1500,compression = 'lzw')
g
dev.off()

######

score_Zhong=score_Zhong[,p$Unnamed..0]
score_Katleen=score_Katleen[,p$Unnamed..0]
Score_EMTAB=data.frame('Zhong et al'=(score_Zhong),
                          'De Preter et al'=(score_Katleen),
                          'DL model'=rt$p)
write.csv(Score_EMTAB,'./3.Alternatives/result/Score_EMTAB.csv')
