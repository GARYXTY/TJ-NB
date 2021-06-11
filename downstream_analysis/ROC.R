rm(list = ls())                                    
options(stringsAsFactors = F)
library(ggplot2)
### ROC 

plotdata_GSE=t(read.csv('./2.Appraisal/GSE49710_ROC.csv'))
plotdata_GSE=as.data.frame(plotdata_GSE)
colnames(plotdata_GSE)=c('x','y')
plotdata_GSE$group='GSE49710'


plotdata_EMTAB=t(read.csv('./2.Appraisal/E_MATAB_ROC.csv'))
plotdata_EMTAB=as.data.frame(plotdata_EMTAB)
colnames(plotdata_EMTAB)=c('x','y')
plotdata_EMTAB$group='EMTAB'


plotdata=rbind(plotdata_GSE,plotdata_EMTAB)


g <- ggplot(plotdata) + 
  geom_path(aes(x = x, y = y,colour = group), size=1) + 
  labs(x="1 - Specificity", y = "Sensitivity") +
  ggpubr::theme_classic2()+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme(axis.text=element_text(size=6,face="bold"),axis.title=element_text(size=14,face="bold"))+
  geom_abline(intercept=0, slope=1,linetype="dashed")+
  scale_colour_hue(name="my legend", labels=c(
    'EMTAB: 0.891',
    'GSE49710: 0.968'))+
  theme(legend.title=element_blank())+theme(legend.text = element_text(size = 12,face = "bold"))+
  theme(legend.position=c(.766,.25))
tiff('./2.Appraisal/result/ROC.tiff',res = 300,width = 1600,height = 1000,compression = 'lzw')
g
dev.off()
######

### tROC for GSE49710
######
rm(list = ls())
library(survival)
library(timeROC)


rt=read.csv('./pre/GSE49710/GSE49710_phenotype.csv')[,c(13,14)]
rt$p=read.csv('./2.Appraisal/GSE49710_ph.csv')$status_p
rt$os_days=rt$os_days/365

roc<-timeROC(T=rt$os_days,
        delta=rt$os_bin,
        marker=rt$p,
        cause=1,weighting="marginal",
        times=c(3,5,10),
        iid=TRUE)

x <- unlist(roc$FP[,3])  ##提取x值
y <- unlist(roc$TP[,3])
plotdata_10 <- data.frame(x,y) 
plotdata_10$group <- "10"

x <- unlist(roc$FP[,1])  ##提取x值
y <- unlist(roc$TP[,1])
plotdata_3 <- data.frame(x,y) 
plotdata_3$group <- "3"

x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_5 <- data.frame(x,y) 
plotdata_5$group <- "5"

plotdata <- rbind(plotdata_3,plotdata_5,plotdata_10)
plotdata$group=factor(plotdata$group,levels = c(3,5,10))

g <- ggplot(plotdata) + 
  geom_path(aes(x = x, y = y,colour = group), size=1) + 
  labs(x="1 - Specificity", y = "Sensitivity") +
  ggpubr::theme_classic2()+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme(axis.text=element_text(size=6,face="bold"),axis.title=element_text(size=14,face="bold"))+
  geom_abline(intercept=0, slope=1,linetype="dashed")+
  scale_colour_hue(name="my legend", labels=c(
    paste("AUC of 3-year-survival: ",round(roc$AUC[1],3)),
    paste("AUC of 5-year-survival: ",round(roc$AUC[2],3)),
    paste("AUC of 10-year-survival: ",round(roc$AUC[3],3))
    )
    )+
  theme(legend.title=element_blank())+theme(legend.text = element_text(size = 12,face = "bold"))+
  theme(legend.position=c(.7,.2))
tiff('./2.Appraisal/result/tROC_GSE.tiff',res = 300,width = 1600,height = 1000,compression = 'lzw')
g
dev.off()
######

### tROC for EMTAB
######
rm(list = ls())
library(survival)
library(timeROC)


rt=read.csv('./pre/E-MTAB/phenotype_EMTAB.csv')[,c(5,7)]
colnames(rt)=c('os_bin','os_days')
rt$p=read.csv('./2.Appraisal/E_MATAB_ph (1).csv')$status_p
roc<-timeROC(T=rt$os_days,
             delta=rt$os_bin,
             marker=rt$p,
             cause=1,weighting="marginal",
             times=c(3,5,10),
             iid=TRUE)

x <- unlist(roc$FP[,3])  ##提取x值
y <- unlist(roc$TP[,3])
plotdata_10 <- data.frame(x,y) 
plotdata_10$group <- "10"

x <- unlist(roc$FP[,1])  ##提取x值
y <- unlist(roc$TP[,1])
plotdata_3 <- data.frame(x,y) 
plotdata_3$group <- "3"

x <- unlist(roc$FP[,2])  ##提取x值
y <- unlist(roc$TP[,2])
plotdata_5 <- data.frame(x,y) 
plotdata_5$group <- "5"

plotdata <- rbind(plotdata_3,plotdata_5,plotdata_10)
plotdata$group=factor(plotdata$group,levels = c(3,5,10))
g <- ggplot(plotdata) + 
  geom_path(aes(x = x, y = y,colour = group), size=1) + 
  labs(x="1 - Specificity", y = "Sensitivity") +
  ggpubr::theme_classic2()+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme(axis.text=element_text(size=6,face="bold"),axis.title=element_text(size=14,face="bold"))+
  geom_abline(intercept=0, slope=1,linetype="dashed")+
  scale_colour_hue(name="my legend", labels=c(
    paste("AUC of 3-year-survival: ",round(roc$AUC[1],3)),
    paste("AUC of 5-year-survival: ",round(roc$AUC[2],3)),
    paste("AUC of 10-year-survival: ",round(roc$AUC[3],3))
    )
    )+
  theme(legend.title=element_blank())+theme(legend.text = element_text(size = 12,face = "bold"))+
  theme(legend.position=c(.7,.2))
tiff('./2.Appraisal/result/tROC_EMTAB.tiff',res = 300,width = 1600,height = 1000,compression = 'lzw')
g
dev.off()
######