rm(list=ls())
options(stringsAsFactors = F)


library(survival)
library(survminer)
library(forestplot)
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
rt$age=factor(ifelse(rt$age>1.5,'older','younger'),levels = c('younger','older'))
rt$MYCN=factor(ifelse(rt$MYCN==0,'Non-amplified','Amplified'),levels = c('Non-amplified','Amplified'))
rt$gender=as.factor(rt$gender)
rt$stage=factor(ifelse(rt$stage==1| rt$stage==2,'1/2','3/4/4S'),levels = c('1/2','3/4/4S'))
rt$Risk=factor(ifelse(rt$Risk==1,'High','Low'),levels = c('Low','High'))
rt$Probability=factor(ifelse(rt$Probability<0.5,'Low','High'),levels = c('Low','High'))


outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  HR_mean=round(coxSummary$conf.int[,"exp(coef)"],3)
  HR.95L=round(coxSummary$conf.int[,"lower .95"],3)
  HR.95H=round(coxSummary$conf.int[,"upper .95"],3)
  outTab=rbind(outTab,
               cbind(characteristics=rownames(coxSummary$coefficients),
                     HR_mean=round(coxSummary$conf.int[,"exp(coef)"],3),
                     HR.95L=round(coxSummary$conf.int[,"lower .95"],3),
                     HR.95H=round(coxSummary$conf.int[,"upper .95"],3),
                     HR_CL=paste(HR_mean,'(',HR.95L,'-',HR.95H,')',sep = ''),
                     pvalue=round(coxSummary$coefficients[,"Pr(>|z|)"],3))
  )
  
}



outTab$characteristics=c('Age (>1.5 vs ≤1.5 years)',
                         'MYCN (Amplified vs Non-amplified)',
                         'Gender (Male vs Female)',
                         'Stage (3/4/4S vs 1/2)',
                         'Risk (High vs Low)',
                         'Probability (High vs Low)')

for (i in colnames(outTab)[c(2,3,4,6)]) {
  outTab[,i]=as.numeric(outTab[,i])
}

outTab$pvalue[c(1,2,4,5,6)]='<0.001'



tabletext <- cbind(c("\nUnivariate analysis",NA, outTab$characteristics),
                   c("Hazard Ratio\n(95% CI)", NA, outTab$HR_CL),
                   c("P-value", NA, outTab$pvalue))


tiff('./5.clinic/result/clin_uni.tiff',res = 300, width = 4000, height = 1200, compression = "lzw")
forestplot(labeltext=tabletext, #图中的文本
           mean=c(NA,1,outTab$HR_mean),#HR
           lower=c(NA,1,outTab$HR.95L), #95%置信区间下限
           upper=c(NA,1,outTab$HR.95H),#95%置信区间上限
           #title="Hazard Ratio",
           graph.pos=3,#图在表中的列位置
           graphwidth = unit(.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           boxsize = 0.5,
           col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           #虚线(可多条)及其横坐标、颜色、线宽
           #xticks = c(round(min(data$mean),1), 1,round((max(data$mean)-min(data$mean)/2),1), round(max(data$mean),1)), #横坐标刻度根据具体情况设置
           lwd.xaxis=2, #X轴线宽
           xlab="Hazard Ratio",#X轴标题
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#第三行顶部加黑线，引号内数字标记行位置
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#加阴影，弱项不建议使用
                           "9" = gpar(lwd=2, col="black")),#最后一行底部加黑线,""中数字为nrow(data)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#各种字体大小设置
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           clip =c(0, 30),
           #is.summary = c(T,rep(F,27)),#首行字体类型设置
           lineheight = unit(1,"cm"),#固定行高
           #align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
           #cex=10, 
           colgap = unit(1,"cm"),#列间隙  
           mar=unit(rep(1.25, times = 4), "cm"),#图形页边距
           new_page = F#是否新页
)
dev.off()


write.csv(outTab,'./5.clinic/result/clin-uni.csv')


rt=rt[,-5]
multiCox <- coxph(Surv(futime, fustat) ~ ., data = rt)
multisum=summary(multiCox)
multisum

HR_mean=round(multisum$conf.int[,"exp(coef)"],3)
HR.95L=round(multisum$conf.int[,"lower .95"],3)
HR.95H=round(multisum$conf.int[,"upper .95"],3)
multiout=cbind(characteristics=c('Age (>1.5 vs ≤1.5 years)',
                                 'MYCN (Amplified vs Non-amplified)',
                                 'Stage (3/4/4S vs 1/2)',
                                 'Risk (High vs Low)',
                                 'Probability (High vs Low) '),
               HR_mean=round(multisum$conf.int[,"exp(coef)"],3),
               HR.95L=round(multisum$conf.int[,"lower .95"],3),
               HR.95H=round(multisum$conf.int[,"upper .95"],3),
               HR_CL=paste(HR_mean,'(',HR.95L,'-',HR.95H,')',sep = ''),
               pvalue=round(multisum$coefficients[,"Pr(>|z|)"],3))
multiout=as.data.frame(multiout)
multiout$pvalue[c(2,4,5)]='<0.001'
for (i in colnames(multiout)[c(2,3,4)]) {
  multiout[,i]=as.numeric(multiout[,i])
}

tabletext <- cbind(c("\nMultivariate analysis",NA, multiout$characteristics),
                   c("Hazard Ratio\n(95% CI)", NA, multiout$HR_CL),
                   c("P-value", NA, multiout$pvalue))


tiff('./5.clinic/result/clin_mul.tiff',res = 300, width = 4000, height = 2400, compression = "lzw")
forestplot(labeltext=tabletext, #图中的文本
           mean=c(NA,1,multiout$HR_mean),#HR
           lower=c(NA,1,multiout$HR.95L), #95%置信区间下限
           upper=c(NA,1,multiout$HR.95H),#95%置信区间上限
           #title="Hazard Ratio",
           graph.pos=3,#图在表中的列位置
           graphwidth = unit(.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           boxsize = 0.5,
           col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           #虚线(可多条)及其横坐标、颜色、线宽
           #xticks = c(round(min(data$mean),1), 1,round((max(data$mean)-min(data$mean)/2),1), round(max(data$mean),1)), #横坐标刻度根据具体情况设置
           lwd.xaxis=2, #X轴线宽
           xlab="Hazard Ratio",#X轴标题
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#第三行顶部加黑线，引号内数字标记行位置
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#加阴影，弱项不建议使用
                           "8" = gpar(lwd=2, col="black")),#最后一行底部加黑线,""中数字为nrow(data)+5,   #####注意修改#######
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#各种字体大小设置
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           clip =c(0, 15),
           #is.summary = c(T,rep(F,27)),#首行字体类型设置
           lineheight = unit(1,"cm"),#固定行高
           #align=c("l","c","c"),#每列文字的对齐方式，偶尔会用到
           #cex=10, 
           colgap = unit(1,"cm"),#列间隙  
           mar=unit(rep(1.25, times = 4), "cm"),#图形页边距
           new_page = F#是否新页
)
dev.off()

write.csv(multiout,'5.clinic/result/clin-multi.csv')