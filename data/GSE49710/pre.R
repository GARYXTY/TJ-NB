rm(list = ls())                                    
options(stringsAsFactors = F)

library(GEOquery)
gset <- getGEO('GSE49710', destdir="./pre/GSE49710/",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)       ## 平台文件
a=gset[[1]] #
dat_GSE49710=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat_GSE49710)#看一下dat_GSE49710这个矩阵的维度
pd_GSE49710=pData(a) #通过查看说明书知道取对象a里的临床信息用pData

gpl<-getGEO('GPL16876', destdir="./pre/GSE49710/",
            AnnotGPL = F,     ## 注释文件
            getGPL = F)       ## 平台文件
ids=Table(gpl)[,c(1,18)]

ids=ids[-which(ids$GeneSymbol==''),]
ids=ids[!duplicated(ids$ID),]
ids=ids[!duplicated(ids$GeneSymbol),]
ids=ids[ids$ID %in% rownames(dat_GSE49710),]
ids$ID=as.character(ids$ID)
dat_GSE49710=dat_GSE49710[ids$ID,] 

rownames(dat_GSE49710)=ids$GeneSymbol

phenotype_GSE49710=pd_GSE49710[,42:50]
supply_pd=read.table('./pre/GSE49710/GSE62564_series_matrix.txt',header = T,sep = '\t',comment.char = '!')
supply_pd=supply_pd[2:7,-1]
rownames(supply_pd)=c('sex','age','efs-days','efs-bin','os-day','os-bin')
phenotype_GSE49710$efs_days=t(supply_pd[3,])
phenotype_GSE49710$efs_bin=t(supply_pd[4,])
phenotype_GSE49710$os_days=t(supply_pd[5,])
phenotype_GSE49710$os_bin=t(supply_pd[6,])

dat_GSE49710=log2(dat_GSE49710+1)


save(dat_GSE49710,phenotype_GSE49710,file = './pre/GSE49710/GSE49710.Rdata')
write.csv(dat_GSE49710,file = './pre/GSE49710/GSE49710_matrix.csv')
write.csv(phenotype_GSE49710,file = './pre/GSE49710/GSE49710_phenotype.csv')
