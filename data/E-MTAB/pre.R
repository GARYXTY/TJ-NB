rm(list = ls())                                    
options(stringsAsFactors = F)

dat_EMTAB=read.table('./pre/E-MTAB/expression_values.txt',header = T,sep = '\t')

library(GEOquery)
gpl<-getGEO('GPL16876', destdir="./pre/GSE49710/",
            AnnotGPL = F,     ## 注释文件
            getGPL = F)       ## 平台文件
ids=Table(gpl)[,c(1,18)]

ids=ids[-which(ids$GeneSymbol==''),]
ids=ids[!duplicated(ids$ID),]
ids=ids[!duplicated(ids$GeneSymbol),]
ids=ids[ids$ID %in% rownames(dat_EMTAB),]
ids$ID=as.character(ids$ID)
dat_EMTAB=dat_EMTAB[ids$ID,] 

rownames(dat_EMTAB)=ids$GeneSymbol
dat_EMTAB=dat_EMTAB[,-1]

pd=read.table('./pre/E-MTAB/E-MTAB-8248.sdrf.txt.txt',header = T,sep = '\t')
pd=pd[order(pd$Source.Name),]
rownames(pd)=pd$Source.Name
pd=pd[,c(4,8,9,10,11,15,19,22,23,24,25)]

dat_EMTAB=log2(dat_EMTAB+1)
# library(preprocessCore)
# dat_EMTAB=normalize.quantiles(as.matrix(dat_EMTAB))
# rownames(dat_EMTAB)=ids$EntrezGeneID
# colnames(dat_EMTAB)=rownames(pd)

save(dat_EMTAB,pd,file = './pre/E-MTAB/E-MTAB.Rdata')
write.csv(dat_EMTAB,file = './pre/E-MTAB/dat_EMTAB.csv')
write.csv(pd,file = './pre/E-MTAB/phenotype_EMTAB.csv')
