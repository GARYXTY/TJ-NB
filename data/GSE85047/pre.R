rm(list = ls())                                    
options(stringsAsFactors = F)

library(GEOquery)
gset <- getGEO('GSE85047', destdir="./pre/GSE85047/",
               AnnotGPL = F,     ## 注释文件
               getGPL = F)       ## 平台文件
a=gset[[1]] #
dat_GSE85047=exprs(a) #a现在是一个对象，取a这个对象通过看说明书知道要用exprs这个函数
dim(dat_GSE85047)#看一下dat_GSE85047这个矩阵的维度
pd_GSE85047=pData(a) #通过查看说明书知道取对象a里的临床信息用pData

library(huex10sttranscriptcluster.db)
ids=toTable(huex10sttranscriptclusterSYMBOL)

ids=ids[ids$symbol != '',]
ids=ids[!duplicated(ids$probe_id),]
ids=ids[!duplicated(ids$symbol),]
ids=ids[ids$probe_id %in%  rownames(dat_GSE85047),]

dat_GSE85047[1:4,1:4]   
dat_GSE85047=dat_GSE85047[ids$probe_id,] 

ids$median=apply(dat_GSE85047,1,median) #ids新建median这一列，列名为median，同时对dat_GSE85047这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat_GSE85047=dat_GSE85047[ids$probe_id,] #新的ids取出probe_id这一列，将dat_GSE85047按照取出的这一列中的每一行组成一个新的dat_GSE85047
rownames(dat_GSE85047)=ids$symbol#把ids的symbol这一列中的每一行给dat_GSE85047作为dat_GSE85047的行名

phenotype_GSE85047=pd_GSE85047[,38:44]

dat_GSE85047=log2(dat_GSE85047+1)


save(dat_GSE85047,phenotype_GSE85047,file = './pre/GSE85047/GSE85047.Rdata')
write.csv(dat_GSE85047,file = './pre/GSE85047/GSE85047_matrix.csv')
write.csv(phenotype_GSE85047,file = './pre/GSE85047/GSE85047_phenotype.csv')
