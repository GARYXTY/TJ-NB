rm(list = ls())
options(stringsAsFactors = F)

library(rjson)

all_files=dir("./result2/")

### AUC
######
result_auc=matrix(nrow = 200,ncol = 4)

for (i in 1:200) {
  json_file_name=paste("./result2/",all_files[i],sep="")
  json_file=fromJSON(file=json_file_name)
  result_auc[i,1]=json_file$train_auc
  result_auc[i,2]=json_file$test_auc
  result_auc[i,3]=json_file$val_auc
  result_auc[i,4]=json_file$dataset2_auc
}

colnames(result_auc)=c("train_auc","test_auc","validation_auc","External_validation")

mean(result_auc[,1]);sd(result_auc[,1])
mean(result_auc[,2]);sd(result_auc[,2])
mean(result_auc[,3]);sd(result_auc[,3])
mean(result_auc[,4]);sd(result_auc[,4])
######

library(survival)
library(timeROC)

### train t-ROC AUC
######
rm(list = ls())
options(stringsAsFactors = F)

all_files=dir("./result2/")
train_t_auc=matrix(nrow = 200,ncol = 3)

for (i in 1:200) {
  json_file_name=paste("./result2/",all_files[i],sep="")
  json_file=fromJSON(file=json_file_name)
  
  rt=read.csv('./pre/GSE49710/GSE49710_phenotype.csv')[,c(1,13,14)]
  rt$os_days=rt$os_days/365
  rownames(rt)=rt$X
  
  x=unlist(json_file$dataset1_result$train_result)
  rt=rt[names(x),]
  rt$p=x
  
  ROC<-timeROC(T=rt$os_days,
               delta=rt$os_bin,
               marker=rt$p,
               cause=1,weighting="marginal",
               times=c(3,5,10),
               iid=TRUE)
  
  tAUC=ROC$AUC
  
  train_t_auc[i,1]=tAUC[1]
  train_t_auc[i,2]=tAUC[2]
  train_t_auc[i,3]=tAUC[3]
}
mean(train_t_auc[,1]);sd(train_t_auc[,1])
mean(train_t_auc[,2]);sd(train_t_auc[,2])
mean(train_t_auc[,3]);sd(train_t_auc[,3])

### test t-ROC AUC
######
rm(list = ls())
options(stringsAsFactors = F)

all_files=dir("./result2/")
test_t_auc=matrix(nrow = 200,ncol = 3)

for (i in 1:200) {
  json_file_name=paste("./result2/",all_files[i],sep="")
  json_file=fromJSON(file=json_file_name)
  
  rt=read.csv('./pre/GSE49710/GSE49710_phenotype.csv')[,c(1,13,14)]
  rt$os_days=rt$os_days/365
  rownames(rt)=rt$X
  
  x=unlist(json_file$dataset1_result$test_result)
  rt=rt[names(x),]
  rt$p=x
  
  ROC<-timeROC(T=rt$os_days,
               delta=rt$os_bin,
               marker=rt$p,
               cause=1,weighting="marginal",
               times=c(3,5,10),
               iid=TRUE)
  
  tAUC=ROC$AUC
  
  test_t_auc[i,1]=tAUC[1]
  test_t_auc[i,2]=tAUC[2]
  test_t_auc[i,3]=tAUC[3]
}
mean(test_t_auc[,1]);sd(test_t_auc[,1])
mean(test_t_auc[,2]);sd(test_t_auc[,2])
mean(test_t_auc[,3]);sd(test_t_auc[,3])

### val t-ROC AUC
######
rm(list = ls())
options(stringsAsFactors = F)

all_files=dir("./result2/")
val_t_auc=matrix(nrow = 200,ncol = 3)

for (i in c(1:69,71:86,88:100,102:121,123:200)) {
  json_file_name=paste("./result2/",all_files[i],sep="")
  json_file=fromJSON(file=json_file_name)
  
  rt=read.csv('./pre/GSE49710/GSE49710_phenotype.csv')[,c(1,13,14)]
  rt$os_days=rt$os_days/365
  rownames(rt)=rt$X
  
  x=unlist(json_file$dataset1_result$val_result)
  rt=rt[names(x),]
  rt$p=x
  
  ROC<-timeROC(T=rt$os_days,
               delta=rt$os_bin,
               marker=rt$p,
               cause=1,weighting="marginal",
               times=c(3,5,10),
               iid=TRUE)
  
  tAUC=ROC$AUC
  
  val_t_auc[i,1]=tAUC[1]
  val_t_auc[i,2]=tAUC[2]
  val_t_auc[i,3]=tAUC[3]
}
val_t_auc=val_t_auc[-c(70,87,101,122),]
val_t_auc=val_t_auc[-c(81,118),]
mean(val_t_auc[,1]);sd(val_t_auc[,1])
mean(val_t_auc[,2]);sd(val_t_auc[,2])
mean(val_t_auc[,3]);sd(val_t_auc[,3])

### EMTAB t-ROC AUC
######
rm(list = ls())
options(stringsAsFactors = F)

all_files=dir("./result2/")
EMTAB_t_auc=matrix(nrow = 200,ncol = 3)

for (i in 1:200) {
  json_file_name=paste("./result2/",all_files[i],sep="")
  json_file=fromJSON(file=json_file_name)
  
  rt=read.csv('./pre/E-MTAB/phenotype_EMTAB.csv')[,c(1,5,7)]
  colnames(rt)=c('X','os_bin','os_days')
  rownames(rt)=rt$X
  
  x=unlist(json_file$dataset2_result)
  rt=rt[gsub(".","-",names(x),fixed = T),]
  rt$p=x
  
  ROC<-timeROC(T=rt$os_days,
               delta=rt$os_bin,
               marker=rt$p,
               cause=1,weighting="marginal",
               times=c(3,5,10),
               iid=TRUE)
  
  tAUC=ROC$AUC
  
  EMTAB_t_auc[i,1]=tAUC[1]
  EMTAB_t_auc[i,2]=tAUC[2]
  EMTAB_t_auc[i,3]=tAUC[3]
}
mean(EMTAB_t_auc[,1]);sd(EMTAB_t_auc[,1])
mean(EMTAB_t_auc[,2]);sd(EMTAB_t_auc[,2])
mean(EMTAB_t_auc[,3]);sd(EMTAB_t_auc[,3])
