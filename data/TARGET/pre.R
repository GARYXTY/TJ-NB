rm(list = ls())
options(stringsAsFactors = F)

# read gdc download files
raw_count_info<-read.table("./pre/TARGET/gdc_sample_sheet.2020-01-29.tsv",sep = "\t",header = T)
raw_count=matrix(nrow =60483,ncol = 161 )
for (i in 1:161) {
  gf=gzfile(paste("./pre/TARGET/gdc/",raw_count_info$File.ID[i],"/",raw_count_info$File.Name[i],sep = ""),"rt")
  data=read.table(gf,sep = "\t",comment.char = "_")
  raw_count[,i]=data[,2]
}
rownames(raw_count)=data$V1
colnames(raw_count)=raw_count_info$Sample.ID
colnames(raw_count)=substr(colnames(raw_count),1,19)
rownames(raw_count)=substr(rownames(raw_count),1,15)

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
raw_count <- apply(raw_count,2,fpkmToTpm)
colSums(raw_count)

library(rtracklayer)
AnnoData = import('D:/Homo_sapiens.GRCh38.99.gtf')
index = which(AnnoData$type == 'gene')

########mRNA
Target = data.frame(Ensembl_ID = AnnoData$gene_id[index], Symbol = AnnoData$gene_name[index], Biotype = AnnoData$gene_biotype[index])
Target=Target[Target$Biotype=='protein_coding',]
# make a intersection between RNA-seq and gencode
common = intersect(Target$Ensembl_ID, rownames(raw_count))

mRNA<- raw_count[common,]
# 19600

# if needed, transform id from Ensemble into Symbol
ids<- Target[Target$Ensembl_ID %in% rownames(mRNA),]

# find the duplicated genes
table(duplicated(ids$Symbol))

ids$median=apply(mRNA,1,median) #ids新建median这一列，列名为median，同时对mRNA这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$Symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$Symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
mRNA=mRNA[rownames(mRNA) %in% ids$Ensembl_ID,] #新的ids取出Ensembl_ID这一列，将mRNA按照取出的这一列中的每一行组成一个新的mRNA
rownames(mRNA)=ids$Symbol#把ids的symbol这一列中的每一行给mRNA作为mRNA的行名
# 19597

save(mRNA,file = './pre/TARGET/TCGA_mRNA.Rdata')
write.csv(mRNA,file = './pre/TARGET/TCGA_mRNA.csv')
#########


########lncRNA
Target = data.frame(Ensembl_ID = AnnoData$gene_id[index], Symbol = AnnoData$gene_name[index], Biotype = AnnoData$gene_biotype[index])
Target=Target[Target$Biotype=='lncRNA',]
# make a intersection between RNA-seq and gencode
common = intersect(Target$Ensembl_ID, rownames(raw_count))

lncRNA<- raw_count[common,]
# 14083

# if needed, transform id from Ensemble into Symbol
ids<- Target[Target$Ensembl_ID %in% rownames(lncRNA),]

# find the duplicated genes
table(duplicated(ids$Symbol))

ids$median=apply(lncRNA,1,median) #ids新建median这一列，列名为median，同时对mRNA这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$Symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$Symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
lncRNA=lncRNA[rownames(lncRNA) %in% ids$Ensembl_ID,] #新的ids取出Ensembl_ID这一列，将mRNA按照取出的这一列中的每一行组成一个新的mRNA
rownames(lncRNA)=ids$Symbol#把ids的symbol这一列中的每一行给mRNA作为mRNA的行名
# 14080

save(lncRNA,file = './pre/TARGET/TCGA_lncRNA.Rdata')
write.csv(mRNA,file = './pre/TARGET/TCGA_lncRNA.csv')
#########

########miRNA
Target = data.frame(Ensembl_ID = AnnoData$gene_id[index], Symbol = AnnoData$gene_name[index], Biotype = AnnoData$gene_biotype[index])
Target=Target[Target$Biotype=='miRNA',]
# make a intersection between RNA-seq and gencode
common = intersect(Target$Ensembl_ID, rownames(raw_count))

miRNA<- raw_count[common,]
# 1448

# if needed, transform id from Ensemble into Symbol
ids<- Target[Target$Ensembl_ID %in% rownames(miRNA),]

# find the duplicated genes
table(duplicated(ids$Symbol))

ids$median=apply(miRNA,1,median) #ids新建median这一列，列名为median，同时对mRNA这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$Symbol,ids$median,decreasing = T),]#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$Symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
miRNA=miRNA[rownames(miRNA) %in% ids$Ensembl_ID,] #新的ids取出Ensembl_ID这一列，将mRNA按照取出的这一列中的每一行组成一个新的mRNA
rownames(miRNA)=ids$Symbol#把ids的symbol这一列中的每一行给mRNA作为mRNA的行名
# 1448

save(lncRNA,file = './pre/TARGET/TCGA_miRNA.Rdata')
write.csv(mRNA,file = './pre/TARGET/TCGA_miRNA.csv')
#########

phenotype=read.table('./pre/TARGET/clinical.tsv',header = T,sep = '\t')


phenotype$days_to_last_follow_up[is.na(phenotype$days_to_last_follow_up)]=0
phenotype$days=phenotype$days_to_last_follow_up
phenotype=data.frame(phenotype$submitter_id,phenotype$gender,
                    phenotype$race,phenotype$vital_status,
                    phenotype$age_at_diagnosis,phenotype$inss_stage,
                    phenotype$tissue_or_organ_of_origin,phenotype$days)
table(duplicated(substr(colnames(mRNA),1,16)))
intersect(substr(colnames(mRNA),1,16),phenotype$phenotype.submitter_id)

phenotype_TARGET=data.frame()

for (i in 1:161) {
  phenotype_TARGET=rbind(phenotype_TARGET,phenotype[phenotype$phenotype.submitter_id==substr(colnames(mRNA),1,16)[i],])
}


save(phenotype,file = './pre/TARGET/phenotype.Rdata')
write.csv(phenotype,file = './pre/TARGET/phenotype.csv')







