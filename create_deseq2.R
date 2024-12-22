################################################################################
# Author: Anne Muldoon
# Title: create_deseq2
# Purpose: R script to create deseq_results.csv from norm_count data
# DE: explant 0hr vs explant 72hr
################################################################################
library(dplyr)
library(stringr)
library(readr)
library(tidyverse)
library(DESeq2)

norm_counts<-read.delim('GSE64403_FPKM_table.txt',sep='\t',header=TRUE)
gene_to_geneid<-norm_counts%>%select(Gene,GeneID)
norm_counts<-as_tibble(norm_counts)%>%
  select(-Coordinates)%>%
  select(Gene,ex0hr_1,ex0hr_2,ex72hr_1,ex72hr_2)
#filter out zero variance
non_zero_var <- apply(norm_counts[-1],MARGIN=1,var)!=0
counts<-norm_counts[non_zero_var,]

#make counts a matrix
counts<-as.data.frame(counts)%>%column_to_rownames(var='Gene')%>%as.matrix()

#create rowdata
rowdata<-rownames(counts)
#filter metadata
metadata<-read.delim('metadata.tsv',header=TRUE,sep='\t')
coldata<-as_tibble(metadata)%>%
  filter(sample_name=='ex0hr_1'|sample_name=='ex0hr_2'|
           sample_name=='ex72hr_1'|sample_name=='ex72hr_2')%>%
  select(sample_name,timepoint)
#make SE object
se<-SummarizedExperiment(assays=list(counts=round(counts)),colData=coldata,rowData=rowdata)
print(se)
#run DESEQ2
dds <- DESeqDataSet(se,design= ~ timepoint)
dds <- DESeq(dds)
res<-results(dds)
res_tbl<-as_tibble(res)%>%mutate(Gene=rownames(res)) %>%arrange(pvalue)



#join with 'GeneID' from metadata
res_tbl<-res_tbl%>%left_join(gene_to_geneid)
res_df<-as.data.frame(res_tbl)

write.csv(res_df,"deseq_results.csv", row.names = FALSE)





