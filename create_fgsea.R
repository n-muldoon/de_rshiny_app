################################################################################
# Author: Anne Muldoon
# Title: create_fgsea
# Purpose: R script to create fgsea_results.csv from deseq_2 data
# DE: explant 0hr vs explant 72hr
# Thresholds: padj<0.05
################################################################################
library(dplyr)
library(readr)
library(tidyverse)
library(biomaRt)
library(fgsea)
#Read in deseq_results.csv, filter based on threshold
deseq_res_csv<-read.csv('deseq_results.csv',header=TRUE)
deseq_res<-as.data.frame(deseq_res_csv)%>%
  arrange(desc(log2FoldChange)) %>%
  filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
  dplyr::select(GeneID, log2FoldChange, padj)%>%
  distinct(GeneID, .keep_all = TRUE)

ranks<-deseq_res
ranks <- deseq_res$log2FoldChange
names(ranks) <- deseq_res$GeneID
ranks <- ranks[is.finite(ranks)]

gene_sets <- gmtPathways("m2.cp.v2024.1.Mm.symbols.gmt")
fgsea_results <- fgsea(pathways = gene_sets,
                       stats = ranks_up,
                       minSize = 15,
                       maxSize = 500)
fgsea_res<-as.data.frame(fgsea_results)


# explant_up<-as.data.frame(up_reg$Gene)
# explant_down<-as.data.frame(down_reg$Gene)
# 
# write_tsv(explant_up, file = "explant_up_genes.csv")
# write_tsv(explant_down, file = "explant_down_genes.csv")
# #After went to DAVID, downloaded cluster data and annotation chart for explant_up and explant_down
# 
# #need to combine into one doc for upload
# up_res<-read.delim('david_up_go_res.tsv',header=TRUE,sep='\t')
# up_res<-as_tibble(up_res)%>%mutate('Regulation'='up')
# 
# down_res<-read.delim('david_down_go_res.tsv',header=TRUE,sep='\t')
# down_res<-as_tibble(down_res)%>%mutate('Regulation'='down')
# fgsea_res<-bind_rows(up_res,down_res)
write_csv(fgsea_res, file = "fgsea_res.csv")
