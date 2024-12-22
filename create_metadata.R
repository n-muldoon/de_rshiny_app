################################################################################
# Author: Anne Muldoon
# Title: create_metadata
# Purpose: R script to create metadata.csv from norm_count data
################################################################################
library(dplyr)
library(stringr)
library(readr)
norm_counts<-read.delim('GSE64403_FPKM_table.txt', sep='\t', header=TRUE)
metadata<-as_tibble(tail(colnames(norm_counts), -3))%>%
  rename('sample_name'=value)%>%
  mutate(description=case_when(str_detect(sample_name,'ESC') ~ "embryonic stem cell",
                               str_detect(sample_name,'MES') ~ "mesoderm",
                               str_detect(sample_name,'CP') ~ "cardiac progenitor",
                               str_detect(sample_name,'CM') ~ "embryonic myocyte",
                               str_starts(sample_name,'ex') ~ "explanted cardiac myocytes",
                               str_starts(sample_name,'iP') ~ "isolated cardiomyocyte from neonatal heart ventrical",
                               str_starts(sample_name,'iD7R') ~ "isolated cardiomyocyte from heart in resected surgery",
                               str_starts(sample_name,'iD7S') ~ "isolated cardiomyocyte from heart in sham surgery",
                               str_starts(sample_name,'vP') ~ "neonatal heart ventrical",
                               str_starts(sample_name,'vAd') ~ "adult heart ventrical",
                               str_detect(sample_name,'R') ~ "heart apex from resected surgery",
                               str_detect(sample_name,'S') ~ "heart apex from sham surgery",
                               ))%>%
  mutate(timepoint=case_when(str_detect(sample_name,'ESC') ~ "day0",
                             str_detect(sample_name,'MES') ~ "day2",
                             str_detect(sample_name,'CP') ~ "day1",
                             str_detect(sample_name,'CM') ~ "day10",
                             str_detect(sample_name,'0hr') ~ "day0",
                             str_detect(sample_name,'24hr') ~ "day1",
                             str_detect(sample_name,'48hr') ~ "day2",
                             str_detect(sample_name,'72hr') ~ "day3",
                             str_detect(sample_name,'P0') ~ "day0",
                             str_detect(sample_name,'P4') ~ "day4",
                             str_detect(sample_name,'P7') ~ "day7",
                             str_detect(sample_name,'D1') ~ "day1",
                             str_detect(sample_name,'D7') ~ "day7",
                             str_detect(sample_name,'Ad') ~ "week8",
                             TRUE ~ NA_character_))%>%
  mutate(replicate=case_when(str_ends(sample_name,'_1') ~ '1',
                             str_ends(sample_name,'_2') ~ '2',
                             str_ends(sample_name,'_3') ~ '3',
                             ))
write_tsv(metadata, file = "metadata.tsv")
# col_names<-col_names%>%mutate(timepoint=case_when(startsWith(sample_name,'vP0') ~ "neonatal heart ventrical",
#                                                     startsWith(sample_name,'i') ~ "isolated_cardiomyocyte",
#                                                     startsWith(sample_name,'v') ~ "vivo"))
#col_names$description<-ifelse(startsWith(sample_name,'i'), 'isolated_cardiomyocyte','vivo')


