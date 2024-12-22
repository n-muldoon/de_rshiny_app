# de_rshiny_app
Differential Expression Analysis Rshiny App

Actual app is the app.R file, raw input file is 'GSE64403_FPKM_table.txt', which is manipulated in create_metadataa.R to create 'metadata.tsv', create_deseq2.R to create 'deseq_results.csv'. 'deseq_results.csv' is alsoo further manipulated in create_fgsea.R to create 'fgsea_res.csv'.

All these input files can be uploaded to the corresping tabs of the app to analyze the data. Differential Expression analysis is performed on the explant0hrs and explant72hrs samples using DESeq2 with the standard threhsolds, and fgsea analysis is conducted using 'CP: CanonicalPathways' gmt files from the GSEA MSigDB mouse collection.
