# Neoantigen searching ##
project 201905-201908


### pan-cancer differential expression among hervs ###
process:   
1.TCGA case herv exprssion(normalized), ***S1-sdt12_cancer_herv_exp.Rdata*** from paper *Endogenous retroviral signatures predict immunotherapy response in clear cell renal cell carcinoma*  
2.TCGA normal RNA-seq raw data, 127normal normalized herv expression saved in ***S2-127normal_23cancer.normalized.csv***  
3.calculate average expression among samples for each herv to get LOG 2 Foldchange 【LOG(ave_tumor/ave_normal,2)】,p-value from t-test. see ***S3-8678case_127normal.statistic.csv***  
4.Volcano Plot with R script,extract Log2 FoldChange >4 and -log10(p-value) > 50, 146 hervs bed file saved in ***S4-146herv.bed***  



