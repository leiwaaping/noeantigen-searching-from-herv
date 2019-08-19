# Neoantigen searching ##
project 201905-201908


### pan-cancer differential expression among hervs ###
process:   
1.TCGA case herv exprssion(normalized), ***S1-sdt12_cancer_herv_exp.Rdata*** from paper *Endogenous retroviral signatures predict immunotherapy response in clear cell renal cell carcinoma*  
2.TCGA normal RNA-seq raw data, 127normal normalized herv expression saved in ***S2-127normal_23cancer.normalized.csv***  
3.calculate average expression among samples for each herv to get LOG 2 Foldchange 【LOG(ave_tumor/ave_normal,2)】,p-value from t-test. see ***S3-8678case_127normal.statistic.csv***  
4.Volcano Plot with R script,extract Log2 FoldChange > 4 and -log10(p-value) > 50
5.extract target hervs bed file, 146 hervs bed file saved in ***S4-146herv.bed*** 
```
$ head herv_146_2test.csv
herv_ID,logFC,p_value,mwu_p_value
herv_4126,5.7443457133411115,78.6939638953669,6.989572662638197
herv_611,5.603246137835455,309.7772383086611,31.01750729774527
herv_4417,5.147029295939092,71.85619582951088,7.203442282391108
……
$cut -d "," -f1 herv_146_2test.csv >herv_id.txt
$cat ~/workshop/07hervQuant_2/ref/hervquant_final_reference.fa |grep ">"|sed "s/>//g"|awk -v FS='_' '{print $1,$0}' >herv_pos.txt
$awk 'NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' herv_pos.txt herv_id.txt >146herv.bed
$cut -d ":" -f2 146herv.bed |awk -v FS='-' '{print $2-$1}' >1.tmp
```

### site-wides differential expression ###  
data : 34 paired TCGA-COAD , data info in ***COAD_pairedsample.csv***
1.extract depth for each base.
```
$ ls COAD/
18625fe4-3c19-45d9-9d7c-a295fbf83f2eAligned.out.filtered.bam  a4d598d7-896b-495c-9e6b-fef43193f9e8Aligned.out.filtered.bam
1bd24aac-6813-42d2-9f98-232f18e39a75Aligned.out.filtered.bam  afb3c381-bf54-453f-ae30-b0c1fc0e55fcAligned.out.filtered.bam
……
$ ls COAD/|sed 's/.bam//g' >file.list
$ head file.list
18625fe4-3c19-45d9-9d7c-a295fbf83f2eAligned.out.filtered
1bd24aac-6813-42d2-9f98-232f18e39a75Aligned.out.filtered
……
$ cat command.sh
#!/bin/bash
path=COAD/
samtools rmdup -sS $path$1.bam $path$1.rmdup.bam &&
samtools sort $path$1.rmdup.bam -o $path$1.rmdup.sorted.bam &&
samtools index $path$1.rmdup.sorted.bam &&
samtools depth -a $path$1.rmdup.sorted.bam >$path$1.txt &&
sort -nk1 -nk2 $path$1.txt >$path$1.txt.sorted &&
rm $path$1.rmdup.bam
rm $path*.bam.*
cut -d "_" -f1 $path$1.txt.sorted >1.txt
cut -f3 $path$1.txt.sorted >2.txt
cut -f 2 $path$1.txt.sorted |paste 1.txt -|paste - 2.txt >$path$1.txt
rm $path*.sorted.bam
rm 1.txt
rm 2.txt
rm *.sorted
$ for i in `cat file.list`;do bash command.sh ${i};done
[bam_rmdupse_core] 167 / 696 = 0.2399 in library '	'
[bam_rmdupse_core] 67400 / 152308 = 0.4425 in library '	'
[bam_rmdupse_core] 38102 / 144502 = 0.2637 in library '	'
[bam_rmdupse_core] 39408 / 139108 = 0.2833 in library '	'
……
```

2.generated depth file ***S5-peptide_depth.csv*** 【R script】  
3.side wides volcano plot【R script】 

### from site-wide base to peptides ###  
