# Neoantigen searching ##
project 201905-201908

### calculating herv expression level ###
```
#running hervQuant with docker ytding/hervQuant
$ docker run -it --rm --name hpli_hervq --cpus 30 -m 50G -u $(id -u) -v /alldata/HPli_data/data_normal:/home/project/input -v /home/hpli/workshop:/home/project/ref -v/alldata/HPli_data/hervquant:/home/project/output ytding/hervquant /bin/bash
$ bash hervquant_docker.sh |tee ./out.log

#get quant.sf data
$ cat 29herv_analysis/extractquant.sh
#!/bin/bash
path=/alldata/HPli_data/hervquant/

ls ${path}/*/quant.sf |cut -d "/" -f6 >file.list
sample=`head -1 file.list`
cat ${path}/${sample}/quant.sf |sort -t "_" -nk1|grep -v "Name" |cut -f 1|tr "\n" ","|sed -e 's/$/\n/'|sed "s/^/Sample_ID,#total_reads,&/"|sed "s/,$/\n/" >quant_normal.csv
for i in `cat file.list`; do nread=`cat ${path}/${i}/*.final.out|head |grep "Number of input reads" |cut -d "|" -f2 |sed 's/\t//g'` ; cat ${path}/${i}/quant.sf |sort -t "_" -nk1 |cut -f 5|grep -v "Num"|tr "\n" ","|sed -e 's/ $/\n/'|sed "s/^/${i},${nread},&/"|sed "s/,$/\n/" >>quant_normal.csv;done
sed -i '/^$/d' quant_normal.csv
```

### pan-cancer differential expression amoung hervs ###  
1.TCGA case herv exprssion(normalized), ***TCGA_tumor_herv_expr.Rdata*** from paper *Endogenous retroviral signatures predict immunotherapy response in clear cell renal cell carcinoma*  
2.TCGA normal RNA-seq raw data, 534 normal normalized herv expression saved in ***normalized_normal_534.csv***  
3.calculate average expression among samples for each herv to get LOG 2 Foldchange 【LOG(ave_tumor/ave_normal,2)】,p-value from t-test. see ***pancancer_8470_534_herv_level_statistic.csv***  
4.Volcano Plot with R script,extract 66hervs (Log2 FoldChange > 4 and -log10(p-value) > 50 ) ***pancancer_highlydiff_66herv.csv*** 
5.extract target hervs bed file, 66 hervs bed file saved in ***pancancer_highlydiff_66herv.csv*** 
6.get bed file of 66 hervs ***highdiff_herv.bed***
```
$ head pancancer_highlydiff_66herv.csv
hervID,avg_tumor,avg_normal,logFC,pvalue,threshold
herv_4227,0.00899683679438682,0,10,1.22159453566286e-64,Up
herv_3786,0.0180878425622008,0.000211703443820225,6.41683179158543,4.72341384204794e-128,Up
herv_2741,0.00645892050894506,0.000270108529962547,4.57968195919543,6.53361720174555e-79,Up
herv_1090,0.0083225867544987,0.000392748558052434,4.40535412017913,2.34426808526767e-68,Up
……
#extract bed files
$ cat ../14readdepth/generateBed.sh
#! /bin/bash

cut -d "," -f1 $1 |sed -e 's/herv_//g' -e '/^$/d' >herv_id.tmp
cat ~/workshop/07hervQuant_2/ref/hervquant_final_reference.fa |grep ">"|sed "s/>//g"|awk -v FS='_' '{print $1,$0}' >herv_pos.tmp
awk 'NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' herv_pos.tmp herv_id.tmp >diffherv.tmp
cut -d ":" -f2 diffherv.tmp |awk -v FS='-' '{print 0","$2-$1}' >1.tmp
paste -d "," diffherv.tmp 1.tmp >highdiff_herv.bed
sed -i -e "s/\s//g" -e "s/,/\t/g" highdiff_herv.bed
rm *.tmp

```
7.statistical analysis and volcano plot(herv_level) in ***R script***

### site-wides differential expression ###  
data : 34 paired TCGA-COAD , data info in ***COAD_pairedsamples.csv***
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
$ cat 14readdepth/getdepth.sh

path=COAD3/
#rm ${path}*rmdup*
#for i in `ls ${path}*.filtered.bam`;do samtools rmdup -sS $i $i.rmdup.bam &&samtools sort $i.rmdup.bam -o $i.rmdup.sorted.bam && samtools index $i.rmdup.sorted.bam;done
#ls $path/*.rmdup.sorted.bam >1.tmp
samtools depth -a -b highdiff_herv.bed -f 1.tmp >unitdepth.txt
#rm 1.tmp
sed -i -e '1i\file_ID' -e '1i\basepos'  -e 's/COAD3\/\///g' -e 's/Aligned.*.bam//g' -e 's/.rmdup.*.bam//g' 1.tmp
cat 1.tmp |xargs >hdiffherv_region_depth.csv
cat unitdepth.txt >>hdiffherv_region_depth.csv
#rm 1.tmp unitdepth.txt
```

#add header to depth files，the unitdepth.txt only calculate the regions in S4-146herv.bed file.

$ bash command.sh
[bam_rmdupse_core] 167 / 696 = 0.2399 in library '	'
[bam_rmdupse_core] 67400 / 152308 = 0.4425 in library '	'
[bam_rmdupse_core] 38102 / 144502 = 0.2637 in library '	'
[bam_rmdupse_core] 39408 / 139108 = 0.2833 in library '	'
……
```

2.generated depth file ***S5-peptide_depth.csv*** 【R script】  
3.side wides volcano plot【R script】 

### from site-wide base to peptides ###  【python3 script】
here we have a human exciting peptides dataset from uniport (unit all possible human protein, covert them back to peptides chains and use sliding widows to split into 7-11 aa chains. then those peetides are supporsed to be the human normally expressed peptide). In this research we will exclude all those human-noramlly-expr site. The dataset info saved in ***S8-exclude_position.pkl***  

