#!/bin/bash
path=/alldata/HPli_data/hervquant/COAD_normal/

ls ${path}/*/quant.sf |cut -d "/" -f7 >file.list
sample=`head -1 file.list`
cat ${path}/${sample}/quant.sf |sort -t "_" -nk1|grep -v "Name" |cut -f 1|tr "\n" ","|sed -e 's/$/\n/'|sed "s/^/Sample_ID,#total_reads,&/"| sed "s/,$/\n/" >quant_normal.csv
for i in `cat file.list`; do nread=`cat ${path}/${i}/*.final.out|head |grep "Number of input reads" |cut -d "|" -f2 |sed 's/\t//g'` ; cat ${path}/${i}/quant.sf |sort -t "_" -nk1 |cut -f 5|grep -v "Num"|tr "\n" ","|sed -e 's/ $/\n/'|sed "s/^/${i},${nread},&/"|sed "s/,$/\n/" >>quant_normal.csv;done
sed -i '/^$/d' quant_normal.csv
