#! /bin/bash

cut -d "," -f1 $1 |sed -e 's/herv_//g' -e '/^$/d' >herv_id.tmp
cat ~/workshop/07hervQuant_2/ref/hervquant_final_reference.fa |grep ">"|sed "s/>//g"|awk -v FS='_' '{print $1,$0}' >herv_pos.tmp
awk 'NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' herv_pos.tmp herv_id.tmp >diffherv.tmp
cut -d ":" -f2 diffherv.tmp |awk -v FS='-' '{print 0","$2-$1}' >1.tmp
paste -d "," diffherv.tmp 1.tmp >highdiff_herv.bed
sed -i -e "s/\s//g" -e "s/,/\t/g" highdiff_herv.bed
rm *.tmp
