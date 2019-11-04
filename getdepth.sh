
path=COAD3/
#rm ${path}*rmdup*
for i in `ls ${path}*.filtered.bam`;do samtools rmdup -sS $i $i.rmdup.bam &&samtools sort $i.rmdup.bam -o $i.rmdup.sorted.bam && samtools index $i.rmdup.sorted.bam;done
ls $path/*.rmdup.sorted.bam >1.tmp
samtools depth -a -b $1 -f 1.tmp >unitdepth.txt
#rm 1.tmp
sed -i -e '1i\file_ID' -e '1i\basepos'  -e 's/COAD3\/\///g' -e 's/Aligned.*.bam//g' -e 's/.rmdup.*.bam//g' 1.tmp
cat 1.tmp |xargs >hdiffherv_region_depth.csv
cat unitdepth.txt >>hdiffherv_region_depth.csv
rm 1.tmp unitdepth.txt
