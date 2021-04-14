#!usr/bin/env bash

name=$1
project=`basename $name`
total=`cat $name | sed '1d' |awk 'BEGIN{sum=0}{sum+=$NF}END{print sum}'`
cat $name | cut -f 1,2,3,4,5,9 | sed '1d'|  \
    awk -v total=$total 'BEGIN{OFS="\t"}{print $2,$3,$4,$1,$6/total*1000000,$5}' | \
    grep -v "chr[0-9XYM]\{1,2\}_" | sort -k1,1 -k2,2n | \
    bedtools intersect -a - -b /mnt/Storage2/home/zengshiyang/DB/refGene/hg38/hg38.repeats.main.transcript.bed -wao | \
    awk 'BEGIN{OFS="\t"}($3==$9 && $2-1==$8){$11=$5;print $0;}' | \
    cut -f 7-12 