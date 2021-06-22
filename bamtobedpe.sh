#!/usr/bin/env bash

ibam=$1
bam_tmp="$1".tmp
samtools sort -@ 10 -n  -o $bam_tmp $ibam 
bamToBed -bedpe -i ${bam_tmp} | awk '{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}' | sort -k1,1 -k2,2n | uniq 

rm ${bam_tmp}
