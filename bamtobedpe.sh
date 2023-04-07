#!/usr/bin/env bash

ibam=$1
threads=$2
bam_tmp="$1".tmp
samtools sort -@ ${threads} -n  $ibam -m 4096M |\
     bamToBed -bedpe -i - |\
      awk '{if($2<$5) print $1"\t"$2"\t"$6; else print $1"\t"$5"\t"$3}' | sort -k1,1 -k2,2n | uniq 
