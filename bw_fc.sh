#!/usr/bin/env bash

bw=$1
chrom=$2
obw=$3

project=`basename -s .bw $bw`
bdg="${project}.bdg"
bdg_tmp="${bdg}.tmp"
if [[ ! -n "${obw}" ]];then
obw="${project}.fc.bw"
fi
ave=`bigWigInfo $1 | grep mean | cut -d ":" -f 2 `
echo $ave
bigWigToBedGraph $bw $bdg
awk '{$4=$4/a;print $0}' a="${ave}" $bdg  | sort -k1,1 -k2,2n > $bdg_tmp
bedGraphToBigWig $bdg_tmp $chrom $obw

rm -f $bdg
rm -f $bdg_tmp
