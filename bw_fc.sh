#!/usr/bin/env bash

bw=$1
chrom=$2
obw=$3

if [[ -n $4 ]];then
    scale=$4;
else
    scale=1;
fi;

project=`basename -s .bw $bw`
bdg="${project}.bdg"
bdg_tmp="${bdg}.tmp"
if [[ ! -n "${obw}" ]];then
obw="${project}.fc.bw"
fi

ave=`bigWigInfo $1 | grep mean | cut -d ":" -f 2 `

echo $ave,$scale
bigWigToBedGraph $bw $bdg
echo "generate bdg.tmp..."
awk '{$4=$4/a*b;print $0}' a="${ave}"  b="${scale}" $bdg  | sort -k1,1 -k2,2n > $bdg_tmp
echo "generate new bigwig .."
bedGraphToBigWig $bdg_tmp $chrom $obw
echo "process over, clean up!"
rm -f $bdg
rm -f $bdg_tmp
