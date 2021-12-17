#!/usr/bin/bash

usage="extractCutSiteToBw.sh <genome> <inputbam> [<obw_forword> <obw_reverse> | <odir> ]"

## if input params is less than 3, print usage
if [ $# -lt 3 ]; then
    echo $usage
    exit 1
fi

genome=$1
ibam=$2
## if $4 is not empty, then $3 is obw_forword and $4 is obw_reverse
## else $3 is odir
if [ -z $4 ]; then
    odir=$3
    project=`basename $ibam .bam`
    obw_forword=$odir/${project}.cutsite.plus.bw
    obw_reverse=$odir/${project}.cutsite.neg.bw
else
    obw_forword=$3
    obw_reverse=$4
fi

date +'%D %T: Parse parms over'

date +'%D %T: Extract cut sites'
# extract the cut site
plusbed=$ibam.cutsite.plus.bed
bamToBed -i $ibam| awk '$6=="+"' | tr ' ' $'\t' | cut -f 1-3 | grep -e "chr[0-9XYM]\{1,2\}\b" | awk 'BEGIN{OFS="\t"}{$2=$2; $3=$2+1;print $0}'   > $plusbed
bedSort  $plusbed $plusbed
negbed=$ibam.cutsite.neg.bed
bamToBed -i $ibam| awk '$6=="-"' | tr ' ' $'\t' | cut -f 1-3 | grep -e "chr[0-9XYM]\{1,2\}\b" | awk 'BEGIN{OFS="\t"}{$3=$3; $2=$3-1;print $0}' | awk '$2>0'  > $negbed
bedSort  $negbed $negbed


date +'%D %T: Convert bed to bdg'
# turn cut site to bdg
plusbdg=$ibam.cutsite.plus.bdg
bedtools genomecov -i $plusbed -g $genome -bga > $plusbdg
negbdg=$ibam.cutsite.neg.bdg
bedtools genomecov -i $negbed -g $genome -bga > $negbdg
bedSort $plusbdg $plusbdg
bedSort $negbdg $negbdg

date +'%D %T: Convert bdg to bw'
# turn bdg to bw
bedGraphToBigWig $plusbdg $genome $obw_forword
bedGraphToBigWig $negbdg $genome $obw_reverse

date +'%D %T: Rm tmps'
# rm tmp files
rm $plusbed $negbed $plusbdg $negbdg

date +'%D %T: Over!'