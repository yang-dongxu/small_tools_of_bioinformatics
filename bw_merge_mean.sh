#!/usr/bin/env bash
usage="usage: $0 out.bw num chromsize [args for bigWigMerge(ucsc tools)]"

outbw=$1
shift
num=$1
shift
chromsize=$1
shift
if [ -n "$1" ]
then
    other_args=$*;
else
    echo $usage
    echo "No enough args input, check bigWigMerge "
    exit;
fi;

project=`basename outbw`
bdg_tmp=${project}.sum.bdg.tmp
bdg_mean=${project}.mean.bdg.tmp 
echo "##### merge and get bdg: `date +'%D %T'` #####"
bigWigMerge $other_args $bdg_tmp 
echo "##### calulate mean for each site: `date +'%D %T'` #####"
cat $bdg_tmp | awk -v num=$num 'BEGIN{OFS="\t"}{$4=$4/num; print $0}' | sort -k1,1 -k2,2n > $bdg_mean
echo "##### turn bdg to bigwig: `date +'%D %T'` #####"
bedGraphToBigWig $bdg_mean $chromsize $outbw
echo "##### rm tmps: `date +'%D %T'` #####"
rm -f $bdg_mean
rm -f $bdg_tmp
echo "##### processed over: `date +'%D %T'` #####"
