#!/usr/bin/env bash

idir=$1
oname=$2

idir=${idir%/}

names=`du -a $idir | grep "bed$" | cut -f 2`
for name in $names
do
project=`basename $name .bed`
cat $name | cut -f 4,5 | sed "1ifeature\t${project}" > "$idir/$project.summary.tmp"
echo "### `date "+%D %T"`  $project "
done

paste ${idir}/*summary.tmp | awk 'BEGIN{OFS="\t"}{printf("%s",$1);for (i=2;i<=NF;i++){if (i%2==0) {printf("\t%s",$i)}} printf("\n")}' > "${oname}"

rm ${idir}/*summary.tmp