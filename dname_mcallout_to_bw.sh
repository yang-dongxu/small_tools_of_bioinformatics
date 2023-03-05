#!/usr/bin/bash


usage="$0  chromsize  threshold oname  <beds...>"

chromsize=$1
shift;
threshold=$1
shift;
oname=$1
shift;
beds="$*"

tobdg() {
	cat $beds |  sed '/#/d' |   awk -v t=$threshold 'BEGIN{OFS="\t"} $5>t{$3=$2+1;print $0}' |\
	cut -f 1-4 | grep -e  "^chr[0-9|XYM^_]\{1,2\}\b" | sort -k1,1 -k2,2n  | \
	bedtools merge -c 4 -o mean | cut -f 1-4  | sort -k1,1 -k2,2n  
}

date +"%D %T: start to proceess ${oname} from input ${beds} with threshold: ${threshold}"

if [  ! -f  $chromsize ] ; then
	tobdg |  wigToBigWig -clip stdin <(fetchChromSizes $chromsize ) $oname
else
	tobdg |  wigToBigWig -clip stdin  $chromsize $oname
fi;

date +"%D %T: ${oname} over!"
