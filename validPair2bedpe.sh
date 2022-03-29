#!/usr/bin/bash

allvalidpairs=$1
quality=$2
awk -v q=$quality 'BEGIN{OFS="\t"}($11>q&& $12>q){print $2,$3,$3+1,$5,$6,$6+1}' $allvalidpairs