#!/usr/bin/env bash

ibam=$1
threads=$2
samtools sort   -@ ${threads} -n $ibam | \
    bamToBed -bedpe -i - -mate1  | sort -k1,1 -k2,2n | uniq | \
    awk 'BEGIN{OFS="\t"}
        {printf("%s\t%d\t%d\t%s/1\t%d\t%s\n",$1,$2,$3,$7,$8,$9)
            ;printf("%s\t%d\t%d\t%s/2\t%d\t%s\n",$4,$5,$6,$7,$8,$9) 
        }'
