#!/usr/bin/env bash
## edit from https://github.com/wwang-chcn/wwang_bioinfo_tools/blob/master/genome_source_build/build_genome_source.sh

USAGE="<genomeVersion>  <threads>"

MY_PATH=`echo $MY_PATH`

genomeVersion=${1}
threads=${2}


rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/bigZips/${genomeVersion}.2bit .
twoBitToFa ${genomeVersion}.2bit ${genomeVersion}.fa

cat ${genomeVersion}.fa | awk 'BEGIN{flag=1} (flag==1 && index($1,">")==0){print $0}(index($1,">")!=0){if (index($1,"_")==0) {flag=1;print $0; } else flag=0} ' > ${genomeVersion}_main.fa
#python ${MY_PATH}/del_fa.py ${genomeVersion}.fa > ${genomeVersion}_main.fa

# mapping index
bowtie2-build --threads ${threads} ${genomeVersion}.fa ${genomeVersion}
hisat2-build -p ${threads} ${genomeVersion}.fa ${genomeVersion}
bwa index ${genomeVersion}.fa

bowtie2-build --threads ${threads} ${genomeVersion}_main.fa ${genomeVersion}_main
hisat2-build -p ${threads} ${genomeVersion}_main.fa ${genomeVersion}_main
bwa index ${genomeVersion}_main.fa

mkdir -p bismark_index/{raw,main}
cd bismark_index/raw && ln -s ../../${genomeVersion}.fa .
cd ../main && ln -s ../../${genomeVersion}_main.fa .
cd ../..
bismark_genome_preparation --parallel ${threads} --single_fasta ./bismark_index/raw
bismark_genome_preparation --parallel ${threads} --single_fasta ./bismark_index/main

# STAR index
mkdir -p star/raw
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir star/raw --genomeFastaFiles ${genomeVersion}.fa 

mkdir -p star/main
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir star/main --genomeFastaFiles ${genomeVersion}_main.fa 




twoBitInfo ${genomeVersion}.2bit ${genomeVersion}.chrom.sizes
grep -v "_" ${genomeVersion}.chrom.sizes > ${genomeVersion}_main.chrom.sizes


# annotation
# refGene
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/refGene.txt.gz ${genomeVersion}.refGene.txt.gz
gunzip ${genomeVersion}.refGene.txt.gz
cut -f 2-11 ${genomeVersion}.refGene.txt > ${genomeVersion}.refGene.genePred
cut -f 2-16 ${genomeVersion}.refGene.txt > ${genomeVersion}.refGene.genePredExt
genePredToGtf -utr file ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.gtf
genePredToBed ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.bed

# refGene with out random chromsomes
cat  ${genomeVersion}.refGene.txt | awk 'index($3,"_")==0{print $0}' | sort -k3,3 -k5,5n  > ${genomeVersion}.refGene.main.txt
cut -f 2-11 ${genomeVersion}.refGene.main.txt > ${genomeVersion}.refGene.main.genePred
cut -f 2-16 ${genomeVersion}.refGene.main.txt > ${genomeVersion}.refGene.main.genePredExt
genePredToGtf -utr file ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.main.gtf
genePredToBed ${genomeVersion}.refGene.main.genePredExt ${genomeVersion}.refGene.main.bed

#get table of transcipt id to gene name
cat ${genomeVersion}.refGene.txt | cut -f 2,13  >  ${genomeVersion}.refGene.transcript_to_gene.tsv

# prepare promoter with 2k of tss
cat  ${genomeVersion}.refGene.bed | awk '{if ($6 == "+") {start=$2-2000;end=$2+2000;} else if ($6 == "-" ) {start=$3-2000;end=$3+2000;}else {start="Start";end="End"};printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,start,end,$4,$5,$6)}' | awk '{if ($5 != "Strand") {if ($2>0) start=$2;else {start=0};} else start="Start" ; printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,start,$3,$4,$5,$6)}' | sort -k1,1 -k2,2n > ${genomeVersion}.refGene.promoter_2k.bed
cat  ${genomeVersion}.refGene.main.bed | awk '{if ($6 == "+") {start=$2-2000;end=$2+2000;} else if ($6 == "-" ) {start=$3-2000;end=$3+2000;}else {start="Start";end="End"};printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,start,end,$4,$5,$6)}' | awk '{if ($5 != "Strand") {if ($2>0) start=$2;else {start=0};} else start="Start" ; printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,start,$3,$4,$5,$6)}' | sort -k1,1 -k2,2n > ${genomeVersion}.refGene.main.promoter_2k.bed


# prepare repeats annotations
url="https://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/rmsk.txt.gz"
wget -c ${url} -O ${genomeVersion}.repeats.txt.gz
gunzip ${genomeVersion}.repeats.txt.gz
cat ${genomeVersion}.repeats.txt | cut -f 6,7,8,10,11 | awk 'BEGIN{OFS="\t"}{$6=$4;$4=$5;$5="0";print $0}' | sed '1d' | sort -k1,1 -k2,2n > ${genomeVersion}.repeats.bed
cat ${genomeVersion}.repeats.txt | cut -f "11-13" | sed '1d' | sort | uniq | sed '1i#feature\tcalss\tfamily' > ${genomeVersion}.repeats.classfication.tsv
cat ${genomeVersion}.repeats.bed | awk '{ if (!($4 in name)) {name[$4]=0;}name[$4]+=$3-$2}END{for (var in name) {print var "\t" name[var]}}' | sort > ${genomeVersion}.repeats.length.tsv

cat ${genomeVersion}.repeats.txt | awk 'index($6,"_")==0' | awk '$12!="Simple_repeat"' |awk '$12 != "Low_complexity"'  > ${genomeVersion}.repeats.main.txt
cat ${genomeVersion}.repeats.main.txt | cut -f 6,7,8,10,11 | awk 'BEGIN{OFS="\t"}{$6=$4;$4=$5;$5="0";print $0}' | sed '1d' | sort -k1,1 -k2,2n > ${genomeVersion}.repeats.main.bed
cat ${genomeVersion}.repeats.main.bed | awk '{ if (!($4 in name)) {name[$4]=0;}name[$4]+=$3-$2}END{for (var in name) {print var "\t" name[var]}}' | sort > ${genomeVersion}.repeats.main.length.tsv

## with NR as suffix in col4
cat ${genomeVersion}.repeats.main.bed | awk 'BEGIN{OFS="\t"}{$4=$4 "_" NR;print $0}' > ${genomeVersion}.repeats.main.transcript.bed

cat ${genomeVersion}.repeats.main.transcript.bed | awk '{ if (!($4 in name)) {name[$4]=0;}name[$4]+=$3-$2}END{for (var in name) {print var "\t" name[var]}}' | sort > ${genomeVersion}.repeats.main.transcript.length.tsv

#转换出transc的saf
cat ${genomeVersion}.repeats.main.transcript.bed | awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3,$6}'  | sed '1iGeneID\tChr\tStart\tEnd\tStrand' > ${genomeVersion}.repeats.main.transcript.saf
cat ${genomeVersion}.repeats.main.bed | awk 'BEGIN{OFS="\t"}{print $4,$1,$2,$3,$6}'  | sed '1iGeneID\tChr\tStart\tEnd\tStrand' > ${genomeVersion}.repeats.main.saf

#转换出gtf
bedToGenePred ${genomeVersion}.repeats.main.transcript.bed  stdout | genePredToGtf file stdin ${genomeVersion}.repeats.main.transcript.gtf
bedToGenePred ${genomeVersion}.repeats.transcript.bed  stdout | genePredToGtf file stdin ${genomeVersion}.repeats.transcript.gtf



# ensGene
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/ensGene.txt.gz ${genomeVersion}.ensGene.txt.gz
if [[ -e ${genomeVersion}.ensGene.txt.gz ]]; then
    gunzip ${genomeVersion}.ensGene.txt.gz
    cut -f 2-11 ${genomeVersion}.ensGene.txt > ${genomeVersion}.ensGene.genePred
    cut -f 2-16 ${genomeVersion}.ensGene.txt > ${genomeVersion}.ensGene.genePredExt
    genePredToGtf -utr file ${genomeVersion}.ensGene.genePred ${genomeVersion}.ensGene.gtf
    python ${MY_PATH}/correctGtfGeneID.py ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.gtf
    genePredToBed ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.bed
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genomeVersion}/database/ensemblToGeneName.txt.gz ${genomeVersion}.ensemblToGeneName.txt.gz
    gunzip ${genomeVersion}.ensemblToGeneName.txt.gz
    python ${MY_PATH}/getEnsGeneToGeneName.py ${genomeVersion}
fi

python ${MY_PATH}/genePredExtToSqlite3.py ${genomeVersion}.refGene.genePredExt ${genomeVersion}.refGene.sq3
if [[ -e ${genomeVersion}.ensGene.txt.gz ]]; then
    python ${MY_PATH}/genePredExtToSqlite3.py ${genomeVersion}.ensGene.genePredExt ${genomeVersion}.ensGene.sq3
fi