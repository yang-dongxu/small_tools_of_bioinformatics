#!/usr/bin/env bash
odir=$1
project_acc=$2
sra_filter=$3

if [ ! -d "$odir"  ];then
	mkdir -p $odir
	echo "making new dirs ${odir} "
fi

now_dir=`pwd`

function ena_meta_download {
    project_accession=${1}
    wget -O ${project_accession}.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${project_accession}&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,library_layout,fastq_md5,fastq_ftp&format=tsv&download=true"
}

project_accession=${project_acc}

cd $odir
ena_meta_download ${project_accession}
mv ${project_accession}.tsv ${project_accession}.tsv.tmp
cat <(head -1 ${project_accession}.tsv.tmp) <(egrep ${sra_filter}  ${project_accession}.tsv.tmp) > ${project_accession}.tsv
# rm ${project_accession}.tsv.tmp
gen_download_from_ebi.py ${project_accession}.tsv ${project_accession}.sh
bash ${project_accession}.sh
md5_of_ENAfastq_generate.py ${project_accession}.tsv ${project_accession}.md5
md5sum -c ${project_accession}.md5 2>&1 | tee -a md5check.log
cd $now_dir 
