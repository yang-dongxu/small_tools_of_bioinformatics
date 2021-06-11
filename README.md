# small_tools_of_bioinformatics

| script  | descriptio |
| --- | --- |
| bed_to_matrix.sh| a script to combine bed6 files to a tsv, with 4th col as index  |
| build_genome_environment.sh| a script to build genome environment  |
| bw_fc.sh | a script to transform bigwig file scale by fold of average signal  |
| bw_summary_bed.py | summary bigwig file by mean in region of given bed.  |
| bw_summary_bed_bins.py| summary bigwig file by mean in the region of given bed and bin nums your input  |
| concat_ceas_out.sh| [decrept] concat ceas out together  |
| del_fa.py| [decrept] delete alt and random chrosomes in given fasta  |
| extract_log.py| extract logs of given format. It's a combination of other logs exctractor here |
| extract_log_bowtie2.py| [decrept] |
| extract_log_star.py| [decrept] |
| extract_readdistribution.py| [decrept] |
| gtf_to_table.py| transfer gtf format to tables like csv or tsv |
| normaliseRNA_featureCounts.py| normalise out put of software [featureCounts](http://subread.sourceforge.net/), but at first you have to turn the last column name to Reads. May decrept soon |
| normalize_RNA_expression.py| normalise RNA expression by given count tsv (col| feature,counts) and feature length tsv(col| feature, length) |
| normalize_homer_counts.sh| normalise software [HOMER](http://homer.ucsd.edu/homer/) output of repeatsAnalyze repeats subcommand. |
| peak_stat.py| a script to analyze frip |
| sort_noheader.sh| [non-complete] |
| tableStack.py | a script to combine tables together, vstack or hstack(merge, not complete) |
| vstack_files.py| [decrept] a script to vstack tabless |
| gtf_to_table.py | a script to transfomr gtf to a table in tsv or csv |
| bwcor.py | a script to calculate bigwigs correlation on given bed |
| sra_downlaod.py | a pipeline to download sra and fastq-dump to given folder, from ncbi |
| region_enrichment.py | a script to calculate region on given bed |
| snpSelection.py | a script to split bam file based on snp, source from Yanghui () |
| get_from_ena.sh | a script to download fastq from ena, source from [wangwen](https://github.com/wwang-chcn/) and [Reemann](https://github.com/Reemann) |
