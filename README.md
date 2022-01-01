# small_tools_of_bioinformatics


## prepare steps

| script | description | 
| -- | -- |
| download_raw_seq.py | download fastq.gz from ena instead of ncbi sra |
| sra_download.py | if the project is absent from ena, then you can use this to download sra from ncbi, then it helps to extract fastq and rm sras |
| build_genome_envrionment.sh | quickly build a well-equiped annotation environment |
| genePredExtToSqlite3.py | build ceas annotation |


## preprocess steps
| script | description | 
| -- | -- |
| bw_fc.sh | normalise bigwig signal based on global | 
| bw_merge_mean.sh | merge many bws by mean to a single one | 
| cal_arrayValue_from_mnase_signal.py | arrayValue describes how well the nucleosome is positioned |
| bamtobedpe.sh | turn pe bam to fragments | 
| bamtoreads.sh | turn pe bam to reads, with pairs in adjecent lines |
| snpSelection.py | a script to split bam file based on snp, source from Yanghui () |
| normalize_RNA_expression.py | normalize outputs from generalized counts files and a width file. **decrepted ** |
| extract_log.py| extract logs of given format. It's a combination of other logs exctractor here |
| extractCutSiteToBw.sh | extract cut sites of DNase-seq, for foot-print analysis | 

## QC steps
| script | description | 
| -- | -- |
| bwcor.py | calculate bw signal correlation within bed regions |
| peak_stat.py | a script to analyze frip | 
| region_enrichment.py | calculate how enriched of bed1 in regions defined in bed2 | 
| ucscHubTojbrowser2Config.py | jbrowser2 is a genome browser similiar with ucsc genome browser, but is quicker. the scripts turn ucsc hub file to jbrowser2 config. |


## usefull tools 
| script | description | 
| -- | -- |
| gtf_to_table.py | turn gtf to a tsv table, usefull for stat stringtie output | 
| get_idr_regions.py | output IDR regions based on alphafold and other process | 
| bed_to_matrix.sh| a script to combine bed6 files to a tsv, with 4th col as index  |
| tableStack.py | a script to combine tables together, vstack or hstack(merge, not complete) |
| 

## others

| script  | description |
| --- | --- |
| bw_summary_bed.py | summary bigwig file by mean in region of given bed. **used bigwigAverageOverBed instead** |
| bw_summary_bed_bins.py| summary bigwig file by mean in the region of given bed and bin nums your input  |
| normaliseRNA_featureCounts.py| normalise out put of software [featureCounts](http://subread.sourceforge.net/), but at first you have to turn the last column name to Reads. May decrept soon |
| normalize_RNA_expression.py| normalise RNA expression by given count tsv (col| feature,counts) and feature length tsv(col| feature, length) |
| normalize_homer_counts.sh| normalise software [HOMER](http://homer.ucsd.edu/homer/) output of repeatsAnalyze repeats subcommand. |
| vstack_files.py| [decrept] a script to vstack tabless |
| get_from_ena.sh | *use download_raw_seq.py instead*,a script to download fastq from ena, source from [wangwen](https://github.com/wwang-chcn/) and [Reemann](https://github.com/Reemann) |
| cellRangeoutv2Tomatrix.py | cellRange output are three files, which is not convient to process by you self. This helps you. | 
