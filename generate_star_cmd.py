#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys


# In[11]:


try:
    usage=f"##python {__file__} project seq1 seq2  index threads outdir"
    if len(sys.argv)<5:
        print(usage)
        sys.exit(1)
    project=sys.argv[1]
    seq1=sys.argv[2]
    seq2=sys.argv[3]
    index=sys.argv[4]
    threads=int(sys.argv[5])
    if len(sys.argv)==7:
        outdir=os.path.abspath(sys.argv[6])
    else:
        outdir=os.getcwd()
except:
    usage="python this_script project seq1 seq2"
    project="GR-1210-5"
    seq1="/mnt/Storage2/home/zengshiyang/Work_YDX/1.PGC/j098_raw_result/CSPoutput/trimed/GR-1210-5_2a_1_val_1.fq.gz"
    seq2="/mnt/Storage2/home/zengshiyang/Work_YDX/1.PGC/j098_raw_result/CSPoutput/trimed/GR-1210-5_2a_2_val_2.fq.gz"
    index="/mnt/Storage2/home/zengshiyang/DB/STAR/hg38"
    threads=10
    outdir="STARout"

assert os.path.exists(seq1)
assert os.path.exists(seq2)
assert os.path.exists(index)

    


# In[12]:
outpath=os.path.join(outdir,project)
if not os.path.exists(outpath):
    os.makedirs(outpath)
outname=os.path.join(outdir,project,project+"_")
cmd=f"STAR --genomeDir {index} --runThreadN {threads} --outSAMattributes NH HI NM MD XS AS  --outFilterMultimapNmax 500 --outFileNamePrefix {outname} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --readFilesIn {seq1} {seq2} \n "


cmd1=f"ls {outpath} | grep ort | grep bam  | awk -v i=\"{outpath}\" '{{print \"mv  \" i\"/\"$1 \"\t{outname}sort.bam\\n\"}}' | bash \n"

cmd2=f"featureCounts -g transcript_id -O --fracOverlap 0.1 -p -B -M --fraction -T 10 -a /mnt/Storage2/home/zengshiyang/DB/refGene/hg38.repeats.gtf -o {outname+'featurecounts'} -R BAM {outname+'sort.bam'} \n"

cmd3=f"samtools view -@ 10 -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x XS {outname+'sort.bam.featureCounts.bam'} | grep -e 'XT:Z' > {outname}featurecounts.sam.txt \n"

cmd4=f"redistribute_multiple_aligned_reads.r -f {outname}featurecounts.sam.txt  -r /mnt/Storage2/home/zengshiyang/DB/refGene/hg38.repeats.transcript_id.saf -n ./sam/{project}_distribution.txt -s 50 -m 1 -p 12"
# In[13]:


print(cmd)
print(cmd1)
print(cmd2)
print(cmd3)
print(cmd4)


# In[ ]:


sys.exit(0)


