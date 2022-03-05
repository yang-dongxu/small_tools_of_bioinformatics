from email.policy import default
import os
import io
import sys
import click
import pathlib
import logging 
import tempfile

import subprocess as sp
import numpy as np
import pandas as pd
from scipy import stats

def preprocess_regions(region:list, noMerge):
    '''region: -r input files. a list\n regions: stdin\n noMerge: whether merge'''
    logging.info("## load other input regions")
    all_regions = list(region)
    all_regions_str = " ".join(all_regions)
    f1 = tempfile.NamedTemporaryFile(mode="w",encoding="utf",delete=False)
    path1 = f1.name
    f1.close()
    if noMerge: #combine all peaks without more operation
        cmd = f'''cat {all_regions_str} | sort -k1,1 -k2,2n | cut -f 1-3  > {path1} '''
    else:
        cmd = f'''cat {all_regions_str} | sort -k1,1 -k2,2n | cut -f 1-3 | bedtools merge | sort -k1,1 -k2,2n > {path1} '''
    os.system(cmd)
    return path1

def count_all_reads(file:str):
    kind="bam"
    if file.endswith("am"):
        logging.info(f"## determine {file} as SAM/BAM")
    else:
        logging.info(f"## determine {file} as BED")
        kind="bed"
    
    if kind=="bam":
        cmd = f"samtools view {file} -t 4 | wc -l"
    else:
        cmd = f"cat {file} | wc -l"
    c = sp.run(cmd, capture_output=True,shell=True,text=True)
    logging.info(f"## {file} with reads {int(c.stdout)}")
    return int(c.stdout)

    

def counts(a,b,regions):
    '''a: sample A bam or bed\n b:sample B bam or bed\nregions:a file as bed3 format'''
    logging.info("## count total number of input bam/beds...")
    numA, numB = count_all_reads(a), count_all_reads(b)

    logging.info("## count reads on regions...")
    cmd = f'''bedtools intersect -a {regions} -b {a} -c | bedtools intersect -a - -b {b} -c '''
    c = sp.run(cmd,capture_output=True,shell=True,text=True)
    df = pd.read_csv(io.StringIO(c.stdout),sep="\t",names=["chr","start","end","A","B"])
    df["total"] = df["A"] + df["B"]
    df["A_total"] = numA
    df["B_total"] = numB  
    os.remove(regions)
    return df

@click.command(help="find significant different peaks between given sample A and B. ")
@click.option("-a",type=click.Path(exists=True),help="Bam or bed of sample A")
@click.option("-b",type=click.Path(exists=True),help="Bam or bed of sample A")
@click.option("-r",'--region',multiple = True,type=click.Path(exists=True), required =False, help="region files. ")
@click.option("-o","--ofile",type=click.File("w",encoding="utf8"), default=sys.stdout, help="output path. Default is stdout")
@click.option("--noMerge",is_flag=True, flag_value=False,help="Not merge region and call diff separately. By efault all peaks are merged by `bedtools merge` at first.")
@click.option("-A","--pa",type=str,default="A",help="alias for sampleA. It will be used as column name")
@click.option("-B","--pb",type=str,default="B",help="alias for sampleB. It will be used as column name")
def diffbind(a,b,region,ofile,nomerge,pa,pb):
    logging.info("## start to process input regions...")
    regions = preprocess_regions(region,nomerge)

    logging.info("## start to count reads on regions...")
    df = counts(a,b,regions)

    logging.info("## start to calculate pvalue...")
    _a = df.head(5).to_csv(sep='\t',index=False).replace('\n','\n## ')
    logging.info(f"## {_a}")
    df["pvalue"] = df.apply(lambda x: stats.binom_test(x["A"], x["total"], x["A_total"]/(x["A_total"]+x["B_total"])),axis=1)
    df["B_larger"] = df["B"] > df["total"]*df["B_total"]/(df["A_total"]+df["B_total"])

    logging.info("## start to output...")
    df.columns = [i.replace("A",pa).replace("B",pb) for i in df.columns]
    df.to_csv(ofile,sep="\t",index=False)
    logging.info("## over! see you~")
    return

if __name__ == "__main__":
    logging.basicConfig(level = logging.INFO)
    diffbind()