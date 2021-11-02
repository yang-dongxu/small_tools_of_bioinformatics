#!/usr/bin/env python
import enum
import os
import re
import io
import logging
import sys
import click
from click.decorators import group

import pandas as pd


@click.command()
@click.option("-i","--ifile",type=click.File("r"), default=sys.stdin,help="your download project info file. Generate by download_raw_seq.py")
@click.option("-o","--ofile",type=click.File("w"),default=sys.stderr,help="out name")
@click.option("-f","--filter",type=str,default=".", help = "filter you input file")
@click.option("-g","--groupby",multiple=True,default=["experiment_accession","sample_title"],help = "how to concentrate rows together")
@click.option('-s',"--sep",default=";",help="seperate cols")
@click.option("-c","--collasp",default=",",help="seperate iterms in single col")
@click.option("-d","--ifs",default="\\t",help="input file delimiter")
@click.option("--idir",default=os.path.join(os.path.abspath(os.getcwd()),"rawdata"),help="the idir col in output file")
@click.option("--run",default="run_accession",help="run accession col name")
@click.option("--suffix",multiple=True, default=["_1.fastq.gz","_2.fastq.gz"], help="the suffix you want to append to run accession")
def generate(ifile,ofile,filter,groupby,sep,collasp,ifs,idir,run,suffix):
    text = ifile.read().strip()

    if len(text) ==0:
        logging.warning("Input file is blank!")
        sys.exit(1)
    
    lines = [i.strip() for i in text.split("\n") if len(i.strip())]
    lines = [i for i in lines if re.search(filter,i)] ## filter lines
    text = "\n".join(lines)

    df = pd.read_csv(io.StringIO(text),sep=ifs)
    print("load input file: ")
    print("#"+df.head(3).to_csv(sep="\t").replace("\n","\n#"))

    groupby = list(groupby)
    dfo = df.groupby(groupby)[run].apply(lambda x: collasp.join(x)).reset_index()
    dfo.columns = groupby+["samples"]

    suffix_titles = []
    for i,s in enumerate(suffix):
        dfo[f"O{i}"]=dfo["samples"].apply(lambda x:x.replace(collasp,f"{s}{collasp}")+s)
        suffix_titles.append(f"O{i}")

    for i in groupby:
        dfo[i]=dfo[i].apply(lambda x: "_".join(x.split()))
    
    dfo["idir"] = idir

    ocols = groupby + ["idir"] + suffix_titles
    ofile.write(dfo[ocols].to_csv(sep=sep,index=False, header=False))
    logging.info("Out put over!")
    return dfo[ocols]

if __name__ == "__main__":
    generate()
    