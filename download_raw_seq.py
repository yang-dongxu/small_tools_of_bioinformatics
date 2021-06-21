import os
import sys
import logging
import pathlib
import argparse
import subprocess
import re

import numpy as np
import pandas as pd


# create logger
logger_name = "download raw data"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.DEBUG)

sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)

fmt = "%(asctime)-15s %(levelname)s %(name)s : pid - %(process)d :  %(message)s"
datefmt = "#%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
sh.setFormatter(formatter)
logger.addHandler(sh)

# '''this script aims to download fastq files from ncbi or ena, default is ena'''

class Filter:
    def __init__(self,*args):
        self.standard=args

    def filter(self,target:str):
        for i in self.standard:
            if re.search(i,target):
                return True
        return False

    def __str__(self):
        return ";".join(self.standard)

def generate_opt(): 
    opt=argparse.ArgumentParser()
    opt.add_argument("-p", "--project", dest="project", 
        action="store", required=True, help="the project you want to download, used to down metafile \n")
    opt.add_argument("-r", "--re",dest="re",action="store_true",required=False,default=False,help="force to use re but not file as filter, if filter happen to be same with a file name \n")
    opt.add_argument("-f", "--filter",dest="filter",action="append",required=False,default=[],help="how to filter out sra you want to download, \nregex and list in file is allowed. defalut is all \n") 
    	
    opt.add_argument("-o", "--outdir",dest="odir",action="store",required=False,default="rawdata",help="where to save downloaded files\n")

    opt.add_argument("-c", "--cmd",dest="cmd",action="store",default="download.sh",help="where you want to output download cmd scripts\n")

    args=opt.parse_args()

    return args


def validate_opt(args:argparse.ArgumentParser):
    logger.info("start to validate opts...")

    logger.info(f"input project are {args.project}")


    fs = []
    for f in args.filter:
        if os.path.isfile(f):
            if not args.re:
                logger.info("treat filter as file input")
                with open(args.filter) as fi:
                    ps = fi.read().split("\n")
                fs += ps
            else:
                logger.warning(f"filter {args.filter} is exist as a file, but forced to treat as string filter")
        fs += [f]
    if len(fs)==0:
        fs=[".*"]
    filterTool = Filter(*fs)
    logger.info(f"choosed filter is {filterTool}")

    logger.info(f"outdir is {args.odir}")
    logger.info(f"output cmd path is {args.cmd}")

    ia = os.path.join(os.getcwd(), f"{args.project}.info.all.tsv")
    ii = os.path.join(os.getcwd(), f"{args.project}.info.tsv")

    logger.info(f"out put all info tsv is {ia}")
    logger.info(f"out put all info tsv is {ii}")

    args.ia=ia
    args.ii=ii
    args.filterTool=filterTool

    return args, filterTool

def get_metainfo(project:str,filterTool:Filter,ia:str,ii:str):
    logger.info(f"start to download meta info to {ia}")

    cmd=f'''wget -O {ia} "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,library_layout,fastq_md5,fastq_ftp&format=tsv&download=true" '''
    code = 0
    code=subprocess.run(cmd,shell=True)
    if code==1:
        logger.error(f"download error! cmd: {cmd}")
        sys.exit(1)

    logger.info(f"start to filter meta info into {ii}")

    ia = pathlib.Path(ia)
    ii = pathlib.Path(ii)
    ia.parent.mkdir(parents=True, exist_ok=True)
    ii.parent.mkdir(parents=True, exist_ok=True)
    with open(ia, 'r') as f1, open(ii, 'w') as f2:
        for line in f1:
            if filterTool.filter(line):
                f2.write(line)
                f2.write("\n")
    logger.info(f"meta info processed! see {ii}")
    df_metainfo = pd.read_csv(str(ii), sep="\t")
    return df_metainfo


def run():
    args=generate_opt()
    args,filterTool=validate_opt(args)
    df_metainfo=get_metainfo(args.project,filterTool,args.ia,args.ii)

if __name__=='__main__':
    run()