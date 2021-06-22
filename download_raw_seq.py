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
    def __init__(self, *args):
        self.standard = args

    def filter(self, target: str):
        for i in self.standard:
            if re.search(i, target):
                return True
        return False

    def __str__(self):
        return ";".join(self.standard)


def generate_opt():
    opt = argparse.ArgumentParser()
    opt.add_argument("-p", "--project", dest="project",
        action="store", required=True, help="the project you want to download, used to down metafile \n")

    opt.add_argument("-r", "--re",dest="re",action="store_true",required=False,default=False,help="force to use re but not file as filter, if filter happen to be same with a file name \n")

    opt.add_argument("-f", "--filter",dest="filter",action="append",required=False,default=[],help="how to filter out sra you want to download, \nregex and list in file is allowed. defalut is all \n") 

    opt.add_argument("-o", "--outdir",dest="odir",action="store",required=False,default="rawdata",help="where to save downloaded files\n")

    opt.add_argument("-c", "--cmd",dest="cmd",action="store",default="download.sh",help="where you want to output download cmd scripts\n")

    opt.add_argument("-m", "--md5", dest="md5", action="store", default="md5.tsv", help="where to save your md5 info")

    args = opt.parse_args()

    return args


def validate_opt(args: argparse.ArgumentParser):
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
    if len(fs) == 0:
        fs = [".*"]
    filterTool = Filter(*fs)
    logger.info(f"choosed filter is {filterTool}")

    logger.info(f"outdir is {args.odir}")
    logger.info(f"output cmd path is {args.cmd}")

    ia = os.path.join(os.getcwd(), f"{args.project}.info.all.tsv")
    ii = os.path.join(os.getcwd(), f"{args.project}.info.tsv")

    logger.info(f"out put all info tsv is {ia}")
    logger.info(f"out put selected info tsv is {ii}")

    args.ia = ia
    args.ii = ii
    args.filterTool = filterTool


    ###odir
    odir=pathlib.Path(args.odir)
    odir.mkdir(parents=True, exist_ok=True)
    logger.info(f"odir is {odir}")

    ###md5
    md5=args.md5
    md5=pathlib.Path(md5)
    md5.resolve().parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"md5 summary file is {md5}")
    return args, filterTool


def get_metainfo(project: str, filterTool: Filter, ia: str, ii: str) -> pd.DataFrame:
    logger.info(f"start to download meta info to {ia}")

    cmd=f'''wget -O {ia} "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={project}&result=read_run&fields=sample_accession,experiment_accession,run_accession,scientific_name,library_layout,fastq_md5,fastq_ftp&format=tsv&download=true" '''
    code = 0
    #code = subprocess.run(cmd, shell=True)
    if code == 1:
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
    df_metainfo = pd.read_csv(ii, sep="\t")
    assert isinstance(df_metainfo, pd.DataFrame)
    return df_metainfo


def get_raw_seqs(df_metainfo, odir: str) -> str:
    def get_raw_seq_byrow(url, odir) -> str:
        urls = url.split(";")
        for i, j in enumerate(urls):
            fastq_file = j.rsplit('/', 1)[1]
            fastq_file_path = j.split('/', 1)[1]
            cmd = f'''cd {odir} && if [ ! -f {fastq_file} ]; then ~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:{fastq_file_path} .; fi '''
            yield cmd

            cmd = f'''cd {odir} && if [ ! -f {fastq_file} ]; then wget ftp://ftp.sra.ebi.ac.uk/{fastq_file_path}; fi '''
            yield cmd

    cmds = []
    assert isinstance(df_metainfo,pd.DataFrame)
    for index, row in df_metainfo.iterrows():
        for cmd in get_raw_seq_byrow(row["fastq_ftp"], odir):
            cmds.append(cmd)

    return cmds


def get_md5_sum(df_metainfo: pd.DataFrame, odir: str, md5_name: str) -> str:

    logger.info("generate md5sum cmd ... ")

    def get_md5_sum_byrow(md5, url):
        md5s = md5.split(";")
        urls = url.split(";")
        fqs = [os.path.join(odir, i.rsplit('/', 1)[1]) for i in urls]
        for i, j in zip(md5s, fqs):
            cmd = f'''{i}\t{j}'''
            yield cmd

    with open(md5_name, 'w') as f:
        for index, row in df_metainfo.iterrows():
            for line in get_md5_sum_byrow(row["fastq_md5"], row["fastq_ftp"]):
                f.write(line)
                f.write("\n")
    cmd = f" md5sum -c {md5_name} "

    return [cmd]


def run():
    args = generate_opt()
    args, filterTool = validate_opt(args)
    df_metainfo = get_metainfo(args.project, filterTool, args.ia, args.ii)
    assert isinstance(df_metainfo, pd.DataFrame)
    #print(df_metainfo)
    cmds = get_raw_seqs(df_metainfo, args.odir)
    cmds += get_md5_sum(df_metainfo, args.odir, args.md5)
    return "\n".join(cmds)


if __name__ == '__main__':
    cmd=run()
    print(cmd)
