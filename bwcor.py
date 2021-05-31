#!/usr/bin/env python 

import os
import sys
import pyBigWig as pbw
import argparse
import logging 
import scipy
import pandas as pd
import seaborn as sns
import numpy as np
from copy import deepcopy

# create logger
logger_name = "bwcor (bigwig correlation caculation)"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.DEBUG)

sh=logging.StreamHandler()
sh.setLevel(logging.DEBUG)

fmt = "%(asctime)-15s %(levelname)s %(name)s : pid - %(process)d :  %(message)s"
datefmt = "#%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
sh.setFormatter(formatter)
logger.addHandler(sh)

def generate_opt() -> argparse.ArgumentParser:
    opt=argparse.ArgumentParser()

    opt.add_argument("-b","--bed",action="append",dest="beds",default=[],required=True)
    opt.add_argument("-w","--bigwig",action="append",dest="bws",required=True)

    opt.add_argument("-p","--project",action="append",dest="projects",default=[])
    opt.add_argument("-na","--navalue",action="store",type=str,dest="na",default="na")

    opt.add_argument("-o","--outname",action="store",type=str,dest="oname",default="std")
    opt.add_argument("-d","--ofs",action="store",type=str,dest="ofs",default="\t")

    return opt


def validate_arg(arg) -> tuple :
    logger.info("start to validate input args")
    beds=",".join(arg.beds).split(",")
    bws=",".join(arg.bws).split(",")
    projects=",".join(arg.projects).split(",")

    logger.info(f"input bed files are :{beds}")
    logger.info(f"input bw files are :{bws}")
    logger.info(f"input  projects are :{projects}")



    if len(beds) != len(bws):
        logger.error(f"beds and bigwigs are different at length! check -b and -w")
        beds=beds
        #sys.exit(1)

    if str(arg.na).lower() in ["na","none","nan"]:
        arg.na=pd.NA
    else:
        try:
            na=float(arg.na)
        except:
            na=arg.na
        finally:
            arg.na=na
    logger.info(f"treat na values as {arg.na}")

    for bed in beds:
        if not os.path.isfile(bed):
            logger.error(f"bed file {bed} has errors! May not exisi.")
            sys.exit(1)
    
    for bw in bws:
        if not  os.path.isfile(bw):
            logger.error(f"bw file {bed} has errors! May not exisi.")
            sys.exit(1)

    if len(projects) != len(bws):
        if len(projects) ==0 :
            logger.warning("No projects info provide! Inferring from bws!")
            projects=[i.split(".")[0] for i in bws]
            logger.warning(f"Inferred projects are {projects}")
        else:
            logger.error(f"Length of project is not equal as bws! Check them!: {projects}")
            sys.exit(1)

    
    if arg.oname != "std":
        oname=os.path.abspath(arg.oname)
        path=os.path.split(oname)[0]
        if not os.path.isdir(path):
            os.makedirs(path)
        arg.oname=oname
    logger.info(f"out name : {arg.oname}")
    
    return beds, bws, projects,arg.na

def get_bw_stat(x:pd.Series,bw) -> float:
    strand="+"
    try:
        x["start"]=max(x["start"],0)
        outs=bw.stats(x["chr"],int(x["start"]),int(x["end"]),nBins=1)[0]
        outs=np.array(outs)
    except Exception as e:
        logger.error(f'Error at (chr, start, end): ({x["chr"],int(x["start"]),int(x["end"])})')
        logger.warning(f'Treat result of (chr, start, end): ({x["chr"],int(x["start"]),int(x["end"])}) as 0')
        logger.warning(e)
        outs=np.nan
    return outs

def process(beds,bws,projects,*args) -> pd.DataFrame:
    dfs=[]
    for bed in beds:
        logger.info(f"start loading bed file {bed}")
        try:
            df=pd.read_csv(bed,sep="\t",comment="#",header=None)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        names=["chr","start","end"]
        if len(df.columns) <3:
            logger.error(f"Not enought columns provided at {bed}! At least three columns!")
            sys.exit(1)
        df=df.iloc[:,:3]
        df.columns=names
        dfs.append(df)
    logger.info("All bed files are loaded. Start to merge")    
    regions=pd.concat(dfs)

    tmp_name=f"bwcor_{abs(np.random.randn()*100)}.bed"
    tmp_name_o=tmp_name+"out.bed"
    regions.to_csv(tmp_name,sep="\t",header=False,index=False)
    cmd=f"cat {tmp_name} | sort -k1,1 -k2,2n | bedtools merge -d -1 -i - > {tmp_name_o} "
    os.system(cmd)
    regions=pd.read_csv(tmp_name_o,sep="\t",names=["chr","start","end"])
    os.remove(tmp_name)
    os.remove(tmp_name_o)

    odf={}
    for bw,project in zip(bws,projects):
        bw_hander=pbw.open(bw)
        logger.info(f"start to process bigwig file {bw}")
        odf[project]=regions.apply(lambda x: get_bw_stat(x,bw_hander), axis=1)
        logger.info(f"end to process bigwig file {bw}")
    
    odf=pd.DataFrame(odf)

    p=odf.head().to_csv(sep="\t")
    logger.debug("stats example: \n## " +p.replace("\n","\n## "))
    p=odf.dropna().copy()
    for z in p.columns:
        p.loc[:,z]=deepcopy(p[z].apply(float))
    df_cor=p.corr(method ='pearson')
    print(df_cor)
    return df_cor

def output(df_cor:pd.DataFrame,oname="std",sep="\t"):

    logger.info("prepare output content ...")
    out=df_cor.to_csv(sep=sep)
    if oname=="std":
        fo=sys.stdout
    else:
        fo=open(oname,'w',encoding='utf8')
    logger.info("outputing ... ")
    fo.write(out)
    fo.close()
    logger.info(f"output over! see {oname}")
    return out

if __name__ == '__main__':
    logger.info("Hello ~")
    opt=generate_opt()
    arg=opt.parse_args()
    beds, bws, projects, na=validate_arg(arg)
    odf=process(beds,bws,projects,na)
    output(odf,arg.oname,arg.ofs)
    logger.info("see you~")
    



