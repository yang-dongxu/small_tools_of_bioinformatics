import os
import sys
import re
import argparse
import logging

import pandas as pd
import numpy as np
import pyBigWig as pbw

from copy import deepcopy


# create logger
logger_name = "bigWigSummary-bins"
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
    parser=argparse.ArgumentParser()
    parser.add_argument('-b',"--bed",action="store",type=str, dest="bed", required=True, help="the bed file your input, must have 6 columns or more, and the 4th col will be used as id in output\n")
    parser.add_argument('-w',"--bigwig",action="store",type=str, dest="bw", required=True, help="the bigwig file you want to scan\n")
    parser.add_argument('-n',"--bins",action="store",type=int,default=1,dest="bins",help="how many data points you want to scan in the each region\n")
    parser.add_argument('-5',"--upstream",action="store",type=int,default=2000,dest="upstream",help="how many bp to extense the bed region up stream (strand specific, 5' direction)\n ")
    parser.add_argument('-3',"--downstream",action="store",type=int,default=2000,dest="downstream",help="how many bp to extense the bed region down stream (strand specific, 3' direction) \n")

    parser.add_argument('--5bins',action="store",type=int,default=100,dest="bins_5",help="how many bins  of 5'direction tail\n")
    parser.add_argument('--3bins',action="store",type=int,default=100,dest="bins_3",help="how many bins  of 3'direction tail\n")
    parser.add_argument('-o','--outname',action="store",type=str,default="std",dest="outname",help="output file name, tsv file. first col is the id in given bed, second is the description of part (upstream, region, and downstream), third is the bin id  in given part, and 4th is the real intensity\n")
    parser.add_argument('-0','--non-directional',action="store_false",default=True,dest="directional",help="output the bins in the order of 5->3, or by the order of site, regardless of strand if choose, not recommend!\n")

    return parser


def validate_opt(args:argparse.ArgumentParser):
    flag=True
    if not  os.path.isfile(args.bed):
        logger.error(f"input bed file {args.bed} is not exist! ")
        flag=False
    elif not os.path.isfile(args.bw):
        logger.error(f"input bigwig file {args.bw} is not exist! ")
        flag=False
    elif args.upstream< args.bins_5 or args.downstream< args.bins_3:
        logger.error("bins num > region length!")
        args.bins_5=args.upstream
        args.bins_3=args.downstream
        logger.warning(f"turn 5bins to {args.bins_5}")
        logger.warning(f"turn 3bins to {args.bins_3}")
    
    if not flag:
        sys.exit(1)
    
    return args

def get_bw_stat(x:pd.Series,bw,nbins=1):
    strand=x["strand"]
    outs=bw.stats(x["chr"],int(x["start"]+1),int(x["end"]),nBins=nbins)
    if strand=="-":
        outs=outs[::-1]
    return outs

def process_each_part(df,bw,nbins,label):
    df_p=df.query("start+1 < end")
    df_result=pd.DataFrame(list(df_p.apply(get_bw_stat,args=(bw,nbins),axis=1)))
    df_result["name"]=df_p["name"]
    df_result=df_result.melt(id_vars="name",var_name="site",value_name="intensity")
    df_result["label"]=label
    return df_result

def split_df_bed(df_bed,upstream,downstream):
    df_upper=deepcopy(df_bed)
    df_upper["start"]=df_bed.apply(lambda x: x["end"] if x["strand"]=="+" else x["start"]-downstream,axis=1)
    df_upper["end"]=df_bed.apply(lambda x: x["end"]+downstream if x["strand"]=="+" else x["start"],axis=1)

    df_down=deepcopy(df_bed)
    df_down["start"]=df_bed.apply(lambda x: x["start"]-upstream if x["strand"]=="+" else x["end"],axis=1)
    df_down["end"]=df_bed.apply(lambda x: x["start"] if x["strand"]=="+" else x["end"]+upstream,axis=1)

    return df_upper,df_bed,df_down

def output(df:pd.DataFrame,name):
    if name=="std":
        df.to_csv(sep="\t",index=False)    
    else:
        df.to_csv(name,sep="\t",index=False)
    return True

def process(args):
    df_bed= pd.read_csv(args.bed,sep="\t",comment="#",usecols=range(0,6),header=None)
    assert isinstance(df_bed,pd.DataFrame)
    df_bed.columns=["chr","start","end","name","score","strand"]
    bw=pbw.open(args.bw)

    labels=["0","1","2"]
    results=[]
    bins=[args.bins_5,args.bins,args.bins_3]
    for df, label, nbin in zip(split_df_bed(df_bed,args.upstream,args.downstream),labels,bins):
        if bins==0:
            logger.warning(f"skip label {label} becaus bins num is 0")
            continue
        new_df=process_each_part(df,bw,nbin,label)
        results.append(new_df)
        new_df==None
        logger.info(f"label {label} part is processed!")

    df_result=pd.concat(results).fillna(0).sort_values(["name","label","site"])
    logger.info("prepare to output")
    output(df_result,args.outname)
    logger.info("output over! see you ~")
    return df_result
    

if __name__=='__main__':
    opt=generate_opt()
    args=opt.parse_args()
    logger.info("validating opts input ...")
    validate_opt(args)
    logger.info("input options pass check!")
    process(args)