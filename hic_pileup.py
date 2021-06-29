import os
import sys
import io
import argparse
import pathlib 
import subprocess

from straw import straw
import logging

import numpy as np
from numpy.core.defchararray import endswith
import pandas as pd


from tqdm import tqdm

logger_name = "hic pileup"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.INFO)

sh=logging.StreamHandler()
sh.setLevel(logging.DEBUG)

fmt = "%(asctime)-15s %(levelname)s %(name)s : pid - %(process)d :  %(message)s"
datefmt = "#%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
sh.setFormatter(formatter)
logger.addHandler(sh)

def generate_opt() -> argparse.ArgumentParser:
    opt=argparse.ArgumentParser("Hic pileuper")

    opt.add_argument("-i","--input",dest="hic",action="store"
        ,required=True,help="input hic file. (fic format)")
    opt.add_argument("-b","--bed",dest="bed",action="store"
        ,required=True,help="input bed centers you want to pileup, just first three columns are used")

    opt.add_argument("-g","--genome",dest="genome",action="store",
        required=True,help="the chrom size file, download from uscs, tab-delimter")
    opt.add_argument("-j","--juicertools",dest="juicer",action="store",
        default="straw", help="where your juicer_tools.jar exist. if not set, will used straw instead.(Recommend)")
    opt.add_argument("-o","--oname",dest="oname",action="store"
        ,required=True,help="where to store your stat file")

    opt.add_argument('-r',"--resolution",dest="resolution",action="store",type=int,
        default=25000,help="resolution you want, corresponds to BP option in juicer_tools dump. must in [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000] ")
    opt.add_argument('-f',"--flank",dest="flank",action="store",type=int,
        default=100000,help="flank distance of center you want to observe, default is 100k bp. center position is the bin where region center is in .")

    opt.add_argument("-d","--sep",dest="ofs",action="store",
        default="\t",help="out put file delimter")
    opt.add_argument("-c","--compression",dest="compression",action="store",
        default="gzip",
        help='how to compression your output file, default is gzip, choices:["infer", "gzip", "bz2", "zip", "xz", "None"] for detail see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html')

    args=opt.parse_args()
    return args

def validate_opts(arg:argparse.ArgumentParser) -> argparse.ArgumentParser:
    logger.info("start to validate opts...")
    ## input file validate
    arg.hic=pathlib.Path(os.path.expandvars(arg.hic))
    if not arg.hic.is_file():
        logger.critical(f"input file {arg.hic.resolve()} is not exist!")
        sys.exit(1)
    logger.info(f"input hic file is {arg.hic}")
    
    arg.bed=pathlib.Path(os.path.expandvars(arg.bed))
    if not arg.bed.is_file():
        logger.critical(f"input file {arg.bed.resolve()} is not exist!")
        sys.exit(1)
    logger.info(f"input bed file is {arg.bed}")

    ## juicertools
    if arg.juicer.strip() == "straw":
        logger.info(f" choose straw as info extractor ...")
        arg.juicer=="straw"
    else:
        arg.juicer=pathlib.Path(os.path.expandvars(arg.juicer))
        if not arg.juicer.is_file():
            logger.critical(f"input file {arg.juicer.resolve()} is not exist!")
            sys.exit(1)
        logger.info(f"input juicer_tools path is {arg.juicer}")

    ## genome
    arg.genome=pathlib.Path(os.path.expandvars(arg.genome))
    if not arg.genome.is_file():
        logger.critical(f"input file {arg.genome.absolute()} is not exist!")
        sys.exit(1)
    logger.info(f"input genome path is {arg.genome}")
    
    ## compression
    if arg.compression in ["infer", "gzip", "bz2", "zip", "xz", "None"]:
        if arg.compression == "None":
            logger.info(f"choose compression method None. (No compression)")
            arg.compression=None
        elif arg.oname in ["-", "std"]:
            logger.error("you can't output comressed out to std, set compression as None")
            arg.compression=None
        else:
            logger.info(f"choose compression method {arg.compression}")
    else:
        logger.error("Not correct compression-format, set as gzip")
        arg.compression="gzip"

    ## create output dir
    if arg.oname == "-" or arg.oname == "std":
        logger.info(f"output stat to stdout")
        arg.fo=sys.stdout
    else:
        arg.oname=pathlib.Path(arg.oname)
        if not arg.oname.parent.exists():
            logger.warning("output dir is not exist, will be created now")
        arg.oname.parent.mkdir(parents=True,exist_ok=True)
        if arg.compression:
            arg.fo=open(arg.oname,'wb')
        else:
            arg.fo=open(arg.oname,'w')
        logger.info(f"output stat file is {arg.oname}")


    ## validate flank and resultion
    if arg.flank <0:
        logger.error(f"input flank is smaller than 0 ,set to 100000 (100k)")
        arg.flank=100000
    if not arg.resolution  in [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000] :
        logger.error(f"input resolution is smaller than 0 ,set to 25000 (25k)")
        arg.resolution=25000
    if arg.flank <= arg.resolution:
        logger.error(f"input flank distance is smaller than resolution, reset!")
        arg.flank=100000
        arg.resolution=25000
    logger.info(f"input flank is {arg.flank}")
    logger.info(f"input resolution is {arg.resolution}")

    logger.info("pass opt validation.")
    return arg

def get_center_bins(bed,genome,resolution) -> pd.DataFrame:
    while True:
        tmp=f"{np.random.rand()}.bed.tmp"
        if not os.path.exists(tmp):
            break
    cmd=f'''bedtools makewindows -g {genome} -w {resolution} > {tmp} && \
        cat {bed} | awk 'BEGIN{{OFS="\t"}}{{center=int(($2+$3)/2);$2=center;$3=$2+1;print $0}}' | cut -f 1-3 | \
        bedtools intersect -c -a {tmp} -b -  \
        '''
    p=subprocess.run(cmd,shell=True,capture_output=True,text=True)
    os.remove(tmp)
    logger.info(cmd)
    df=pd.read_csv(io.StringIO(p.stdout),sep="\t",names=["chr","start","end","count"])
    df=df.groupby(["chr","start","end"])["count"].sum().reset_index().query("count > 0").copy()
    return df 

def get_dump(juicertool,hic,chr,start,end,flank,resolution,weight) -> pd.DataFrame:
    while True:
        oname=f"{np.random.rand()}.tsv.tmp"
        if not os.path.exists(oname):
            break

    start, end = int(start), int(end)
    
    if start-flank>0:
        start_p=start - flank
    else:
        start_p=0
    end_p=end+flank
    region=f"{chr}:{start_p}:{end_p}"
    cmd=f'''java -jar {juicertool} dump observed KR {hic} {region} {region} BP {resolution} {oname} >> /dev/null 2>&1'''
    logger.debug(f"cmd: {cmd}")
    subprocess.run(cmd,shell=True)
    df=pd.read_csv(f"{oname}",sep="\t",names=["region1","region2","intensity"])
    os.remove(oname)
    df["region1"]=df["region1"]-start
    df["region1"]=(df["region1"]/resolution).astype(int)
    df["region2"]=df["region2"]-start
    df["region2"]=(df["region2"]/resolution).astype(int)

    dfs=[]
    for i in range(weight):
        dfs.append(df.copy())
    try:
        df=pd.concat(dfs)
    except:
        df=pd.DataFrame(columns=["region1","region2","intensity"])
    
    logger.debug(f"{df.to_csv()}".replace("\n","\n##").replace(",","\t"))
    return df

def get_dump_straw(straw_matrix,start,end,flank,resolution,weight):
    if start-flank>0:
        start_p=start - flank
    else:
        start_p=0
    end_p=end+flank
    results=straw_matrix.query(f"region1 < {end_p} and region1 >= {start_p}  and region2 < {end_p} and region2 >= {start_p}")
    df=results.copy()
    df["region1"]=df["region1"]-start
    df["region1"]=(df["region1"]/resolution).astype(int)
    df["region2"]=df["region2"]-start
    df["region2"]=(df["region2"]/resolution).astype(int)

    dfs=[]
    for i in range(weight):
        dfs.append(df.copy())
    try:
        df=pd.concat(dfs)
    except:
        df=pd.DataFrame(columns=["region1","region2","intensity"])
    logger.debug(f"{df.to_csv()}".replace("\n","\n##").replace(",","\t"))
    return df

def process(arg:argparse.ArgumentParser) -> pd.DataFrame:
    logger.info("start to generate center region...")
    df_center=get_center_bins(arg.bed,arg.genome,arg.resolution)

    logger.info("start to dump info by juicer_tools...")
    df=pd.DataFrame(columns=["chr","base","bin","region1","region2","intensity"])
    pbar=tqdm(range(len(df_center)))
    for (_,row),i in zip(df_center.iterrows(),pbar):
        chrom=row["chr"]
        start=row["start"]
        end=row["end"]
        dft=get_dump(arg.juicer,arg.hic,chrom,start,end,arg.flank,arg.resolution,row["count"])
        dft["chr"]=chrom
        dft["base"]=start
        dft["bin"]=arg.resolution
        df=df.append(dft,ignore_index=True)
        pbar.set_description("processing status ... ")
    
    return df

def spread_center(df:pd.DataFrame,flank:int,resolution:int):
    dfs=[]
    ## build a common dataframe
    rs=np.arange(-flank,flank+resolution,step=resolution)
    bins=np.arange(-(flank//resolution),flank//resolution+1,1)
    df1=pd.DataFrame({"start2":rs,"region2":bins})
    for r,b in zip(rs,bins):
        dfp=df1.copy()
        dfp["start1"]=r
        dfp["region1"]=b
        dfs.append(dfp)
    df_common=pd.concat(dfs)

    dfs=[]
    for _,row in df.iterrows():
        r=row["start"]
        for c in range(row["count"]):
            dfp=df_common.copy()
            dfp["start1"]=r+dfp["start1"]
            dfp["start2"]=r+dfp["start2"]
            dfs.append(dfp)
    dfo=pd.concat(dfs)
    return dfo


def process_straw(arg:argparse.ArgumentParser) -> pd.DataFrame:
    logger.info("start to generate center region by straw...")
    df_center=get_center_bins(arg.bed,arg.genome,arg.resolution)
    logger.info("start to dump info by straw...")

    dfs=[]
    for chrom, df in df_center.groupby(["chr"]):

        logger.info(f"process chrom {chrom} now ... ")
        logger.info(f"start to spread region on {chrom}")
        dft=spread_center(df,arg.flank,arg.resolution)
        logger.info("spreaded_region header: ")
        logger.info("\n##"+dft.head().to_csv(sep="\t").replace("\n","\n##"))

        logger.info(f"spread over, start to load hic file...")

        try:
            matrix=straw("KR",str(arg.hic),chrom, chrom, 'BP', arg.resolution)
        except Exception as e:
            logger.error(e)
            logger.error(f"{chrom} may not exist in the hic file")
            continue
        df_matrix=pd.DataFrame({"start1":matrix[0],"start2":matrix[1],"intensity":matrix[2]}).sort_values(["start1","start2"])


        logger.info(f"start to extrach info of  chrom {chrom} now ... ")
        dft=pd.merge(df_matrix,dft,on=["start1","start2"],how="inner")
        dft["chr"]=chrom
        dft["bin"]=arg.resolution

        dft=dft[["chr","start1","start2","region1","region2","bin","intensity"]].sort_values(["chr","start1","start2"])
        dfs.append(dft)
    dfo=pd.concat(dfs)

    return dfo

def run(arg:argparse.ArgumentParser):
    logger.info("start to process...")
    if arg.juicer=="straw":
        df=process_straw(arg)
    else:
        df=process(arg)
    logger.info("dump over, start to output...")
    df.to_csv(arg.fo,sep=arg.ofs,index=False,compression=arg.compression)
    logger.info("process over!")
    return df

if __name__=='__main__':
    arg=generate_opt()
    arg=validate_opts(arg)
    run(arg)


