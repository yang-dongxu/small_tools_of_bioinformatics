#!/usr/bin/env python
import os
import sys
import re

import argparse
import logging
from typing import Pattern

import pandas as pd
from decorator import decorator

logger_name = "log extractor"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.INFO)

sh=logging.StreamHandler()
sh.setLevel(logging.DEBUG)

fmt = "%(asctime)-15s %(levelname)s %(name)s : pid - %(process)d :  %(message)s"
datefmt = "#%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
sh.setFormatter(formatter)
logger.addHandler(sh)


#@decorator
#def process(func, )

def generate_opt():
    opt=argparse.ArgumentParser()

    #chipseq related
    opt.add_argument('-b','--bowtie2',action="store",type=str,dest="bowtie2",default=None,help="bowtie2 logs name, with comma ',' as delim\n")
    opt.add_argument('-m','--macs',action="store",type=str,dest="macs",default=None,help="macs logs name, with comma ',' as delim\n")
    opt.add_argument('-l','--peak',action="store",type=str,dest="peak",default=None,help="peak_stat.py logs name, with comma ',' as delim\n")

    #rnaseq related
    opt.add_argument('-s','--star',action="store",type=str,dest="star",default=None,help="STAR logs name, with comma ',' as delim\n")
    opt.add_argument('-r','--read-distribution',action="store",type=str,dest="readsDistribution",default=None,help="peak_stat.py logs name, with comma ',' as delim\n")


    opt.add_argument('-p','--project',action="store",type=str,dest="project",default="",help="projects, if not provided, script will try to infer from logs input',' as delim\n")

    opt.add_argument('-o',"--outname",action="store",type=str,dest="oname",default="std",help="output file name, if std will print out(default) \n")
    opt.add_argument("--ofs",action="store",type=str,dest="ofs",default="\t",help="out file delimiter, default is \t as tsv, comma may also a choice to form csv \n")

    opt.add_argument("--debug",action="store_true",dest="debug",default=False,help="debug infos \n")

    return opt

def validate_args(args:argparse.ArgumentParser):
    ifiles={}

    flag=False
    if args.bowtie2 !=None and len(args.bowtie2.strip()):
        ifiles["bowtie2"]=args.bowtie2.split(",")
    if args.macs !=None and len(args.macs.strip()):
        ifiles["macs"]=args.macs.split(",")
    if args.star !=None and len(args.star.strip()):
        ifiles["star"]=args.star.split(",")
    if args.peak  !=None and len(args.peak.strip()):
        ifiles["peak"]=args.peak.split(",")
    if args.readsDistribution  !=None and len(args.readsDistribution.strip()):
        ifiles["readsDistribution"]=args.readsDistribution.split(",")

    if len(ifiles)==0:
        logger.critical("No input files! broken!")
        flag=True

    lengths=[]
    for key,value in ifiles.items():
        for file in value:
            if not os.path.isfile(file):
                logger.error(f"input file {file} in param {key} is not exist!")
                flag=True
        lengths.append(len(value))
        ## infer project
        projects=[os.path.split(i)[-1].split(".")[0] for i in value]
    
    if len(set(lengths))!=1:
        logger.error(f"input files num are not same, or files not paired! ")
        flag=True
    
    if  len(args.project.strip())==0:
        args.project=",".join(projects)
        logger.warning(f"set project info from input {key}, cause no info provided")
    projects=args.project.split(",")

    if len(lengths):
        if len(projects)!=lengths[0]:
            logger.error("Error in project name info! check your -p, or other params!")
            flag=True

    
    args.ifiles=ifiles
    if flag and not args.debug:
        sys.exit(1)
    return ifiles,projects

def not_definded(name):
    logger.error(f"You have input a not defined type of log: {name}")
    return {}

def parse_log(name,header,indexs,pattern):
    f=open(name)
    txt=f.read()
    f.close()
    nums=re.findall(pattern,txt)
    logger.debug(f"file: {name}\t nums:{nums}\theader:{header}\tindex:{indexs}")
    result={}
    for head,index in zip(header,indexs):
        result[head]=nums[index]
    return result    

def parser_bowtie2(name):
    pattern="\s*(\d+)\s[\(r]"
    header=["total","uniq","multi"]
    indexs=[0,3,4]
    result=parse_log(name,header,indexs,pattern)
    return result

def parser_macs(name):
    # pattern="total\sfragments\sin\streatment:\s(\d+)}|after\sfiltering\sin\streatment:\s(\d+)|Redundant\srate\sof\streatment:\s([0-9\.]+)"
    p1 = "total.* in treatment:\s(\d+)"
    p2 = "after\sfiltering\sin\streatment:\s(\d+)"
    p3 = "Redundant\srate\sof\streatment:\s([0-9\.]+)"
    pattern = "|".join([p1,p2,p3])

    try:
        header=["total","filtered","redundant"]
        indexs=[0,1,2]
        result=parse_log(name,header,indexs,pattern)
        for head,index in zip(header,indexs):
            result[head]=result[head][index]
    except Exception as e:
        logger.error(f"Error in parse {name} with {e}")
        logger.warning("only grab total reads info")
        header=["total"]
        indexs=[0]
        result=parse_log(name,header,indexs,p1)

    return result

def parser_star(name):
    p1="Number of input reads\s\|\s(\d+)\n"
    p2="Uniquely mapped reads number\s\|\s(\d+)\n"
    p3="Number of reads mapped to multiple loci\s\|\s(\d+)\n"
    pattern="|".join([p1,p2,p3])
    #logger.debug(pattern)
    header=["total","uniq","multi"]
    indexs=[0,1,2]
    result=parse_log(name,header,indexs,pattern)
    for head,index in zip(header,indexs):
        result[head]=result[head][index]
    return result    

def parser_peak(name):
    pattern="nums:\s(\d*)"
    header=["peaks","reads_in_peak"]
    indexs=[0,1]
    result=parse_log(name,header,indexs,pattern)
    return result

def extract_num(line):
    return [int(i) for i in re.findall("(\d\d+)",line)]

def parser_readsDistribution(name):
    with open(name) as f:
        data=f.read().splitlines()

    header=["total_tags","Extron_tags","Intron_tags","exon_ratio","intron_ratio"]
    total,cds,utr_5,utr_3,intron= data[1],data[5],data[6],data[7],data[8]
    values=[extract_num(i) for i in (total,cds,utr_5,utr_3,intron)]
    total_tags=values[0][0]
    extron_tags=values[1][1]+values[2][1]+values[3][1]
    intron_tags=values[4][1]

    exon_ratio=extron_tags/total_tags
    intron_ratio=intron_tags/total_tags
    results=[total_tags,extron_tags,intron_tags,exon_ratio,intron_ratio]
    result_dict={key:value for key,value in zip(header,results)}
    return result_dict

functions={
    "bowtie2":parser_bowtie2,
    "macs":parser_macs,
    "star":parser_star,
    "peak":parser_peak,
    "readsDistribution":parser_readsDistribution
    }


def run(ifiles={},projects=[]):
    results={i:{} for i in projects}
    for key,files in ifiles.items():
        func=functions.get(key,not_definded)
        for file, project in zip(files,projects):
            try:
                results[project][key]=func(file)
            except:
                results[project][key]={}
                logger.error(f"{project} {key} has wrong! see {file}")
    return results

def output(result:dict,arg:argparse.ArgumentParser,**kwargs):
    '''result is a dict, with {project:{part:{attribute:value}}} '''
    logger.debug(f"final result:\n{result}")
    lines=[]
    for project, info in result.items():
        line={"project":project}
        for part, attribute in info.items():
            for key,value in attribute.items():
                line[part+"_"+key]=value
        lines.append(line)
    logger.debug(f"out lines:\t{lines}")
    df=pd.DataFrame(lines)
    if args.oname=="std":
        print(df.to_csv(sep=args.ofs,index=False))
    else:
        df.to_csv(args.oname,sep=args.ofs,index=False)
    return 


if __name__=='__main__':
    opt=generate_opt()
    args=opt.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    ifiles,projects=validate_args(args)
    logger.debug(f"projects: {projects}")
    result=run(ifiles,projects)
    output(result,args)
    logger.info("process over!")

    
    