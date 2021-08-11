import os
import sys
import argparse
import subprocess
import logging
from numpy.lib.arraysetops import intersect1d
import pandas as pd
import numpy as np


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

def generate_opt():
    opt=argparse.ArgumentParser()

    opt.add_argument("-a","--querys",dest="querys",type=str,action="append",required=True,help= '''query regions. Multi input files can add by comma seperated or multi -t options ''')
    opt.add_argument("-b","--databases",dest="dbs",type=str,action="append",required=True,help= ''' databases regions. Multi input files can add by comma seperated or multi -t options ''' )

    opt.add_argument("--query_tag",dest="query_tags",default=[],type=str,action="append",required=False,help= '''query region tags. Multi input files can add by comma seperated or multi -t options ''')
    opt.add_argument("--db_tag",dest="db_tags",type=str,default=[],action="append",required=False,help= '''query region tags. Multi input files can add by comma seperated or multi -t options ''')

    opt.add_argument("-s","--strand",dest="strand",default=False,action="store_true",help='''strand specific ''')
    opt.add_argument('-g','--genome',dest="genome",default="mm",help= ''' input genome size, can be mm for mouse, hg for human, or exact number. To scale output enrichment ''')
    opt.add_argument("-o","--output",dest="oname",default="enrichment.table",help='''out put file name''')
    opt.add_argument("-d","--ofs",dest="ofs",default="\t",help="output delim. defalut is \\t")


    args=opt.parse_args()

    return args

def select_genome(genome):
    genomes={
        "mm":1.87e9, ##mouse
        "hs":2.7e9, ## human
        "hg":2.7e9, ## human alias
        "ce": 9e7, ## C.elegans
        "dm": 1.2e8 ## fruitfly
    }

    try:
        a=float(genome)
    except ValueError as e:
        if genome in genomes:
            a=genomes[genome]
        else:
            logger.error("Non known genome! Input exact number. eg. 1e9, 2.7e9")
            sys.exit(255)
    return a


def validate_opt(opt:argparse.ArgumentParser):
    logger.info("validate options...")
    querys=[j for i in opt.querys for j  in i.split(",")]
    dbs=[j for i in opt.dbs for j  in i.split(",")]
    query_tags=[j for i in opt.query_tags for j  in i.split(",")]
    db_tags=[j for i in opt.db_tags for j  in i.split(",")]
    logger.info(f"input querys are {querys}")
    logger.info(f"input dbs are {dbs}")
    logger.info(f"input query tags are {query_tags}")
    logger.info(f"input db tags are {dbs}")

    opt.genome=select_genome(opt.genome)
    logger.info(f"input genome: {opt.genome}")


    if opt.strand:
        logger.info("Choose strand specific mode")
    else:
        logger.info("Choose strand ignorant mode")

    if len(querys) == 0:
        logger.error("No input querys! exit")
        sys.exit(255)
    if len(dbs) == 0:
        logger.error("No input querys! exit")
        sys.exit(255)

    if len(querys) != len(query_tags):
        logger.warning("query_tags are not pair to querys. Infer tags from input names!")
        query_tags=[".".join(os.path.split(i)[1].split(".")[:-1])  if "." in os.path.split(i)[1] else os.path.split(i)[1] for i in querys]
        logger.warning(f"query tags are inferred as {query_tags}")
    
    if len(dbs) != len(db_tags):
        logger.warning("db_tags are not pair to querys. Infer tags from input names!")
        db_tags=[".".join(os.path.split(i)[1].split(".")[:-1])  if "." in os.path.split(i)[1] else os.path.split(i)[1]  for i in dbs]
        logger.warning(f"query tags are inferred as {db_tags}")

    opt.querys=querys
    opt.query_tags=query_tags
    opt.dbs=dbs
    opt.db_tags=db_tags

    logger.info("validation over!")
    return opt

def output(data:pd.DataFrame,oname="std",sep="\t"):
    logger.info("start to output...")
    if oname=="std":
        fo=sys.stdout
    else:
        p=os.path.abspath(oname)
        p=os.path.split(p)[0]
        if not os.path.exists(p):
            os.makedirs(p)
        fo=open(oname,'w')
    fo.write(data.to_csv(sep=sep,index=False))
    logger.info(f"output over! see {oname}")
    return 0

def process(key,strand=False,genome_size=1):
    q_tag,query,db_tag,db =key
    tmp1=f".{q_tag}.{np.random.randn()}.{np.random.randn()}.tmp"
    tmp2=f".{db_tag}.{np.random.randn()}.{np.random.randn()}.tmp"
    cmd=f"cat {query} | bedtools merge -i - > {tmp1}"
    cmd+=f"&& cat {db} | bedtools merge -i - > {tmp2} "
    if strand:
        cmd+=f" \n bedtools intersect -a {tmp1} -b {tmp2} -s | "
    else:
        cmd+=f" \n bedtools intersect -a {tmp1} -b {tmp2}  | "
    cmd+=f" awk 'BEGIN{{s=0}}{{b=$3-$2;s+=b}}END{{print s}}' "
    overlap=subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout
    
    cmd=f"cat {tmp1} | awk 'BEGIN{{s=0}}{{b=$3-$2;s+=b}}END{{print s}}' "
    a_len=subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout

    cmd=f"cat {tmp2} | awk 'BEGIN{{s=0}}{{b=$3-$2;s+=b}}END{{print s}}' "
    b_len=subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout

    cmd=f"rm {tmp1} {tmp2}"
    subprocess.run(cmd,shell=True)
    overlap, a_len, b_len = int(overlap), int(a_len), int(b_len)
    logger.info(f"overlap, {q_tag} len, {db_tag} len: ( {overlap}, {a_len}, {b_len})")
    enrichment=(overlap/a_len)/(b_len/genome_size)
    return enrichment

def get_len(file):
    cmd=f"wc -l {file} | cut -f 1 | sed 's;^ *;;' | cut -f 1 -d ' ' "
    p=subprocess.run(cmd,capture_output=True, text=True)
    return int(p.stdout)

def get_intersect_len(a,b):
    cmd=f"bedtools intersect -a {a} -b {b} -wa -u | wc -l "
    p=subprocess.run(cmd,capture_output=True, text=True)
    a_intersected=int(p.std)

    cmd=f"bedtools intersect -a {b} -b {a} -wa -u | wc -l "
    p=subprocess.run(cmd,capture_output=True, text=True)
    b_intersected=int(p.std)
    return a_intersected, b_intersected


def run(args) ->pd.DataFrame:
    '''input a class with attribute querys, db, query_tags, db_tags, and output a pd.DataFrame '''
    logger.info("start to process...")
    querys=args.querys
    dbs=args.dbs
    query_tags=args.query_tags
    db_tags=args.db_tags
    datas=[]
    q_len={q_tag:get_len(q) for q_tag,q in zip(querys,query_tags) }
    db_len={db_tag:get_len(db) for db_tag,db in zip(dbs,db_tags) }
    for query,q_tag in zip(querys,query_tags):
        for db, db_tag in zip(dbs,db_tags):
            key=(q_tag,query,db_tag,db)
            logger.info(f"Process (q_tag,query,db_tag,db) : {(q_tag,query,db_tag,db)}.")
            enrichment=process(key,args.strand,args.genome)
            intersect_nums=get_intersect_len(query,db)
            data={"query_tag":q_tag,"db_tag":db_tag,"query":query,"db":db,"enrichment":enrichment,"query_len":q_len[q_tag],"db_len":db_len[db_tag],"query_overlap":intersect_nums[0],"query_overlap":intersect_nums[1]}
            logger.info(f"Process (q_tag,query,db_tag,db) : {(q_tag,query,db_tag,db)}. Enrichment: {enrichment:.3f}")
            datas.append(data)
    result=pd.DataFrame(datas)
    logger.info("process over")
    return result

if __name__=="__main__":
    opt=generate_opt()
    args=validate_opt(opt)
    data=run(args)
    output(data,args.oname,args.ofs)
    logger.info("process over! see you~")




