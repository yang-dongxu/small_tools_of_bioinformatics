#!/usr/bin/env python

import os
import sys
import argparse
import logging

import pandas as pd
import numpy as np


logger_name = "log extracter"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.DEBUG)

sh=logging.StreamHandler(sys.stderr)
sh.setLevel(logging.DEBUG)

fmt = "%(asctime)-15s %(levelname)s %(name)s : pid - %(process)d :  %(message)s"
datefmt = "#%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
sh.setFormatter(formatter)
logger.addHandler(sh)

def generate_opt():
    opt=argparse.ArgumentParser(prog="tableStack.py")
    opt.add_argument

    subs=opt.add_subparsers(help='sub-command help')

    generate_opt_vstack(subs)

    return opt

def generate_opt_vstack(sub:argparse.ArgumentParser):
    opt=sub.add_parser("vstack",help="stack tables vertical")
    assert isinstance(opt,argparse.ArgumentParser)

    opt.add_argument("-i","--idir",dest="idir",type=str,action="store",help="if your tables are all in a dir, then this dir only contain your tables, this may helps.")
    opt.add_argument("-t","--table",dest="tables",default=[],type=str,action="append")
    opt.add_argument("-p","--project",dest="projects",default=[],type=str,action="append")
    opt.add_argument('-o','--oname',dest="oname",default="-",type=str,action="store",help="out file name, '-' means stdout")

    opt.add_argument('-n','--header',dest="header",default=0,type=int,action="store",help="0: first row is header, -1: no header in table")
    opt.add_argument('--names',dest="names",type=str,default="",action="store",help="provide header by hand, use comma(,) as delim ")
    opt.add_argument("--comment",dest="comment",default="#",action="store",help="ignore comment lines start with # ")
    opt.add_argument("--skip",dest="skip",default=0,type=int,action="store",help="skip first n columns ")
    opt.add_argument('-d',"--sep",dest="sep",default="\t",type=str,action="store",help="input file delim. (for tsv, set like this : -d $'\t'. for csv, use comma(,) )")

    opt.add_argument('-k','--key',dest="key",default="project",type=str,action="store",help="the column name of projects")

    opt.add_argument("--out-header",dest="oheader",default=True,action="store_true",help="whether output table with header")
    opt.add_argument("-s","--sortby",dest="sortby",default=[],type=list,action="append",help="out table will sort by col names in given order")
    opt.add_argument("--ascending",dest="ascending",default="T",type=str,action="append",help="ascending or not. set only one value will apply to all keys, otherwise it have to be same with --sortby terms, can be add many times, by comma or --ascending again. T: ascending. F:descending")
    opt.add_argument("--ofs",dest="ofs",default="\t",type=str,action="store",help="output file delim")


    opt.set_defaults(func=Vstack)
    return opt


class Stack:
    def __init__(self,args) -> None:
        self.args=args
        #logger.info("choose vstack sub-command")

        self.tables=[]
        self.projects=[]
        self.dfs=[]
        self.df=pd.DataFrame()

        self.validate_opt()
        self.get_tables()

    def validate_opt(self):
        logger.info("start to validate options")
        tables=self.args.tables
        idir=self.args.idir
        projects=self.args.projects

        ## process input tables 
        if len(tables)==0:
            logger.warning("your have not input tables by -t ")
            if not os.path.isdir(idir):
                logger.error(f"your input dir {idir} does not exist! please checks it ! exit.")
                sys.exit(1)
            tables=os.listdir(idir)
            if len(tables)==0:
                logger.error(f"your input dir {idir} is empty! please checks it ! exit.")
                sys.exit(1)
        else:
            logger.info(f"input tables: {tables}")
            tables=self.parse_list(tables)
        self.tables=tables
        logger.info(f"determined input tables: {tables}")

        error_files=[]
        for i in tables:
            if not os.path.isfile(i):
                error_files.append(i)
        if len(error_files):
            logger.error(f"input-files not exist: ERROR FILES: {error_files}")
            sys.exit(1)
        else:
            del error_files


        ## process input projects
        if len(projects)==0:
            logger.warning("your have not input projects by -p, so infer from tables information")
            projects=[os.path.split(i)[-1].split(".")[0] for i in tables]
        else:
            logger.info(f"input projects: {projects}")
            projects=self.parse_list(projects)
        if len(set(projects)) != len(projects):
            logger.warning("you have duplicated projects!")
        self.projects=projects
        logger.info(f"determined projects: {projects}")
        if len(projects) != len(tables):
            logger.error("tables and projects are NOT paired! EXIT!")
            sys.exit(1)
        
        logger.info("options validate over, start to read tables ...")

        ## determine comment type
        if len(arg.comment) == 1:
            logger.info(f"your comment is {arg.comment}")
        elif len(arg.comment) == 2 and arg.comment[0] =="\\":
            logger.info(f"your comment is {arg.comment}")
        elif arg.comment == "none":
            arg.comment=None
            logger.info(f"your comment is {arg.comment}")
        elif len(arg.comment) >=2 :
            logger.warning(f"your input comment {arg.comment} is too longer! Only single character or backsplash once is support! Turn comment to NA")
            arg.comment=None
        return True

    def get_tables(self):
        logger.info("start to load tables")
        dfs=[]
        header=self.args.header
        names=[ i.strip()  for  i in self.args.names.split(",") if len(i.strip())>0]
        comment=self.args.comment
        skip=self.args.skip
        sep=self.args.sep
        key=self.args.key
        if header <0:
            header=None
            if len(names) ==0:
                logger.info("set none header")
            else:
                logger.warning(f"set header as: {names}")
        else:
            logger.info(f"set header as line {header}")
        cols=[]
        for table,project in zip(self.tables,self.projects):
            try:
                logger.info(f"try to load project {project} with table {table}, sep as {sep}")
                sep="%s"%sep
                if header !=None:
                    df=pd.read_csv(table,sep=sep,comment=comment,skip_blank_lines=True,skiprows=skip,header=header)
                else:
                    df=pd.read_csv(table,sep=sep,comment=comment,skip_blank_lines=True,skiprows=skip,header=None)
                    if len(names):
                        if len(names)==len(df.columns):
                            df.columns=names
                        elif len(names)<=len(df.columns):
                            logger.error(f"input names length < table cols, so name won't be set!")
                        else:
                            new_names=names[0:len(df.columns)]
                            logger.error(f"input names length > table cols, so names are choosed first {len(df.columns)}: {new_names}!")
                            df.columns=new_names

                if key in df.columns:
                    logger.warning(f"your project col name {key} is duplicate with exist columns!")
                df[key]=project
                cols.append(set(df.columns))
                dfs.append(df)
                logger.debug(f"project {project} table {table}:\n{df.head().to_string()}")
            except:
                logger.error(f"parse table {table} of project {project} faild, check it!")
        if len(dfs)==0:
            logger.error(f"no table loaded sucessfully, check them!")
            sys.exit(1)
        intersect_cols=set.intersection(*cols)
        total_cols=set.union(*cols)
        logger.info(f"common cols: {intersect_cols}")
        logger.info(f"all cols together: {total_cols}")
        if len(intersect_cols)<=1:
            logger.warning(f"too little common cols! common: {len(intersect_cols)}, total: {len(total_cols)}")
        logger.info("load tables finined, start to process")
        self.dfs=dfs
        return True
        
    def process(self):
        logger.warning("no process method in base class!")

    def parse_list(self,term):
        term_expand=[i.split(",") for i in term] ## 展开tables
        terms=[j for i in term_expand for j in i]
        return terms

    def sort(self):
        sortby=self.args.sortby
        ascending=self.args.ascending
        df=self.df
        if isinstance(ascending,str):
            ascending= ascending=="T"
        else:
            ascending=self.parse_list(ascending)
            if len(ascending) >=len(sortby):
                logger.warning(f"length of ascend is larger than sortby keys! cut first {len(sortby)}")
                ascending=ascending[0:len(sortby)]
                logger.warning(f"choosed ascending: {ascending}")
                ascending=[i =="T" for i in ascending ]
            else:
                ascending=ascending[0]*len(sortby)
                logger.error(f"length of ascend is smaller than sortby keys! choose first one: {ascending[0]}!")
        logger.info(f"determine ascending: {ascending}") 
        if len(sortby):
            incongruous_col_names=[]
            right_cols=[]
            orders=[]
            for col,order in zip(sortby,ascending):
                if col not in df.columns:
                    incongruous_col_names.append(col)
                else:
                    right_cols.append(col)
                    orders.append(order)
            if len(incongruous_col_names):
                logger.error(f"{incongruous_col_names} are not appear in col names of given tables")
            if len(right_cols):
                logger.info(f"output table will be sorted by the order {right_cols}")
                df.sort_values(by=right_cols,axis=0,ascending=orders)
        self.df=df
        return df
    
    def output(self):
        ofs=self.args.ofs
        if ofs=="\\t":
            ofs=="\t"
        oheader=self.args.oheader
        df=self.sort()
        if self.args.oname=="-":
            f=sys.stdout
        else:
            name=self.args.oname
            path=os.path.split(os.path.abspath(name))[0]
            if not os.path.isdir(path):
                os.makedirs(path)
                logger.warning(f"output path {path} you given is not exisit, the script create it")
            elif os.path.isfile(name):
                logger.warning(f"output name {name} is exist, this output will overwrite it.")
            f=open(name,'w')
        sys.stdout=f

        out=self.df.to_csv(sep=ofs,header=oheader,index=False,encoding="utf8")
        try:
            sys.stdout.write(out)
        except BrokenPipeError:
            logger.error("Errno 32, Broken pipe. output stop.")
        f.close()
        logger.info(f"out name: {self.args.oname}")

class Vstack(Stack):
    def __init__(self, args) -> None:
        super().__init__(args)
        logger.info("choose vstack sub-command")
        self.process()

    def process(self):
        logger.info("vstack process!")
        if len(self.dfs)>=2:
            self.df=pd.concat(self.dfs)
        else:
            self.df=self.dfs[0]
        self.output()


if __name__=='__main__':
    opt=generate_opt()
    arg=opt.parse_args()
    arg.func(arg)
    logger.info("process over, see you ~")