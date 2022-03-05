#!/usr/bin/env python

from email.policy import default
import os
import sys
import click 

import pandas as pd
import logging



@click.command()
@click.option('-i','--input',type=click.File('r'),default=sys.stdin,help="input seqpair table")
@click.option('-o','--odir',default="rawdata_refined",help="path to save output fastqs")
@click.option('--run/--no-run',is_flag=True,default=False,help="Run or not. It not set, default is output shell command. Otherwise, run the shell command.")
@click.option("-d",'--ifs',default=";",help="input file separator")
@click.option('-s','--sample_id_col',default=1,help="sample_id column,0-based")
@click.option('-p','--path_col',default=2,help="path column,0-based")
@click.option('-1','--fq1_col',default=3,help="fq1 column,0-based")
@click.option('-2','--fq2_col',default=4,help="fq2 column,0-based")
def process(input,odir,run,ifs,sample_id_col,path_col,fq1_col,fq2_col):
    """
    This script is used to prepare the input for HiCPro.
    """

    logging.info("Start to prepare HiCPro input")
    logging.info("Input seqpair table: {}".format(input.name))

    odir = os.path.abspath(odir) 

    logging.info("Output directory: {}".format(odir))
    if not os.path.exists(odir):
        os.makedirs(odir)

    df = pd.read_csv(input,sep=ifs,header=None).dropna()

    cmds = []
    for _,row in df.iterrows():
        sample_id = row[sample_id_col]
        sample_odir = os.path.join(odir,sample_id)
        if not os.path.exists(sample_odir):
            os.makedirs(sample_odir)
            
        path = row[path_col]
        fq1s = row[fq1_col]
        fq2s = row[fq2_col]
        oseq1 = os.path.join(sample_odir,"{}_R1.fastq".format(sample_id))
        oseq2 = os.path.join(sample_odir,"{}_R2.fastq".format(sample_id))
        

        if not os.path.exists(path):
            logging.error("Path {} does not exist".format(path))
            continue
        
        for (a,b) in zip(fq1s.split(","),fq2s.split(",")):
            fq1 = os.path.join(path,a)
            fq2 = os.path.join(path,b)
            if not os.path.exists(fq1):
                logging.error("Fq1 {} does not exist".format(fq1))
                continue
            if not os.path.exists(fq2):
                logging.error("Fq2 {} does not exist".format(fq2))
                continue
        cmd = f'''cd {path}; cat {fq1s.replace(',',' ')} > {oseq1}; cat {fq2s.replace(',',' ')} > {oseq2}'''
        cmds.append(cmd)
    
    cmds.append(f'''cd {os.getcwd()} ''')
    
    for cmd in cmds:
        print(cmd)

    if run:
        for cmd in cmds:
            logging.info("Run shell command: {}".format(cmd))
            os.system(cmd)

        

if __name__ == '__main__':
    process()