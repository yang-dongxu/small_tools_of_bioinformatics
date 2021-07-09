#!usr/bin/env python

import os 
import sys
import argparse

import numpy as np
import pandas as pd


def read_homer_output(ifile) ->pd.DataFrame :
    df=pd.read_csv(ifile,sep="\t",skiprows=1,names=["repeats_id","chr","start","end","strand","length","copys","annotation","count"])
    df=df[["repeats_id","chr","start","end","strand","length","count"]]
    df["start"]=df["start"]-1 ## homer output is 1base
    #print(df.head())
    #df["count"]=df["count"]*1000000 it's already normalize to 1m 

    return df

def process_insertion(ifile, annotation)->pd.DataFrame:
    df=read_homer_output(ifile)
    df_gtf=pd.read_csv(annotation,sep="\t",names=["chr","start","end","name","score","strand"])
    dfo=df_gtf.merge(df,on=["chr","start","end","strand"])
    dfo=dfo[["chr","start","end","name","count","strand"]]
    return dfo

def process_subfamily(ifile)->pd.DataFrame: 
    df=read_homer_output(ifile)
    df["name"]=df["repeats_id"].str.split("\|",expand=True)[0]
    return df[["chr","start","end","name","count","strand"]]

usage=f"## python {__file__} [insertion|subfamily] ifile [annotation] "
if len(sys.argv) < 3:
    pring(usage)
    sys.exit()
kind=sys.argv[1]
ifile=sys.argv[2]

if kind=="insertion":
    if len(sys.argv)<4:
        print(usage)
        print("insertion mode has proviede a annotation bed!")
    annotation=sys.argv[3]
    dfo=process_insertion(ifile,annotation)
else:
    dfo=process_subfamily(ifile)

dfo.columns=["chr","start","end","name","score","strand"]
print(dfo.to_csv(sep="\t",header=True,index=False))

