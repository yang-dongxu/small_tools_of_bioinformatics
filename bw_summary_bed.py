#!/usr/bin/env python

from logging import warning
import numpy as np
import pyBigWig as pbw
import pandas as pd

import os
import sys

usage=f"## usage: python {__file__} bed6 bw [oname]"

try:
    assert len(sys.argv)>3
except:
    print(usage)
    sys.exit(1)

bed=sys.argv[1]
bigwig=sys.argv[2]
try:
    oname=sys.argv[3]
except:
    oname="std"

df_bed=pd.read_csv(bed,sep="\t",names=["#chr","start","end","name","score","strand"])
bw = pbw.open(bigwig)

def get_score(site:pd.Series,bw=bw):
    try :
        out=bw.stats(site["#chr"],site["start"],site["end"])[0]
        if out==None:
            out=0
        return out
    except:
        return 0

df_bed["score"]=df_bed.apply(get_score,axis=1)
#df_bed["score"]=df_bed["score"].astype(float)
if oname=="std":
    df_bed.to_csv(sep="\t",index=False)
else:
    df_bed.to_csv(oname,sep="\t",index=False)
sys.exit(0)