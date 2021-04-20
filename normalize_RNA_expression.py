import sys
import os
import re

import pandas as pd
import numpy as np

usage=f"python {__file__} inputname widthfile  [oname] "
'''if oname == std, then print out result as default '''

if len(sys.argv)>=3:
    pass
else:
    print(usage)
    sys.exit(1)


try:
    oname=sys.argv[3]
except:
    oname="std"


def output(oname,content):
    if oname=="std":
        print(content)
    else:
        with open(oname,'w') as f:
            f.write(content)
    return 0

inputname = str(sys.argv[1]).strip()
widthfile=str(sys.argv[2]).strip()

assert os.path.exists(inputname)
assert os.path.exists(widthfile)

df_width=pd.read_csv(widthfile,sep="\t",comment="#",names=["feature","width"],dtype={'feature': str, 'width': np.int}) # with col name of GeneID and length
df_reads=pd.read_csv(inputname,sep="\t",skip_blank_lines=True,comment="#",names=["feature","intensity"],dtype={'feature': str, 'intensity': np.float}) ## first column as GeneID, and second as intensity
df=pd.merge(df_reads,df_width,on="feature")

df["TPM"]=1e6*(df["intensity"]/df["width"])/np.sum(df["intensity"]/df["width"])
df["FPKM"]=1e9*df["intensity"]/np.sum(df["intensity"])/df["width"]

if oname=="std":
    print(df.to_csv(sep="\t",index=False))
else:
    df.to_csv(oname,sep="\t",index=False)
sys.exit(0)