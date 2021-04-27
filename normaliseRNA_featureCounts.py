import os
import sys

import pandas as pd
import numpy as np


usage="##usage: python {__file__} kind inputfile oname group_key(optional,if not provided, will take as a single sample) \n## kind: FPKM,TPM; \n##inputfile: with header feature	Chr	Start	End	Strand	Length	Reads [group_key]\n## group_key: the column name of group_key"

featureCounts_Reads_name="Reads"
out_Column_name="intensity"

if len(sys.argv) <=4 or len(sys.argv)>5:
    print(usage)
if len(sys.argv)<4 or len(sys.argv)>5:
    sys.exit(0)

kind=sys.argv[1]
ifile=sys.argv[2]
oname=sys.argv[3]
assert (kind in ["TPM", "FPKM"])
assert os.path.exists(ifile)

if len(sys.argv)>=5:
    group_key=sys.argv[4]
    need_group=True
else:
    group_key="ADD_TO_COMPATIBLE"
    need_group=False

df=pd.read_table(ifile,sep="\t",skip_blank_lines=True,comment="#")
if need_group==False:
    df[group_key]=group_key

def cal_FPKM(idf:pd.DataFrame,group_key):
    df=idf
    df_t=df.groupby(group_key).agg({featureCounts_Reads_name:"sum"}).rename(columns={featureCounts_Reads_name:"total"})
    df=df.merge(df_t,how="left",on=group_key)
    df["FPKM"]=(df[featureCounts_Reads_name]/df["total"])/df["Length"]*1e9
    name=df.columns[0]
    df_out=df[[name,"FPKM",group_key]]
    return df_out.rename(columns={"FPKM":out_Column_name})

def cal_TPM(idf:pd.DataFrame,group_key):
    df=idf
    df["RPK"]=1e3*df[featureCounts_Reads_name]/df["Length"]
    df_t=df.groupby(group_key).agg({"RPK":"sum"}).rename(columns={"RPK":"total"})
    df=df.merge(df_t,how="left",on=group_key)
    df["TPM"]=1e6*df["RPK"]/df["total"]
    name=df.columns[0]
    df_out=df[[name,"TPM",group_key]]
    return df_out.rename(columns={"TPM":out_Column_name})

def process(idf:pd.DataFrame,need_group,group_key:str,kind):
    if kind=="TPM":
        df= cal_TPM(idf,group_key)
    elif kind =="FPKM":
        df=cal_FPKM(idf,group_key)
    else:
        raise TypeError("Not supported normalized kind : {kind}")
    if need_group:
        return df
    else:
        return df.drop(axis=1,columns=group_key)

result=process(df,need_group,group_key,kind)
if oname=="std":
    print(result.to_csv(sep="\t",index=False))
else:
    result.to_csv(oname,sep="\t",index=False)

sys.exit()

