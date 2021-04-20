#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
from copy import deepcopy

import numpy as np

usage=f"##usage: python {__file__} ifile1,ifile2,ifile3...  [project1,project2,project3...] [outname] [idname]\n## if outname is std, out will print in the tsv format, else write tsv in the given file, default is std \n ## this scripts vstack the table together and add a column describe the part comes from, and skip comment lines"

if len(sys.argv) <=2 or len(sys.argv)>5:
    print(usage)
if len(sys.argv)==1 or len(sys.argv)>5:
    sys.exit(0)


ifiles=[i for i in sys.argv[1].split(",") if len(i)] 
default_projects=[os.path.split(i)[1].split(".")[0].strip() for i in ifiles] ##split project name from filename

idname="project"

if len(sys.argv) ==2:
    projects=default_projects
    outname="std"
elif len(sys.argv)==3:
    term=sys.argv[2]
    if "," in term:
        projects = [i for i in sys.argv[2].split(",") if len(i)]
        outname="std"
    else:
        projects=default_projects
        outname=term
elif len(sys.argv)==4: 
    term=sys.argv[2]
    if "," in term:
        projects=[i for i in sys.argv[2].split(",") if len(i)]
        outname=sys.argv[3]
    else: ## 省略了projects这一项
        outname=sys.argv[2]
        idname=sys.argv[3]
        projects=default_projects
else:
    projects = [i for i in sys.argv[2].split(",") if len(i)]
    outname=sys.argv[3]
    idname=sys.argv[4]



dfs=[]

for i in range(len(ifiles)):
    filename=ifiles[i]
    project=projects[i]
    df=pd.read_table(filename,comment="#",skip_blank_lines=True)
    df[idname]=project
    dfs.append(df)

df_final=deepcopy(dfs[0])
for df in dfs[1:]:
    df_final=df_final.append(df,ignore_index=True)
df_final.fillna(value=0,inplace=True)

if outname=="std":
    print(df_final.to_csv(sep="\t",index=False))
else:
    df_final.to_csv(outname,index=False,sep="\t")
sys.exit()

