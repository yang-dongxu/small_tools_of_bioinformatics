#!/usr/bin/env python
# coding: utf-8

import os
import sys

try:
    inputdir=sys.argv[1]
except:
    inputdir="/mnt/Storage2/home/zengshiyang/Work_YDX/1.PGC/j098_raw_result/CSPoutput/trimed"

def extract_project(name):
    project=name.split("_")[0]
    return project

def extract_suffix(name):
    suffix=name.split(".")[-1]
    return suffix

def determine_pair1(seq1,seq2):
    assert len(seq1)==len(seq2)
    for i in range(len(seq1)):
        if seq1[-i]==seq2[-i]:
            continue
        try:
            na=int(seq1[-i])
            nb=int(seq2[-i])
            if na<nb:
                ta,tb=seq1,seq2
            else:
                tb,ta=seq1,seq2
            return ta,tb
        except:
            continue
    return seq1,seq2




names={}
##names={project:{suffix1:[name1,name2]}}
for root,dirs,files in os.walk(inputdir):
    for file in files:
        project=extract_project(file)
        suffix=extract_suffix(file)
        pathname=os.path.join(os.path.abspath(root),file)
        if project in names:
            if suffix in names[project]:
                names[project][suffix].append(pathname)
            else:
                names[project][suffix]=[pathname]
        else:
            names[project]={suffix:[pathname]}



for project, suffix_files in names.items():
    for suffix, files in suffix_files.items():
        if len(files)==2:
            seqa=files[0]
            seqb=files[1]
            ta,tb=determine_pair1(seqa,seqb)
            seqpair=f"{project}\t{ta}\t{tb}"
            print(seqpair)
            



sys.exit(0)






