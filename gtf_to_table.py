import os
import sys
import argparse
import pandas as pd

## this script aims to transfer gtf format to table like tsv or csv, for downstream analysis.

opt=argparse.ArgumentParser()
opt.add_argument("ifile",default="-",action="store",type=str,help="input gtf file, with - means stdin. comment with # start will be skipped, and col name will appear at last line as comment.")
opt.add_argument("-o","--oname",default="-",dest="oname",action="store",type=str,help="output name, default is '-' as stdout, other path will be auto created")
opt.add_argument("-l","--lines",default=100,type=int,dest="lines",action="store",help="how many lines used to infer col names")
opt.add_argument("-f","--feature",default="transcript",dest="feature",action="store",type=str, help="which kind of data to keep: transcript, gene, or exon. dot(.) means all. ")
opt.add_argument("-d","--ofs",default="\t",dest="ofs",type=str,action="store",help="out put delim, default \\t ")

arg=opt.parse_args()

## start to validata arg
if arg.ifile=="-":
    fin=sys.stdin
else:
    if not os.path.isfile(arg.ifile):
        print("Error! Input file IS NOT exist! : {arg.ifile} ")
        sys.exit(1)
    fin=open(arg.ifile)

if arg.oname=="-":
    fout=sys.stdout
else:
    path=os.path.abspath(os.path.split(arg.oname)[0])
    if not os.path.exists(path):
        os.makedirs(path)
    fout=open(arg.oname,'w')

## start to process
#### load data
header_basic=["#chrom","source","feature","start","end","score","strand","phrase"] 
header_extend=[] ## attributes col name
lineID=0
data=[] ## lineId : { col: value}
for line in fin:
    lineID+=1
    lineSplit=line.split("\t") ## [chrom, source, feature, start, end, score, strand, phrase, attributes]
    
    if "#" in lineSplit[0]: ## is comment
        lineID-=1
        continue
    if arg.feature!=".":
        if arg.feature != lineSplit[2]:
            lineID-=1
            continue
    lineData={i:j for i,j in zip(header_basic,lineSplit)}
    attributes=[i.strip().replace("\"","") for i in lineSplit[-1].split(";") if len(i.strip().strip('"'))]

    if lineID<=arg.lines: ## infer extra headers
        for i in attributes:
            fs=i.split()
            #print(fs)
            key=fs[0]
            value=fs[1]
            if key not in header_extend:
                header_extend.append(key)
    for i in attributes:
        key=fs[0]
        value=fs[1]
        fs=i.strip().split()
        lineData[key]=value
    data.append(lineData)

#### output
header=header_basic+header_extend
data_df = pd.DataFrame(data)[header]
data_df.to_csv(fout, sep=arg.ofs, index=False) ###
fout.write("##"+arg.ofs.join(header) + '\n') ###

try:
    fin.close()
    fout.close()
except:
    pass
sys.exit(0)





