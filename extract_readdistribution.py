import os
import sys
import re
from copy import deepcopy

usage="##usage: python {__file__} log1,log2,log3...  [project1,project2,project3...] [outname]\n## if outname is std, out will print in the tsv format, else write tsv in the given file, default is std"

if len(sys.argv) <=2 or len(sys.argv)>4:
    print(usage)
if len(sys.argv)==1 or len(sys.argv)>4:
    sys.exit(0)

## process input files
logs=[i for i in sys.argv[1].split(",") if len(i)] 
default_projects=[os.path.split(i)[1].split(".")[0].strip() for i in logs] ##split project name from filename

if len(sys.argv) ==2:
    projects=default_projects
    outname="std"
elif len(sys.argv)==3:
    term=sys.argv[2]
    if "," in term:
        projects = term.split(",")
        outname="std"
    else:
        projects=default_projects
        outname=term
else:
    projects = [i for i in sys.argv[2].split(",") if len(i)]
    outname=sys.argv[3]



# def used funtions
def extract_num(line):
    return [int(i) for i in re.findall("(\d\d+)",line)]

def process(name,project,header):
    with open(name) as f:
        data=f.read().splitlines()
    total,cds,utr_5,utr_3,intron= data[1],data[5],data[6],data[7],data[8]
    values=[extract_num(i) for i in (total,cds,utr_5,utr_3,intron)]
    total_tags=values[0][0]
    extron_tags=values[1][1]+values[2][1]+values[3][1]
    intron_tags=values[4][1]

    exon_ratio=extron_tags/total_tags
    intron_ratio=intron_tags/total_tags
    results=[project,total_tags,extron_tags,intron_tags,exon_ratio,intron_ratio]
    result_dict={key:value for key,value in zip(header,results)}
    return deepcopy(result_dict)


# start process
header=["project","total_tags","Extron_tags","Intron_tags","exon_ratio","intron_ratio"]

result=[process(name,project,header) for name,project in zip(logs,projects)]

line="\t".join(["{"+i+"}" for i in header])+"\n"

format_result="\t".join(header)+"\n"
for item in result:
    format_result+=line.format_map(item)

if outname=="std":
    print(format_result)
else:
    with open(outname,'w',encoding='utf8') as f:
        f.write(format_result)
sys.exit()
    
