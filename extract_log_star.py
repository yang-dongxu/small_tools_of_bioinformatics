import os
import sys
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
    return line.split("|")[-1].strip()

def process(name,project):
    with open(name) as f:
        data=f.read().splitlines()
    total_line,unique_line, multi_line= data[5],data[8],data[23]
    total=extract_num(total_line)
    unique=extract_num(unique_line)
    multi=extract_num(unique_line)
    result={"project":project,"total_reads":total,"uniq_mapping_reads":unique,"multi_mapping_reads":multi}
    return deepcopy(result)


# start process
result=[process(name,project) for name,project in zip(logs,projects)]

header=["project","total_reads","uniq_mapping_reads","multi_mapping_reads"]
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
    
