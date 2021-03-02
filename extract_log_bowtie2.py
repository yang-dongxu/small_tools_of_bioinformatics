# %%
import sys
import os
import re
# %%

usage=f"python {__file__} logname outname [mapping_part_name] [peakcalling_part_name]"

if len(sys.argv) <3:
    print(usage)
    sys.exit(1)

log_name=sys.argv[1]
outname=sys.argv[2]

mapping_part=sys.argv[3] if len(sys.argv) >=4 else "bowtie2"
peakcalling_part=sys.argv[4] if len(sys.argv) >=5 else "peakcalling"

f=open(log_name)
text=f.read()
f.close()

f=open(outname,'w')
columns=["project","total_reads","uniq_mapping_reads","multi_mapping_reads","filtered_reads","mapping_ratio"]
f.write(",".join(columns)+"\n")

# %%
aligment_pattern=f"#{{5}}\s*start\s*(\S+)\s*{mapping_part}\s*#{{5}}[\s\S]+?(\d+ reads[\s\S]+?alignment rate)\s#{{5}}\s*stop (\S+) {mapping_part}\s*#{{5}}"
macs_pattern="#{{5}}\s*start\s*{project}\s*{peakcalling_part}\s*#{{5}}([\s\S]+?)#{{5}}\s*stop\s*{project}\s*{peakcalling_part}\s*#{{5}}"
# %%
all_mapping_results=re.findall(aligment_pattern,text)
for p1,body,p2 in all_mapping_results:
    if p1.strip() !=p2.strip():
        continue
    print(p1)

    ### get filtered reads num from macs2 log
    macs_pattern_project=macs_pattern.format(project=p1,peakcalling_part=peakcalling_part)
    macs_log=re.findall(macs_pattern_project,text)
    assert len(macs_log)==1
    extract_peak_reads_pattern="fragments after filtering in treatment: (\d+)"
    filtered_reads=re.findall(extract_peak_reads_pattern,macs_log[0])[0]


    ##### get aligment reads num from bowtie2 log
    parts=[i.strip() for i in body.split('----') ]
    assert len(parts)>= 2
    total_reads=parts[0].split('\n')[0].strip().split()[0]
    uniq_reads=parts[0].split('\n')[3].strip().split()[0]
    dup_reads=parts[0].split('\n')[4].strip().split()[0]
    mapping_ratio=(int(uniq_reads)+int(dup_reads))/int(total_reads)
    result="{project},{total},{unique},{dup},{filter_unique},{mapping_ratio}\n".format(project=p1,total=total_reads,unique=uniq_reads,dup=dup_reads,filter_unique=filtered_reads,mapping_ratio=mapping_ratio)
    f.write(result)

f.close()
print("##alignment stats complete!")
sys.exit()



# %%
