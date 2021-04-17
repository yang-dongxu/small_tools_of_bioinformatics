import sys
import os
import datetime

usage=f"python {__file__} peak.bed bam [oname]"
if len(sys.argv)<3:
    print("## not enough params!",file=sys.stderr)
    print(f"## usage {usage}",file=sys.stderr)
    sys.exit(1)

peak=os.path.abspath(sys.argv[1])
bam=os.path.abspath(sys.argv[2])

if not os.path.isfile(peak):
    print(f"## ERROR! {peak} is not exist!",file=sys.stderr)
    sys.exit(1)
if not os.path.isfile(bam):
    print(f"## ERROR! {bam} is not exist!",file=sys.stderr)
    sys.exit(1)


try:
    oname=sys.argv[3]
    fout=open(oname,'w')
except:
    oname="std"
    fout=sys.stdout

def now():
    a=datetime.datetime.now()
    return a.strftime("%Y-%m-%d %H:%M:%S")

print(f"## input peak file: {peak}")


print(f"## start to cal peak nums at {now()}",file=fout)
cmd_peak=f"cat {peak} | wc -l"
peaks_num=int(os.popen(cmd_peak).read())
print(f"peak nums: {peaks_num}",file=fout)

print(f"## start to read nums in peaks at {now()}",file=fout)
cmd_frip=f"bedtools intersect -wa -u -a {bam} -b {peak}  | wc -l "
reads_in_peaks=int(os.popen(cmd_frip).read())
print(f"read in peak nums: {reads_in_peaks}",file=fout)
print(f"## end process at {now()}",file=fout)



