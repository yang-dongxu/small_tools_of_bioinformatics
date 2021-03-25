import os
import sys

usage=f"python {__file__} inputfile"
'''
del the random and alt chromsomes in the given fasta file, and print result line by line
'''
if len(sys.argv) !=2:
    print(usage)
    sys.exit(1)

name=sys.argv[1]

f=open(name)

flag=True
for line in f:
    if ">" in line:
        if "_" in line:
            flag=False
        else:
            flag=True
    if flag:
        print(line,end="")

f.close()