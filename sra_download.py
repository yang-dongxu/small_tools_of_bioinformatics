import subprocess
import sys

with open(sys.argv[1]) as f:
    p=f.read().split("\n")

outdir=sys.argv[2]
sra_numbers=[i.strip() for i in p if len(i.strip())>1 ] 
print(sra_numbers)

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder input by argv[2]'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = f"fastq-dump --outdir {outdir}  --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip ~/ncbi/public/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
