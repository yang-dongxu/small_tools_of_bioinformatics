
import os
import sys

import click
import logging

logger_name = "fragments to pseudo reads"
logger = logging.getLogger(logger_name)
logger.setLevel(logging.INFO)

sh=logging.StreamHandler()
sh.setLevel(logging.DEBUG)


fmt = "%(asctime)-15s %(levelname)s %(name)s : pid - %(process)d :  %(message)s"
datefmt = "#%a %d %b %Y %H:%M:%S"
formatter = logging.Formatter(fmt, datefmt)
sh.setFormatter(formatter)
logger.addHandler(sh)

@click.command()
@click.option('-i','--ifile',type=click.File("r"), default = sys.stdin, help="input fragments file, stdin supported, at least first 3 columns are required")
@click.option('-o','--ofile',type=click.File("w"), default = sys.stdout,help="output pseudo reads file, stdout supported")
@click.option('-l','--length',type=int, default = 150,help="length of pseudo reads")
@click.option("-p",''"--prefix",type=str, default = "pseudo_",help="prefix of pseudo reads")
def run(ifile,ofile,length,prefix):
    for lid, line in enumerate(ifile):
        if line.startswith("#"):
            continue
        line=line.strip()
        if len(line)==0:
            continue
        items=line.split("\t")
        if len(items)<3 :
            logger.error("invalid line: %s",line)
            logger.error("at least 3 columns are required")
            sys.exit(1)
        chrom = items[0]
        start = int(items[1])
        end = int(items[2])
        name = items[3] if len(items) > 3 else "%d" % (lid)
        name = f"{prefix}{name}"

        if start > end:
            logger.error("invalid line: %s",line)
            logger.error("start should be less than end")
            sys.exit(1)
        if end - start < length:
            logger.warning("pseudo reads length is longer than fragment length, line: %s",line)
            logger.warning(f"fragments {name} will be droped")

        # generate pseudo reads
        ofile.write(f"{chrom}\t{start}\t{start+length}\t{name} /1\n")
        ofile.write(f"{chrom}\t{end-length}\t{end}\t{name} /2\n")
    
if __name__ == "__main__":
    run()
        
