import os
import sys 

import click
import pandas as pd



@click.command()
@click.option("-i","--ifile",type=click.File("r"), default = sys.stdin,help="input table")
@click.option("-o","--ofile",type=click.File("w"), default = sys.stdout,help="output table")
@click.option("--ifs", default=",",help = "input table seperator, default is comma")
@click.option("--ofs", default=",",help = "output table seperator, default is comma")
def transpose(ifile, ofile, ifs, ofs):
    table = {}
    for line in ifile:
        line = line.strip()
        if line:
            line = line.split(ifs)
            table[line[0]] = line[1:]
    df = pd.DataFrame(table)
    df.to_csv(ofile, sep=ofs, index=False)

if __name__ == '__main__':
    transpose()