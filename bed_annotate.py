#!/usr/env/bin python
import click
from click.termui import prompt
import pandas
import sys
#import pybedtools

@click.group()
def annotate():
    click.echo("Welcome to use annotate peak functions!")
    pass

def echo(*args):
    for i in args:
        click.echo(f"{i}", err=True)

@click.command("annotateByGenepredExt")
@click.option("-a","--annotation",required=True)
@click.option('-b',"--bed",required=True)
@click.option("-m","--mode",type=click.Choice(["full","matrix"]),default="full", help="output format. Full means output not only overlap matrix,\
     but also gene list. matrix means only binding matrix. ")
@click.option("-o","--oname",default="annotate_results",help="where to store your output. for full mode, it determines the output dir name. for matrix mode, matrix name.")
@click.option("-p","--project",default="",help="When full mode choosed, output files will be named as project.gene.list or so on. Ignored on matrix mode")
def annotateByGenepredEXT(annotation, bed, mode, oname, project):
    sys.stderr.write("## start to annotate peak!")

annotate.add_command(annotateByGenepredEXT)
if __name__ == '__main__':
    annotate()
    