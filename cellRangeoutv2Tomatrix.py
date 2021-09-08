#!/usr/bin/env python
import os
import io
import sys
import click
import gzip
from numpy.core.defchararray import index
import pandas as pd




def validate_opts(*args,**kwargs):
    print("validate")
    print(args)
    print(kwargs)
    pass

def gzip_wrap(*args,**kwargs):
    outs=[]
    for i in args:
        print(type(i).__bases__)
        assert isinstance(i,io.IOBase) ## has to be a file reader or writer
        if i.name.endswith("gz"):
            ig=gzip.open(i)
            outs.append(ig)
        else:
            outs.append(i)
    return outs

def get_odf(df,transpose,sparse,symbol):
    if symbol:
        gene="gene_symbol"
    else:
        gene="gene_id"
    if sparse:
        odf=df[["gene_id","gene_symbol","cell_barcode","umi"]]
    elif transpose:
        odf = df[[gene,"cell_barcode","umi"]].set_index([gene,"cell_barcode"]).unstack()["umi"].reset_index()
    else:
        odf = df[[gene,"cell_barcode","umi"]].set_index(["cell_barcode",gene]).unstack()["umi"].reset_index()
    return odf

@click.command()
@click.option("-g","--genes",type=click.File('rb'),required=True, help="cell range v2 out gene features, tsv format \n")
@click.option("-b","--barcodes",type=click.File('rb'),required=True, help="cell range v2 out barcode, a txt list, each line with a barcode \n")
@click.option("-m","--matrix",type=click.File('rb'),required=True,help="cell range v2 out matrix, a special tsv without considering first three lines. col 1: gene id. col2: barcode id. col 3: umi count\n")
@click.option("-o","--output",default="-",type=click.File('w'),help="output name\n, default is tsv format, with cols are genes and rows are cells. To transpose it, use -T params. To output sparse form, use -S flags. -S is higher at prioity \nDefault is - for stdout")
@click.option("-T","--transpose",default=False,is_flag=True,help="transpose output table, with rows as genes and cols as cells")
@click.option("-S","--sparse",default=False,is_flag=True,help="sparse form output table, tsv, 1st col: gene id, 2nd col: gene symbol, 3nd col: cell barcode/cell id, 4th col: umi count")
@click.option('-s',"--symbol",default=False,is_flag=True, help= "use gene symbol but not gene id as names")
def main(genes,barcodes,matrix,output,transpose,sparse,symbol,*args,**kwargs):
    genes,barcodes,matrix= gzip_wrap(genes,barcodes,matrix)

    barcodes_list=[i.strip() for i in barcodes.read().decode().split("\n") if len(i.strip())]
    genes_list=[ (i.split()[0], i.split()[1]) for i in genes.read().decode().split("\n") if len(i.strip()) ]
    matrix_list=[i.split() for i in  matrix.readlines()[3:]]
    
    #df_exp_matrix_sparse=pd.read_csv(matrix,skip_blank_lines=True,sep="\t",skiprows=3,names=["gene","cell","umi"])
    df_exp_matrix_sparse=pd.DataFrame(matrix_list,columns=["gene","cell","umi"])
    df_exp_matrix_sparse["gene"]=df_exp_matrix_sparse["gene"].astype(int)
    df_exp_matrix_sparse["cell"]=df_exp_matrix_sparse["cell"].astype(int)
    df_exp_matrix_sparse["umi"]=df_exp_matrix_sparse["umi"].astype(int)

    df_exp_matrix_sparse["cell_barcode"]=df_exp_matrix_sparse["cell"].apply(lambda x: barcodes_list[x-1])
    df_exp_matrix_sparse["gene_id"]=df_exp_matrix_sparse["gene"].apply(lambda x: genes_list[x-1][0])
    df_exp_matrix_sparse["gene_symbol"]=df_exp_matrix_sparse["gene"].apply(lambda x: genes_list[x-1][1])

    odf=get_odf(df_exp_matrix_sparse,transpose,sparse,symbol)
    output.write(odf.to_csv(sep="\t",index=False))

    ## clean up
    genes.close()
    barcodes.close()
    matrix.close()
    output.close()
    pass

if __name__ == '__main__':
    main()