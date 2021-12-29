import os
import sys
import json
import re

import itertools
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio import PDB
from Bio.SeqUtils import seq1

import requests
import click
import logging

'''
#### Source from Shuang Hou #### 
Description: AlphaFold pLDDT<70, then use mobidb method to merge short streches and filter length
Usage: python get_idr_regions.py specie uniprot_id/uniprot_id_filename output_filename
Warning: can only calculate proteins in AlphaFold databse (length<2700 aa)
Example: python get_idr_regions.py human P23771 out.csv
#### Source from Shuang Hou #### 
'''

        
def intervals_extract(iterable):
    '''Convert list of sequential number into intervals'''
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]
        

def intervals_format(intervals):
    if len(intervals)==0:
        formated_intervals = ''
    else:
        formated_intervals = [str(i)+'..'+str(j) for i,j in intervals]
        formated_intervals = ','.join(formated_intervals)
    return formated_intervals


def get_seq_and_name(file):
    for record in SeqIO.parse(file, "pdb-seqres"):
        name =record.description.split(' ')[-1].split('_')[0]
    return str(record.seq), len(record.seq), name 


def getIdr(structure,threshold):
    '''position is 1-based'''
    idr = list()
    num = 0
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                num += 1
                if residue.has_id("CA"):
                    ca = residue["CA"]
                    if ca.get_bfactor() < threshold:
                        idr.append(num)
    return idr


def matMorphology(seq, rmax=3):
    # One, two and three residues-long ID stretches flanked on both sides by one, 
    # two and three residue-long ordered stretches are converted to order and vice- versa.
    
    # Disorder expansion
    seq = rmax*"D"+seq+rmax*"D"

    for r in range(1,rmax +1):
        pattern = 'D'*r + 'S'*r + 'D'*r
        new_pattern = 'D'*r + 'D'*r + 'D'*r

        for i in range(0, r +1):
            seq = seq.replace(pattern,new_pattern)

    # Disorder contraction
    seq = rmax*"S"+seq[rmax:-rmax]+rmax*"S"

    for r in range(1,rmax +1):
        pattern = 'S'*r + 'D'*r + 'S'*r
        new_pattern = 'S'*r + 'S'*r + 'S'*r

        for i in range(0, r +1):
            seq = seq.replace(pattern,new_pattern)

    return seq[rmax:-rmax]


def mobidb_run(seq,idr,length):
    state = ['D' if i in idr else 'S' for i in range(1,length+1)]
    state = ''.join(state)

    consensus = matMorphology(state,3)
    
    #Structured stretches of up to 10 consecutive residues are then converted to ID 
    #if they are flanked by two disordered regions of at least 20 residues
    flag=True
    while flag:
        m = re.search("D{21,}S{1,10}D{21,}", consensus)
        if m:
            matchLength = m.end(0) - m.start(0)
            consensus = consensus[:m.start(0)] + "D" * matchLength + consensus[m.end(0):]
        else: 
            flag = False
            
    position = [i for i in np.arange(1,length+1) if consensus[i-1]=='D']
    
    idr_intervals = list(intervals_extract(position))
    idr_pos = list()
    idr_seq = list()
    for i,j in idr_intervals:
        if j-i>=19:
            idr_seq.append(seq[(i-1):j])
            for pos in range(i,j+1):
                idr_pos.append(pos)
    
    idr_binary = np.zeros(length)
    idr_binary[list(np.array(idr_pos)-1)]=1
    idr_binary = [str(int(i)) for i in idr_binary]
    idr_intervals = list(intervals_extract(idr_pos))
    
    return ''.join(idr_binary),','.join(idr_seq),intervals_format(idr_intervals)

def download(uid,database,url):
    '''Download pdb file from alpha-fold database'''
    r = requests.get(url%uid)
    if r.status_code != 200:
        logging.error('%s not found in %s'%(uid,database))
    content = r.content.decode('utf-8')
    with open(os.path.join(database,uid+'.pdb'),'w') as f:
        f.write(content)
    return os.path.join(database,uid+'.pdb')

@click.command()
@click.option('-u','--uniprot_id', type=str, multiple=True, default = ["P49711"], help='Uniprot ID you want to query')
@click.option('-f','--file',is_flag=True, flag_value = True, default = False, help = "-u input is a file contain uids")
@click.option('-o','--output', type=click.File(mode = 'w'), default = sys.stdout, help="where to output")
@click.option('-r',"--remote", is_flag=True, flag_value = True, default = False, help="download from remote, use base url in url option. Result will be saved to local path in -d. If not choose, the script will try to find the file in local path")
@click.option('-l','--url', type=str, default = 'https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v2.pdb', help="url of alphaFold database")
@click.option("-d","--database", type=str, default = "alphaFoldDataBase", help="local database path")
@click.option("-t","--threshold", type=float, default = 70, help="threshold for IDR")
def run(uniprot_id,file,output,remote,url,database,threshold):
    if file:
        with open(uniprot_id[0]) as f:
            uniprot_id = [line.strip() for line in f]
    else:
        uniprot_id = uniprot_id
    logging.warning(f"input uids: {uniprot_id}")

    ## if remote, download them
    files = []
    if remote:
        for uid in uniprot_id:
            files.append(download(uid,database,url))
    else:
        for uid in uniprot_id:
            files.append(os.path.join(database,uid+'.pdb'))
    
    ## check exist
    for file in files:
        if not os.path.exists(file):
            logging.error('%s not found'%file)
    files = [file for file in files if os.path.exists(file)]
    uniprot_id = [uid for uid in uniprot_id if os.path.exists(os.path.join(database,uid+'.pdb'))]

    ## process part
    rows = []
    for file,uid in zip(files,uniprot_id):
        seq, length, name = get_seq_and_name(file)
        idr = getIdr(PDB.PDBParser(QUIET=True).get_structure("_",file),threshold)
        idr_binary, idr_seq, idr_positions = mobidb_run(seq, idr,length)
        idr_length = idr_binary.count('1')
        idr_percentage = round(idr_length/length,3)
        rows.append([name,uid,seq,length,idr_length,idr_percentage,idr_binary,idr_seq,idr_positions])
    df = pd.DataFrame(rows, columns=['protein_name', 'uid', 'sequence','length', 'idr_length','idr_percentage','idr_binary','idr_seq','idr_positions'])
    otable = df.to_csv(index=False)
    output.write(otable)
    return 0    


if __name__ == "__main__":
    run()
