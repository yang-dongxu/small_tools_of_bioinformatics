#!/usr/bin/env python
# -*- coding: utf-8 -*-

# The script is used to calculate the array value from mnase signal.
# Array value is defined in the paper: doi:10/f5q8qx, describes how well-orgainze the nucleosome is
# Denpendency:
#   - python3
#   - numpy
#   - pandas
#   - numba
#   - bedtools  (bedtools makewindows)
#   - UCSCtools (bigWigAverageOverBed, bedGraphPack, bedGraphToBigWig)


import os
import sys
import pathlib
import click

import numpy as np
import numba as nb
import pandas as pd


# use bedtools, generate the bin_bed from genome and resolution
def generate_genome_window_bed(genome, resolution, bin_bed):
    if not os.path.isfile(genome):
        raise FileNotFoundError("The genome file is not found: {}".format(genome))
    if not os.path.isfile(bin_bed):
        cmd = f'''bedtools makewindows -g {genome} -w {resolution} | awk 'BEGIN{{OFS="\t"}}{{$4=$1 "," $2 "," $3; print $0}}'> {bin_bed}'''
        os.system(cmd)
    else:
        print("The windows file is already generated: {}".format(bin_bed))
    return bin_bed

# summarize the mnase signal from ibw in the bin_bed by bigWigAverageOverBed, to fix_step_signal_bed
def summarize_mnase_signal(ibw, bin_bed, fix_step_signal_bed):
    if not os.path.isfile(ibw):
        raise FileNotFoundError("The mnase signal file is not found: {}".format(ibw))
    if not os.path.isfile(bin_bed):
        raise FileNotFoundError("The windows file is not found: {}".format(bin_bed))
    if not os.path.isfile(fix_step_signal_bed):
        cmd = f"bigWigAverageOverBed {ibw} {bin_bed} {fix_step_signal_bed}"
        os.system(cmd)
    else:
        print("The fix step signal file is already generated: {}".format(fix_step_signal_bed))
    return fix_step_signal_bed

# turn output of bigWigAverageOverBed to a fixed step bedgraph  
def fix_step_signal_bedgraph_from_bed(fix_step_signal_bed, fix_step_signal_bedgraph, nafill=True):
    if not os.path.isfile(fix_step_signal_bed):
        raise FileNotFoundError("The fix step signal file is not found: {}".format(fix_step_signal_bed))
    if not os.path.isfile(fix_step_signal_bedgraph):
        if nafill:
            cmd = f'''cat  {fix_step_signal_bed} | cut -f 1,5 | tr $'\t' ',' > {fix_step_signal_bedgraph}'''
        else:
            cmd = f'''cat  {fix_step_signal_bed} | cut -f 1,6 | tr $'\t' ',' > {fix_step_signal_bedgraph}'''
        os.system(cmd)
    else:
        print("The fix step signal bedgraph file is already generated: {}".format(fix_step_signal_bedgraph))
    return fix_step_signal_bedgraph

# turn input ibw to a fixed step bedgraph
def ibw_to_fixed_step_bedgraph(ibw,resolution, genome, project, odir, nafill=True):
    ## steps: generate the bin_bed, summarize the mnase signal, generate the fix_step_signal_bed, turn it to a fixed step bedgraph

    ## project is the stem of ibw filename
    ## species is the stem of the genome file name
    species = os.path.basename(genome).split(".")[0]
    bin_bed = os.path.join(odir, species + "_" + str(resolution) + "bp.bed")
    fix_step_signal_bed = os.path.join(odir, project + "_" + str(resolution) + "bp_fix_step_signal.bed")
    fix_step_signal_bedgraph = os.path.join(odir, project + "_" + str(resolution) + "bp_fix_step_signal.bedgraph")

    # if fix_step_signal_bedgraph exist, skip the step
    if os.path.isfile(fix_step_signal_bedgraph):
        print("The fix step signal bedgraph file is already generated: {}".format(fix_step_signal_bedgraph))
    else:
        ## generate the bin_bed
        generate_genome_window_bed(genome, resolution, bin_bed)

        ## summarize the mnase signal
        summarize_mnase_signal(ibw, bin_bed, fix_step_signal_bed)

        # turn output of bigWigAverageOverBed to a fixed step bedgraph
        fix_step_signal_bedgraph_from_bed(fix_step_signal_bed, fix_step_signal_bedgraph, nafill=nafill)

        ## rm fix_step_signal_bed
        os.remove(fix_step_signal_bed)
    return fix_step_signal_bedgraph

# get gaussian smooth kernal
def get_gaussian_smooth_kernal(win, bandwidth):
    a = (np.arange(1,win*2)-win)/bandwidth
    kernal = np.exp(-0.5*a**2)
    kernal = kernal/np.sum(kernal)
    return kernal

# find localMax in given signal arrays
@nb.jit(nopython=True)
def find_local_max(signal:np.array, win:int):
    local_maxs = np.zeros_like(signal)
    for i in np.arange(win,len(signal)-win):
        this = signal[i]
        nearby = np.max(signal[(i-win):(i+win+1)])
        if this == nearby:
            local_maxs[i] = True
    return local_maxs

# get array_value from signal and local_maxs
@nb.jit(nopython=True)
def __get_array_value(signal:np.array, localMax:np.array):
    array_values = np.zeros_like(signal)
    max_index = np.where(localMax)[0]
    for a, b in zip(max_index[:-1],max_index[1:]):
        ya, yb = signal[a], signal[b]
        scores = np.arange(0,b-a)
        scores = scores*(yb-ya)/(b-a)+ya
        array_values[a:b] = scores
    return array_values

# load the fix_step_signal_bedgraph to a pandas dataframe, and calculate the array value
def get_array_value(bedgraph,resolution, window, smooth_band):
    # steps: load bedgraph, smooth, get diff and abs ,get local_maxs, get array_values

    # win is the window size // resolution
    win = window // resolution
    band = smooth_band // resolution

    if win==0:
        win = window
    if band==0:
        band = smooth_band

    ## get smooth kernal
    kernal = get_gaussian_smooth_kernal(win, band)
    print("The gaussian smooth kernal is: {}".format(kernal))
    
    print("Loading bedgraph file: {}".format(bedgraph))
    df = pd.read_csv(bedgraph, sep=",", header=None, names=["chr", "start", "end", "value"])
    df["value"] = df["value"].astype(float)

    print("Smoothing the signal...")
    df["smooth"] = np.convolve(df["value"], kernal, mode="same")
    df["diff"] = np.abs(df["smooth"].diff().fillna(0))

    print("Get array_value...")
    df["local_max"] = find_local_max(np.array(df["diff"]), win)
    df["array_value"] = __get_array_value(np.array(df["smooth"]), np.array(df["local_max"]))
    return df

# get the array value from the fix_step_signal_bedgraph, and output as a bedgraph, then pack it, finally output as a bigwig
def get_array_value_bedgraph(fix_step_signal_bedgraph, resolution, window, smooth_band, project, odir,genome):
    # steps: get array value, output a tmp bedGraph, then pack it, output a bigwig, rm tmp bedGraph

    df = get_array_value(fix_step_signal_bedgraph, resolution, window, smooth_band)
    tmp_bedgraph = os.path.join(odir, project + "_tmp_array_value.bedgraph")
    print("Outputing the array value bedgraph file: {}".format(tmp_bedgraph))
    df[["chr","start","end","array_value"]].to_csv(tmp_bedgraph, sep="\t", header=False, index=False)

    # pack the tmp bedgraph
    print("Packing the array value bedgraph file...")
    tmp_bedgraph_packed = os.path.join(odir, project + "_tmp_array_value.packed.bedgraph")
    cmd = f'''bedGraphPack {tmp_bedgraph} {tmp_bedgraph_packed}'''
    os.system(cmd)

    # sort bedgraph
    print("Sorting the array value bedgraph file...")
    tmp_bedgraph_sorted = os.path.join(odir, project + "_tmp_array_value.sorted.bedgraph")
    cmd = f'''sort -k1,1 -k2,2n {tmp_bedgraph_packed} > {tmp_bedgraph_sorted}'''
    os.system(cmd)


    ## output the bigwig
    print("Outputing the bigwig file...")
    bw = os.path.join(odir, project + f".{resolution}bp" + ".array_value.bw")
    cmd = f'''bedGraphToBigWig {tmp_bedgraph_sorted} {genome} {bw}'''
    os.system(cmd)

    ## rm tmps
    print("Removing the tmp files...")
    os.remove(tmp_bedgraph)
    os.remove(tmp_bedgraph_sorted)
    os.remove(tmp_bedgraph_packed)
    return bw





## get params from command line
@click.command(help="The script is used to calculate the array value from mnase signal. See doi:10/f5q8qx")
@click.option('-i','--ibw',type=click.Path(exists=True), help='The mnase signal file in bigwig format')
@click.option('-r','--resolution',type=int, help='The resolution of the windows')
@click.option('-g','--genome',type=click.Path(), help='The genome file in fasta format')
@click.option('-o','--odir',type=click.Path(), default=os.getcwd(),help='The output directory')
@click.option('-p','--project',default="",help='The project name')
@click.option('-w','--window',type=int, default=73,help='The half window size. Default is 73bp, half of DNA len in histone')
@click.option('-b','--bandwidth',type=int, default=30,help='The bandwidth of the gaussian smooth')
@click.option('-n','--nonafill',is_flag=True,flag_value=False, help='Whether to fill the NA value with 0', default=True)
def main(ibw, resolution, genome, odir, project, window, bandwidth, nonafill):
    if len(project.strip())==0:
        project = os.path.basename(ibw).split(".")[0]
    else:
        project = project.strip()
    bdg = ibw_to_fixed_step_bedgraph(ibw, resolution, genome, project, odir, nafill=not(nonafill))
    get_array_value_bedgraph(bdg, resolution, window, bandwidth, project, odir, genome)
    print("Done!")

if __name__ == "__main__":
    main()

