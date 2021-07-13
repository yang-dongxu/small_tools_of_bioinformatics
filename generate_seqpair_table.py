import os
import sys
import re
import pandas as pd
import numpy as np


import argparse


def generate_opt():
    opt = argparse.ArgumentParser(
        usage="helps to build a seqpair table for totalAnalysisPipline")

    opt.add_argument(
        "-i", "--idir", required=True, dest="idir", help="where you fastq files/ input files are in")
    opt.add_argument(
        "-f", "--filter", dest="filter",type=str, default="fq.gz", help="the re pattern to filter")

    arg = opt.parse_args()
    return arg


def get_score(a: str, b: str):
    score = 0
    for i, j in zip(a, b):
        if i == j:
            score += 1
        else:
            break
    return score


def get_commen(*args):
    end = ""
    if len(args) == 2:
        score = get_score(args[0], args[1])
        end = args[0][:score]
    elif len(args) == 1:
        end = args[0].split(".")[0]
    return end


def get_files(idir: str, pattern: str):
    files = list(os.listdir(idir))
    files = sorted([i for i in files if re.search(pattern, i)])
    return files


def get_pair_score(files):
    a = len(files)
    scores = np.zeros((a, a))
    for i in range(a):
        for j in range(i + 1, a):
            scores[i, j] = get_score(files[i], files[j])
    pairs = get_pair(scores, [])
    out_pairs = [[files[i[0]], files[i[1]]] for i in pairs]
    all_files = set([j for i in out_pairs for j in i])
    left_files = set(files) - all_files
    out_pairs += [[i] for i in list(left_files)]

    return out_pairs


def get_pair(scores: np.array, results=[]):
    if np.sum(scores) <= 0:
        return results
    else:
        a = np.max(scores, axis=1)
        i = np.argmax(a)
        j = np.argmax(scores[i, :])
        scores[i, :] = 0
        scores[j, :] = 0
        scores[:, i] = 0
        scores[:, j] = 0
        results.append([i, j])
        #print(f"pair: {(i,j)}")
        return get_pair(scores, results=results)
    pass


def out(result, idir):
    for i in result:
        line = ";".join(i)
        project = get_commen(*i)
        line = f'{project};{idir};{line}'
        print(line)
    return


def run(arg):
    idir = arg.idir
    idir = os.path.abspath(idir)
    files = get_files(idir, arg.filter)
    results = get_pair_score(files)
    out(results, idir)
    return idir, files


if __name__ == '__main__':
    arg = generate_opt()
    run(arg)
