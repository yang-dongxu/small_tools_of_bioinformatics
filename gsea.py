import os 
import sys 
import numpy as np 
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests



def  get_ES(rank:pd.Series, genelist:list, signalAsWeight:bool = True):
    '''
    rank: pd.Series, index is gene name, value is rank
    genelist: list, gene list
    signalAsWeight: bool, if True, use the rank as weight, otherwise use 1
    '''
    # Sort genes by score in decreasing order
    rank = rank.copy().sort_values(ascending=False)

    # signalAsWeight to int
    signalAsWeight = 1 if signalAsWeight else 0

    # calculate Nr to normalize the rank

    Nr = np.sum([
     abs(rank.iloc[i])**signalAsWeight 
        for i in range(len(rank)) if rank.index[i] in genelist
     ])

    # get hits and non-hits
    hits = [i for i in rank.index if i in genelist]
    hits_idx = [i for i in range(len(rank)) if rank.index[i] in genelist]
    non_hits = [i for i in rank.index if i not in genelist]

    # calculate running ES
    running_es = [0]
    for i,v in rank.items():
        weight = abs(v) ** signalAsWeight / Nr
        if i in hits:
            new_es = running_es[-1] + weight
        else:
            new_es = running_es[-1] - 1.0 / len(non_hits)
        running_es.append(new_es)

    maximum, minimum = max(running_es), min(running_es)
    ES_score = maximum if maximum+minimum>0 else minimum
    
    return ES_score, hits, hits_idx, running_es

def permutated(rank):
    '''
    rank: pd.Series, index is gene name, value is rank
    '''
    # Sort genes by score in decreasing order
    rank = rank.copy().sort_values(ascending=False)
    # create a random permutation of the ranks
    values_perm = np.random.permutation(rank.values)
    # create a new rank series
    rank_perm = pd.Series(values_perm, index=rank.index)
    return rank_perm

def get_NES(rank, genelist, signalAsWeight, permutations = 1000, random_state = 0):
    '''
    rank: pd.Series, index is gene name, value is rank
    genelist: list, gene list
    signalAsWeight: bool, if True, use the rank as weight, otherwise use 1
    permutations: int, number of permutations
    '''
    # calculate ES
    np.random.seed(random_state)
    ES_score, hit, hits_idx, running_es = get_ES(rank, genelist, signalAsWeight)

    # calculate null distribution
    null_distribution = []
    for i in range(permutations):
        rank_perm = permutated(rank)
        ES_perm, hit_perm,_, running_es_perm = get_ES(rank_perm, genelist, signalAsWeight)
        null_distribution.append(ES_perm)

    # calculate NES
    NES_score = ES_score / np.mean(null_distribution)
    return NES_score, ES_score, hit, hits_idx, running_es, null_distribution


def get_pvalue(ES_score, null_distribution):
    '''
    ES_score: float, enrichment score
    null_distribution: list, null distribution of ES
    '''

    # calculate pvalue
    es_norm = (ES_score - np.mean(null_distribution)) / np.std(null_distribution)
    pvalue = (1 - norm.cdf(es_norm))
    
    return pvalue


def get_gsea(rank, genelist, signalAsWeight, permutations = 1000, method = 'fdr_bh'):
    '''
    rank: pd.Series, index is gene name, value is rank
    genelist: list, gene list
    signalAsWeight: bool, if True, use the rank as weight, otherwise use 1
    permutations: int, number of permutations
    method: str, method for FDR
    '''
    rank = rank.copy().sort_values(ascending=False)
    # calculate NES
    NES_score, ES_score, hit, hits_idx, running_es, null_distribution = get_NES(rank, genelist, signalAsWeight, permutations = permutations)

    # calculate pvalue
    pvalue = get_pvalue(ES_score, null_distribution)


    return NES_score, ES_score, hit, hits_idx, running_es, null_distribution, pvalue
