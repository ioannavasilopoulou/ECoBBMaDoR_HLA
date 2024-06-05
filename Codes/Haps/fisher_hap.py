import argparse, itertools
from collections import Counter, defaultdict
from parse_haplotypes import haplomat_output_parce
from functools import reduce
import pandas as pd
import numpy as np
import sys
from tqdm import tqdm
import scipy.stats as stats
from itertools import combinations



def compare_haplotypes(pop_1,         #haplomat output for one population
                       pop_2,      #haplomat output for another population
                       num_1,     #number of samples in one population
                       num_2,      #number of samples in another population
                       merged_haplotypes,
                       no_diagram=None, #returns venn diagramm
                       pvals=None       #returns venn p-values
                      ):

    '''
    The statistic comparison bewtween the haplotypes will be estimated with Fisher's exact test,
    due to the small number of counts for the haplotypes. Otherwise, chi_square would be used,
    like in the allele frequency comparison.
    '''
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
    pop_1_haplotypes = haplomat_output_parce(pop_1, num_1)
    pop_2_haplotypes = haplomat_output_parce(pop_2, num_2)
    #union_haplotypes = list(set(pop_1_haplotypes.keys()) | set(pop_2_haplotypes.keys()))
    print('Files parced!')
    common_haplotypes = list(set(pop_1_haplotypes.keys()) & set(pop_2_haplotypes.keys()))
    #print(common_haplotypes)
    print(f'We have {len(common_haplotypes)} common haplotypes.')
    #comparison_list_counts = [[x,0,0,0,0] for x in common_haplotypes]
    comparison_list_counts = [[x,0,0,0,0] for x in merged_haplotypes]
    #print(comparison_list_counts)
    for x in comparison_list_counts:


        x[1] += round(pop_1_haplotypes[x[0]]*2*num_1) if x[0] in pop_1_haplotypes.keys() else 0
        x[2] += round(pop_2_haplotypes[x[0]]*2*num_2) if x[0] in pop_2_haplotypes.keys() else 0
        x[3] += pop_1_haplotypes[x[0]] if x[0] in pop_1_haplotypes.keys() else 0
        x[4] += pop_2_haplotypes[x[0]] if x[0] in pop_2_haplotypes.keys() else 0



    #print(comparison_list_counts)
    #comparison_list_counts = sorted(comparison_list_counts, key=lambda x:x[1], reverse=True)
    fisher_inputs = [[[2*num_1 - x[1], x[1]], [2*num_2 - x[2], x[2]]] for x in comparison_list_counts]
    p_values = defaultdict(lambda:defaultdict(lambda:0))
    for x, hap in zip(fisher_inputs, comparison_list_counts):
        if hap[0] in common_haplotypes:
            oddsratio, pvalue = stats.fisher_exact(x)
            p_values[hap[0]][f'{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}'] = pvalue
        else:
            p_values[hap[0]][f'{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}'] = 0


    #hapl_pval = dict(zip([x[0] for x in comparison_list_counts], p_values))
    output = defaultdict(lambda:defaultdict(lambda:[]))
    #print(comparison_list_counts)

    for x in comparison_list_counts:
        output[x[0]][pop_1.split('htf_')[1]] = [x[1], x[3]]
        output[x[0]][pop_2.split('htf_')[1]] = [x[2], x[4]]
        output[x[0]][f'Fisher_{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}'] = p_values[x[0]][f'{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}']

    #print(output)
    #print(len(output))
    #sys.exit(-1)
    #hapl_pval = pd.DataFrame(
    #                        [[x[0] for x in comparison_list_counts],
    #                         [x[1] for x in comparison_list_counts],
    #                         [x[2] for x in comparison_list_counts],
    #                         [x[3] for x in comparison_list_counts],
    #                         [x[4] for x in comparison_list_counts],
    #                         p_values]
    #                        )
    if pvals:
        return output,pop_1_haplotypes,pop_2_haplotypes
        #return hapl_pval
