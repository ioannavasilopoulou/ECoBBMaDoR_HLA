import argparse, itertools
from collections import Counter, defaultdict
from functools import reduce
import pandas as pd
import numpy as np
import sys
from tqdm import tqdm
import scipy.stats as stats
###########################################################################
##                    Haplotype Frequencies                              ##
###########################################################################
def haplomat_output_parce(dat_file, #haplomat output
                          n, #number of samples
                         ):

    print(f'Read file: {dat_file}')
    path = '/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes'
    #path = '/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/Haplotype_Analyses/Haplotype_data/'
    with open(f'{path}/{dat_file}.dat') as f:
        haplotypes_initial = [x.replace('\n', '').replace('g', '').split() if 'g' in x else x.replace('\n', '').split() for x in f.readlines()]
    haplotypes_initial = [( '-'.join(sorted(x[0].split('~'))  ), float(x[1].replace(',', '.'))) for x in haplotypes_initial]
    haplotypes = {x[0]:x[1] for x in haplotypes_initial}
    haplotypes = {k:v for k,v in haplotypes.items() if v >= 1/(2*n)}

    #df = pd.DataFrame(list(zip(list(haplotypes.keys()),
    #                        [int(np.round(x*n*2)) for x in list(haplotypes.values())],
    #                        [round(x,4) for x in haplotypes.values()])),
    #                        columns = ['Haplotypes', 'N', 'HF'])


    #df2 = pd.DataFrame(
    #                {'Haplotypes': [f'Total: {len(haplotypes_initial)} - After: {len(haplotypes)}'],
    #                    'N': [''], 'HF': ['']}
    #                  )
    #display(df2)
    #df = pd.concat([df,df2], ignore_index = True)
    #df.to_excel(f'haplotypes_{file.split("htf_")[1]}.xlsx', index = False)
    return haplotypes
#haplomat_output_parce('RunX_htf_TOTAL_REGIONS_864_5loci', 864)
