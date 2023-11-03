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
    with open(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/Haplotype_Analyses/Haplotype_data/{dat_file}.dat') as f:
        haplotypes_initial = [x.replace('\n', '').replace('g', '').split() if 'g' in x else x.replace('\n', '').split() for x in f.readlines()]
    haplotypes_initial = [( '-'.join(sorted(x[0].split('~'))  ), float(x[1])) for x in haplotypes_initial]
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




'''
files = [
         f'RunX_htf_THESSALY_148_3loci',f'RunX_htf_ST_ELLADA_19_3loci', f'RunX_htf_PELOPONNISOS_34_3loci',
         f'RunX_htf_N_AIGAIO_9_3loci', f'RunX_htf_KRITI_370_3loci', f'RunX_htf_KEN_MAKED_1130_3loci',
         f'RunX_htf_IONIA_7_3loci',f'RunX_htf_HPEIROS_32_3loci', f'RunX_htf_DYT_MAKED_136_3loci',
         f'RunX_htf_DYT_ELLADA_57_3loci', f'RunX_htf_B_AIGAIO_6_3loci', f'RunX_htf_ATTIKI_409_3loci',
         f'RunX_htf_AN_MAKED_243_3loci','RunX_htf_3loci_2843_Cyprus_CBUs_g_resolution_161023',
         'RunX_htf_Greece_10000_181022_3loci',
         'RunX_htf_THESSALYvsTOTAL_2452_3loci', 'RunX_htf_ST_ELLADAvsTOTAL_2581_3loci', 'RunX_htf_PELOPONNvsTOTAL_2566_3loci',
         'RunX_htf_N_AIGAIOvsTOTAL_2591_3loci','RunX_htf_KRITIvsTOTAL_2230_3loci', 'RunX_htf_KEN_MAKEDvsTOTAL_1470_3loci',
         'RunX_htf_IONIAvsTOTAL_2593_3loci', 'RunX_htf_HPEIROSvsTOTAL_2568_3loci', 'RunX_htf_DYT_MAKEDvsTOTAL_2464_3loci',
         'RunX_htf_DYT_ELLADAvsTOTAL_2543_3loci', 'RunX_htf_B_AIGAIOvsTOTAL_2594_3loci','RunX_htf_ATTIKIvsTOTAL_2191_3loci',
         'RunX_htf_AN_MAKEDvsTOTAL_2357_3loci',
         'RunX_htf_TOTAL_REGIONS_2600_3loci'
        ]

files = [
         f'RunX_htf_THESSALY_20_5loci',f'RunX_htf_ST_ELLADA_14_5loci', f'RunX_htf_PELOPONNISOS_28_5loci',
         f'RunX_htf_N_AIGAIO_4_5loci', f'RunX_htf_KRITI_364_5loci', f'RunX_htf_KEN_MAKED_42_5loci',
         f'RunX_htf_IONIA_1_5loci',f'RunX_htf_HPEIROS_14_5loci', f'RunX_htf_DYT_MAKED_3_5loci',
         f'RunX_htf_DYT_ELLADA_34_5loci', f'RunX_htf_B_AIGAIO_4_5loci', f'RunX_htf_ATTIKI_325_5loci',
         f'RunX_htf_AN_MAKED_11_5loci','RunX_htf_5loci_2702_Cyprus_CBUs_g_reolution_161023',
         'RunX_htf_Greece_10000_181022_5loci',
         'RunX_htf_THESSALYvsTOTAL_844_5loci', 'RunX_htf_ST_ELLADAvsTOTAL_850_5loci', 'RunX_htf_PELOPONNvsTOTAL_836_5loci',
        'RunX_htf_N_AIGAIOvsTOTAL_860_5loci','RunX_htf_KRITIvsTOTAL_500_5loci', 'RunX_htf_KEN_MAKEDvsTOTAL_822_5loci',
        'RunX_htf_IONIAvsTOTAL_863_5loci', 'RunX_htf_HPEIROSvsTOTAL_850_5loci', 'RunX_htf_DYT_MAKEDvsTOTAL_861_5loci',
        'RunX_htf_DYT_ELLADAvsTOTAL_830_5loci', 'RunX_htf_B_AIGAIOvsTOTAL_860_5loci','RunX_htf_ATTIKIvsTOTAL_539_5loci',
        'RunX_htf_AN_MAKEDvsTOTAL_853_5loci',
        'RunX_htf_TOTAL_REGIONS_864_5loci'
        ]

for file in files:
    if 'Cyprus' in file:
        haplomat_output_parce(file, 2843)
    elif 'Greece' in file:
        haplomat_output_parce(file, 10000)
    else:
        haplomat_output_parce(file, int(file.split('_')[-2]))
'''
#Total: 1813 - After: 1671
#Total: 956 - After: 934
#

#result = haplomat_output_parce('/Users/vasou/Downloads/RunX_htf_kentrik_maked_58_5loci.dat', 58)
#result.to_excel('hap_freqs_kentrik_maked_58_5loci.xlsx', index = False)
