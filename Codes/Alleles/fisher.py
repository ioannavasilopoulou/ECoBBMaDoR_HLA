
import pandas as pd
import numpy as np
import argparse, itertools
from collections import defaultdict
from scipy.stats import fisher_exact
from utils import read_data, gather_alleles, create_sim_data
from alllele_frequencies import allele_freq
from Alleles_from_Haps_manos import af_estimation_from_HaploMat_g_group


def fisher_anal( banks,
            target_gene,
            all_alleles = None,
            comb_iterlist = None,
            regions = False,
            loci = None,
          ):

    results = defaultdict(lambda:[])
    if comb_iterlist:
        iterlist = itertools.combinations(banks, 2)

    else:
        iterlist = [(info,banks[0]) for info in banks[1:]]

    if not all_alleles:
        all_alleles = gather_alleles(banks, target_gene, loci = loci)#TODO: CHECK FOR cyprus)
        i=0
    for a,b in iterlist:
        i+=1
        if i>=40:#Stop if you have done 40 pairs
            break
        print(a[0],b[0])
        #in case where we want to compare with DKMS alleles (alleles from haplotypes)
        if 'dat' in a[0] and 'dat' in b[0]:
            bank_1_info = af_estimation_from_HaploMat_g_group(a[1], target_gene, int(a[1].split('_')[-2]))
            total_bank_1 = sum([bank_1_info[key]['N'] for key in bank_1_info])
            bank_2_info = af_estimation_from_HaploMat_g_group(b[1], target_gene, int(b[1].split('_')[-2]))
            total_bank_2 = sum([bank_2_info[key]['N'] for key in bank_2_info])
        elif 'dat' in a[0]:
            bank_1_info = af_estimation_from_HaploMat_g_group(a[1], target_gene, int(a[1].split('_')[-2]))
            total_bank_1 = sum([bank_2_info[key]['N'] for key in bank_2_info])
            bank_2_info, total_bank_2 = allele_freq(b[1], target_gene, output_for_fisher = True)
        elif 'dat' in b[0]:
            bank_1_info, total_bank_1 = allele_freq(a[1], target_gene, output_for_fisher = True)
            bank_2_info = af_estimation_from_HaploMat_g_group(b[1], target_gene, int(b[1].split('_')[-2]))
            total_bank_2 = sum([bank_2_info[key]['N'] for key in bank_2_info])
            #print(bank_1_info, bank_2_info)
        else:

            bank_1_info, total_bank_1 = allele_freq(a[1], target_gene, output_for_fisher = True)
            bank_2_info, total_bank_2 = allele_freq(b[1], target_gene, output_for_fisher = True)

        #print(bank_1_info)
        #print(a[0], b[0])
        #all_alleles = sorted(set(all_alleles)|(set(bank_1_info.keys())|set(bank_2_info.keys())))
        common_alleles = sorted(set(bank_1_info.keys())&set(bank_2_info.keys()))
        print(f'all_alleles:{len(all_alleles)}')
        print(f'common_alleles:{len(common_alleles)}')
        p_values = []
        for allele in all_alleles:
            if allele in common_alleles:

                m = np.array(
                            [
                                [bank_1_info[allele]['N'], total_bank_1-bank_1_info[allele]['N']],
                                [bank_2_info[allele]['N'], total_bank_2-bank_2_info[allele]['N']]
                            ]
                            )

                oddsratio, pvalue = fisher_exact(m)
                #if allele == '01:01':
                #    print(m)
                #    print(oddsratio, pvalue)
                #print(f'------>{pvalue}<------------')
                p_values.append(pvalue)
            else:
                p_values.append(0)
        print(f'p_values:{len(p_values)}')
        if not regions:
            results['Alleles'] = all_alleles
            results[f'N_{a[0]}'] = [bank_1_info[allele]['N'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'N_{b[0]}'] = [bank_2_info[allele]['N'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'AF_{a[0]}'] = [bank_1_info[allele]['AF'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'AF_{b[0]}'] = [bank_2_info[allele]['AF'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'Fisher_p_value_{a[0]}/{b[0]}'] = p_values

        else:
            results['Alleles'] = all_alleles
            results[f'N_{a[0]}'] = [bank_1_info[allele]['N'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'N_{b[0]}'] = [bank_2_info[allele]['N'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'AF_{a[0]}'] = [bank_1_info[allele]['AF'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'AF_{b[0]}'] = [bank_2_info[allele]['AF'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'Fisher_p_value_{a[0]}/{b[0]}'] = p_values

    return results
