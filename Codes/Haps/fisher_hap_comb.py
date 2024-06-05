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
from fisher_hap import compare_haplotypes

parser = argparse.ArgumentParser(description='Terminal arguments')
parser.add_argument("--onetoall", action="store_true", dest="onetoall", default=False,
					help="Compare one population with others.")
args = parser.parse_args()

with pd.ExcelWriter(f"test.xlsx", engine = "openpyxl", mode = "w") as writer:
    for loci in ['3loci','5loci']:
        if loci == '3loci':
            #files = [
            #         'Greece/RunX_htf_Greek_CBUs_2fields_3012_3loci','Greece/RunX_htf_Greek_BMDs_2fields_75599_3loci',
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_2fields_2843_3loci', 'Cyprus/RunX_htf_cyprus_donors_2fields_96155_3loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields__98998_3loci'
            #        ]
            #files = [
            #         'Banks/RunX_htf_DHTOB_428_3loci','Banks/RunX_htf_BRFAA_592_3loci','Banks/RunX_htf_PAP_1992_3loci',
            #         'Greece/RunX_htf_Greek_CBUs_2fields_3012_3loci','Greece/RunX_htf_Greek_BMDs_2fields_75599_3loci',
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_2fields_2843_3loci', 'Cyprus/RunX_htf_cyprus_donors_2fields_96155_3loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields__98998_3loci'
            #        ]
            #files = [
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields__98998_3loci',
            #         'DKMS/RunX_htf_Austria_10000_3loci', 'DKMS/RunX_htf_Bosnia_10000_3loci','DKMS/RunX_htf_Croatia_10000_3loci',
            #         'DKMS/RunX_htf_Dutch_10000_3loci', 'DKMS/RunX_htf_France_10000_3loci','DKMS/RunX_htf_Greece_10000_3loci',
            #         'DKMS/RunX_htf_Italy_10000_3loci', 'DKMS/RunX_htf_Portugal_10000_3loci','DKMS/RunX_htf_Romania_10000_3loci',
            #         'DKMS/RunX_htf_Spain_10000_3loci', 'DKMS/RunX_htf_Turkey_10000_3loci',
            #        ]
            #files = [
            #         'Banks/RunX_htf_DHTOB_428_3loci','Banks/RunX_htf_BRFAA_592_3loci','Banks/RunX_htf_PAP_1992_3loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields__98998_3loci',
            #         'DKMS/RunX_htf_Austria_10000_3loci', 'DKMS/RunX_htf_Bosnia_10000_3loci','DKMS/RunX_htf_Croatia_10000_3loci',
            #         'DKMS/RunX_htf_Dutch_10000_3loci', 'DKMS/RunX_htf_France_10000_3loci','DKMS/RunX_htf_Greece_10000_3loci',
            #         'DKMS/RunX_htf_Italy_10000_3loci', 'DKMS/RunX_htf_Portugal_10000_3loci','DKMS/RunX_htf_Romania_10000_3loci',
            #         'DKMS/RunX_htf_Spain_10000_3loci', 'DKMS/RunX_htf_Turkey_10000_3loci',
            #        ]
            #files = [
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci',
            #         'RunX_htf_Group_1_156155_3loci','RunX_htf_Group_2_50000_3loci','RunX_htf_12_Countries_206155_3loci',
            #        ]
            files = [
                     'Greece/RunX_htf_Greek_CBUs_2fields_3012_3loci',
                     'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci',
                     'RunX_htf_Group_1_156155_3loci','RunX_htf_Group_2_50000_3loci','RunX_htf_12_Countries_206155_3loci',
                    ]
        else:
            #files = [
            #         'Greece/RunX_htf_Greek_CBUs_2fields_1092_5loci','Greece/RunX_htf_Greek_BMDs_2fields_69130_5loci',
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_2fields_2702_5loci', 'Cyprus/RunX_htf_cyprus_donors_2fields_73004_5loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields_75706_5loci'
            #        ]
            #files = [
            #         'Banks/RunX_htf_DHTOB_428_5loci','Banks/RunX_htf_BRFAA_564_5loci',
            #         'Greece/RunX_htf_Greek_CBUs_2fields_1092_5loci','Greece/RunX_htf_Greek_BMDs_2fields_69130_5loci',
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_2fields_2702_5loci', 'Cyprus/RunX_htf_cyprus_donors_2fields_73004_5loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields_75706_5loci'
            #        ]
            #files = [
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields_75706_5loci',
            #         'DKMS/RunX_htf_Austria_10000_5loci', 'DKMS/RunX_htf_Bosnia_10000_5loci','DKMS/RunX_htf_Croatia_10000_5loci',
            #         'DKMS/RunX_htf_Dutch_10000_5loci', 'DKMS/RunX_htf_France_10000_5loci','DKMS/RunX_htf_Greece_10000_5loci',
            #         'DKMS/RunX_htf_Italy_10000_5loci', 'DKMS/RunX_htf_Portugal_10000_5loci','DKMS/RunX_htf_Romania_10000_5loci',
            #         'DKMS/RunX_htf_Spain_10000_5loci', 'DKMS/RunX_htf_Turkey_10000_5loci',
            #        ]
            #files = [
            #         'Banks/RunX_htf_DHTOB_428_5loci','Banks/RunX_htf_BRFAA_564_5loci',
            #         'Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields_75706_5loci',
            #         'DKMS/RunX_htf_Austria_10000_5loci', 'DKMS/RunX_htf_Bosnia_10000_5loci','DKMS/RunX_htf_Croatia_10000_5loci',
            #         'DKMS/RunX_htf_Dutch_10000_5loci', 'DKMS/RunX_htf_France_10000_5loci','DKMS/RunX_htf_Greece_10000_5loci',
            #         'DKMS/RunX_htf_Italy_10000_5loci', 'DKMS/RunX_htf_Portugal_10000_5loci','DKMS/RunX_htf_Romania_10000_5loci',
            #         'DKMS/RunX_htf_Spain_10000_5loci', 'DKMS/RunX_htf_Turkey_10000_5loci',
            #        ]
            #files = [
            #         'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci',
            #         'RunX_htf_Group_1_133004_5loci','RunX_htf_Group_2_50000_5loci','RunX_htf_12_Countries_183004_5loci',
            #        ]
            files = [
                     'Greece/RunX_htf_Greek_CBUs_2fields_1092_5loci',
                     'Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci',
                     'RunX_htf_Group_1_133004_5loci','RunX_htf_Group_2_50000_5loci','RunX_htf_12_Countries_183004_5loci',
                    ]

        if not args.onetoall:
            pops = [haplomat_output_parce(file, int(file.split('_')[-2])) for file in files]

            df = pd.DataFrame()
            iterlist = combinations(files, 2)
            all_haplotypes = [list(pop.keys()) for pop in pops]
            all_haplotypes = sum(all_haplotypes, [])
            all_haplotypes = sorted(set(all_haplotypes))
            df['Haplotypes'] = all_haplotypes
            i=0
            for file_1, file_2 in list(iterlist):
                i+=1
                if i>=37 and loci == '3loci':
                    break
                elif i>=37 and loci == '5loci':
                    break
                print(file_1, file_2)
                results,_,_ = compare_haplotypes(
                                            pop_1 = file_1,        #haplomat output for one population
                                            pop_2 = file_2,      #haplomat output for another population
                                            num_1 = int(file_1.split('_')[-2]),     #number of samples in one population
                                            num_2 = int(file_2.split('_')[-2]), #number of samples in another population
                                            merged_haplotypes = all_haplotypes,
                                            no_diagram = None, #returns venn diagramm
                                            pvals = True,
                                            )

                df[f'counts_{file_1.split("htf_")[1]}'] = [results[haplotype][f'{file_1.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                df[f'counts_{file_2.split("htf_")[1]}'] = [results[haplotype][f'{file_2.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                df[f'HF_{file_1.split("htf_")[1]}'] = [results[haplotype][f'{file_1.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                df[f'HF_{file_2.split("htf_")[1]}'] = [results[haplotype][f'{file_2.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                df[f'Fisher_Pvalue_{file_1.split("htf_")[1]}-{file_2.split("htf_")[1]}'] = [results[haplotype][f'Fisher_{file_1.split("htf_")[1]}-{file_2.split("htf_")[1]}'] for haplotype in all_haplotypes]

                name = f'{loci}'

                df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)


        else:
            pops = [haplomat_output_parce(file, int(file.split('_')[-2])) for file in files]
            df = pd.DataFrame()
            all_haplotypes = [list(pop.keys()) for pop in pops]
            all_haplotypes = sum(all_haplotypes, [])
            all_haplotypes = sorted(set(all_haplotypes))
            df['Haplotypes'] = all_haplotypes
            for file in files[1:]:
                results,_,_ = compare_haplotypes(
                                            pop_1 = file,        #haplomat output for one population
                                            pop_2 = files[0],      #haplomat output for another population
                                            num_1 = int(file.split('_')[-2]),     #number of samples in one population
                                            num_2 = int(files[0].split('_')[-2]), #number of samples in another population
                                            merged_haplotypes = all_haplotypes,
                                            no_diagram = None, #returns venn diagramm
                                            pvals = True,
                                            )
                df[f'counts_{files[0].split("htf_")[1]}'] = [results[haplotype][f'{files[0].split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                df[f'HF_{files[0].split("htf_")[1]}'] = [results[haplotype][f'{files[0].split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                df[f'counts_{file.split("htf_")[1]}'] = [results[haplotype][f'{file.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                df[f'HF_{file.split("htf_")[1]}'] = [results[haplotype][f'{file.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                df[f'Fisher_Pvalue_{file.split("htf_")[1]}-{files[0].split("htf_")[1]}'] = [results[haplotype][f'Fisher_{file.split("htf_")[1]}-{files[0].split("htf_")[1]}'] for haplotype in all_haplotypes]

                name = f'{loci}'
                #name = f'3loci'

                df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)
