import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from alllele_frequencies import allele_freq
from Alleles_from_Haps_manos import af_estimation_from_HaploMat_g_group

def read_data(
                file,
                sheet = None
             ):
    path = '/Users/vasou/Documents/HLA/Matching Coverage/Data/'
    if sheet == 'natives':
        data = pd.read_excel(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/DATA_regions/{file}.xlsx', sheet_name = None)
        #print(list(data.keys())[1])
        name = list(data.keys())[1]
        native_data = data[name]
        native_data = native_data.replace(r'\?.*','not found', regex = True)
        native_data = native_data.replace(np.nan,'not found')
        print(f'Done: Read File {list(data.keys())[1]}!,{name}')
        return native_data, name
    elif sheet == 'mixed_1':
        data = pd.read_excel(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/DATA_regions/{file}.xlsx', sheet_name = None)
        #print(list(data.keys())[1])
        name = list(data.keys())[2]
        mixed_1_data = data[name]
        mixed_1_data = mixed_1_data.replace(r'\?.*','not found', regex = True)
        mixed_1_data = mixed_1_data.replace(np.nan,'not found')
        print(f'Done: Read File {name}!')
        return mixed_1_data, name
    elif sheet == 'mixed_2':
        data = pd.read_excel(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/DATA_regions/{file}.xlsx', sheet_name = None)
        #print(list(data.keys())[1])
        name = list(data.keys())[3]
        mixed_2_data = data[name]
        mixed_2_data = mixed_2_data.replace(r'\?.*','not found', regex = True)
        mixed_2_data = mixed_2_data.replace(np.nan,'not found')
        print(f'Done: Read File {name}!')
        return mixed_2_data, name

    elif sheet == 'region':
        data = pd.read_excel(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/DATA_regions/{file}.xlsx', sheet_name = None)
        name = list(data.keys())[0]
        all_region = data[name]
        all_region = all_region.replace(r'\?.*','not found', regex = True)
        all_region = all_region.replace(np.nan,'not found')
        print(f'Done: Read File {list(data.keys())[0]}!')
        return all_region, name
    elif sheet == 'cyprus':
        data = pd.read_excel(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/DATA/{file}.xlsx')
        #data = pd.read_excel(f'~/HLA/data/{file}.xlsx')
        data = data.replace(r'\?.*','not found', regex = True)
        data = data.replace(np.nan,'not found')
        print('Done: Read File!')
        return data, '_'.join(file.split('_')[:-3])#'100%_Cyprus_4581'

    elif sheet:
        data = pd.read_excel(f'{path}/{file}.xlsx', sheet_name = sheet)#sheet_name = sheet for banks
        #data = pd.read_excel(f'~/HLA/data/{file}.xlsx')
        data = data.replace(r'\?.*','not found', regex = True)
        data = data.replace(np.nan,'not found')
        data = data[~data.apply(lambda row: row.astype(str).str.contains('not found')).any(axis=1)]
        print(f'Done: File {file} Readed!')
        print(f'The shape of data: {data.shape}')
        return data

    else:
        data = pd.read_excel(f'{path}/{file}.xlsx')
        #data = pd.read_excel(f'~/HLA/data/{file}.xlsx')
        data = data.replace(r'\?.*','not found', regex = True)
        data = data.replace(np.nan,'not found')
        data = data[~data.apply(lambda row: row.astype(str).str.contains('not found')).any(axis=1)]
        print(f'Done: File {file} Readed!')
        print(f'The shape of data: {data.shape}')
        return data

def gather_alleles(banks, target_gene, loci=None, cyprus = None, region = None,):
    if cyprus and not region:
        bank_4_info, _ = allele_freq(banks[3][1], target_gene, output_for_fisher = True)
        bank_5_info, _ = allele_freq(banks[4][1], target_gene, output_for_fisher = True)
        all_alleles = sorted(set(bank_4_info.keys())|set(bank_5_info.keys()))
        return all_alleles
    elif cyprus and region:
        bank_1_info, _ = allele_freq(banks[0], target_gene, output_for_fisher = True)
        bank_2_info, _ = allele_freq(banks[1], target_gene, output_for_fisher = True)
        all_alleles = sorted(set(bank_1_info.keys())|set(bank_2_info.keys()))
        return all_alleles

    elif region:
        all_region, _ = allele_freq(banks, target_gene, output_for_fisher = True)
        all_alleles = sorted(set(all_region.keys()))
        return all_alleles
    else:
        all_alleles = set()
        for bank in banks:
            if 'DKMS' in bank[1] or 'dat' in bank[1]:
                bank_info = af_estimation_from_HaploMat_g_group(bank[1], target_gene, int(bank[1].split('_')[-2]))
            else:
                bank_info, _ = allele_freq(bank[1], target_gene, output_for_fisher = True, loci = loci)
            all_alleles = all_alleles|set(bank_info.keys())
        all_alleles = sorted(all_alleles)

        return all_alleles


def create_sim_data(
                    data_1,
                    data_2,
                   ):
    if data_1[1].shape[0]<data_2[1].shape[0]:

        sim_data = data_2[1].sample(frac = (data_1[1].shape[0]/data_2[1].shape[0]),
                                    replace = True
                                   )


        return data_1, (data_2[0], sim_data)

    elif data_1[1].shape[0]>data_2[1].shape[0]:
        sim_data = data_1[1].sample(frac = (data_2[1].shape[0]/data_1[1].shape[0]),
                                    replace = True
                                   )


        return (data_1[0], sim_data), data_2
