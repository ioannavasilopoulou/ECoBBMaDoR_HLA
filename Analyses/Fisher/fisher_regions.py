
import pandas as pd
import numpy as np
import argparse, itertools
from collections import defaultdict
from scipy.stats import fisher_exact
from utils import read_data,gather_alleles
from alllele_frequencies import allele_freq
from fisher_banks import fisher

parser = argparse.ArgumentParser(description='Terminal arguments')
parser.add_argument("--cyprus", action="store_true", dest="cyprus", default=False,
					help="If you want to compare the greek regions with cyprus.")
parser.add_argument("--GreekNatives", action="store_true", dest="GreekNatives", default=False,
					help="If you want to compare the greek, 1739 Greek CBUs, with cyprus.")

if __name__ == "__main__":
    args = parser.parse_args()

    if args.GreekNatives:
        data_all = [
                    'ANATOLIKH_MAKEDONIA_142_100%and200_50%samples_271023', 'KENTRIKH_MAKEDONIA_866_100%and574_50%samples_271023',
                    'DYTIKH_MAKEDONIA_81_100%and116_50%samples_271023', 'HPEIROS_12_100%and49_50%samples_271023', 'THESSALY_69_100%and155_50%samples_271023',
                    'ST_ELLADA_5_100%and37_50%samples_271023', 'IONIA_1_100%and12_50%samples_271023', 'DYT_ELLADA_13_100%and79_50%samples_271023',
                    'PELOPONNISOS_9_100%and65_50%samples_271023', 'ATTIKI_238_100%and374_50%samples_271023', 'B_AIGAIO_2_100%and17_50%samples_271023',
                    'N_AIGAIO_1_100%and11_50%samples_271023', 'KRITI_300_100%and142_50%samples_271023',
                    'CYPRUS_cbus_4581_13_10_23',
                   ]

        all_greek_data = read_data('GreekCBUS_2921')
        all_cyprus_data = read_data('CYPRUS_cbus_4581_13_10_23')
        regions_list = []
        for region in data_all[:-1]:
            data_region, SheetName = read_data(region, sheet = 'natives')
            regions_list.append((region.split('and')[0], data_region))

        data_region, SheetName = read_data(data_all[-1], sheet = 'cyprus')
        regions_list.append((SheetName, data_region))

        all_greeks = pd.DataFrame()
        for region in regions_list[:-1]:
            all_greeks = pd.concat([all_greeks, region[1]], axis = 0)

        all_greeks = ('100%_GREEKS_1739', all_greeks)

        with pd.ExcelWriter(f'fisher_per_1739_GREEKS_100%_and_cyprus.xlsx', engine = "openpyxl", mode = "w") as writer:
            for gene in ['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']:#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1'],['A', 'B','DRB1',]
                print(gene)
                all_region_alleles = gather_alleles((all_greek_data, all_cyprus_data), gene, cyprus = True, region = True)
                print('All alleles DONE!!!')
                results = fisher(banks = [all_greeks,regions_list[-1] ],
                                target_gene = gene,
                                all_alleles = all_region_alleles,
                                regions = True,
                                cyprus = True,
                                )



                #print(results)
                df = pd.DataFrame(results)
                df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)


    else:



        data_all = [
                    'ANATOLIKH_MAKEDONIA_142_100%and200_50%samples_271023', 'KENTRIKH_MAKEDONIA_866_100%and574_50%samples_271023',
                    'DYTIKH_MAKEDONIA_81_100%and116_50%samples_271023', 'HPEIROS_12_100%and49_50%samples_271023', 'THESSALY_69_100%and155_50%samples_271023',
                    'ST_ELLADA_5_100%and37_50%samples_271023', 'IONIA_1_100%and12_50%samples_271023', 'DYT_ELLADA_13_100%and79_50%samples_271023',
                    'PELOPONNISOS_9_100%and65_50%samples_271023', 'ATTIKI_238_100%and374_50%samples_271023', 'B_AIGAIO_2_100%and17_50%samples_271023',
                    'N_AIGAIO_1_100%and11_50%samples_271023', 'KRITI_300_100%and142_50%samples_271023',
                    'NOT_FOUND_58_100%and27_50%samples_271023', 'ABROAD_43_100%and125_50%samples_271023', 'GREEK_IMM_9_100%and179_50%samples_271023',
                   ]
        data_big_regions = [
                            'ANATOLIKH_MAKEDONIA_142_100%and200_50%samples_271023', 'KENTRIKH_MAKEDONIA_866_100%and574_50%samples_271023',
                            'DYTIKH_MAKEDONIA_81_100%and116_50%samples_271023', 'THESSALY_69_100%and155_50%samples_271023',
                            'ATTIKI_238_100%and374_50%samples_271023', 'KRITI_300_100%and142_50%samples_271023',
                            'NOT_FOUND_58_100%and27_50%samples_271023', 'ABROAD_43_100%and125_50%samples_271023',
                           ]

        data_regions_cyprus = [
                                'KENTRIKH_MAKEDONIA_866_100%and574_50%samples_271023',
                                'ATTIKI_238_100%and374_50%samples_271023', 'KRITI_300_100%and142_50%samples_271023',
                                'CYPRUS_cbus_4581_13_10_23'
                              ]

        data =  data_regions_cyprus#data_big_regions
        regions_list = []

        if not args.cyprus:
            all_data = read_data('GreekCBUS_2921')
            for region in data:
                data_region, SheetName = read_data(region, sheet = 'natives')
                regions_list.append((region.split('and')[0], data_region))
        else:
            all_greek_data = read_data('GreekCBUS_2921')
            all_cyprus_data = read_data('CYPRUS_cbus_4581_13_10_23')

            for region in data[:-1]:
                data_region, SheetName = read_data(region, sheet = 'natives')
                regions_list.append((region.split('and')[0], data_region))

            data_region, SheetName = read_data(data[-1], sheet = 'cyprus')
            regions_list.append((SheetName, data_region))

        '''
        TODO!!!!!!!!!!!!: CHANGE THE NAME OF THE FILE
        '''
        with pd.ExcelWriter(f'fisher_per_greek_regions_and_cyprus.xlsx', engine = "openpyxl", mode = "w") as writer:
            for gene in ['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']:#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1'],['A', 'B','DRB1',]
                print(gene)
                if not args.cyprus:
                    all_region_alleles = gather_alleles(all_data, gene, region = True)
                    results = fisher(banks = regions_list,
                                    target_gene = gene,
                                    all_alleles = all_region_alleles,
                                    regions = True,
                                    )
                else:
                    all_region_alleles = gather_alleles((all_greek_data, all_cyprus_data), gene, cyprus = True, region = True)
                    print('All alleles DONE!!!')
                    results = fisher(banks = regions_list,
                                    target_gene = gene,
                                    all_alleles = all_region_alleles,
                                    regions = True,
                                    cyprus = True,
                                    )



                #print(results)
                df = pd.DataFrame(results)
                #for key in results:
                #    df[key] = results[key]
                #print(df)

                df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)
