import pandas as pd
import numpy as np
import argparse, itertools
from collections import defaultdict
from scipy.stats import fisher_exact
from utils import read_data,gather_alleles
from alllele_frequencies import allele_freq



parser = argparse.ArgumentParser(description='Terminal arguments')
#parser.add_argument("--cyprus", action="store_true", dest="cyprus", default=False,
#					help="If you want to compare the greek banks with cyprus.")



if __name__ == "__main__":
    args = parser.parse_args()
    all_data = read_data('GreekCBUS_2921')
    data = [
            'ANATOLIKH_MAKEDONIA_142_100%and200_50%samples_271023', 'KENTRIKH_MAKEDONIA_866_100%and574_50%samples_271023',
            'DYTIKH_MAKEDONIA_81_100%and116_50%samples_271023', 'HPEIROS_12_100%and49_50%samples_271023', 'THESSALY_69_100%and155_50%samples_271023',
            'ST_ELLADA_5_100%and37_50%samples_271023', 'IONIA_1_100%and12_50%samples_271023', 'DYT_ELLADA_13_100%and79_50%samples_271023',
            'PELOPONNISOS_9_100%and65_50%samples_271023', 'ATTIKI_238_100%and374_50%samples_271023', 'B_AIGAIO_2_100%and17_50%samples_271023',
            'N_AIGAIO_1_100%and11_50%samples_271023', 'KRITI_300_100%and142_50%samples_271023',
            'NOT_FOUND_58_100%and27_50%samples_271023', 'ABROAD_43_100%and125_50%samples_271023', 'GREEK_IMM_9_100%and179_50%samples_271023',
           ]

    with pd.ExcelWriter(f'Allele_freqs_regions_banks.xlsx', engine = "openpyxl", mode = "w") as writer:
        for gene in ['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']:#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1'],['A', 'B','DRB1',]
            print(gene)
            all_region_alleles = gather_alleles(all_data, gene, region = True)
            df = pd.DataFrame(all_region_alleles)
            df.columns = [f'Allele_{gene}']
            for region in data:
                for target_sheet in ['natives', 'mixed_1', 'mixed_2']:
                    data_region, SheetName = read_data(region, sheet = target_sheet)
                    for bank in ['CRETE', 'BRFAA', 'PAP']:
                        data_bank = data_region[data_region['BANK'] == bank]

                        result,_ = allele_freq(
                                                data = data_bank,
                                                target_gene = gene,
                                                Grouped = True, #if the alleles are grouped: Grouped = True
                                                output_for_fisher = True,
                                              )

                        tmp = pd.DataFrame([
                                            #all_region_alleles,
                                            [result[key]['N'] if key in list(result.keys()) else 0 for key in all_region_alleles],
                                            [result[key]['AF'] if key in list(result.keys()) else 0 for key in all_region_alleles]
                                          ]).T
                        tmp.columns = [f'counts_{bank}_{SheetName}', f'AF_{bank}_{SheetName}']
                        df = pd.concat([df,tmp], axis = 1)

            #df = df.drop(df[df[f'Allele {gene}'] == 'not found'].index)
            df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)
