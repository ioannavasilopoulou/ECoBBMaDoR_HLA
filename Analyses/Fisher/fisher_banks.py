
import pandas as pd
import numpy as np
import argparse, itertools
from collections import defaultdict
from scipy.stats import fisher_exact
from utils import read_data, gather_alleles, create_sim_data
from alllele_frequencies import allele_freq


def fisher( banks,
            target_gene,
            all_alleles = None,
            cyprus = None,
            regions = False,
          ):

    results = defaultdict(lambda:[])
    if not cyprus:
        iterlist = itertools.combinations(banks, 2)
    else:
        iterlist = [(info,banks[-1]) for info in banks[:-1]]

    if not all_alleles:
        all_alleles = gather_alleles(banks, target_gene, cyprus)
    for a,b in iterlist:

        bank_1_info, total_bank_1 = allele_freq(a[1], target_gene, output_for_fisher = True)
        bank_2_info, total_bank_2 = allele_freq(b[1], target_gene, output_for_fisher = True)

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
                #print(f'------>{pvalue}<------------')
                p_values.append(pvalue)
            else:
                p_values.append(0)
        print(f'p_values:{len(p_values)}')
        if not regions:
            results['Alleles'] = all_alleles
            results[f'N_{a[0].split("_")[0]}'] = [bank_1_info[allele]['N'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'N_{b[0].split("_")[0]}'] = [bank_2_info[allele]['N'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'AF_{a[0].split("_")[0]}'] = [bank_1_info[allele]['AF'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'AF_{b[0].split("_")[0]}'] = [bank_2_info[allele]['AF'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'Fisher_p_value_{a[0].split("_")[0]}/{b[0].split("_")[0]}'] = p_values

        else:
            results['Alleles'] = all_alleles
            results[f'N_{a[0]}'] = [bank_1_info[allele]['N'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'N_{b[0]}'] = [bank_2_info[allele]['N'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'AF_{a[0]}'] = [bank_1_info[allele]['AF'] if allele in bank_1_info else 0 for allele in all_alleles]
            results[f'AF_{b[0]}'] = [bank_2_info[allele]['AF'] if allele in bank_2_info else 0 for allele in all_alleles]
            results[f'Fisher_p_value_{a[0]}/{b[0]}'] = p_values

    return results

parser = argparse.ArgumentParser(description='Terminal arguments')
parser.add_argument("--cyprus", action="store_true", dest="cyprus", default=False,
					help="If you want to compare the greek banks with cyprus.")


if __name__ == "__main__":
    args = parser.parse_args()
    with pd.ExcelWriter('Fisher_between_greek_banks_and_cyprus.xlsx', engine = "openpyxl", mode = "w") as writer:
        data_bank_1 = read_data('BRFAA_samples781_210623_FINAL')
        data_bank_2 = read_data('Papanikolaou_samples1896_21062023_FINAL')
        data_bank_3 = read_data('Crete_428CBUs_260723')
        data_bank_4 = read_data('GreekCBUS_2921')
        data_bank_5 = read_data('CYPRUS_cbus_4581_13_10_23')

        banks = [
                ('BRFAA_samples781_210623_FINAL',data_bank_1), ('Papanikolaou_samples1896_21062023_FINAL',data_bank_2), ('Crete_428CBUs_260723',data_bank_3),
                ('GreekCBUS_2921',data_bank_4), ('CYPRUS_cbus_4581_13_10_23',data_bank_5)
                ]

        for gene in ['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']:#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1'],['A', 'B','DRB1',]
            print(gene)
            result = fisher(banks,
                            gene,
                            cyprus = args.cyprus
                            )
            #print(result)
            if not args.cyprus:
                df = pd.DataFrame([result['Alleles'], result[f'N_{banks[0][0].split("_")[0]}'], result[f'AF_{banks[0][0].split("_")[0]}'],
                result[f'N_{banks[1][0].split("_")[0]}'], result[f'AF_{banks[1][0].split("_")[0]}'],result[f'Fisher_p_value_{banks[0][0].split("_")[0]}/{banks[1][0].split("_")[0]}'],
                result[f'N_{banks[2][0].split("_")[0]}'], result[f'AF_{banks[2][0].split("_")[0]}'],result[f'Fisher_p_value_{banks[0][0].split("_")[0]}/{banks[2][0].split("_")[0]}'],
                result[f'Fisher_p_value_{banks[1][0].split("_")[0]}/{banks[2][0].split("_")[0]}']]).T

                df.columns = [f'Allele_{gene}', f'counts_{banks[0][0].split("_")[0]}', f'AF_{banks[0][0].split("_")[0]}',
                            f'N_{banks[1][0].split("_")[0]}', f'AF_{banks[1][0].split("_")[0]}',f'Fisher_p_value_{banks[0][0].split("_")[0]}/{banks[1][0].split("_")[0]}',
                            f'N_{banks[2][0].split("_")[0]}', f'AF_{banks[2][0].split("_")[0]}',f'Fisher_p_value_{banks[0][0].split("_")[0]}/{banks[2][0].split("_")[0]}',
                            f'Fisher_p_value_{banks[1][0].split("_")[0]}/{banks[2][0].split("_")[0]}'
                            ]
            else:
                df = pd.DataFrame(
                                [result['Alleles'],
                                result[f'N_{banks[-1][0].split("_")[0]}'], result[f'AF_{banks[-1][0].split("_")[0]}'],
                                result[f'N_{banks[0][0].split("_")[0]}'], result[f'AF_{banks[0][0].split("_")[0]}'],
                                result[f'Fisher_p_value_{banks[0][0].split("_")[0]}/{banks[-1][0].split("_")[0]}'],
                                result[f'N_{banks[1][0].split("_")[0]}'], result[f'AF_{banks[1][0].split("_")[0]}'],result[f'Fisher_p_value_{banks[1][0].split("_")[0]}/{banks[-1][0].split("_")[0]}'],
                                result[f'N_{banks[2][0].split("_")[0]}'], result[f'AF_{banks[2][0].split("_")[0]}'],result[f'Fisher_p_value_{banks[2][0].split("_")[0]}/{banks[-1][0].split("_")[0]}'],
                                result[f'N_{banks[3][0].split("_")[0]}'], result[f'AF_{banks[3][0].split("_")[0]}'],result[f'Fisher_p_value_{banks[3][0].split("_")[0]}/{banks[-1][0].split("_")[0]}'],
                                ]).T

                df.columns = [f'Allele_{gene}',
                            f'counts_{banks[-1][0].split("_")[0]}', f'AF_{banks[-1][0].split("_")[0]}',
                            f'N_{banks[0][0].split("_")[0]}', f'AF_{banks[0][0].split("_")[0]}',f'Fisher_p_value_{banks[0][0].split("_")[0]}/{banks[-1][0].split("_")[0]}',
                            f'N_{banks[1][0].split("_")[0]}', f'AF_{banks[1][0].split("_")[0]}',f'Fisher_p_value_{banks[1][0].split("_")[0]}/{banks[-1][0].split("_")[0]}',
                            f'N_{banks[2][0].split("_")[0]}', f'AF_{banks[2][0].split("_")[0]}',f'Fisher_p_value_{banks[2][0].split("_")[0]}/{banks[-1][0].split("_")[0]}',
                            f'N_{banks[3][0].split("_")[0]}', f'AF_{banks[3][0].split("_")[0]}',f'Fisher_p_value_{banks[3][0].split("_")[0]}/{banks[-1][0].split("_")[0]}',
                            ]

            #df = df.drop(df[df[f'Allele {gene}'] == 'not found'].index)
            df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)
