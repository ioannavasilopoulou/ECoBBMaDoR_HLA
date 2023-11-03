import pandas as pd
import numpy as np
import sys
from tqdm import tqdm
import argparse, itertools
from collections import defaultdict
from scipy.stats import fisher_exact
from utils import read_data, gather_alleles, create_sim_data
from alllele_frequencies import allele_freq

def sim_fisher(region_one, region_two, target_gene, all_alleles):
    results = defaultdict(lambda:defaultdict(lambda:[]))

    for round in range(1000):
        print(f'Number of samples before: {region_one[1].shape[0], region_two[1].shape[0]}')
        region_1, region_2 = create_sim_data(region_one, region_two)
        print(f'Number of samples after simulation: {region_1[1].shape[0], region_2[1].shape[0]}')
        region_1_info, total_region_1 = allele_freq(region_1[1], target_gene, output_for_fisher = True)
        region_2_info, total_region_2 = allele_freq(region_2[1], target_gene, output_for_fisher = True)
        #print(region_2_info, total_region_2)
        common_alleles = sorted(set(region_1_info.keys())&set(region_2_info.keys()))
        #print(f'all_alleles:{len(all_alleles)}')
        #print(f'common_alleles:{len(common_alleles)}')
        p_values = defaultdict(lambda:[])
        for allele in all_alleles:
            if allele in common_alleles:
                m = np.array(
                            [
                                [region_1_info[allele]['N'], total_region_1-region_1_info[allele]['N']],
                                [region_2_info[allele]['N'], total_region_2-region_2_info[allele]['N']]
                            ]
                            )
                oddsratio, pvalue = fisher_exact(m)
                #print(f'------>{pvalue}<------------')
                p_values[allele].append(pvalue)
                #print(p_values[allele])
            else:
                p_values[allele].append(0)
            #print(region_2_info[allele]['N'] if allele in region_2_info else 0, type(region_2_info[allele]['N'] if allele in region_2_info else 0))
            results[allele][f'N_{region_1[0]}'].append(region_1_info[allele]['N'] if allele in region_1_info else 0)
            results[allele][f'N_{region_2[0]}'].append(region_2_info[allele]['N'] if allele in region_2_info else 0)
            results[allele][f'AF_{region_1[0]}'].append(region_1_info[allele]['AF'] if allele in region_1_info else 0)
            results[allele][f'AF_{region_2[0]}'].append(region_2_info[allele]['AF'] if allele in region_2_info else 0)
            results[allele][f'Fisher_p_value_{region_1[0]}/{region_2[0]}'].append(p_values[allele][0])
    #print(results)
    #print(f'p_values:{len(p_values)}')
    average_results =defaultdict(lambda:[])
    #print(results[allele][f'N_{region_1[0]}'])

    average_results['Alleles'] = [key for key in results]
    average_results[f'N_{region_1[0]}'] = [sum(results[key][f'N_{region_1[0]}'])/1000 for key in results]
    average_results[f'N_{region_2[0]}'] = [sum(results[key][f'N_{region_2[0]}'])/1000 for key in results]
    average_results[f'AF_{region_1[0]}'] = [sum(results[key][f'AF_{region_1[0]}'])/1000 for key in results]
    average_results[f'AF_{region_2[0]}'] = [sum(results[key][f'AF_{region_2[0]}'])/1000 for key in results]
    #print(results[key][f'Fisher_p_value_{region_1[0]}/{region_2[0]}'], len(results[key][f'Fisher_p_value_{region_1[0]}/{region_2[0]}']))
    average_results[f'Fisher_p_value_{region_1[0]}/{region_2[0]}'] = [sum(results[key][f'Fisher_p_value_{region_1[0]}/{region_2[0]}'])/1000 for key in results]
    #print(results)
    #print(average_results[f'Fisher_p_value_{region_1[0]}/{region_2[0]}'], len(average_results[f'Fisher_p_value_{region_1[0]}/{region_2[0]}']))
    #sys.exit(-1)

    return average_results



parser = argparse.ArgumentParser(description='Terminal arguments')
parser.add_argument("--cyprus", action="store_true", dest="cyprus", default=False,
					help="If you want to compare the greek regions with cyprus.")


if __name__ == "__main__":
    args = parser.parse_args()
    if not args.cyprus:
        all_data = read_data('GreekCBUS_2921')
    else:
        all_greek_data = read_data('GreekCBUS_2921')
        all_cyprus_data = read_data('CYPRUS_cbus_4581_13_10_23')

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

    data =  data_regions_cyprus
    regions_list = []
    if args.cyprus:
        for region in data[:-1]:
            data_region, SheetName = read_data(region, sheet = 'natives')
            regions_list.append((region.split('and')[0], data_region))

        data_region, SheetName = read_data(data[-1], sheet = 'cyprus')
        regions_list.append((SheetName, data_region))
    else:
        for region in data:
            data_region, SheetName = read_data(region, sheet = 'natives')
            regions_list.append((region.split('and')[0], data_region))


    if not args.cyprus:
        iterlist = itertools.combinations(regions_list, 2)
    else:
        iterlist = [(info,regions_list[-1]) for info in regions_list[:-1]]

    #print(iterlist)
    for a,b in iterlist:
        print(a[0], b[0])
        #print(a,b)
        with pd.ExcelWriter(f'sim_fisher_per_regions_100%_{a[0]}-{b[0]}.xlsx', engine = "openpyxl", mode = "w") as writer:
            for gene in tqdm(['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']):#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1'],['A', 'B','DRB1',]
                print(gene)
                if not args.cyprus:
                    all_region_alleles = gather_alleles(all_data, gene, region = True)
                else:
                    all_region_alleles = gather_alleles((all_greek_data, all_cyprus_data), gene, cyprus = True, region = True)
                    print('All alleles DONE!!!')


                results = sim_fisher(a, b, gene, all_region_alleles)

                df = pd.DataFrame(results)


                df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)
