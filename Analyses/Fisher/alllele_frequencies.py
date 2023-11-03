import pandas as pd
import numpy as np
from collections import defaultdict

###############################################################
##      Calculates the frequency of each allele              ##
###############################################################
def allele_freq(data,
                target_gene,
                output_for_fisher,
                Grouped = True, #if the alleles are grouped: Grouped = True
               ):

    '''
    The Cretan  HLA profile should include all the HLA genes (A,B,C,DRB1,DQB1,DPB1).
    So we keep all the individuals who are fully genotyped on these 6 HLA genes.
    '''

    case = {'A':['A1_GROUPED', 'A2_GROUPED'], #indexes
            'B':['B1_GROUPED', 'B2_GROUPED'],
            'C':['C1_GROUPED', 'C2_GROUPED'],
            'DRB1':['DRB1_1_GROUPED', 'DRB1_2_GROUPED'],
            'DQB1':['DQB1_1_GROUPED', 'DQB1_2_GROUPED'],
            'DPB1':['DPB1_1_GROUPED', 'DPB1_2_GROUPED'],
            }

    alleles = [x for i in zip(data[case[target_gene][0]], data[case[target_gene][1]]) for x in i]
    if not output_for_fisher:
        counter_dict = defaultdict(lambda: 0)
        results = defaultdict(lambda: 0)
        for value in alleles:
            if value != 'not found':
                counter_dict[value] += 1

        freqs = {k:v/(sum(counter_dict.values())) for k,v in counter_dict.items()}
        results['Alleles'] = counter_dict.keys()
        results[f'N_{target_gene}'] = counter_dict.values()
        results[f'AF_{target_gene}'] = freqs.values()
        return results

    else:
        send_to_fisher = defaultdict(lambda: defaultdict(lambda: 0))
        for value in alleles:
            if value != 'not found':
                send_to_fisher[value]['N'] += 1
        sum_n = sum([send_to_fisher[value]['N'] for value in send_to_fisher.keys()])
        for value in alleles:
            if value != 'not found':
                try:
                    send_to_fisher[value]['AF'] = send_to_fisher[value]['N']/sum_n
                except:
                    send_to_fisher[value]['AF'] = 0


        return send_to_fisher,sum_n

if __name__ == "__main__":
    #BRFAA_samples781_210623_FINAL, Papanikolaou_samples1896_21062023_FINAL, Crete_428CBUs_260723, GreekCBUS_2921, CYPRUS_cbus_4581_13_10_23
    read_data = pd.read_excel(f'/Users/vasou/Documents/GitHub/HLA.github.io/Analyses/DATA/cbus_cyprus_4581_13_10_23.xlsx')
    read_data = read_data.replace(r'\?.*','not found', regex = True)
    read_data = read_data.replace(np.nan,'not found')
    print('Done: Read File!')

    #result = allele_freq(data = read_data,
    #            target_gene = 'A',
    #            output_for_fisher = True
    #            )
    #print(result)

    with pd.ExcelWriter('Allele_frequencies_CYPRUS.xlsx', engine = "openpyxl", mode = "w") as writer:
        for gene in ['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']:#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1'],['A', 'B','DRB1',]
            print(gene)
            result = allele_freq(data = read_data,
                            target_gene = gene,
                            Grouped = False, #if the alleles are grouped: Grouped = True
                            output_for_fisher = False,
                            )

            df = pd.DataFrame([result['Alleles'], result[f'N_{gene}'], result[f'AF_{gene}']]).T
            df.columns = [f'Allele_{gene}', f'counts_{gene}', f'AF_{gene}']
            #df = df.drop(df[df[f'Allele {gene}'] == 'not found'].index)
            df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)
