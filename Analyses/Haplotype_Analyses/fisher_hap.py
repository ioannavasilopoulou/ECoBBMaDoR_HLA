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



def compare_haplotypes(pop_1,         #haplomat output for one population
                       pop_2,      #haplomat output for another population
                       num_1,     #number of samples in one population
                       num_2,      #number of samples in another population
                       merged_haplotypes,
                       no_diagram=None, #returns venn diagramm
                       pvals=None       #returns venn p-values
                      ):

    '''
    The statistic comparison bewtween the haplotypes will be estimated with Fisher's exact test,
    due to the small number of counts for the haplotypes. Otherwise, chi_square would be used,
    like in the allele frequency comparison.
    '''
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
    pop_1_haplotypes = haplomat_output_parce(pop_1, num_1)
    pop_2_haplotypes = haplomat_output_parce(pop_2, num_2)
    #union_haplotypes = list(set(pop_1_haplotypes.keys()) | set(pop_2_haplotypes.keys()))
    print('Files parced!')
    common_haplotypes = list(set(pop_1_haplotypes.keys()) & set(pop_2_haplotypes.keys()))
    #print(common_haplotypes)
    print(f'We have {len(common_haplotypes)} common haplotypes.')
    #comparison_list_counts = [[x,0,0,0,0] for x in common_haplotypes]
    comparison_list_counts = [[x,0,0,0,0] for x in merged_haplotypes]
    #print(comparison_list_counts)
    for x in comparison_list_counts:


        x[1] += round(pop_1_haplotypes[x[0]]*2*num_1) if x[0] in pop_1_haplotypes.keys() else 0
        x[2] += round(pop_2_haplotypes[x[0]]*2*num_2) if x[0] in pop_2_haplotypes.keys() else 0
        x[3] += pop_1_haplotypes[x[0]] if x[0] in pop_1_haplotypes.keys() else 0
        x[4] += pop_2_haplotypes[x[0]] if x[0] in pop_2_haplotypes.keys() else 0



    #print(comparison_list_counts)
    #comparison_list_counts = sorted(comparison_list_counts, key=lambda x:x[1], reverse=True)
    fisher_inputs = [[[2*num_1 - x[1], x[1]], [2*num_2 - x[2], x[2]]] for x in comparison_list_counts]
    p_values = defaultdict(lambda:defaultdict(lambda:0))
    for x, hap in zip(fisher_inputs, comparison_list_counts):
        if hap[0] in common_haplotypes:
            oddsratio, pvalue = stats.fisher_exact(x)
            p_values[hap[0]][f'{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}'] = pvalue
        else:
            p_values[hap[0]][f'{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}'] = 0


    #hapl_pval = dict(zip([x[0] for x in comparison_list_counts], p_values))
    output = defaultdict(lambda:defaultdict(lambda:[]))
    #print(comparison_list_counts)

    for x in comparison_list_counts:
        output[x[0]][pop_1.split('htf_')[1]] = [x[1], x[3]]
        output[x[0]][pop_2.split('htf_')[1]] = [x[2], x[4]]
        output[x[0]][f'Fisher_{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}'] = p_values[x[0]][f'{pop_1.split("htf_")[1]}-{pop_2.split("htf_")[1]}']

    #print(output)
    #print(len(output))
    #sys.exit(-1)
    #hapl_pval = pd.DataFrame(
    #                        [[x[0] for x in comparison_list_counts],
    #                         [x[1] for x in comparison_list_counts],
    #                         [x[2] for x in comparison_list_counts],
    #                         [x[3] for x in comparison_list_counts],
    #                         [x[4] for x in comparison_list_counts],
    #                         p_values]
    #                        )
    if pvals:
        return output,pop_1_haplotypes,pop_2_haplotypes
        #return hapl_pval

    sig_hapl_pval = {k:v for k,v in hapl_pval.items() if v < 0.05}

    if no_diagram:
        print(f'Total Cretan haplotypes: {len(pop_1_haplotypes)}')
        print(f"Total others haplotypes: {len(pop_2_haplotypes)}")
        print(f'Common haplotypes: {len(common_haplotypes)}')
        print(f'Statistcal Significant: {len(sig_hapl_pval)}\n')
        return

    white_patch = mpatches.Patch(color='none', label=f'{len(sig_hapl_pval)} / {len(common_haplotypes)}\nstatistically\nsignificant\ndifference')
    leg = plt.legend(handles=[white_patch], loc='upper right', frameon=False)
    leg.set_alpha(0.0)

    venn2(subsets = (len(others_haplotypes)-len(common_haplotypes), len(cretan_haplotypes)-len(common_haplotypes),
                     len(common_haplotypes)), set_labels = (f'OTHERS\n({num_others})', f'Crete\n({num_Cretans})'),
                     set_colors=('skyblue', 'red'), alpha = 0.7)
#     venn2(subsets = (len(DKMS_haplotypes)-len(common_haplotypes), len(cretan_haplotypes)-len(common_haplotypes), len(common_haplotypes)), set_labels = (f'DKMS\n({num_DKMS})', f'Crete\n({num_Cretans})'), set_colors=('darkgrey', 'slategrey'), alpha = 0.7)
    plt.title('Haplotype comparison between\nCretans and others', fontsize=15)

#     plt.savefig(f'Vehn with {len(common_haplotypes)} common haplotypes.png', dpi=300)

    plt.show()



parser = argparse.ArgumentParser(description='Terminal arguments')
parser.add_argument("--cyprus", action="store_true", dest="cyprus", default=False,
					help="Compare greek regions with cyprus.")

#parser.add_argument("--dkms", action="store_true", dest="dkms", default=False,
#					help="Compare greek regions with DKMS.")

#parser.add_argument("--minus", action="store_true", dest="minus", default=False,
#					help="Compare greek regions with total greek regions minus.")

parser.add_argument("--region_pairs", action="store_true", dest="region_pairs", default=False,
					help="Between regions")

parser.add_argument("--loci5", action="store_true", dest="loci5", default=False,
					help="5loci Analyses")
parser.add_argument("--loci3", action="store_true", dest="loci3", default=False,
					help="3loci Analyses")

#parser.add_argument("--num_1", type = int, dest="num_1", required = True,
#					help="The number of samples in the first population.")
#parser.add_argument("--cyprus", type = str, dest="cyprus", required = True,
#					help="Compare greek regions with cyprus")


if __name__ == "__main__":
    args = parser.parse_args()



    if args.region_pairs and args.loci3:


        files = [
                 f'RunX_htf_AM_142_100%_141_3loci',f'RunX_htf_KM_866_100%_842_3loci', f'RunX_htf_DM_81_100%_78_3loci',
                 f'RunX_htf_THESSALY_69_100%_67_3loci', f'RunX_htf_ATTIKI_238_100%_238_3loci', f'RunX_htf_KRITI_300_100%_300_3loci',
                 f'RunX_htf_NOTFOUND_58_100%_57_3loci', f'RunX_htf_ABROAD_43_100%_43_3loci',
                ]
        pops = [haplomat_output_parce(file, int(file.split('_')[-2])) for file in files]
        #print(len(set(all_haplotypes)))
        #print(all_haplotypes)

        with pd.ExcelWriter(f"Haplotype_fisher_between_regions_{args.loci3}.xlsx", engine = "openpyxl", mode = "w") as writer:
            for i in range(len(files) - 1):
                df = pd.DataFrame()
                iterlist = [(files[i],file) for file in files[i+1:]]
                all_haplotypes = [list(pop.keys()) for pop in pops[i:]]
                all_haplotypes = sum(all_haplotypes, [])
                all_haplotypes = sorted(set(all_haplotypes))
                #all_comb = combinations(files, 2)
                df['Haplotypes'] = all_haplotypes
                for file_1, file_2 in list(iterlist):
                    results,_,_ = compare_haplotypes(
                                                pop_1 = file_1,        #haplomat output for one population
                                                pop_2 = file_2,      #haplomat output for another population
                                                num_1 = int(file_1.split('_')[-2]),     #number of samples in one population
                                                num_2 = int(file_2.split('_')[-2]), #number of samples in another population
                                                merged_haplotypes = all_haplotypes,
                                                no_diagram = None, #returns venn diagramm
                                                pvals = True,
                                                )
                    #print(type(results))
                    #print(len(results)
                    df[f'counts_{file_1.split("htf_")[1]}'] = [results[haplotype][f'{file_1.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                    df[f'counts_{file_2.split("htf_")[1]}'] = [results[haplotype][f'{file_2.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                    df[f'HF_{file_1.split("htf_")[1]}'] = [results[haplotype][f'{file_1.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                    df[f'HF_{file_2.split("htf_")[1]}'] = [results[haplotype][f'{file_2.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                    df[f'Fisher_Pvalue_{file_1.split("htf_")[1]}-{file_2.split("htf_")[1]}'] = [results[haplotype][f'Fisher_{file_1.split("htf_")[1]}-{file_2.split("htf_")[1]}'] for haplotype in all_haplotypes]

                    #print(df)
                    #results = results.T
                    #df.columns = ['Haplotype',f'counts_{file_1.split("htf_")[1]}', f'counts_{file_2.split("htf_")[1]}', f'freq_{file_1.split("htf_")[1]}', f'freq_{file_2.split("htf_")[1]}', 'fisher_p_value']
                name = f'{files[i].split("htf_")[1]}'

                df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)

        #df = pd.DataFrame(results,index = [0]).T.reset_index()#, columns=['Haplotype', 'fisher_p_value'])
        #df = df.T
        #df.columns = ['Haplotype',f'counts_{args.pop_1.split("htf_")[1]}', f'counts_{args.pop_2.split("htf_")[1]}', f'freq_{args.pop_1.split("htf_")[1]}', f'freq_{args.pop_2.split("htf_")[1]}', 'fisher_p_value']
        #df.to_excel(f'{args.pop_1.split("htf_")[1]} - {args.pop_2.split("htf_")[1]}.xlsx', index = False)
    elif args.region_pairs and args.loci5:
        files = [
                 f'RunX_htf_ATTIKI_238_100%_227_5loci', f'RunX_htf_KRITI_300_100%_299_5loci',
                 f'RunX_htf_NOTFOUND_58_100%_53_5loci', f'RunX_htf_ABROAD_43_100%_33_5loci',
                ]
        pops = [haplomat_output_parce(file, int(file.split('_')[-2])) for file in files]
        #print(len(set(all_haplotypes)))
        #print(all_haplotypes)

        with pd.ExcelWriter(f"Haplotype_fisher_between_regions_5loci.xlsx", engine = "openpyxl", mode = "w") as writer:
            for i in range(len(files) - 1):
                df = pd.DataFrame()
                iterlist = [(files[i],file) for file in files[i+1:]]
                all_haplotypes = [list(pop.keys()) for pop in pops[i:]]
                all_haplotypes = sum(all_haplotypes, [])
                all_haplotypes = sorted(set(all_haplotypes))
                #all_comb = combinations(files, 2)
                df['Haplotypes'] = all_haplotypes
                for file_1, file_2 in list(iterlist):
                    results,_,_ = compare_haplotypes(
                                                pop_1 = file_1,        #haplomat output for one population
                                                pop_2 = file_2,      #haplomat output for another population
                                                num_1 = int(file_1.split('_')[-2]),     #number of samples in one population
                                                num_2 = int(file_2.split('_')[-2]), #number of samples in another population
                                                merged_haplotypes = all_haplotypes,
                                                no_diagram = None, #returns venn diagramm
                                                pvals = True,
                                                )
                    #print(type(results))
                    #print(len(results)
                    df[f'counts_{file_1.split("htf_")[1]}'] = [results[haplotype][f'{file_1.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                    df[f'counts_{file_2.split("htf_")[1]}'] = [results[haplotype][f'{file_2.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                    df[f'HF_{file_1.split("htf_")[1]}'] = [results[haplotype][f'{file_1.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                    df[f'HF_{file_2.split("htf_")[1]}'] = [results[haplotype][f'{file_2.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                    df[f'Fisher_Pvalue_{file_1.split("htf_")[1]}-{file_2.split("htf_")[1]}'] = [results[haplotype][f'Fisher_{file_1.split("htf_")[1]}-{file_2.split("htf_")[1]}'] for haplotype in all_haplotypes]

                    #print(df)
                    #results = results.T
                    #df.columns = ['Haplotype',f'counts_{file_1.split("htf_")[1]}', f'counts_{file_2.split("htf_")[1]}', f'freq_{file_1.split("htf_")[1]}', f'freq_{file_2.split("htf_")[1]}', 'fisher_p_value']
                name = f'{files[i].split("htf_")[1]}'

                df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)

    elif args.cyprus:
        files = [
                 'RunX_htf_KM_866_100%_842_3loci','RunX_htf_ATTIKI_238_100%_238_3loci', 'RunX_htf_KRITI_300_100%_300_3loci',
                 'RunX_htf_Cyprus_2843_3loci',
                ]
        files = [
                 'RunX_htf_ATTIKI_238_100%_227_5loci', 'RunX_htf_KRITI_300_100%_299_5loci',
                 'RunX_htf_Cyprus_2702_5loci',
                ]

        pops = [haplomat_output_parce(file, int(file.split('_')[-2])) for file in files]

        with pd.ExcelWriter(f"Haplotype_fisher_regions_VS_CYPRUS_5loci.xlsx", engine = "openpyxl", mode = "w") as writer:
            df = pd.DataFrame()
            all_haplotypes = [list(pop.keys()) for pop in pops]
            all_haplotypes = sum(all_haplotypes, [])
            all_haplotypes = sorted(set(all_haplotypes))
            df['Haplotypes'] = all_haplotypes
            for file in files[:-1]:
                results,_,_ = compare_haplotypes(
                                            pop_1 = file,        #haplomat output for one population
                                            pop_2 = files[-1],      #haplomat output for another population
                                            num_1 = int(file.split('_')[-2]),     #number of samples in one population
                                            num_2 = int(files[-1].split('_')[-2]), #number of samples in another population
                                            merged_haplotypes = all_haplotypes,
                                            no_diagram = None, #returns venn diagramm
                                            pvals = True,
                                            )

                df[f'counts_{files[-1].split("htf_")[1]}'] = [results[haplotype][f'{files[-1].split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                df[f'counts_{file.split("htf_")[1]}'] = [results[haplotype][f'{file.split("htf_")[1]}'][0] for haplotype in all_haplotypes]
                df[f'HF_{files[-1].split("htf_")[1]}'] = [results[haplotype][f'{files[-1].split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                df[f'HF_{file.split("htf_")[1]}'] = [results[haplotype][f'{file.split("htf_")[1]}'][1] for haplotype in all_haplotypes]
                df[f'Fisher_Pvalue_{file.split("htf_")[1]}-{files[-1].split("htf_")[1]}'] = [results[haplotype][f'Fisher_{file.split("htf_")[1]}-{files[-1].split("htf_")[1]}'] for haplotype in all_haplotypes]

            name = f'regionsvsCYPRUS'

            df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)


    else:
        files = [
                'RunX_htf_TOTAL_1739_100%_1707_3loci', 'RunX_htf_TOTAL_1072_1054_3loci'
                ]
        files = [
                'RunX_htf_TOTAL_1739_100%_574_5loci', 'RunX_htf_TOTAL_1072_331_5loci'
                ]
        files = [
                'RunX_htf_TOTAL_2921_2870_3loci', 'RunX_htf_Cyprus_2843_3loci'
                ]
        files = [
                'RunX_htf_TOTAL_2921_992_5loci', 'RunX_htf_Cyprus_2702_5loci'
                ]
        files = [
                'RunX_htf_TOTAL_1739_100%_1707_3loci', 'RunX_htf_Cyprus_2843_3loci'
                ]
        files = [
                'RunX_htf_TOTAL_1739_100%_574_5loci', 'RunX_htf_Cyprus_2702_5loci'
                ]
        files = [
                'RunX_htf_TOTAL_1072_1054_3loci', 'RunX_htf_Cyprus_2843_3loci'
                ]
        files = [
                'RunX_htf_TOTAL_1072_331_5loci', 'RunX_htf_Cyprus_2702_5loci'
                ]

        #ena = '1739'
        #duo = '1072'

        #ena = '1072'


        #ena = '2921'
        duo = 'CYPRUS'

        loci = '5loci'

        pops = [haplomat_output_parce(file, int(file.split('_')[-2])) for file in files]
        #print(len(set(all_haplotypes)))
        #print(all_haplotypes)

        with pd.ExcelWriter(f"Haplotype_fisher_{ena}_VS_{duo}_{loci}.xlsx", engine = "openpyxl", mode = "w") as writer:
            df = pd.DataFrame()
            all_haplotypes = [list(pop.keys()) for pop in pops]
            all_haplotypes = sum(all_haplotypes, [])
            all_haplotypes = sorted(set(all_haplotypes))
            df['Haplotypes'] = all_haplotypes

            results,_,_ = compare_haplotypes(
                                        pop_1 = files[0],        #haplomat output for one population
                                        pop_2 = files[1],      #haplomat output for another population
                                        num_1 = int(files[0].split('_')[-2]),     #number of samples in one population
                                        num_2 = int(files[1].split('_')[-2]), #number of samples in another population
                                        merged_haplotypes = all_haplotypes,
                                        no_diagram = None, #returns venn diagramm
                                        pvals = True,
                                        )

            df[f'counts_{files[0].split("htf_")[1]}'] = [results[haplotype][f'{files[0].split("htf_")[1]}'][0] for haplotype in all_haplotypes]
            df[f'counts_{files[1].split("htf_")[1]}'] = [results[haplotype][f'{files[1].split("htf_")[1]}'][0] for haplotype in all_haplotypes]
            df[f'HF_{files[0].split("htf_")[1]}'] = [results[haplotype][f'{files[0].split("htf_")[1]}'][1] for haplotype in all_haplotypes]
            df[f'HF_{files[1].split("htf_")[1]}'] = [results[haplotype][f'{files[1].split("htf_")[1]}'][1] for haplotype in all_haplotypes]
            df[f'Fisher_Pvalue_{files[0].split("htf_")[1]}-{files[1].split("htf_")[1]}'] = [results[haplotype][f'Fisher_{files[0].split("htf_")[1]}-{files[1].split("htf_")[1]}'] for haplotype in all_haplotypes]

            name = f'{ena}vs{duo}'

            df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)







'''
5 loci me g RESOLUTION
3 me g RESOLUTION
to 5 loci de ginetai G kefaleio gt yparxoun Null alllele


----------------->  BETWEEN REGIONS 3 LOCI  <-------------------: DONE

----------------->  BETWEEN REGIONS and TOTAL MINUS REGIONS 3 LOCI  <-------------------: DONE

----------------->  BETWEEN REGIONS AND DKMS 3loci  <-------------------: DONE

----------------->  BETWEEN REGIONS AND cyprus 3loci  <-------------------: DONE

----------------->  BETWEEN cyprus AND DKMS 3loci  <-------------------: DONE

----------------->  BETWEEN cyprus AND total greeks 3loci  <-------------------: DONE

----------------->  BETWEEN DKMS AND total greeks  3loci  <-------------------: DONE



----------------->  BETWEEN REGIONS 5LOCI  <-------------------: DONE

----------------->  BETWEEN REGIONS and TOTAL MINUS REGIONS 5LOCI  <-------------------:

----------------->  BETWEEN REGIONS AND DKMS 5LOCI  <-------------------: DONE

----------------->  BETWEEN REGIONS AND cyprus 5LOCI  <-------------------: DONE

----------------->  BETWEEN cyprus AND DKMS 5LOCI  <-------------------: DONE

----------------->  BETWEEN cyprus AND total greeks 5LOCI  <-------------------: DONE

----------------->  BETWEEN DKMS AND total greeks  5LOCI  <-------------------: DONE






'''
'''
elif args.dkms:

    files = [
             'RunX_htf_THESSALY_148_3loci','RunX_htf_ST_ELLADA_19_3loci', 'RunX_htf_PELOPONNISOS_34_3loci',
             'RunX_htf_N_AIGAIO_9_3loci', 'RunX_htf_KRITI_370_3loci', 'RunX_htf_KEN_MAKED_1130_3loci',
             'RunX_htf_IONIA_7_3loci','RunX_htf_HPEIROS_32_3loci', 'RunX_htf_DYT_MAKED_136_3loci',
             'RunX_htf_DYT_ELLADA_57_3loci', 'RunX_htf_B_AIGAIO_6_3loci', 'RunX_htf_ATTIKI_409_3loci',
             'RunX_htf_AN_MAKED_243_3loci', 'RunX_htf_Greece_10000_181022_3loci'
            ]

    files = [
             f'RunX_htf_THESSALY_20_5loci',f'RunX_htf_ST_ELLADA_14_5loci', f'RunX_htf_PELOPONNISOS_28_5loci',
             f'RunX_htf_N_AIGAIO_4_5loci', f'RunX_htf_KRITI_364_5loci', f'RunX_htf_KEN_MAKED_42_5loci',
             f'RunX_htf_IONIA_1_5loci',f'RunX_htf_HPEIROS_14_5loci', f'RunX_htf_DYT_MAKED_3_5loci',
             f'RunX_htf_DYT_ELLADA_34_5loci', f'RunX_htf_B_AIGAIO_4_5loci', f'RunX_htf_ATTIKI_325_5loci',
             f'RunX_htf_AN_MAKED_11_5loci', 'RunX_htf_Greece_10000_181022_5loci'
            ]


    with pd.ExcelWriter(f"Haplotype_fisher_regions_VS_DKMS_{loci}.xlsx", engine = "openpyxl", mode = "w") as writer:
        for file in files[:-1]:
            df = compare_haplotypes(
                                    pop_1 = file,        #haplomat output for one population
                                    pop_2 = files[-1],      #haplomat output for another population
                                    num_1 = int(file.split('_')[-2]),     #number of samples in one population
                                    num_2 = 10000,      #number of samples in another population
                                    no_diagram = None, #returns venn diagramm
                                    pvals = True       #returns venn p-values
                                   )
            df = df.T
            df.columns = ['Haplotype',f'counts_{file.split("htf_")[1]}', f'counts_DKMS_{loci}', f'freq_{file.split("htf_")[1]}', f'freq_DKMS_{loci}', 'fisher_p_value']
            name = f'{file.split("htf_")[1].split(f"_{loci}")[0]}-DKMS'

            df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)
elif args.minus:

    files_regions = [
             'RunX_htf_THESSALY_148_3loci','RunX_htf_ST_ELLADA_19_3loci', 'RunX_htf_PELOPONNISOS_34_3loci',
             'RunX_htf_N_AIGAIO_9_3loci', 'RunX_htf_KRITI_370_3loci', 'RunX_htf_KEN_MAKED_1130_3loci',
             'RunX_htf_IONIA_7_3loci','RunX_htf_HPEIROS_32_3loci', 'RunX_htf_DYT_MAKED_136_3loci',
             'RunX_htf_DYT_ELLADA_57_3loci', 'RunX_htf_B_AIGAIO_6_3loci', 'RunX_htf_ATTIKI_409_3loci',
             'RunX_htf_AN_MAKED_243_3loci',
            ]

    files_minus = ['RunX_htf_THESSALYvsTOTAL_2452_3loci', 'RunX_htf_ST_ELLADAvsTOTAL_2581_3loci', 'RunX_htf_PELOPONNvsTOTAL_2566_3loci',
                   'RunX_htf_N_AIGAIOvsTOTAL_2591_3loci','RunX_htf_KRITIvsTOTAL_2230_3loci', 'RunX_htf_KEN_MAKEDvsTOTAL_1470_3loci',
                   'RunX_htf_IONIAvsTOTAL_2593_3loci', 'RunX_htf_HPEIROSvsTOTAL_2568_3loci', 'RunX_htf_DYT_MAKEDvsTOTAL_2464_3loci',
                   'RunX_htf_DYT_ELLADAvsTOTAL_2543_3loci', 'RunX_htf_B_AIGAIOvsTOTAL_2594_3loci','RunX_htf_ATTIKIvsTOTAL_2191_3loci',
                   'RunX_htf_AN_MAKEDvsTOTAL_2357_3loci',
                 ]

    files_regions = [
             f'RunX_htf_THESSALY_20_5loci',f'RunX_htf_ST_ELLADA_14_5loci', f'RunX_htf_PELOPONNISOS_28_5loci',
             f'RunX_htf_N_AIGAIO_4_5loci', f'RunX_htf_KRITI_364_5loci', f'RunX_htf_KEN_MAKED_42_5loci',
             f'RunX_htf_IONIA_1_5loci',f'RunX_htf_HPEIROS_14_5loci', f'RunX_htf_DYT_MAKED_3_5loci',
             f'RunX_htf_DYT_ELLADA_34_5loci', f'RunX_htf_B_AIGAIO_4_5loci', f'RunX_htf_ATTIKI_325_5loci',
             f'RunX_htf_AN_MAKED_11_5loci',
            ]
    files_minus = ['RunX_htf_THESSALYvsTOTAL_844_5loci', 'RunX_htf_ST_ELLADAvsTOTAL_850_5loci', 'RunX_htf_PELOPONNvsTOTAL_836_5loci',
                   'RunX_htf_N_AIGAIOvsTOTAL_860_5loci','RunX_htf_KRITIvsTOTAL_500_5loci', 'RunX_htf_KEN_MAKEDvsTOTAL_822_5loci',
                   'RunX_htf_IONIAvsTOTAL_863_5loci', 'RunX_htf_HPEIROSvsTOTAL_850_5loci', 'RunX_htf_DYT_MAKEDvsTOTAL_861_5loci',
                   'RunX_htf_DYT_ELLADAvsTOTAL_830_5loci', 'RunX_htf_B_AIGAIOvsTOTAL_860_5loci','RunX_htf_ATTIKIvsTOTAL_539_5loci',
                   'RunX_htf_AN_MAKEDvsTOTAL_853_5loci',
                 ]




    with pd.ExcelWriter(f"Haplotype_fisher_regions_VS_TOTAL_MINUS_{loci}.xlsx", engine = "openpyxl", mode = "w") as writer:
        for file_1,file_2 in zip(files_regions,files_minus):
            df = compare_haplotypes(
                                    pop_1 = file_1,        #haplomat output for one population
                                    pop_2 = file_2,      #haplomat output for another population
                                    num_1 = int(file_1.split('_')[-2]),     #number of samples in one population
                                    num_2 = int(file_2.split('_')[-2]),      #number of samples in another population
                                    no_diagram = None, #returns venn diagramm
                                    pvals = True       #returns venn p-values
                                   )
            df = df.T
            df.columns = ['Haplotype',f'counts_{file_1.split("htf_")[1]}', f'counts_TOTAL_MINUS', f'freq_{file_1.split("htf_")[1]}', f'freq_TOTAL_MINUS', 'fisher_p_value']
            name = f'{file_1.split("htf_")[1].split(f"_{loci}")[0]}-TOTAL_MINUS'

            df.to_excel(writer, sheet_name = name, engine = "openpyxl", index = False)
'''
