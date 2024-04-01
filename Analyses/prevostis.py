import pandas as pd
import numpy as np
import seaborn as sns
import warnings, argparse
import matplotlib.pyplot as plt
from collections import defaultdict

user_path = '/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/'
def haplomat_output_parce(hfs_dat, n):
    with open(hfs_dat) as f:
        cretan_haplotypes = [x.replace('\n', '').replace('g', '').split() for x in f.readlines()]
    cretan_haplotypes = [( '-'.join(sorted(x[0].split('~'))  ), float(x[1])) for x in cretan_haplotypes]
    cretan_haplotypes = {x[0]:x[1] for x in cretan_haplotypes}
    cretan_haplotypes = {k:v for k,v in cretan_haplotypes.items() if v >= 1/(2*n)}
    return cretan_haplotypes


def calculate_distances(population_list, population_names, euclidean=None, prevosti=None, chi=None):
    if sum([x for x in [euclidean,prevosti] if x])==0 or sum([x for x in [euclidean,prevosti] if x])==2:
        raise Exception("Choose either Euclidean or Prevosti's distance.")
    arr = np.zeros((len(population_names),len(population_names)))
    for x in range(len(population_list)):
        for y in range(len(population_names)):
            if population_list[x] != population_list[y]:

                all_haplotypes = list(set(list(population_list[x].keys())+list(population_list[y].keys())))
                a = [population_list[x][z] if z in population_list[x].keys() else 0 for z in all_haplotypes]
                b = [population_list[y][z] if z in population_list[y].keys() else 0 for z in all_haplotypes]
                if euclidean:
                    arr[x][y] += np.linalg.norm(np.array(a)-np.array(b))
                elif prevosti: # similar to Manhattan distance
                    arr[x][y] += np.sum(np.abs(np.array(a)-np.array(b)))*0.5
    df = pd.DataFrame(data=arr, index=population_names, columns=population_names)
    return df

Greece = f'{user_path}/Greece/RunX_htf_Greek_BMDs_2fields_70077_5loci.dat'
crete = f'{user_path}/Greece/RunX_htf_CreteAdults_1582_5loci.dat'
Cyprus = f'{user_path}/Cyprus/RunX_htf_cyprus_donors_2fields_73004_5loci.dat'
Austria = f'{user_path}/DKMS/RunX_htf_Austria_10000_5loci.dat'
Bosnia = f'{user_path}/DKMS/RunX_htf_Bosnia_10000_5loci.dat'
Croatia = f'{user_path}/DKMS/RunX_htf_Croatia_10000_5loci.dat'
Dutch = f'{user_path}/DKMS/RunX_htf_Dutch_10000_5loci.dat'
France = f'{user_path}/DKMS/RunX_htf_France_10000_5loci.dat'
Greece_DKMS = f'{user_path}/DKMS/RunX_htf_Greece_10000_5loci.dat'
Italy = f'{user_path}/DKMS/RunX_htf_Italy_10000_5loci.dat'
Portugal = f'{user_path}/DKMS/RunX_htf_Portugal_10000_5loci.dat'
Romania = f'{user_path}/DKMS/RunX_htf_Romania_10000_5loci.dat'
Spain = f'{user_path}/DKMS/RunX_htf_Spain_10000_5loci.dat'
Turkey = f'{user_path}/DKMS/RunX_htf_Turkey_10000_5loci.dat'

countries = [Greece, crete, Cyprus, Austria, Bosnia, Croatia, Dutch, France, Greece_DKMS, Italy, Portugal, Romania, Spain, Turkey]
ready_country = [haplomat_output_parce(country,int(country.split('_')[-2])) for country in countries]
country_names = [
                'Greece', 'Crete', 'Cyprus', 'Austria','Bosnia', 'Croatia', 'Dutch','France',
                'Greece_DKMS', 'Italy', 'Portugal', 'Romania', 'Spain','Turkey'
                ]

parser = argparse.ArgumentParser(description='',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metric', action='store_true')
args = parser.parse_args()

'''
Pass
warning: https://github.com/widdowquinn/pyani/issues/73
'''
warnings.filterwarnings("ignore")

print('Calculating Prevostis distances.')
prev_distances = calculate_distances(
                                    population_list = ready_country, #haplotype info
                                    population_names = country_names, #name of countris
                                    prevosti = args.metric, #metric
                                    )
print('Done!')

g = sns.clustermap(prev_distances,)
plt.tight_layout()
g.savefig(f"ClusterMap.png", dpi=300)
print('Saved png: ClusterMap.png')
