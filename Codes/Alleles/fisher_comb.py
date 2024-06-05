
import pandas as pd
import numpy as np
import argparse, itertools
from collections import defaultdict
from scipy.stats import fisher_exact
from utils import read_data,gather_alleles
from alllele_frequencies import allele_freq
from fisher import fisher_anal
from itertools import combinations

parser = argparse.ArgumentParser(description='Terminal arguments')
parser.add_argument("--alltoone", action="store_true", dest="alltoone", default=False,
					help="If you want to compare multiple variables with a specific one.")
parser.add_argument("--alltoall", action="store_true", dest="alltoall", default=False,
					help="If you want to compare multiple variables by two.")

if __name__ == "__main__":
	args = parser.parse_args()

	if args.alltoone:
		pass
	elif args.alltoall:
		#data_files = ['Greece/Greek_CBUS_2fields_3loci', 'Greece/Greek_BMDs_2fields_3loci_75599', 'Greece/Greek_CBUs_BMDs_2fields_3loci_78611',
		#              'Cyprus/cyprus_cbus_2_fields_3loci', 'Cyprus/cyprus_donors_2_fields_3loci', 'Cyprus/cyprus_cbus_plus_donors_2_fields_3loci',
		#             ]
		#data_files = ['Greece/Greek_CBUS_2fields_5loci', 'Greece/Greek_BMDs_2fields_5loci_69130', 'Greece/Greek_CBUs_BMDs_2fields_5loci_70222',
		#              'Cyprus/cyprus_cbus_2_fields_5loci', 'Cyprus/cyprus_donors_2_fields_5loci', 'Cyprus/cyprus_cbus_plus_donors_2_fields_5loci',
		#             ]
		#data_files = ['Other/DHTOB_for_HAPLOMAT', 'Other/IIBEAA_for_HAPlOMAT', 'Other/PAP_for_Haplomat',
		#              'Greece/Greek_BMDs_2fields_3loci_75599',
		#              'Cyprus/cyprus_cbus_2_fields_3loci', 'Cyprus/cyprus_donors_2_fields_3loci',
		#              'Cyprus/cyprus_cbus_plus_donors_2_fields_3loci',
		#             ]

		#data_files = ['Other/DHTOB_for_HAPLOMAT', 'Other/IIBEAA_for_HAPlOMAT',
		#              'Greece/Greek_BMDs_2fields_5loci_69130',
		#              'Cyprus/cyprus_cbus_2_fields_5loci', 'Cyprus/cyprus_donors_2_fields_5loci',
		#              'Cyprus/cyprus_cbus_plus_donors_2_fields_5loci',
		#             ]

		data_files = [
					  #'Greece/Greek_BMDs_2fields_3loci_75599','Cyprus/cyprus_donors_2_fields_3loci',
					 ]


		dkms_files = [
					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_CBUs_2fields_1092_5loci.dat',
					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci.dat',

					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_CBUs_2fields_3012_3loci.dat',
					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci.dat',

					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_DHTOB_428_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_BRFAA_592_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_PAP_1992_3loci.dat',
					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_DHTOB_428_5loci.dat',
					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_BRFAA_564_5loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields__98998_3loci.dat',
					  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Cyprus/RunX_htf_cyprus_cbus_plus_donors_2fields_75706_5loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Austria_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Bosnia_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Croatia_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Dutch_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_France_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Greece_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Italy_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Portugal_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Romania_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Spain_10000_3loci.dat',
					  f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/DKMS/RunX_htf_Turkey_10000_3loci.dat',
					  ]



		'''
		dkms_files = [  #f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_CBUs_2fields_1092_5loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_BMDs_CBUs_2fields_70222_5loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_DHTOB_428_5loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_BRFAA_564_5loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/RunX_htf_Group_1_133004_5loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/RunX_htf_Group_2_50000_5loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/RunX_htf_12_Countries_183004_5loci.dat',
						#'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_CBUs_2fields_3012_3loci.dat',
						#'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Greece/RunX_htf_Greek_BMDs_CBUs_2fields_78611_3loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_DHTOB_428_3loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_BRFAA_592_3loci.dat',
						#f'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/Banks/RunX_htf_PAP_1992_3loci.dat',
						#'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/RunX_htf_Group_1_156155_3loci.dat',
						#'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/RunX_htf_Group_2_50000_3loci.dat',
						#'/Users/vasou/Documents/HLA/Matching Coverage/Estimated_Haplotypes/RunX_htf_12_Countries_206155_3loci.dat',

					 ]
		'''



		#dkms_files = []
		banks = []
		for file in data_files:
			banks.append((file.split('/')[1], read_data(file, sheet = '3loci') if 'for' in file  else read_data(file,)))#sheet = '3loci' for banks

		for dkms in dkms_files:
			banks.append((dkms.split('RunX_htf_')[1], dkms))
		#print(banks)
		with pd.ExcelWriter('testttt.xlsx', engine = "openpyxl", mode = "w") as writer:
			for gene in ['A', 'B', 'DRB1',]:#['A', 'B', 'C','DRB1', 'DQB1', 'DPB1']
				#df = pd.DataFrame()
				print(gene)
				result = fisher_anal(banks = banks,
								target_gene = gene,
								comb_iterlist = True,
								loci = '3loci'
								)
				#df = pd.concat([df, pd.DataFrame(result)], axis = 1)
				df = pd.DataFrame(result)

				df.to_excel(writer, sheet_name = f'{gene}', engine = "openpyxl", index = False)
