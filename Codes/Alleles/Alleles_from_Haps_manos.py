### In order to compare the allele frequencies of the Cretan population with the DKMS minorities we must
### group the alleles with g resolution as the DKMS did. We will estimate first the haplotype frequencies
### with the Hapl-o-Mat software, because it can perform the g grouping and then, from the haplotype
### frequencies we will calculate the grouped allele frequencies and compare them with the DKMS minorities'.
### Arlequin cannot perform the grouping, so that is why it is not a suitable program for this task, although
### it does not drop the genotypes with at least one missing value. If we do not drop these lines, the
### estimated allele frequency would be closer to the real frequency of the sample.

### After we run the Hapl-o-Mat with the g grouping without missing data we calculate the allele frequencies

def af_estimation_from_HaploMat_g_group(hfs_dat, hla_gene, num_of_people, alleles=None, counts=None):

    '''
    !!!README!!!

    This function calculates the allele frequency of the HLA genes A, B, C, DRB1, DQB1 and DPB1
    if the haplotypes have them.

    Hapl-o-Mat orders the genes in the haplotypes alphabetically.

    The haplotypes with 4 genes should have the HLA genes below in this order only:
    A~B~C~DRB1

    The haplotypes with 5 genes should have the HLA genes below in this order only:
    A~B~C~DQB1~DRB1

    The haplotypes with 6 genes should have the HLA genes below in this order only:
    A~B~C~DPB1~DQB1~DRB1

    Other combinations of the HLA genes can not be used as input.

    This function returns the estimated allele frequencies/counts and alleles of the sample from Hapl-o-Mat.
    '''

    if hla_gene not in ('A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1'):
        raise Exception("The hla_gene can only be 'A', 'B', 'C', 'DRB1', 'DQB1' or 'DPB1'")

    if counts and alleles:
        raise Exception('Only 1 output can come from this function.')

    with open(hfs_dat, 'r') as f:
        f = f.readlines()

    f = [x.strip('\n').split('\t') for x in f]
    allele_freq = [(x[0].split('~'), float(x[1].replace(',','.'))) for x in f]
    allele_freq = [([y.strip('g') for y in x[0]],x[1]) for x in allele_freq]

    '''
    allele_freq = (['A*02:01', 'B*18:01', 'C*07:01', 'DQB1*03:01', 'DRB1*11:04'], 0.02880986622617),
    '''

    if hla_gene == 'A':
        n = 0
    elif hla_gene == 'B':
        n = 1
    elif hla_gene == 'C':
        n = 2

    elif hla_gene == 'DRB1': # Hapl-o-Mat orders the genes in the haplotypes alphabetically.
        if len(allele_freq[0][0]) == 3:
            n = 2
        elif len(allele_freq[0][0]) == 4:
            n = 3
        elif len(allele_freq[0][0]) == 5:
            n = 4
        elif len(allele_freq[0][0]) == 6:
            n = 5
        else:
            raise Exception('The haplotypes should have 4, 5 or 6 genes in a specific order. Read the documentation.')

    elif hla_gene == 'DQB1':
        if len(allele_freq[0][0]) == 4:
            raise Exception('DQB1 should not be in a haplotype with only 4 genes. Read the documentation.')
        elif len(allele_freq[0][0]) == 5:
            n = 3
        elif len(allele_freq[0][0]) == 6:
            n = 4
        else:
            raise Exception('The haplotypes should have 4, 5 or 6 genes in a specific order. Read the documentation.')

    elif hla_gene == 'DPB1':
        if len(allele_freq[0][0]) == 4:
            raise Exception('DPB1 should not be in a haplotype with only 4 genes. Read the documentation.')
        elif len(allele_freq[0][0]) == 5:
            raise Exception('DPB1 should not be in a haplotype with only 5 genes. Read the documentation.')
        elif len(allele_freq[0][0]) == 6:
            n = 3
        else:
            raise Exception('The haplotypes should have 4, 5 or 6 genes in a specific order. Read the documentation.')

    HLA_gene = {el:{'N':0, 'AF':0} for el in [x[0][n].split('*')[1] for x in allele_freq]}
    for x in allele_freq:
        HLA_gene[x[0][n].split('*')[1]]['AF'] += x[1]
    for k,v in HLA_gene.items():
        #print(v)
        HLA_gene[k]['N'] = round(v['AF']*2*num_of_people)
    #HLA_gene = {k:v for k, v in sorted(HLA_gene.items(), key=lambda item: item[1], reverse=True)}

    if alleles:
        return list(HLA_gene.keys())

    elif counts:
        return {k:round(v*num_of_people) for k,v in HLA_gene.items()}
        # No need to multiply by 2, because of the 2 genotypes on MAC input

    else:
        result = {k:v for k,v in HLA_gene.items() if v['AF'] > round(1/(2*num_of_people), 6)}
        return result
