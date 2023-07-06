# Calculation of Fixation Index (Fst)

This repository contains a Python code for performing various calculations and analyses related to allele frequencies, heterozygosity, and F-statistics. The code is organized into several functions that can be used in a pipeline to process genetic data.

## Usage
To use the code, follow the steps below:

* Import the required functions from the code file.


* Call the __fst_pipeline()__ function, passing in your data as an argument:

    ```javascript
    result = fst_pipeline(data)
    ```
    The data parameter represents the input data required for the analysis. It's a      dictionary with multiple keys (chania, rethymno, heraklio, lasithi) and each value is a dataframe that contains, for each sample, the ID and the typing for the genes we are interested (A, B, C, DRB1, DQB1, DPB1)

* The fst_pipeline function will return a tuple containing various results:
    ```javascript
    (alleles_all_regions, allele_count, 
    obs_gen_count, allele_freq,counts_HWE,
    excess_sub,excess_total,obs_het,expected_het,
    allele_bar, Het_indices, Fstatistics) = fst_pipeline(data)
    ```
    You can use the returned results for further analysis or visualization as needed.

* Also, you can call each step of this pipeline separately.

## Explanation of each function

* __unique_and_sorted_allele(data)__:
    
    * This function is called to extract the unique alleles (for each gene) from the given data across all regions.
    * The alleles are sorted and stored in the __alleles_all_regions__ dictionary.


* __allele_genotype_counts (data, type_count)__

    * This function ireturns the counts of alleles or observed genotypes for each region. Just specify in __type_count__ variable the calculation you want (Allele or Genotype)
    * The results are stored in the __allele_count__ or __obs_gen_count__ variables, respectively.

* __allele_frequencies(data)__

    * This function is called to compute the allele frequencies for each region.
    * The results are stored in the __allele_freq__ dictionary.

* __HWE(data)__

    * The Hardy-Weinberg Equilibrium (HWE) function is called to determine the expected genotypic counts under Hardy-Weinberg equilibrium.
    * The results are stored in the __counts_HWE__ dictionary.

* __deficiencies_sub(data)__

    * This function calculates excess or dificiencies of observed homozygotes in each subpopulation relative to __HWE__ genotypic counts.
    * The results are stored in the __excess_sub__ dictionary.

* __deficiencies_total(data)__
   
    * This function calculates excess or dificiencies of observed homozygotes in the total population
    * The results are stored in the __excess_total__ dictionary.

* __obs_heterozygosities(data)__

    * This function calculates the local observer heterozygosities of each subpopulation.
    * The results are stored in the __obs_het__ dictionary.

* __exp_heterozygosities(data)__

    * This function calculates the expected heterozygosities of each subpopulation.
    * The results are stored in the __expected_het__ dictionary.

* __p_allele_bar(data)__

    * Calculates the global frequency of each allele across all subpopulations and weighted by each subpopulation size. Also, calculates the global frequency of each allele for each pair of regions.
    * The results are stored in the __expected_het__ dictionary.

* __heterozygosity_indices(data)__
    * Calculates heterozygosity indices:
        * Hi=weighted average of observed heterozygosities
        * Hs=weighted average of expected heterozygosities
        * Ht=expected hetrozygosty based on the global gene frequencies
    * The results are stored in the __Het_indices__ dictionary.

* __global_Fstatistics(data)__
    * Calculates global Fstatistics:
        * FIS_observed heterozygosities
        * FST_expected heterozygosities
        * FIT_expected heterozygosities
    * The results are stored in the __Fstatistics__ dictionary.



