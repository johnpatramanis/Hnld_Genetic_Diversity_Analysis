import os
import os.path
from os import listdir
from os.path import isfile, join
import sys
import random
import numpy as np
from itertools import combinations
import statistics
import matplotlib.pyplot as plt
import itertools


##################################################################################################################################################################################################################################################################################################
##Set Up


GENOTYPE_FILE=open(sys.argv[1],'r')
SAMPLE_NAME=sys.argv[1].split('.gntp')[0]

OUTPUT=open(F'{SAMPLE_NAME}.metric','w')

PROTEIN_COVERAGE_FILE=open('Protein_Coverage.txt','r')



##### What is the test?

N_of_Ancient_Protein_Samples=20
N_of_APSamples_Polymorphisms=0




########### Data Structure


GENOTYPES=[]   ## List of lists, of lists - paired genotypes for each locus
FLAT_GENOTYPES=[] ## List of lists, genotypes for each locus

EXPECTED_HETEROZ=[] ## List of floats


###############################
################# Prepare data



###############################
###### Protein Data, get total number of AA covered by this analysis and use it to calculate the length of the underlying sequence
Number_of_AA=[]

for line in PROTEIN_COVERAGE_FILE:
    line=line.strip().split()
    if len(line)>1:
        AA=line[1].split(',')
        for K in AA:
            Number_of_AA.append(K)

Number_of_AA=len(Number_of_AA)
Length_of_Sequence=3*Number_of_AA






###############################
###### Prep Genotypes


count=1
for LINE in GENOTYPE_FILE:
    LINE=LINE.strip().split('\t')
    
    ##### N of dataset
    NUMBER_OF_SAMPLES=len(LINE) ## How many samples

    ##### Get Genotypes, check which format is used by file 
    #### When format is 0/0
    if ('0/0' or '0/1' or '1/0' or '1/1') in LINE:
        GENOTYPES_HERE=[ x.split('/') for x in LINE ]
        
    ##### When format is 0|0
    if ('0|0' or '0|1' or '1|0' or '1|1') in LINE:
        GENOTYPES_HERE=[ x.split('|') for x in LINE ]
        

    FLAT_GENOTYPES_HERE=[ l for x in GENOTYPES_HERE for l in x]
    
    FREQ_GENOTYPES_HERE=[int(x) for x in FLAT_GENOTYPES_HERE if x!='.']
    FREQ=statistics.mean(FREQ_GENOTYPES_HERE)
    print(F'\nVariant {count} with frequency: {FREQ}\n')
    count+=1
    
    GENOTYPES.append(GENOTYPES_HERE) ## List of lists of pairs of genotypes
    FLAT_GENOTYPES.append(FLAT_GENOTYPES_HERE) ## All Genotypes in one list





################################################################################################################################################### #################################################################################################################################################
######## Calculate PIE - Expected Heterozygosity    
    
    ALLELE_NUMBER=len([x for x in FLAT_GENOTYPES_HERE if x!='.'])
    Q_FREQ=FLAT_GENOTYPES_HERE.count('0')/ALLELE_NUMBER
    P_FREQ=FLAT_GENOTYPES_HERE.count('1')/ALLELE_NUMBER
    
    EXPECTED_HETEROZ.append(2*Q_FREQ*P_FREQ)
    
    
    
EXPECTED_HETEROZ=sum(EXPECTED_HETEROZ)
    
    
    
print(F'\nFound {len(GENOTYPES)} potentially polymorphic sites for {NUMBER_OF_SAMPLES} individuals\n')
    
    
    
    
    
    
    
    
    
    
    

    
################################################################################################################################################# #################################################################################################################################################
##### Calculate PIE - AVG Heterozygosity




OBSERVED_HETEROZ=[]


for IND in range(0,NUMBER_OF_SAMPLES):

    IND_HETEROZ=[]
    
    for GENOTYPES_HERE in GENOTYPES:
        SAMPLE=GENOTYPES_HERE[IND]

        if (SAMPLE.count('0')==1) and (SAMPLE.count('1')==1):
            IND_HETEROZ.append(1)
        if (SAMPLE.count('0')==2) or (SAMPLE.count('1')==2):
            IND_HETEROZ.append(0)
            
            
    OBSERVED_HETEROZ.append(sum(IND_HETEROZ))

OBSERVED_HETEROZ=sum(OBSERVED_HETEROZ)/len(OBSERVED_HETEROZ)



















##################################################################################################################################################################################################################################################################################################
##### Calculate Watterssons estimate


##### Check that sites are indeed polymorphic in population
SEG_SITES=0
FIXED_SITES=0

for J in range(0,len(FLAT_GENOTYPES)):
    FLAT_GENOTYPES_HERE=[int(x) for x in FLAT_GENOTYPES[J] if x != '.']
    FLAT_GENOTYPES_HERE=list(set(FLAT_GENOTYPES_HERE))

    if len(FLAT_GENOTYPES_HERE)>1:
        SEG_SITES+=1
    if len(FLAT_GENOTYPES_HERE)==1:
        FIXED_SITES+=1
 


HARMONIC=sum([ 1/x for x in range(1,SEG_SITES*2)]) #### Harmonic for number of segregating sites * 2 for diploidy
if HARMONIC!=0:
    WTRSNS_E = SEG_SITES/HARMONIC
else:
    WTRSNS_E=0
    
    
    
print(F'\nNumber of Segregating sites: {SEG_SITES}\n')
print(F'Number of Fixed sites: {FIXED_SITES}\n')
print(F'Harmonic Number of Samples: {HARMONIC}\n')
print(F'\nThis results to a Watterssons Estimate of: {WTRSNS_E}\n')






















 
###############################################################################################################################################################################################################################################################################################################################
###### Random Sampling Loop

###
### How many times to perform test
MAX_LOOP=1000

### N_of_Ancient_Protein_Samples == number of available samples from proteins


##### If enough samples in vcf , then sample without replacement
if NUMBER_OF_SAMPLES >= N_of_Ancient_Protein_Samples*MAX_LOOP: ### Sample without replacement
    SAMPLINGS=random.sample(range(NUMBER_OF_SAMPLES), N_of_Ancient_Protein_Samples*MAX_LOOP)
    random.shuffle(SAMPLINGS)
    

    
##### If not enough samples in vcf , then sample with replacement   
if NUMBER_OF_SAMPLES < N_of_Ancient_Protein_Samples*MAX_LOOP: ### Sample with replacement
    SAMPLINGS=[]
    
    for k in range(0,MAX_LOOP):
        CHOICES=random.choices(range(NUMBER_OF_SAMPLES), k=N_of_Ancient_Protein_Samples)
        for j in CHOICES:
            SAMPLINGS.append(j)
            
            
    random.shuffle(SAMPLINGS)
### SAMPLINGS contains the sampled samples







##### These will be used to get the sampling metrics across the quartets

SAMPLING_SUCCESS_OR_NOT=[]
HOMOZYGOTE_SUCCESS_OR_NOT=[]
TOTAL_VARIANT_SUCCESS_OR_NOT=[]
TWO_ALTERNATIVES_SUCCESS_OR_NOT=[]



##### These will be used to get the average diveristy metrics for the quartets

WTRSNS_E_QUARTET=[]
TAJIMAS_THETA_QUARTET=[]
OBSERVED_HETEROZ_QUARTET=[]
EXPECTED_HETEROZ_QUARTET=[]




#### Sampling Loop, sample 'N_of_Ancient_Protein_Samples' individuals and check all segregating sites, do you find varaition in any of them?
for LOOP in range(0,MAX_LOOP):
    
    
    #### These will be used to get the sampling metrics for this quartet
    VARIANT_SPOTTED=0
    TOTAL_VARIANT_SPOTTED=0
    HOMOZYG_SPOTTED=0
    TWO_ALTERNATIVES_SPOTTED=0
    
    SAMPLING_NUMBER=SAMPLINGS[LOOP*N_of_Ancient_Protein_Samples:LOOP*N_of_Ancient_Protein_Samples+N_of_Ancient_Protein_Samples]
    
    ### This will be used for the Diversity Metricsf for this  quartet
    SEG_SITES_LOCAL=0
    WTRSNS_E_LOCAL=0
    TAJIMAS_THETA_LOCAL=[]
    OBSERVED_HETEROZ_LOCAL=[]
    EXPECTED_HETEROZ_LOCAL=[]
    
    
    
    
    
    
    

    

    
    #### For each site
    for SNP in range(0,len(GENOTYPES)):
    
        #### Load site data from matrix
        GENOTYPES_HERE=GENOTYPES[SNP]
        
        ##### Sample using method 1 # Sample diploid individuals
        SAMPLING=[ GENOTYPES_HERE[x] for x in SAMPLING_NUMBER] #### Get X individuals (diploid) using flat genotypes
        SAMPLING_DIPLOID=SAMPLING #### Keep X individuals (diploid) using flat genotypes
        SAMPLING=[l for x in SAMPLING for l in x] ### clean up data to be a flat list of 2*X alleles
        

        
        if '.' in SAMPLING:
            SAMPLING.remove('.')
        
        SAMPLING_UNIQ=list(set(SAMPLING)) #### either a (0) a (1) or (1,0)


        if len(SAMPLING_UNIQ)>1: ##### Check if any variation exists
            VARIANT_SPOTTED=1        ##### Count sucessful test
            TOTAL_VARIANT_SPOTTED+=1 ##### Count how many times Variant spotted has been triggered within a quartet!
            SEG_SITES_LOCAL+=1
            
            
        if (['1','1'] in SAMPLING_DIPLOID) and (['0','0'] in SAMPLING_DIPLOID): #### Check if homozygous individuals for the variant exist
            HOMOZYG_SPOTTED=1
        
        if (SAMPLING.count('1'))>=2: ### Check if more than 1 alternative allele exists
            TWO_ALTERNATIVES_SPOTTED=1
        
        
        
        
        
        
        #### Calc Expected Heterozygosity for this SNP in Sample
        
        ALLELE_NUMBER_LOCAL_SNP=len(SAMPLING) #number of non missing alleles
        Q_FREQ=SAMPLING.count('0')/ALLELE_NUMBER_LOCAL_SNP #Freq of allele 1
        P_FREQ=SAMPLING.count('1')/ALLELE_NUMBER_LOCAL_SNP #Freq of allele 2
        
        EXPECTED_HETEROZ_LOCAL.append(2*Q_FREQ*P_FREQ)
        
        
        
        
        
        
        
        
        
        
        
        #### Calculate Tajima's D for SNP, among all individuals of sampling    
        ### Compute all possible pairs of genotypes for 
        
        ALL_POSSIBLE_PAIRS=[ X for X in itertools.combinations(SAMPLING_DIPLOID,2) ]
        DIFS=[]
        for PAIR in ALL_POSSIBLE_PAIRS:
            
            #### First individual
            GENOTYPE_1=PAIR[0]
            while len(GENOTYPE_1)<2: ### fix if alleles are missing, 
                GENOTYPE_1.append(random.choices([0,1],[Q_FREQ,P_FREQ])[0]) ### here be careful random.choices returns a list of one item
            ###### fix missing genotype by picking one random 0 or 1, depending on the frequency of the alleles
            for x in range(0,len(GENOTYPE_1)):
                if (GENOTYPE_1[x] == '.' ) or (GENOTYPE_1[x] == 'X') or (GENOTYPE_1[x] == '-') or (GENOTYPE_1[x] == '?'):
                    GENOTYPE_1[x]=random.choices([0,1],[Q_FREQ,P_FREQ])[0]
            
            GENOTYPE_1=int(GENOTYPE_1[0]) + int(GENOTYPE_1[1])
            
            
            
            #### Second Individual
            GENOTYPE_2=PAIR[1]
            while len(GENOTYPE_2)<2: ### fix if alleles are missing, 
                GENOTYPE_2.append(random.choices([0,1],[Q_FREQ,P_FREQ])[0]) ### here random.choices returns a list of one item
            ####### fix missing genotype by picking one random 0 or 1, depending on the frequency of the alleles
            for x in range(0,len(GENOTYPE_2)):
                if (GENOTYPE_2[x] == '.' ) or (GENOTYPE_2[x] == 'X') or (GENOTYPE_2[x] == '-') or (GENOTYPE_2[x] == '?'):
                    GENOTYPE_2[x]=random.choices([0,1],[Q_FREQ,P_FREQ])[0]
            
            GENOTYPE_2=int(GENOTYPE_2[0]) + int(GENOTYPE_2[1])
            
            
            
            
            
            ##### Get difference between these 2 individuals on this SNP
            DIFS.append(abs(GENOTYPE_2-GENOTYPE_1))


        TAJIMAS_THETA_LOCAL.append(DIFS)    
        
        
        
        
        
        
 #############################################################################################################################################################################################        
        
        
        





    ##### This is done for every sampling (out of a 1000)


    ###############
    ### Append list with a success or not (1 or 0)
    SAMPLING_SUCCESS_OR_NOT.append(VARIANT_SPOTTED)
    TOTAL_VARIANT_SUCCESS_OR_NOT.append(TOTAL_VARIANT_SPOTTED)
    TWO_ALTERNATIVES_SUCCESS_OR_NOT.append(TWO_ALTERNATIVES_SPOTTED)
    HOMOZYGOTE_SUCCESS_OR_NOT.append(HOMOZYG_SPOTTED)
    
    
    
    ###############
    ### Calc Expected Heterozygosity for the sample set 
    ## Sum from every SNP 
    EXPECTED_HETEROZ_LOCAL=sum(EXPECTED_HETEROZ_LOCAL)
    
    
    
    
    
    
    
    
    ###############
    ### Calc Watterssons E for Sample Set
    
    HARMONIC_LOCAL=sum([ 1/x for x in range(1,SEG_SITES_LOCAL*2)]) #### Harmonic for number of segregating sites * 2 for diploidy
    
    if HARMONIC_LOCAL==0:
        WTRSNS_E_LOCAL=0
    if HARMONIC_LOCAL!=0:
        WTRSNS_E_LOCAL = SEG_SITES_LOCAL/HARMONIC_LOCAL
    
    
    
    
    
    
    
    
    ###############
    ### Calc Tajima's Theta for Sample Set

    AVG_DIF=[]
    TEMP_DIF=[]
    ### Cycle through all pairwise differences for every SNP, len(AVG_DIF) should be == len (all pairwise difs)
    for COMB in range(0,len(TAJIMAS_THETA_LOCAL[0])):
        for SNP in range(0,len(TAJIMAS_THETA_LOCAL)):            
            TEMP_DIF.append(TAJIMAS_THETA_LOCAL[SNP][COMB])
        
        AVG_DIF.append(sum(TEMP_DIF))
        TEMP_DIF=[]
    
    TAJIMAS_THETA_LOCAL=np.mean(AVG_DIF)
    TAJIMAS_THETA_QUARTET.append(TAJIMAS_THETA_LOCAL)













    #### Calc Observed Heterozygosity for sample
    

    for IND in SAMPLING_NUMBER:

        IND_HETEROZ=[]
        
        for GENOTYPES_HERE in GENOTYPES:
            SAMPLE=GENOTYPES_HERE[IND]

            if (SAMPLE.count('0')==1) and (SAMPLE.count('1')==1):
                IND_HETEROZ.append(1)
            if (SAMPLE.count('0')==2) or (SAMPLE.count('1')==2):
                IND_HETEROZ.append(0)
                
                
        OBSERVED_HETEROZ_LOCAL.append(sum(IND_HETEROZ))

    OBSERVED_HETEROZ_LOCAL=sum(OBSERVED_HETEROZ_LOCAL)/len(OBSERVED_HETEROZ_LOCAL)


    WTRSNS_E_QUARTET.append(WTRSNS_E_LOCAL)
    OBSERVED_HETEROZ_QUARTET.append(EXPECTED_HETEROZ_LOCAL)
    EXPECTED_HETEROZ_QUARTET.append(OBSERVED_HETEROZ_LOCAL)

    
    
    
    
    









###### Sampling metrics averages across loop
AT_LEAST_ONE_VARIANT=sum(SAMPLING_SUCCESS_OR_NOT)/len(SAMPLING_SUCCESS_OR_NOT)
ONE_OR_MORE_VARIANT=sum(TOTAL_VARIANT_SUCCESS_OR_NOT)/len(TOTAL_VARIANT_SUCCESS_OR_NOT)
AT_LEAST_ONE_HOMOZ=sum(HOMOZYGOTE_SUCCESS_OR_NOT)/len(HOMOZYGOTE_SUCCESS_OR_NOT)
AT_LEAST_TWO_ALTERNATIVES=sum(TWO_ALTERNATIVES_SUCCESS_OR_NOT)/len(TWO_ALTERNATIVES_SUCCESS_OR_NOT)



###### Diversity Metrics average across quartets

### Print out Wattersons Est, the full list
OUTPUT_WATTERSON=open(F'{SAMPLE_NAME}.wtrrsn','w')
for k in WTRSNS_E_QUARTET:
    OUTPUT_WATTERSON.write(str(k)+'\n')

### Print out Tajimas Theta, the full list
OUTPUT_TAJIMAS_THETA=open(F'{SAMPLE_NAME}.tjmstheta','w')
for j in TAJIMAS_THETA_QUARTET:
    OUTPUT_TAJIMAS_THETA.write(str(j)+'\n')

WTRSNS_E_QUARTET=sum(WTRSNS_E_QUARTET)/len(WTRSNS_E_QUARTET)
OBSERVED_HETEROZ_QUARTET=sum(OBSERVED_HETEROZ_QUARTET)/len(OBSERVED_HETEROZ_QUARTET)
EXPECTED_HETEROZ_QUARTET=sum(EXPECTED_HETEROZ_QUARTET)/len(EXPECTED_HETEROZ_QUARTET)



print(F'Expected Heterozygosity: {EXPECTED_HETEROZ}\nObserved Heterozygosity: {OBSERVED_HETEROZ}\nWattersons Estimator {WTRSNS_E}\n\n\n')
print(F'Average Expected Heterozygosity per Sample Set: {EXPECTED_HETEROZ_QUARTET}\nAverage Observed Heterozygosity per Sample Set: {OBSERVED_HETEROZ_QUARTET}\nAverage Wattersons Estimator per Sample Set {WTRSNS_E_QUARTET}\n\n\n')
print(F'Average Tajimas Theta of Sample Set: {np.mean(TAJIMAS_THETA_QUARTET)}\n')



print(F'Probability of at least one variant: {AT_LEAST_ONE_VARIANT}\n')
print(F'Average number of successes per individual test: {ONE_OR_MORE_VARIANT}\n')
print(F'Probability of at least one homozygous individual: {AT_LEAST_ONE_HOMOZ}\n')
print(F'Probability of at least 2 alternative alleles: {AT_LEAST_TWO_ALTERNATIVES}\n\n\n')


print(F'Total number of Amino Acids: {Number_of_AA}\n')
print(F'Total length of underlying sequence: {Length_of_Sequence}\n')


plt.hist(TOTAL_VARIANT_SUCCESS_OR_NOT)
plt.savefig('Variants_Histogram.pdf',format='pdf')


OUTPUT.write(F'Expected Heterozygosity: {EXPECTED_HETEROZ}\nObserved Heterozygosity: {OBSERVED_HETEROZ}\nWattersons Estimator: {WTRSNS_E}\n')
OUTPUT.write(F'Average Expected Heterozygosity per Sample Set: {EXPECTED_HETEROZ_QUARTET}\nAverage Observed Heterozygosity per Sample Set: {OBSERVED_HETEROZ_QUARTET}\nAverage Wattersons Estimator per Sample Set: {WTRSNS_E_QUARTET}\n')

OUTPUT.write(F'Probability of at least one variant: {AT_LEAST_ONE_VARIANT}\n')
OUTPUT.write(F'Average number of successes per individual test: {ONE_OR_MORE_VARIANT}\n')
OUTPUT.write(F'Probability of at least one homozygous individual: {AT_LEAST_ONE_HOMOZ}\n')
OUTPUT.write(F'Probability of at least 2 alternative alleles: {AT_LEAST_TWO_ALTERNATIVES}\n')


OUTPUT.write(F'Total number of Amino Acids: {Number_of_AA}\n')
OUTPUT.write(F'Total length of underlying sequence: {Length_of_Sequence}\n')
