import os
import os.path
from os import listdir
from os.path import isfile, join
import sys
from Bio import SeqIO



FILES_IN_FOLDER = [f for f in os.listdir('.') if os.path.isfile(f)]
FASTA_FILE_LIST= [f for f in FILES_IN_FOLDER if '.fa' in f]


OUTPUT=open('Protein_Coverage.txt','w')

SAMPLES_TO_INCLUDE=['U.W.101-020','U.W.101-886','U.W.101-511','U.W.110-2','U.W.101-1846','U.W.101-1686','U.W.101-824','U.W.101-516','U.W.102c-589','U.W.101-809','U.W.101-Neo','U.W.101-1396','U.W.101-010','U.W.101-2175','U.W.101-2021','U.W.101-1915','U.W.101-1463','U.W.101-525','U.W.101-583','U.W.101-006']
NUMBER_OF_SAMPLES=len(SAMPLES_TO_INCLUDE)

HIGH_COVERAGE_SAMPLES=[]



####### Get sites covered by each sample
COVERAGE={} ## Dictionary with each sample as key and each protein_name as value (dictionary) each list has the number of sites that are covered

for FILE in FASTA_FILE_LIST:
    fasta_sequences = SeqIO.parse(open(FILE),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        name=name.split('/')[0]
        name=name.split('_')
        sample=name[1]
        protein=name[len(name)-1]
        protein=protein.replace('.fasta','')
        
        if sample in SAMPLES_TO_INCLUDE:
            if sample not in COVERAGE.keys():
                COVERAGE[sample]={}
            
            counter=0          
            for POS in sequence:
                counter+=1
                if ((POS!='?') and (POS!='-') and (POS!='X')):
                    if protein not in COVERAGE[sample].keys():
                        COVERAGE[sample][protein]=[]
                    COVERAGE[sample][protein].append(counter)

            # print(sample,protein)

# for KEY in COVERAGE.keys():
    # print('$$$$$$$$$$')
    # print(KEY)
    # SAMPLE=COVERAGE[KEY]
    # for PROTEIN,SITES in SAMPLE.items():
        # print(PROTEIN,SITES)
        # print('##########\n')








    
TOTAL_COVERAGE={}

for SMPL in COVERAGE.keys():
    for PRTN in COVERAGE[SMPL].keys():
        if PRTN not in TOTAL_COVERAGE.keys():
            TOTAL_COVERAGE[PRTN]=[]
        for PSTN in COVERAGE[SMPL][PRTN]:
            TOTAL_COVERAGE[PRTN].append(PSTN)
            




####### For getting positions covered by ALL 4 P.rob samples
for PRTN in TOTAL_COVERAGE.keys():
    TOTAL_COVERAGE[PRTN]=sorted(TOTAL_COVERAGE[PRTN])
    TOTAL_COVERAGE_UNIQUE=list(set(TOTAL_COVERAGE[PRTN])) ### Sites only once, so we can loop through them

    
    ### Only select positions that are counted 4 times
    TOTAL_COVERAGE[PRTN]=[ str(SITE) for SITE in TOTAL_COVERAGE_UNIQUE if TOTAL_COVERAGE[PRTN].count(SITE)>=NUMBER_OF_SAMPLES ]
    print(PRTN,TOTAL_COVERAGE[PRTN],len(TOTAL_COVERAGE[PRTN]))
    print('#####################\n')
    POSITIONS=','.join(TOTAL_COVERAGE[PRTN])
    OUTPUT.write(F'{PRTN}\t{POSITIONS}\n')







####### For getting positions covered by at least one P.rob sample!
# for PRTN in TOTAL_COVERAGE.keys():
    # TOTAL_COVERAGE[PRTN]=sorted(list(set(TOTAL_COVERAGE[PRTN])))
    # TOTAL_COVERAGE[PRTN]=[str(x) for x in TOTAL_COVERAGE[PRTN]]
    # print(PRTN,TOTAL_COVERAGE[PRTN])
    # POSITIONS=','.join(TOTAL_COVERAGE[PRTN])
    # OUTPUT.write(F'{PRTN}\t{POSITIONS}\n')
    
# print(TOTAL_COVERAGE)