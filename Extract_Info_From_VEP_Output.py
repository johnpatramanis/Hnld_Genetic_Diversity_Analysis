import os
import os.path
from os import listdir
from os.path import isfile, join
import sys

################################################################################################################################################################################################################################################################################################################################################
################################################################### Extract and Filter variants from VEP Output


PROTEIN_POSITIONS_COVERAGE=open('Protein_Coverage.txt','r')

VEP_OUTPUT=open(sys.argv[1],'r')


SAMPLE_NAME=sys.argv[1].split('_VEP.VEP')[0]
OUTPUT=open(F'{SAMPLE_NAME}_Processed_Variants.PV','w')


##### Get which parts of the protein are covered by our ancient samples
COVERAGE={}

for line in  PROTEIN_POSITIONS_COVERAGE:
    line=line.strip().split()
    Protein_Name=line[0]
    
    if len(line)>1:
        Positions=line[1].split(',')
    else:
        Positions=[]
    
    COVERAGE[Protein_Name]=[int(x) for x in Positions]



print(COVERAGE)


### Go through VEP output
FILTERED=[]

for LINE in VEP_OUTPUT:
    if LINE[0]=='#'and LINE[1]!='#':
        LABELS=LINE.strip().split()

    if LINE[0]!='#':
        LINE=LINE.strip().split()
        
        Uploaded_variation=LINE[LABELS.index('#Uploaded_variation')]
        Location=LINE[LABELS.index('Location')]
        
        Gene=LINE[LABELS.index('Gene')] ### Ensemble ID
        Feature=LINE[LABELS.index('Feature')] ### Ensemble ID
        
        Feature_type=LINE[LABELS.index('Feature_type')] ## e.g. 'Transcript'
        VARIANT_CLASS=LINE[LABELS.index('VARIANT_CLASS')] ## e.g. SNV or Insertion
        SYMBOL=LINE[LABELS.index('SYMBOL')] ### Which Protein
        Consequence=LINE[LABELS.index('Consequence')] ## Important! can be multiple things including: synonimous_variant, upstream_variant, missense_variant, splice_donor_variant(?), frameshift_variant
        IMPACT=LINE[LABELS.index('IMPACT')] #### HIGH,LOW other?
        CANONICAL=LINE[LABELS.index('CANONICAL')]
        Protein_position=LINE[LABELS.index('Protein_position')]
        
        
        print(Uploaded_variation,Location,Feature_type,VARIANT_CLASS,SYMBOL,Consequence,IMPACT,Protein_position,CANONICAL)
        
        ###
        ### Filter based on criteria here
        ### FILTERED list contains 
        
        if (VARIANT_CLASS=='SNV') and (Consequence=='missense_variant') and (CANONICAL=='YES'):
            FILTERED.append([SYMBOL,Location,Protein_position])
        

print(LABELS)


for SNP in FILTERED:
    Loc=SNP[1]
    Loc=Loc.split(':')
    Loc='\t'.join(Loc)
    
    Prot=SNP[0]
    Prot_Pos=int(SNP[2])
    
    if Prot in COVERAGE.keys():
        ##COVERAGE[Prot] ## all positions covered by all 4 samples
        print(Loc,Prot,Prot_Pos,COVERAGE[Prot])
        if Prot_Pos in COVERAGE[Prot]:
            OUTPUT.write(Loc+'\n')
