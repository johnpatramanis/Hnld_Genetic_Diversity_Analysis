###### CONDA ENVIRONMENT: conda install -c conda-forge -c bioconda  openssl=1.1 samtools bcftools biopython snakemake bioconductor-genomeinfodbdata ensembl-vep=108.0
###### Requires GRChr38 cache for VEP version 108


import os
import os.path
from os import listdir
from os.path import isfile, join

################################################################################################################################################################################################################################################################################################################################################
################################################################### USEFULL FUNCTIONS




################################################################################################################ SET UP ################################################################################################################################################################################################################################

SAMPLES=[]

SAMPLE_LIST_FILE=open('Vcf_File_List.txt','r')
# GENE_LOCATIONS_FILE=open('Gene_loc_pure.txt','r')


for line in SAMPLE_LIST_FILE:
    line=line.strip().split('.vcf')[0]
    print(line)
    SAMPLES.append(line)










############################################################################################################################################################################################################################################################

####### FINAL SAMPLE SET UP
SAMPLES=list(set(SAMPLES))
SAMPLES.sort()


###### One Rule to Rule them ALL

rule all:
    input:
        'Gene_locs_pure.txt', ### A file with the location of each gene of interest
        'Protein_Coverage.txt', ### A file with all sites covered by the 4 Paranthropus
        expand("{sample}_Gene_Filtered.vcf", sample=SAMPLES), ############ Initial Filter VCF
        expand("{sample}_VEP.VEP", sample=SAMPLES), ######### Get VEP for each Filtered VCF
        expand("{sample}_Processed_Variants.PV", sample=SAMPLES), ######### Second Round of Filtering, keep only SNPs that have an effect on the CANONICAL isoform of the protein and within the coverage of the ancient samples
        expand("{sample}_Second_Filtering.vcf", sample=SAMPLES), ######### Second Round of Filtering, done on the VCF file
        expand("{sample}.gntp", sample=SAMPLES), ######### Output Genotypes for filtered SNPs
        expand("{sample}.metric", sample=SAMPLES), ######### Output Metric for filtered SNPs
        expand("{sample}.windowed.pi",sample=SAMPLES), ###### Output of Genome-wide PIs for VCF, windows of 100000
        expand("{sample}.het",sample=SAMPLES) ###### Output of Genome-wide Heterozygosity for VCF














######################################################################################################################
#### Step 0 - Index VCF if .csi doesn't exits


#### Step 0 #### Index input VCF file
rule Index_VCF_File:
    input:
        VCF_FILE='{sample}.vcf.gz' ############ Will need change name here to work with unzipped VCFs
    output:
        INDEXED_VCF_FILE='{sample}.vcf.gz.csi'
    threads:8
    run:
        shell(F"bcftools index {input.VCF_FILE} -f --threads {threads}")

#### Step 0.5 ### Generate amino acid range, covered by at least X samples (check Get_Protein_Coverage.py to see which and how many samples are considered)
rule Get_Protein_Coverage:
    input:
        VCF_FILES_PREPED=expand('{sample}.vcf.gz.csi',sample=SAMPLES)
    output:
        Protein_Coverage_File='Protein_Coverage.txt'
    threads:1
    run:
        shell(F"python3 Get_Protein_Coverage.py")









######################################################################################################################
#### Step 1 - Filter VCF for SNPs within the boundaries of Genes of interest



#### Step 1 ### 
rule Isolate_Genes_From_Each_VCF_File:
    input:
        GENE_LOCATIONS='Gene_locs_pure.txt',
        INDEXED_VCF_FILE='{sample}.vcf.gz.csi',
        VCF_FILE='{sample}.vcf.gz' ############ Will need change name here to work with unzipped VCFs
    output:
        GENE_FILTERED_VCF='{sample}_Gene_Filtered.vcf'
    threads:8
    run:
        shell(F"bcftools view {input.VCF_FILE} -R {input.GENE_LOCATIONS} --threads {threads} -O v -o {output.GENE_FILTERED_VCF}")




#### Step 1.5 ### 
rule Get_GenomeWide_PI_from_VCF:
    input:
        INDEXED_VCF_FILE='{sample}.vcf.gz.csi',
        VCF_FILE='{sample}.vcf.gz' ############ Will need change name here to work with unzipped VCFs
    output:
        GENE_FILTERED_VCF='{sample}.windowed.pi'
    threads:8
    run:
        file_name=wildcards.sample
        shell(F"vcftools --gzvcf {input.VCF_FILE} --window-pi 100000 --out {file_name}")



#### Step 1.75 ### 
rule Get_GenomeWide_heterozyg_from_VCF:
    input:
        INDEXED_VCF_FILE='{sample}.vcf.gz.csi',
        VCF_FILE='{sample}.vcf.gz' ############ Will need change name here to work with unzipped VCFs
    output:
        GENE_FILTERED_VCF='{sample}.het'
    threads:8
    run:
        file_name=wildcards.sample
        shell(F"vcftools --gzvcf {input.VCF_FILE} --het --out {file_name}")






######################################################################################################################
#### Step 2 - (Second) Filter VCF for SNPs within the boundaries of Genes of interest



#### Step 2 - Run VEP on Filtered VCF
rule RunVEP_on_Filtered_VCF_File:
    input:
        GENE_LOCATIONS='Gene_locs_pure.txt',
        GENE_FILTERED_VCF='{sample}_Gene_Filtered.vcf'
    output:
        VEP_OUTPUT='{sample}_VEP.VEP'
    threads:8
    run:
        shell(F"vep --i {input.GENE_FILTERED_VCF} --tab --species homo_sapiens --offline --dir_cache VEP_Cache/ --output_file {output.VEP_OUTPUT} --force_overwrite --everything")









######################################################################################################################
#### Step 3 - Filter VEP variants for SNPs that: - are within the boundaries of the proteins segments covered by the ancient samples.
####                                             - Occur within exons
####                                             - Are missense variants
####                                             - Have an effect on the Canonical isoform




#### Step 3 - Filter VEP output
rule Second_Round_of_Filtering_Output_Positions:
    input:
        VEP_OUTPUT='{sample}_VEP.VEP',
        PROTEIN_COVERAGE_FILE='Protein_Coverage.txt'
    output:
        PROCESSED_VARIANTS='{sample}_Processed_Variants.PV'
    threads:1
    run:
        shell(F"python3 Extract_Info_From_VEP_Output.py {input.VEP_OUTPUT}")









######################################################################################################################
#### Step 4 - Prepare for second filtering



#### Step 4 - Prepare for second filtering
rule Prepare_for_second_filtering:
    input:
        PROCESSED_VARIANTS='{sample}_Processed_Variants.PV',
        GENE_FILTERED_VCF='{sample}_Gene_Filtered.vcf'
    output:
        GENE_FILTERED_GZVCF='{sample}_Gene_Filtered.vcf.gz',
    threads:8
    run:
        shell(F"bgzip -i -k -f --threads {threads} {input.GENE_FILTERED_VCF}")
        shell(F"bcftools index {output.GENE_FILTERED_GZVCF} -f --threads {threads}")






######################################################################################################################
#### Step 5 - Second Filtering of variants 



#### Step 5 - Filter VEP selected variants
rule Filter_Again_Using_VEP_Variants:
    input:
        PROCESSED_VARIANTS='{sample}_Processed_Variants.PV',
        GENE_FILTERED_GZVCF='{sample}_Gene_Filtered.vcf.gz',
    output:
        SECOND_FILTERING_VCF='{sample}_Second_Filtering.vcf'
    threads:1
    run:
        shell(F"bcftools view {input.GENE_FILTERED_GZVCF} -R {input.PROCESSED_VARIANTS} --threads {threads} -O v -o {output.SECOND_FILTERING_VCF}")









######################################################################################################################
#### Step 6 - Output Genotypes from VCF to a file





#### Step 6 - Output Genotypes from VCF to a file
rule Get_Genotypes:
    input:
        SECOND_FILTERING_VCF='{sample}_Second_Filtering.vcf'
    output:
        GENOTYPES='{sample}.gntp',
        SNP_LOCATIONS='{sample}.lctns'
    threads:1
    run:
        shell(F"bcftools query -f '[%GT\t]\n' {input.SECOND_FILTERING_VCF} > {output.GENOTYPES}")
        shell(F"bcftools query -f '%CHROM %POS %ID %REF %ALT\n' {input.SECOND_FILTERING_VCF} > {output.SNP_LOCATIONS}")






######################################################################################################################
#### Step 7 - Load genotypes into python and do random non-replacement sampling





#### Step 7- Get Metric from Genotypes
rule Get_Metric_From_Genotypes:
    input:
        GENOTYPES='{sample}.gntp'
    output:
        Metric='{sample}.metric'
    threads:1
    run:
        shell(F"python3 Do_Random_Sampling.py {input.GENOTYPES}")

