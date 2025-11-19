# Code for Genetic Diversity Sampling

Here you will find a step-bystep description on how to re-run the 'genetic variation' analysis, as described in Madupe et al 2025.
The results and details of the analysis are described in the supplementary of the same paper.

<br/><br/>

## Requirements
A computer, personal or server, that is running Linux with python3 with biopython, bcftools, vcftools and ensembl's VEP installed.
All required tools can be installed through conda, and the exact environment used is also provided as .yml file (see below).

<br/><br/>
<br/><br/>
<br/><br/>

## Download and Instalation

The first thing we need to do is to Clone this repository to your local machine (your coputer or server). If you ahve Git installed:


```bash
git clone https://github.com/johnpatramanis/Code_for_Genetic_Diversity_Sampling.git
```

Enter the repo and install the required conda environment (Requires conda to be installed: ) which contains all of the necessary prerequisites:

```bash
cd Code_for_Genetic_Diversity_Sampling
conda env create -f Conda_Env_HmNld.yml
```

Once the installation is complete activate the environment:

```bash
conda activate HNaledi
```

<br/><br/>
<br/><br/>
<br/><br/>

## File preparation and set up


The we need some VCF files for this analysis. If you want to run the analysis for the same dataset as we did in our manuscript, you must use the VCF files from here:
https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg
You only need to download the relevant chromosomes. In our case these would be 3,4,1-,11,17,X and Y.

Once you've downloaded them you need to merge them together using:

```bash
ls gnomad.genomes.v3.1.2.hgdp_tgp.chr*.vcf.bgz > GNOMAD_CHROMOSOME_FILES
bcftools concat -f GNOMAD_CHROMOSOME_FILES -o gnomad.genomes_3_4_10_11_17_X_Y.vcf.bgz
bcftools index gnomad.genomes_3_4_10_11_17_X_Y.vcf.bgz
```

Finally the correct version of the VEP cache needs to be downloaded and placed in a folder path "VEP_Cache/homo_sapiens/" within the main directory of this workflow.
Download Cache from: https://ftp.ensembl.org/pub/release-109/variation/vep/homo_sapiens_vep_109_GRCh38.tar.gz using:

```bash
mkdir -p VEP_Cache/hmo_sapiens
cd VEP_Cache/hmo_sapiens
wget https://ftp.ensembl.org/pub/release-108/variation/vep/homo_sapiens_vep_108_GRCh38.tar.gz
tar –xvzf  homo_sapiens_vep_108_GRCh38.tar.gz
```

Now all we need to do is set up which VCF sample to use by adding it in this file: Vcf_File_List.txt
Simply add 'gnomad.genomes_3_4_10_11_17_X_Y.vcf.bgz' in the first line of the file
If you want any additional VCF files to go through the analysis, add the as new lines below that one.

<br/><br/>
<br/><br/>
<br/><br/>

## Execution of the workflow

Finally type:
```bash
snakemake -j4 -F
```

And the analysis will take place!
