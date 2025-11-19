# Code for Genetic Diversity Sampling

Here you will find a step-bystep description on how to re-run the 'genetic variation' analysis, as described in Madupe et al 2025.
The results and details of the analysis are described in the supplementary of the same paper.

<br/><br/>

## Requirements
A computer, personal or server, that is running Linux with python3 with [biopython](https://biopython.org/), [bcftools](https://samtools.github.io/bcftools/bcftools.html), [vcftools](https://vcftools.sourceforge.net/) and ensembl's [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) installed.
All required tools can be installed through [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html), and the exact environment used is also provided as [YAML](https://en.wikipedia.org/wiki/YAML) file (see below).
Additionally around 

<br/><br/>
<br/><br/>
<br/><br/>

## Download and Instalation

The first thing we need to do is to Clone this repository to your local machine (your coputer or server). If you have Git installed simply:

```bash
git clone https://github.com/johnpatramanis/Code_for_Genetic_Diversity_Sampling.git
```


Once the repository is copied, enter the folder and with conda installed and activated install the required **conda environment**. 

```bash
cd Code_for_Genetic_Diversity_Sampling
conda env create -f Conda_Env_HmNld.yml
```
Alternatively, if the above method fails, you can also install the environment manually with the following command:

```
conda create -n HNaledi -c conda-forge -c bioconda  openssl=1.1 samtools bcftools biopython snakemake bioconductor-genomeinfodbdata ensembl-vep=108.0
```


Either of the above commands will create a new conda environment in your computer, named "HNaledi", which contains all of the necessary prerequisites for this analysis. Now activate this environment with:

```bash
conda activate HNaledi
```

<br/><br/>
<br/><br/>
<br/><br/>

## Preparing the 1000 Genomes Data


Now we need some VCF files for this analysis. If you want to run the analysis for the same dataset as we did in our manuscript, you must use the VCF files from the 1000 genomes project. First we need to download the 1000 Genomes VCF file, which contains information on all populations and store it in a folder named **1000_Genomes_Data**:
```
cd 1000_Genomes_Data
bash Download_Data.sh
```

This will download 23 VCF files (and their 23 .tbi file), one for each chromosome. Now we need to merge them together into one:

```
ls *.vcf.gz > To_Merge.txt
bcftools concat -f To_Merge.txt -O z -o 1000_Genomes_Merged_All.vcf.gz
bcftools index 1000_Genomes_Merged_All.vcf.gz
```

Then we want to split this folder into individual-populations folders (1 VCF containing all individuals belonging to 1 population). For this first we need to find how many and which populations there are, and to which population each sample belongs to. Luckily for us here, this has already been done (using a python script named **Extract_Specific_Populations_From_1000_Genomes.py** in combination with a .tsv file with all the necessary information named **igsr_samples.tsv**). So inside the **1000_Genomes_Data** folder, you will find one **XXX_samples.txt** file for each population, which contains the ID of all samples belonging to that population.

```
bash SubSample.sh
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
snakemake -j8 -F
```

And the analysis will take place!


## Plotting the results
