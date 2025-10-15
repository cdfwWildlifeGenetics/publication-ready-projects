# SNP dataset (RADseq)

## 1. Ipyrad
**load IPYRAD**
```
source /share/cdfwwildlife/hallas_dedicated/Miniconda/etc/profile.d/conda.sh

conda create -n py36
conda activate py36
conda install ipyrad -c conda-forge -c bioconda
```
___

**Parameters**
```
------- ipyrad params file (v.0.9.54)-------------------------------------------
nurtiaAug22                    ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
/data/gpfs/assoc/denovo/jhallas/nutria/ipyrad ## [1] [project_dir]: Project dir (made in curdir if not present)
                               ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                               ## [3] [barcodes_path]: Location of barcodes file
/data/gpfs/assoc/denovo/jhallas/nutria/fastq/*fastq.gz                               ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
ddrad                          ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
44                             ## [21] [min_samples_locus]: Min # samples per locus for output
0.2                            ## [22] [max_SNPs_locus]: Max # SNPs per locus
8                              ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
*                              ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file
                               ## [29] [reference_as_filter]: Reads mapped to this reference are removed in step 3
```
___

**Slurm Submission**
```
#!/bin/bash

#SBATCH --job-name=ipyrad_Aug22
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joshua.hallas@gmail.com
#SBATCH --error=/data/gpfs/assoc/denovo/jhallas/nutria/scripts/out_ipyrad/%x_%A.err
#SBATCH --output=/data/gpfs/assoc/denovo/jhallas/nutria/scripts/out_ipyrad/%x_%A.out

####### LOAD SOFTWARE #######
source activate py36

####### FILE DIRECTORY PATHS #######
PARAMS=/data/gpfs/assoc/denovo/jhallas/nutria/scripts/params-nurtiaAug8.txt 

####### MANAGMENT #######
echo ${HOSTNAME}
start_time=$(date +%s)

###### COMMAND ######
ipyrad -p ${PARAMS} -s 1234567 -c ${SLURM_CPUS_PER_TASK}

####### DURATION CALCULATION ######
end_time=$(date +%s)
elapsed=$((end_time - start_time))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
```
___

**filter stats**
```
                            total_filters  applied_order  retained_loci
total_prefiltered_loci                  0              0         335941
filtered_by_rm_duplicates            3424           3424         332517
filtered_by_max_indels                151            151         332366
filtered_by_max_SNPs                    2              0         332366
filtered_by_max_shared_het            512            510         331856
filtered_by_min_sample             274766         274766          57090
total_filtered_loci                278855         278851          57090


## Alignment matrix statistics:
snps matrix size: (57, 85382), 18.68% missing sites.
sequence matrix size: (57, 7274362), 17.26% missing sites.

```
___

## 2. RAxML
**load RAxML**
```
source /share/cdfwwildlife/hallas_dedicated/Miniconda/etc/profile.d/conda.sh

conda activate py36
conda install -c bioconda raxml
```

**Submission File**
```
#!/bin/bash

#SBATCH --job-name=raxml_Aug23
#SBATCH --account=cpu-s1-bionres-0
#SBATCH --partition=cpu-s1-bionres-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joshua.hallas@gmail.com
#SBATCH --error=/data/gpfs/assoc/denovo/jhallas/nutria/scripts/out_raxml/%x_%A.err
#SBATCH --output=/data/gpfs/assoc/denovo/jhallas/nutria/scripts/out_raxml/%x_%A.out

####### LOAD SOFTWARE #######
source activate py36

####### FILE DIRECTORY PATHS #######
PHYFILE=/data/gpfs/assoc/denovo/jhallas/nutria/ipyrad/nurtiaAug22_outfiles/nurtiaAug22.phy
OUTGROUP=NuMD17
OUTDIR=/data/gpfs/assoc/denovo/jhallas/nutria/raxml
OUTFILE=nutria_23aug.raxml.txt

####### MAKE DIRECTORY #######
mkdir -p ${OUTDIR}

cd ${OUTDIR}

####### MANAGMENT #######
echo ${HOSTNAME}
start_time=$(date +%s)

###### COMMAND ######
raxmlHPC-PTHREADS -f a -m GTRGAMMA -# autoMRE -x 12345 -p 12345 -n ${OUTFILE} -s ${PHYFILE} -o ${OUTGROUP} -T ${SLURM_CPUS_PER_TASK}

####### DURATION CALCULATION ######
end_time=$(date +%s)
elapsed=$((end_time - start_time))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
