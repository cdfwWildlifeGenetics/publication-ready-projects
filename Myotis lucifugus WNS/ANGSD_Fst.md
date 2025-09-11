# Calculate smoothed F<sub>ST</sub> in 10 Kbp windows with 10% overlap
## Estimate global SFS and calculate thetas
```
#!/bin/bash

#SBATCH --job-name=angsd_theta
#SBATCH --nodes=1
#SBATCH --time=05-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-8
#SBATCH --output=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/ANGSD/%x_%A_%a.out
#SBATCH --error=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/ANGSD/%x_%A_%a.err

dir="/share/cdfwwildlife/MYLU_NovaSeq"
angsd="$dir/05_AnalysisOutput/ANGSD"
vcf="gatk.snp.super_filtered_autosomes"
pop=`sed "${SLURM_ARRAY_TASK_ID}q;d" $dir/Sample_Lists/pops.txt`

aklog
source /share/cdfwwildlife/Capel_Dedicated/miniconda/etc/profile.d/conda.sh
conda activate base
module load bcftools

echo $HOSTNAME
echo $pop

## edit vcf header for ANGSD compatability and reduce to a single population "pop"
INFO='##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"\">'
tail -n +2 $dir/Sample_Lists/Population_Map.txt | cut -f 1,2 | grep $pop | cut -f 1 > $dir/Sample_Lists/$pop.ind
bcftools view -S $dir/Sample_Lists/$pop.ind $dir/04_GATKvcfs/$vcf.vcf.gz | awk -v infohdr=$INFO '/^#CHROM/ {print infohdr"\n"$0} !/^#CHROM/' | bcftools view -Oz -o $angsd/$vcf/$pop.vcf.gz

## estimate global SFS
angsd -vcf-pl $angsd/$vcf/$pop.vcf.gz -doSaf 1 -anc $dir/REF/HiC_Genome/mMyoLuc1.fa -P 1 -out $angsd/$vcf/$pop
realSFS check $angsd/$vcf/$pop.saf.idx

## calculate thetas
realSFS $angsd/$vcf/$pop.saf.idx -fold 1 -P 1 > $angsd/$vcf/$pop.sfs  # calculate folded SFS due to uknown ancestral states
realSFS saf2theta $angsd/$vcf/$pop.saf.idx -sfs $angsd/$vcf/$pop.sfs -outname $angsd/$vcf/$pop
thetaStat print $angsd/$vcf/$pop.thetas.idx  # view index file
thetaStat do_stat $angsd/$vcf/$pop.thetas.idx -win 10000 -step 9000
cat $angsd/$vcf/$pop.thetas.idx.pestPG

(exit) && echo success
```

## Calculate F<sub>ST</sub> for all pairwise comparisons
```
#!/bin/bash

#SBATCH --job-name=angsd_Fst
#SBATCH --nodes=1
#SBATCH --time=14-00:00:00
#SBATCH --mem=20G
#SBATCH --array=1-8
#SBATCH --ntasks-per-node=16
#SBATCH --output=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/ANGSD/%x_%A_%a.out
#SBATCH --error=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/ANGSD/%x_%A_%a.err

dir="/share/cdfwwildlife/MYLU_NovaSeq"
angsd="$dir/05_AnalysisOutput/ANGSD"
vcf="gatk.snp.qual_hard_filtered_autosomes"
pop=`sed "${SLURM_ARRAY_TASK_ID}q;d" $dir/Sample_Lists/pop_comps.tsv`
pop1=$(echo $pop | awk '{print $1}')
pop2=$(echo $pop | awk '{print $2}')

aklog
source /share/cdfwwildlife/Capel_Dedicated/miniconda/etc/profile.d/conda.sh
conda activate base
module load bcftools

echo $HOSTNAME
echo $vcf
echo $pop

realSFS $angsd/$vcf/$pop1.saf.idx $angsd/$vcf/$pop2.saf.idx -fold 0 -P 16 > $angsd/$vcf/Fst/$pop1.$pop2.ml   # calculate 2d sfs for each pair
realSFS fst index $angsd/$vcf/$pop1.saf.idx $angsd/$vcf/$pop2.saf.idx -sfs $angsd/$vcf/Fst/$pop1.$pop2.ml -fstout $angsd/$vcf/Fst/$pop1.$pop2 -P 1  # prepare Fst for easy window analysis
# realSFS fst print $angsd/$vcf/Fst/$pop1.$pop2.fst.idx | wc -l
# realSFS fst stats $angsd/$vcf/Fst/$pop1.$pop2.fst.idx  # get the global estimate
realSFS fst stats2 $angsd/$vcf/Fst/$pop1.$pop2.fst.idx -win 10000 -step 9000 -P16 > $angsd/$vcf/Fst/10Kb/$pop1.$pop2.win.fst
cat $angsd/$vcf/Fst/10Kb/$pop1.$pop2.win.fst | sed -E 's/Nsites/Nsites\tFst/' >temp && mv temp $angsd/$vcf/Fst/10Kb/$pop1.$pop2.win.fst
realSFS fst stats2 $angsd/$vcf/Fst/$pop1.$pop2.fst.idx -win 5000 -step 4500 -P 16 > $angsd/$vcf/Fst/5Kb/$pop1.$pop2.win.fst
cat $angsd/$vcf/Fst/5Kb/$pop1.$pop2.win.fst | sed -E 's/Nsites/Nsites\tFst/' > temp && mv temp $angsd/$vcf/Fst/5Kb/$pop1.$pop2.win.fs

(exit) && echo success
```
