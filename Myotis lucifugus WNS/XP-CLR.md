```
#!/bin/bash

aklog
source /share/cdfwwildlife/Capel_Dedicated/miniconda/etc/profile.d/conda.sh
conda activate XP-CLR
module load bcftools

# NOTE: index vcf using bcftools prior to executing 

dir="/share/cdfwwildlife/MYLU_NovaSeq"  # root directory
vcf="${dir}/04_GATKvcfs/gatk.snp.super_filtered_autosomes_thin.vcf.gz"  # vcf to be analyzed
pop1="${dir}/Sample_Lists/NYPRE.ind"  # text file with one sample name per line
pop2="${dir}/Sample_Lists/NYPOST.ind"
pop3="${dir}/Sample_Lists/PAPRE.ind"
pop4="${dir}/Sample_Lists/PAPOST.ind"
outdir="${dir}/05_AnalysisOutput/XP-CLR"

c=0
ar1=($pop1 $pop2 $pop3)  # the first n-1 population list variables
ar2=($pop2 $pop3 $pop4)  # the last n-1 population list variables

tabix $vcf

for i in {0..2}
do
    for j in $(seq $c 2)
    do
        bcftools index -s $vcf | cut -f 1 | while read scaff
        do
            p1=${ar1[$i]}
            p1=$(echo ${p1##*/})
            p1=$(echo ${p1%.*})
            p2=${ar2[$i]}
            p2=$(echo ${p2##*/})
            p2=$(echo ${p2%.*})
            echo "RUNNING XP-CLR ON $scaff for $p1 vs $p2..."
            xpclr -F vcf -I $vcf -C $scaff -Sa ${ar1[$i]} -Sb ${ar2[$j]} --size 50000 --step 45000 --minsnps 2 -O $outdir/${scaff}_${p1}_${p2}.out -V 10
        done
    done
    c=$((c+1))
done
```
