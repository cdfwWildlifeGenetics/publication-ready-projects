#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --nodes=1

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
lib=$(cat $4 | grep -w ${sample} | head -n 1 | awk '{print $3}')

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"
echo "${lib}"

module load samtools
echo "Merging sample reads across lanes..."
samtools merge -f -o $2/${lib}/${sample}/${sample}.sorted.rg.md.bam $2/${lib}/${sample}/${sample}*$3*.sorted.rg.md.bam
rm $2/${lib}/${sample}/*$3*.bam
samtools index $2/${lib}/${sample}/${sample}.sorted.rg.md.bam

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime seconds"
