#!/bin/bash

#SBATCH --job-name=align_stats
#SBATCH --nodes=1
#SBATCH --ntasks=7

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
lib=$(cat $6 | grep -w -m 1 ${sample} | awk '{print $3}')
[[ -d $2/${lib}/${sample} ]] || mkdir -p $2/${lib}/${sample}
ref=$3
refpre=$(echo ${ref%.*})

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"
echo "${lib}"

module load samtools
module load bedtools2

echo "Calculating samtools stats:"
echo "  flagstat..."
samtools flagstat $4/${lib}/${sample}/${sample}.sorted.rg.md.bam > $2/${lib}/${sample}/${sample}.dragmap_flagstat.txt
echo "  stats..."
samtools stats -c 1,10000,100 $4/${lib}/${sample}/${sample}.sorted.rg.md.bam > $2/${lib}/${sample}/${sample}.dragmap_stats.txt
echo "  idxstats..."
samtools idxstats $4/${lib}/${sample}/${sample}.sorted.rg.md.bam > $2/${lib}/${sample}/${sample}.dragmap_IDstats.tsv
echo "  depth..."
samtools depth $4/${lib}/${sample}/${sample}.sorted.rg.md.bam | awk '{sum+=$3; cnt++}END{print "no_positions\tmean_depth"; print sum"\t"sum/cnt}' > $2/${lib}/${sample}/${sample}.dragmap_depth.txt
echo "  coverage..."
samtools coverage $4/${lib}/${sample}/${sample}.sorted.rg.md.bam > $2/${lib}/${sample}/${sample}.dragmap_coverage.tsv

echo "Calculating GATK metrics..."
echo "checking mapping quality"
$5/gatk CollectAlignmentSummaryMetrics \
    -R $3 \
    -I $4/${lib}/${sample}/${sample}.sorted.rg.md.bam \
    -O $2/${lib}/${sample}/${sample}.dragmap.gatk.stat.txt

echo "Calculating coverage for 50Kb windows..."
echo -e "scaff\twin_start\twin_end\tmean_depth" > $2/${lib}/${sample}/${sample}.dragmap_50kbwin_meancov.tsv
bedtools coverage -mean -sorted -a ${refpre}_50kb_win.bed -b $4/${lib}/${sample}/${sample}.sorted.rg.md.bam >> $2/${lib}/${sample}/${sample}.dragmap_50kbwin_meancov.tsv

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime seconds"
