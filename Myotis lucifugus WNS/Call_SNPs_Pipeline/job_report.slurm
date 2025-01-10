#!/bin/bash

#SBATCH --job-name=jobrep
#SBATCH --nodes=1

j1=$(sed -n 1p $2/job_IDs.txt)
j2=$(sed -n 2p $2/job_IDs.txt)
j3=$(sed -n 3p $2/job_IDs.txt)
j4=$(sed -n 4p $2/job_IDs.txt)
j5=$(sed -n 5p $2/job_IDs.txt)
j6=$(sed -n 6p $2/job_IDs.txt)
j7=$(sed -n 7p $2/job_IDs.txt)
j8=$(sed -n 8p $2/job_IDs.txt)
j9=$(sed -n 9p $2/job_IDs.txt)
j10=$(sed -n 10p $2/job_IDs.txt)

echo "# JOB SUCCESS REPORT #"

htst=$(ls $1/01_htstream/01_htstream_$j1_*.out | wc -l)
htss=$(cat $1/01_htstream/01_htstream_$j1_*.out | grep success | wc -l)
echo -e "$htst\t$htss" | awk '{if ($2/$1 == 1) print "HTStream jobs (HTS_preproc.slurm) -- 100% of jobs successful"; else if ($2/$1 == 0) print "HTStream jobs (HTS_preproc.slurm) -- FAILED"; else print "HTStream jobs (HTS_preproc.slurm) -- "1-$2/$1*100"% of jobs FAILED; "$2/$1*100"% of jobs successful"}'

hts=$(cat $1/dragmap_hash_$j2.out | grep success)
echo $hts | awk '{if ($0 == "") print "build DRAGMAP hash table (hashDRAGMAP.slurm) -- FAILED"; else print "build DRAGMAP hash table (hashDRAGMAP.slurm) -- successful"}'

dmt=$(ls $1/02_DRAGMAP/02_DRAGMAP_$j3_*.out | wc -l)
dms=$(cat $1/02_DRAGMAP/02_DRAGMAP_$j3_*.out | grep success | wc -l)
echo -e "$dmt\t$dms" | awk '{if ($2/$1 == 1) print "DRAGMAP alignment jobs (alignDRAGMAP.slurm) -- 100% of jobs successful"; else if ($2/$1 == 0) print "DRAGMAP alignment jobs (alignDRAGMAP.slurm) -- FAILED"; else print "DRAGMAP alignment jobs (alignDRAGMAP.slurm) -- "1-$2/$1*100"% of jobs FAILED; "$2/$1*100"% of jobs successful"}'

if [ ${j10} ]
then
    smt=$(ls $1/03_merge/03_merge_$j10_*.out | wc -l)
    sms=$(cat $1/03_merge/03_merge_$j10_*.out | grep success | wc -l)
echo -e "$smt\t$sms" | awk '{if ($2/$1 == 1) print "samtools merge jobs (samtools_merge.slurm) -- 100% of jobs successful"; else if ($2/$1 == 0) print "samtools merge jobs (samtools_merge.slurm) -- FAILED"; else print "samtools merge jobs (samtools_merge.slurm) -- "1-$2/$1*100"% of jobs FAILED; "$2/$1*100"% of jobs successful"}'
fi

gws=$(cat $1/genome_win_$j4.out | grep success)
echo $gws | awk '{if ($0 == "") print "break reference genome into 50kb windows (genome_wins.slurm) -- FAILED"; else print "break reference genome into 50kb windows (genome_wins.slurm) -- successful"}'

ast=$(ls $1/04_stats/04_stats_$j5_*.out | wc -l)
ass=$(cat $1/04_stats/04_stats_$j5_*.out | grep success | wc -l)
echo -e "$ast\t$ass" | awk '{if ($2/$1 == 1) print "alignment statistics jobs (align_stats.slurm) -- 100% of jobs successful"; else if ($2/$1 == 0) print "alignment statistics jobs (align_stats.slurm) -- FAILED"; else print "alignment statistics jobs (align_stats.slurm) -- "1-$2/$1*100"% of jobs FAILED; "$2/$1*100"% of jobs successful"}'

strs=$(cat $1/str_table_$j6.out | grep success)
echo $strs | awk '{if ($0 == "") print "GATK compose STR table (STRtable.slurm) -- FAILED"; else print "GATK compose STR table (STRtable.slurm) -- successful"}'

hct=$(ls $1/05_hapcall/05_hapcall_$j7_*.out | wc -l)
hcs=$(cat $1/05_hapcall/05_hapcall_$j7_*.out | grep success | wc -l)
echo -e "$hct\t$hcs" | awk '{if ($2/$1 == 1) print "GATK HaplotypeCaller jobs (bam_to_gvcf.sbatch) -- 100% of jobs successful"; else if ($2/$1 == 0) print "GATK HaplotypeCaller jobs (bam_to_gvcf.sbatch) -- FAILED"; else print "GATK HaplotypeCaller jobs (bam_to_gvcf.sbatch) -- "1-$2/$1*100"% of jobs FAILED; "$2/$1*100"% of jobs successful"}'

ggt=$(ls $1/06_genotype/06_genotype_$j8_*.out | wc -l)
ggs=$(cat $1/06_genotype/06_genotype_$j8_*.out | grep success | wc -l)
echo -e "$ggt\t$ggs" | awk '{if ($2/$1 == 1) print "GATK GenotypeGVCFs jobs (gvcf_to_vcf_scaff.sbatch) -- 100% of jobs successful"; else if ($2/$1 == 0) print "GATK GenotypeGVCFs jobs (gvcf_to_vcf_scaff.sbatch) -- FAILED"; else print "GATK GenotypeGVCFs jobs (gvcf_to_vcf_scaff.sbatch) -- "1-$2/$1*100"% of jobs FAILED; "$2/$1*100"% of jobs successful"}'

cfs=$(cat $1/07_concat_vcfs_filter_SNPs_$j9.err | grep "ERROR\|Failed\|failed")
echo $cfs | awk '{if ($0 == "") print "finalize VCF job (vcf_scaff_to_snp.vcf.slurm) -- successful"; else print "finalize VCF job (vcf_scaff_to_snp.vcf.slurm) -- FAILED"}'

if [ $3 == "T" ]
then
    rm -r $2
fi
