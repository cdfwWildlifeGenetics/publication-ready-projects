# Run rehh on polarized and BEAGLE-phased data
Scripts:
- **polarize_alleles.slurm** - polarize alleles based on PRE major allele
- **runREHH.slurm** - masster script; submit beagle & rehh scripts with dependencies to slurm
- **runBEAGLE.sbatch** - phase genotypes with beagle
- **REHH_iHH.R** - calculate iHH per population per scaffold
- **REHH_iHH.sbatch** - execute **REHH_iHH.R** as an array of jobs
- **REHH_iHS.R** - combine iHH across scaffolds for each population, calculate iHS, and find per pop candidate regions
- **REHH_iHS.sbatch** - execute **REHH_iHS.R**
- **REHH_Rsb.R** - read in iHS for each population and calculate pariwise Rsb & XP-EHH
- **REHH_Rsb.sbatch** - execute **REHH_Rsb.R**
- **REHH_Rsb_filt.sbatch** - filter 10Kb smoothed Rsb & XP-CLR output to 10% overlapping windows
- **REHH_Rsb_outliers.R** - find candidate regions for Rsb & XP-CLR

## Polarize alleles based on PRE major allele - **polarize_alleles.slurm**
```bash
#!/bin/bash

#SBATCH --job-name=polarize
#SBATCH --nodes=1
#SBATCH -t 10-00:00:00
#SBATCH --mem=50GB
#SBATCH --output=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/polarize_%A.out
#SBATCH --error=/share/cdfwwildlife/MYLU_NovaSeq/SCRIPTS/slurmout/polarize_%A.err

start=`date +%s`
echo $HOSTNAME
dir="/share/cdfwwildlife/MYLU_NovaSeq"
cnd="/share/cdfwwildlife/Capel_Dedicated/miniconda"

aklog
source $cnd/etc/profile.d/conda.sh
conda activate osjdk17
module load bcftools
module load vcftools

## find major allele for PRE pops
bcftools view -S $dir/Sample_Lists/PRE.ind -O z gatk.snp.qual_hard_filtered_autosomes.vcf.gz | \
    vcftools --gzvcf - --freq --out $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE
tail -n +2 $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE.frq | \
    tr ":" "\t" | cut -f 1,2,5- | \
    awk '{if ($4 < 0.5) print $1"\t"$2"\t"$5"\t"$3; else print $1"\t"$2"\t"$3"\t"$5}' | \
    awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | bgzip > $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz
tabix $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz

## polarize REF allele in VCF
echo "Adding AA INFO field to vcf..."
bcftools view $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.vcf.gz | \
    vcf-annotate -a $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz \
                 -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
                 -c CHROM,FROM,TO,INFO/AA | \
    bcftools view -O z -o $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.vcf.gz

echo "Replacing REF and ALT alleles with AA..."
$cnd/envs/osjdk17/bin/java -jar $dir/SCRIPTS/JVARKIT/jvarkit.jar vcffilterjdk \
                           -f $dir/SCRIPTS/test_ancestral_allele/script.js $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.vcf.gz | \
    bgzip > $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.polarized.vcf.gz
rm $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.vcf.gz

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
```

## Phase genotypes by scaffold & population - **runBEAGLE.sbatch**
Some notes on running beagle v5.4:
- Increasing the number of iterations and the number of phase states increases accuracy, but also memory usage
- Reducing the window size reduces the memory usage - reduce to be able to increase # iterations & phase states
- Using more threads will decrease run time, but can also cause a Java memory usage error (heap space runs out) - limit # threads
```bash
#SBATCH --job-name=beagle
#SBATCH --mem=10G
#SBATCH --array=1-88
#SBATCH --ntasks-per-node=8

dir="/share/cdfwwildlife/MYLU_NovaSeq"
ar=`sed "${SLURM_ARRAY_TASK_ID}q;d" $dir/05_AnalysisOutput/beagle/pop_scaff.tsv`
pop=$(echo $ar | awk '{print $1}')
scaff=$(echo $ar | awk '{print $2}')

module load bcftools/1.9
bcftools view -r $scaff -S $dir/Sample_Lists/${pop}.ind $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.vcf.gz | \
  perl -pe "s/\s\.:/\t.\/.:/g" | \
  bcftools view -O z -o $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.${pop}_${scaff}.vcf.gz

module load beagle/5.4
beagle gt=$dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.${pop}_${scaff}.vcf.gz \
       out=$dir/05_AnalysisOutput/beagle/${pop}_${scaff} \
       impute=false \
       window=0.5 \
       overlap=0.025 \
       iterations=40 \
       phase-states=1000 \
       nthreads=16
rm $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.${pop}_${scaff}.vcf.gz

module load vcftools
bcftools view $ind_$scaff.polarized.vcf.gz | \
  bcftools +fill-tags | \
  vcf-annotate -a $dir/04_GATKvcfs/gatk.snp.qual_hard_filtered_autosomes.PRE_AA.bed.gz \
               -d key=INFO,ID=AA,Number=1,Type=String,Description='Ancestral Allele' \
               -c CHROM,-,FROM,INFO/AA | \
    bcftools view -O z -o $ind_$scaff.polarized.AA.phased.vcf.gz
rm $pop_$scaff.polarized.vcf.gz
```
## Caclulate iHS, Rsb, & XP-EHH
### Calculate iHH values for individual scaffolds for each population - **REHH_iHH.R**
```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

.libPaths("/share/cdfwwildlife/Capel_Dedicated/R_lib/4.3.1")
library(rehh)
library(stringr)

# Read in data and generate a master list of dataframes containing per-pop iHH, iES, and inES for given scaff
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/beagle")
file <- list.files(path = ".", pattern = paste0("^",args[2],"_",args[1],".polarized.AA.phased.vcf.gz$"))
chr_num <- str_extract(string = args[1], pattern = "[0-9]+")
hh <- data2haplohh(hap_file = file,
                   polarize_vcf = T,
                   vcf_reader = "data.table",
                   min_perc_geno.mrk = 80)
scan.chrm <- scan_hh(hh,
                     maxgap = 1056,  # calculated as mean+5SD of distribution of all SNP distances
                     discard_integration_at_border = T)
scan.chrm$CHR <- chr_num

# Write output to file
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHH_per_pop_scaff")
write.table(scan.chrm, file = paste0(args[2],"_",args[1],".scan.chrm.tsv"),
            quote = F, sep = "\t", col.names = T, row.names = F)
```

### Combine iHH across scaffolds for each population and calculate iHS - **REHH_iHS.R**
```R
#!/usr/bin/env Rscript
.libPaths("/share/cdfwwildlife/Capel_Dedicated/R_lib/4.3.1")

library(rehh)
library(stringr)

pops <- c("NYPRE", "NYPOST", "PAPRE", "PAPOST", "PRE", "POST")
scan.chrms <- list()
wgscan.iHH <- list()
wgscan.iHH.pops <- list()

setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHH_per_pop_scaff")

# Read in iHH data
cat("Reading in data...\n")
for (i in 1:length(pops)) {
  print(pops[i])
  files <- list.files(".", pattern = paste0("^",pops[i]))
  for (j in 1:length(files)) {
  scaf <- sub(".scan.chrm.tsv", "", sub(paste0(pops[i],"_"), "", files[j]))
  print(scaf)
    scan.chrms[[j]] <- read.csv(files[j], header = T, sep = "\t")
    if (j == 1) {
      wgscan.iHH <- scan.chrms[[j]]
    } else {
      wgscan.iHH <- rbind(wgscan.iHH, scan.chrms[[j]])
    }
  }
  wgscan.iHH.pops[[i]] <- wgscan.iHH
}

# Calculate individual population iHS
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHS_per_pop")
wgscan.iHS.pops <- list()
cr.iHS.pops <- list()
cat("Calculating iHS...\n")
for (i in 1:length(wgscan.iHH.pops)) {
  print(pops[i])
  wgscan.iHS.pops[[i]] <- ihh2ihs(wgscan.iHH.pops[[i]],
                                  freqbin = 0.01, # or 92
                                  min_maf = 0.03333333)
  write.table(wgscan.iHS.pops[[i]], file = paste(pops[i],"_wgscan_iHS.tsv",sep=""), quote = F,
              sep = "\t", col.names = T, row.names = F)
  cr.iHS.pops[[i]] <- calc_candidate_regions(wgscan.iHS.pops[[i]],
                                             threshold = 3,
                                             pval = T,
                                             window_size = 10000,
                                             overlap = 1000,
                                             min_n_extr_mrk = 2)
  write.table(cr.iHS.pops[[i]], file = paste(pops[i],"_cr_iHS.tsv",sep=""), quote = F,
              sep = "\t", col.names = T, row.names = F)
}

cat("Generating iHS plots...\n")
chr.names <- unique(as.numeric(cr.iHS.pops[[1]]$CHR))
pdf("iHS_diagnostic_plots.pdf")
for (i in 1:length(wgscan.iHS.pops)) {
  print(freqbinplot(wgscan.iHS.pops[[i]],
              main = paste("uniHS w/i freq. bins", pops[i])))
  print(distribplot(wgscan.iHS.pops[[i]]$ihs$IHS,
              xlab = "iHS",
              main = paste("Genome-wide distribution", pops[i])))
  print(distribplot(wgscan.iHS.pops[[i]]$ihs$IHS,
              qqplot = T,
              xlab = "iHS"))
  print(manhattanplot(wgscan.iHS.pops[[i]],
                pval = T,
                threshold = 3,  # either 3 or 4* to be comparable to 5SD
                chr.name = chr.names,
                cr = cr.iHS.pops[[i]],
                main = paste("iHS w/p-value cutoff & candidates", pops[i])))
}
dev.off()
```

### Combine iHH across scaffolds for each population and calculate Rsb & XP-EHH - **REHH_Rsb.R**
```R
#!/usr/bin/env Rscript
.libPaths("/share/cdfwwildlife/Capel_Dedicated/R_lib/4.3.1")

library(rehh)
library(stringr)

setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh")
pops <- c("NYPRE", "NYPOST", "PAPRE", "PAPOST", "PRE", "POST")
pop.comps <- read.csv("pop_comps.tsv", header = F, sep = "\t")
pop.comps.n <- data.frame(pop1 = c(1,1,3,3,5), pop2 = c(2,4,2,4,6))
scan.chrms <- list()
wgscan.iHH <- list()
scan.pops <- list()

# Read in iHH data
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/iHH_per_pop_scaff")
cat("Reading in data...\n")
for (i in 1:length(pops)) {
  print(pops[i])
  files <- list.files(".", pattern = paste0("^",pops[i]))
  for (j in 1:length(files)) {
    scaf <- sub(".scan.chrm.tsv", "", sub(paste0(pops[i],"_"), "", files[j]))
    print(scaf)
    scan.chrms[[j]] <- read.csv(files[j], header = T, sep = "\t")
  }
  scan.pops[[i]] <- scan.chrms
}

# Calculate pairwise Rsb
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/Rsb_per_comp")
cat("Calculating pairwise Rsb...\n")
for (i in 1:nrow(pop.comps)) {
  print(paste(pop.comps[i,1], pop.comps[i,2]))
  for (j in 1:length(scan.pops[[i]])) {  # iterate throug scaffs
    print(paste("SUPER__",scan.pops[[i]][[j]][1,1],sep=""))
    Rsb.chrom.left <- ines2rsb(scan_pop1 = scan.pops[[pop.comps.n[i,1]]][[j]],
                               scan_pop2 = scan.pops[[pop.comps.n[i,2]]][[j]],
                               popname1 = pop.comps[i,1],
                               popname2 = pop.comps[i,2],
                               p.side = "left")
    if (j == 1) {
      wgscan.Rsb.left <- Rsb.chrom.left
    } else {
      wgscan.Rsb.left <- rbind(wgscan.Rsb.left, Rsb.chrom.left)
    }
  }
  cat("  Writing file...\n")
  write.table(wgscan.Rsb.left, file = paste("Rsb_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	       quote = F, sep = "\t", col.names = T, row.names = F)

  cat("  Calculating candidate regions...\n")
  cr.Rsb.left <- calc_candidate_regions(wgscan.Rsb.left,
				threshold = -100,
				window_size = 10000,
				overlap = 1000,
				join_neighbors = F,
				min_n_extr_mrk = 0)
  cr.Rsb.left$Pop1 <- pop.comps[i,1]
  cr.Rsb.left$Pop2 <- pop.comps[i,2]
  cat("  Writing file...\n")
  write.table(cr.Rsb.left, file = paste("Rsb_10Kbwin_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	      quote = F, sep = "\t", col.names = T, row.names = F)
}

# Calculate pairwise XP-EHH
setwd("/share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/XP-EHH_per_comp")
cat("Calculating pairwise XP-EHH...\n")
for (i in 1:nrow(pop.comps)) {
  print(paste(pop.comps[i,1], pop.comps[i,2]))
  for (j in 1:length(scan.pops[[i]])) {  # iterate throug scaffs
    print(paste("SUPER__",scan.pops[[i]][[j]][1,1],sep=""))
    xpehh.chrom.left <- ies2xpehh(scan_pop1 = scan.pops[[pop.comps.n[i,1]]][[j]],
                               scan_pop2 = scan.pops[[pop.comps.n[i,2]]][[j]],
                               popname1 = pop.comps[i,1],
                               popname2 = pop.comps[i,2],
                               p.side = "left")
    if (j == 1) {
      wgscan.xpehh.left <- xpehh.chrom.left
    } else {
      wgscan.xpehh.left <- rbind(wgscan.xpehh.left, xpehh.chrom.left)
    }
  }
  cat("  Writing file...\n")
  write.table(wgscan.xpehh.left, file = paste("XP-EHH_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	      quote = F, sep = "\t", col.names = T, row.names = F)

  cat("  Calculating candidate regions...\n")
  cr.xpehh.left <- calc_candidate_regions(wgscan.xpehh.left,
					  threshold = -100,
					  window_size = 10000,
					  overlap = 1000,
					  join_neighbors = F,
					  min_n_extr_mrk = 0)
  cr.xpehh.left$Pop1 <- pop.comps[i,1]
  cr.xpehh.left$Pop2 <- pop.comps[i,2]
  cat("  Writing file...\n")
  write.table(cr.xpehh.left, file = paste("XP-EHH_10Kbwin_",pop.comps[i,1],"-",pop.comps[i,2],".tsv", sep = ""),
   	      quote = F, sep = "\t", col.names = T, row.names = F)
}
```

### Subset smoothed Rsb & XP-CLR down to 10% overlap (90% overlap generated) - **REHH_Rsb_filt.sbatch**
```shell
cd /share/cdfwwildlife/MYLU_NovaSeq/05_AnalysisOutput/rehh/Rsb_per_comp
for file in `ls *10Kb*`
do
    echo -e "Filtering $file ..."
    pre=$(echo $file | sed -E 's/.tsv//')
    head -n 1 $file > ${pre}_10perc_overlap.tsv
    tail -n +2 $file | cut -f 1 | uniq | sort -g | while read chrom
    do
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | head -n 1
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | tail -n +2 | awk 'NR % 9 == 0'
    done >> ${pre}_10perc_overlap.tsv
done

cd ../XP-EHH_per_comp
for file in `ls *10Kb*`
do
    echo -e "Filtering $file ..."
    pre=$(echo $file | sed -E 's/.tsv//')
    head -n 1 $file > ${pre}_10perc_overlap.tsv
    tail -n +2 $file | cut -f 1 | uniq | sort -g | while read chrom
    do
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | head -n 1
        cat $file | awk -v chrom=$chrom '{if ($1 == chrom) print $0}' | tail -n +2 | awk 'NR % 9 == 0'
    done >> ${pre}_10perc_overlap.tsv
done
```

## Find candidate windows of Rsb & XP-CLR using binned 5SD method - **REHH_Rsb_outliers.R**
```R
setwd("C:/Users/SCapel/OneDrive - California Department of Fish and Wildlife/Research/MYLU WGR/R data/rehh/")

library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(grid)

################
# READ IN DATA #
################

files <- list.files(pattern = "Rsb")
data <- data.frame()
for (i in files) {
  print(i)
  comp <- gsub("_10perc_overlap.tsv", replacement = "", x = i)
  comp <- gsub("Rsb_10Kbwin_", replacement = "", x = comp)
  d <- read.csv(i, header = T, sep = "\t")
  d$comp <- comp
  data <- rbind(data, d)
}

files <- list.files(pattern = "XP-EHH")
data <- data.frame()
for (i in files) {
  print(i)
  comp <- gsub("_10perc_overlap.tsv", replacement = "", x = i)
  comp <- gsub("XP-EHH_10Kbwin_", replacement = "", x = comp)
  d <- read.csv(i, header = T, sep = "\t")
  d$comp <- comp
  data <- rbind(data, d)
}


data$pos_center <- data$START + (data$END - data$START)/2

chroms <- unique(data$CHR)
out <- data.frame()
sum <- 0
# calculate cumulative center position for plotting
for (i in 1:length(chroms)) {
  print(chroms[i])
  if (i == 1) {
    len <- 0
  } else  {
    len <- max(na.omit(out[(data$CHR == chroms[i-1]),13]))
  }
  print(len)
  sum <- len + sum
  print(sum)
  d <- data[(data$CHR == chroms[i]),]
  d$center_cum <- d$pos_center + sum
  out <- na.omit(rbind(out, d))
}
data <- out
comps <- unique(data$comp)

###############
# COUNT PEAKS #
###############
wgmin <- 25
datar <- data[(data$N_MRK >= wgmin),]  # adjust for dataset
keep <- c(expression(data$comp == comps[5]),  ## | data$comp == comps[5]),  # NY-PA & PRE-POST
          expression(data$comp != comps[5]))  # & data$comp != comps[5]))  # All others
keepn <- c("NY-PA & PRE-POST", "All Others")
data$SNPbins <- cut(data$N_MRK, breaks = append(seq(0,400,25),max(data$N_MRK)),
                      labels = c("1-25","26-50","51-75","76-100","101-125","126-150",
                                 "151-175","176-200","201-225","226-250","251-275",
                                 "276-300","301-325","326-350","351-375","376-400",
                                 paste("401-",max(data$N_MRK),sep="")))
bins <- levels(data$SNPbins)

#                                   #
## Create binning diagnostic plots ##
#                                   #
td <- data.frame(Comparisons = keepn, vals = c(1,2))
td$Comparisons <- factor(td$Comparisons, levels = keepn)
p <- ggplot(td, aes(x = vals, fill = Comparisons)) +
  geom_bar() +
  scale_fill_manual(values = c("black","grey40")) +
  theme(legend.position = "top")
l1 <- get_plot_component(p, 'guide-box-top', return_all = T)
hist <- ggplot() +
  geom_bar(data[(eval(keep[2])),], mapping = aes(x = SNPbins), fill = "grey40") +
  geom_bar(data[(eval(keep[1])),], mapping = aes(x = SNPbins), fill = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dp1 <- plot_grid(l1, hist, ncol = 1, rel_heights = c(0.1, 0.9))
p <- ggplot(td, aes(x = vals, y = vals, group = Comparisons, linetype = Comparisons)) +
  geom_line() +
  scale_linetype_manual(values = c("dashed","solid")) +
  theme(legend.position = "top", legend.key=element_rect(fill= NA))
l2 <- get_plot_component(p, 'guide-box-top', return_all = T)
dens <- ggplot() +
  geom_density(data[(eval(keep[1])),], mapping = aes(x = MEAN_MRK, group = SNPbins, color = SNPbins)) +
  geom_density(data[(eval(keep[2])),], mapping = aes(x = MEAN_MRK, group = SNPbins, color = SNPbins), linetype = "dashed") +
  theme_classic()
dp2 <- plot_grid(l2, dens, ncol = 1, rel_heights = c(0.1, 0.9))
pdf(paste("../../../Figures/Rsb_qual_filt_",ds,"_binning_diagnostics.pdf",sep=""), width = 10, height = 5)
plot_grid(dp1, dp2, rel_widths = c(0.47,0.53))
dev.off()
cat(paste("Mean # SNPs/window:\n", 
          keepn[1]," ", round(mean(data[(eval(keep[1])),4]),2),"\n",
          keepn[2]," ", round(mean(data[(eval(keep[2])),4]),2)))

#                                    #
## Report number of outlier windows ##
#                                    #
wgdf <- data.frame()
d5sdb <- data.frame()
d5sdwg <- data.frame()
for (i in 1:length(keep)) {
  print(paste("#####",keepn[i],"#####"))
  dfbin <- data.frame()
  dfpop <- data.frame()
  datr <- data[(eval(keep[i]) & data$N_MRK >= wgmin),]
  t <- nrow(datr)
  v5 <- tail(head(datr[order(datr[,5], decreasing=T),5],nrow(datr)*.05),1)
  nv5 <- nrow(datr[(datr$MEAN_MRK >= v5),])
  v1 <- tail(head(datr[order(datr[,5], decreasing=T),5],nrow(datr)*.01),1)
  nv1 <- nrow(datr[(datr$MEAN_MRK >= v1),])
  sd5 <- mean(datr$MEAN_MRK) + sd(datr$MEAN_MRK)*2.5
  nsd5 <- nrow(datr[(datr$MEAN_MRK >= sd5),])
  wg <- data.frame("SNPs/win" = "Whole genome*", "total # win" = t, "top 5% Rsb cutoff" = round(v5,4),
                   "# win in top 5%" = nv5, "top 1% Rsb cutoff" = round(v1,4), "% win in top 1%" = nv1, 
                   "5 SD Rsb cutoff" = round(sd5,4), "# win > 5 SD" = nsd5)
  compsr <- unique(datr$comp)
  for (j in 1:length(compsr)) {
    if (j == 1) {
    }
    v5n <- nrow(datr[datr$comp == compsr[j] & datr$MEAN_MRK >= v5,])
    v1n <- nrow(datr[datr$comp == compsr[j] & datr$MEAN_MRK >= v1,])
    sd5n <- nrow(datr[datr$comp == compsr[j] & datr$MEAN_MRK >= sd5,])
    wgpop <- data.frame(Comp = compsr[j], "WG top 5%" = v5n, "WG top 1%" = v1n, "WG > 5 SD" = sd5n)
    dfpop <- rbind(dfpop, wgpop)
  }
  datar <- data[(eval(keep[i])),]
  pop_sum <- data.frame()
  data5sd <- data.frame()
  for (j in 1:length(bins)) {
    if (j == 1) {
    }
    dat <- datar[(datar$SNPbins == bins[j]),]
    nr <- nrow(dat)
    m <- mean(dat[,5])
    sd5n <- m + sd(dat[,5])*2.5
    nrsd5n <- nrow(dat[(dat$MEAN_MRK > sd5n),])
    t5 <- round(nr*.05,0)
    t1 <- round(nr*.01,0)
    v5 <- tail(head(dat[order(dat[,5], decreasing=T),5],t5),1)
    v1 <- tail(head(dat[order(dat[,5], decreasing=T),5],t1),1)
    bin <- data.frame("SNPs/win" = bins[j], "total # win" = nr, "top 5% Rsb cutoff" = round(v5,4),
                      "# win in top 5%" = t5, "top 1% Rsb cutoff" = round(v1,4), "% win in top 1%" = t1, 
                      "5 SD Rsb cutoff" = round(sd5n,4), "# win > 5 SD" = nrsd5n)
    dfbin <- rbind(dfbin, bin)
    for (k in 1:length(compsr)) {
      d <- dat[(dat$comp == compsr[k]),]
      cv5 <- nrow(d[(d$MEAN_MRK >= v5),])
      cv1 <- nrow(d[(d$MEAN_MRK >= v1),])
      csd5n <- nrow(d[(d$MEAN_MRK > sd5n),])
      df <- data.frame("comp" = compsr[k], "bin" = bins[j], "top5" = cv5, "top1" = cv1, "5SD" = csd5n)
      pop_sum <- rbind(df, pop_sum)
    }
    if (bins[j] != "1-25") {
      datp <- dat[(dat$MEAN_MRK > sd5n),]
      data5sd <- rbind(datp, data5sd)
    }
  }
  binsum <- data.frame("SNPs/win" = "Sum", "total # win" = sum(dfbin$total...win[2:length(dfbin$total...win)]), "top 5% Rsb cutoff" = "",
                       "# win in top 5%" = sum(dfbin$X..win.in.top.5.[2:length(dfbin$total...win)]), "top 1% Rsb cutoff" = "", 
                       "% win in top 1%" = sum(dfbin$X..win.in.top.1.[2:length(dfbin$total...win)]), "5 SD Rsb cutoff" = "", 
                       "# win > 5 SD" = sum(dfbin$X..win...5.SD[2:length(dfbin$total...win)]))
  dfbin <- rbind(dfbin, binsum)
  dfbin <- rbind(dfbin, wg)
  print(dfbin, row.names = F)
  pop_sum <- data.frame(pop_sum[(pop_sum$bin != "1-25"),] %>% group_by(comp) %>% summarise(top5 = sum(top5), top1 = sum(top1), X5SD = sum(X5SD)))
  names(pop_sum) <- c("Comp", "bin top 5%", "bin top 1%", "bin > 5 SD")
  dfpop <- cbind(dfpop, pop_sum[c(2,3,4)])
  print(dfpop, row.names = F)
  wg <- cbind(data.frame("Dataset" = keepn[i]), wg[1,c(3,5,7)])
  wgdf <- rbind(wgdf, wg)
  d5sdb <- rbind(d5sdb, data5sd)
  d5sdwg <- rbind(d5sdwg, datr[(datr$MEAN_MRK > wg[1,4]),])
}
```
