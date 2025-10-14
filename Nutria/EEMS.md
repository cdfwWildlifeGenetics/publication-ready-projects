# Estimated Effective Migration Surface
**EEMS [manual](https://github.com/dipetkov/eems)**
___

## Download
```
git clone https://github.com/dipetkov/eems.git
```

___

## Input Files
*NOTES: there are three files to run EEMS. 1. distance matrix; 2. sample coordinates; 3. habitat coordinates.  Double check manual for file structure.*

**Create boundry file**

In arcgis: once polygon is created
+ ArcToolbox > Data Management > Features > Feature Vertices to points 
+ ArcToolbox > Data Management > Features > Add XY Coordinates
+ ArcToolbox > Conversion Tools > Excel > Table to Excel

*NOTES: make sure that extension is .outer*

## Params File
You need to create at least three for each individual run.

*params-chain1.ini*
```
datapath = /share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/eems_all/nutria_all
mcmcpath = /share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/eems_all/nutria_all_jan9_Indiv267-nSites6809-EEMS-nDemes600-chain1
nIndiv = 267
nSites = 6809
nDemes = 600
diploid = true
numMCMCIter = 10000000
numBurnIter = 1000000
numThinIter = 5000
qEffctProposalS2  = 0.020000
qSeedsProposalS2  = 0.080000
mEffctProposalS2  = 1.150000
mrateMuProposalS2 = 0.001000
mSeedsProposalS2  = 0.007499
```

*NOTES: Proposals will need to fine tuned.  Prelim runs are needed before all three chains are run*

## Slurm File

```
#!/bin/bash

#SBATCH --job-name=eems_nutria_all
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500
#SBATCH --time=14-00:00:00
#SBATCH --array=1-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joshua.hallas@gmail.com
#SBATCH --error=/share/cdfwwildlife/Nutria_LandGen/all_lanes/scripts/out_eems/%x_%A_%a.err
#SBATCH --output=/share/cdfwwildlife/Nutria_LandGen/all_lanes/scripts/out_eems/%x_%A_%a.out

####### LOAD SOFTWARE #######
module load eems/c1849ea

####### CREATES JOB ARRAY #######
REP=${SLURM_ARRAY_TASK_ID} # use job array number as replicate number
SEED=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/seed.list)

####### FILE DIRECTORY PATHS #######
PARAMS=/share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/eems_all

####### MANAGMENT #######
start_time=$(date +%s)

###### COMMAND ######
runeems_snps --params ${PARAMS}/params-chain${REP}.ini --seed ${SEED}

####### DURATION CALCULATION ######
end_time=$(date +%s)
elapsed=$((end_time - start_time))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
```
## Visiualization

*Transfer output directories to laptop*
```
scp -r barbera:/share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/nutria_nov27_Indiv84-nSites-EEMS-nDemes600-chain1 /Users/jhallas/Documents/cdfw/nutria/popgen_analysis 
scp -r barbera:/share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/nutria_nov27_Indiv84-nSites-EEMS-nDemes600-chain2 /Users/jhallas/Documents/cdfw/nutria/popgen_analysis 
scp -r barbera:/share/cdfwwildlife/Nutria_LandGen/all_lanes/landgen/eems/nutria_nov27_Indiv84-nSites-EEMS-nDemes600-chain3 /Users/jhallas/Documents/cdfw/nutria/popgen_analysis 
```

In R
```
library("devtools")

install_github("dipetkov/reemsplots2")

library("reemsplots2")
library("ggplot2")      # Modify the default plots
library("rworldmap")    # Add a geographic map
library("broom")        # Required for the map
library("mapproj")      # Change the projection
library("RColorBrewer") # Change the color scheme

####### selects all three chains
mcmcpath_nutria  <- c("/Users/jhallas/Documents/cdfw/nutria/popgen_analysis/nutria_all_jan9_Indiv267-nSites6809-EEMS-nDemes600-chain1",
                      "/Users/jhallas/Documents/cdfw/nutria/popgen_analysis/nutria_all_jan9_Indiv267-nSites6809-EEMS-nDemes600-chain2",
                      "/Users/jhallas/Documents/cdfw/nutria/popgen_analysis/nutria_all_jan9_Indiv267-nSites6809-EEMS-nDemes600-chain3")

####### generates effective migration surfaces
plots_nutria<-make_eems_plots(mcmcpath_nutria, longlat = F, dpi = 300, add_grid = F,
                            col_grid = "#BBBBBB", add_demes = T, col_demes = "#000000",
                            add_outline = T, col_outline = "black", eems_colors = NULL,
                            prob_level = 0.9, m_colscale = NULL, q_colscale = NULL,
                            add_abline = FALSE)

####### make plots
#names(plots_atratus)
plots_nutria$mrates01 #generates migration rates figure
plots_nutria$qrates01 #generates diversity rates figure
plots_nutria$mrates02
plots_nutria$qrates02
plots_nutria$pilogl01

####### save plots
plots_nutria$mrates01 #generates migration rates figure

dev.print(pdf, "/Users/jhallas/Documents/cdfw/nutria/popgen_analysis/figures/eems_all_jan9_m.pdf")

plots_atratus$qrates01 #generates diversity rates figure

dev.print(pdf, "/Users/jhallas/Documents/cdfw/nutria/popgen_analysis/figures/eems_all_jan1_q.pdf")
```
