ds <- "10Kb"
wgmin <- 20  # 5Kb = 10; 10Kb = 20
binex <- "1-25"
maxsnps <- 500  # 5Kb = 250; 10Kb = 500

setwd(paste("C:/Users/SCapel/OneDrive - California Department of Fish and Wildlife/Research/MYLU WGR/R data/xplcr/qual filter ",ds,"/",sep = ""))

library(ggplot2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(grid)

################
# READ IN DATA #
################

files <- list.files()
data <- data.frame()
for (i in files) {
  print(i)
  comp <- gsub("SUPER__[0-9]+_", replacement = "", x = i)
  comp <- gsub(".out", replacement = "", x = comp)
  d <- read.csv(i, header = T, sep = "\t")
  d$comp <- comp
  data <- rbind(data, d)
}
data <- na.omit(data)
data$chrom <- as.numeric(gsub("SUPER__", replacement = "", data$chrom))
data <- data[order(data$chrom),]
data$pos_center <- data$pos_start + (data$pos_stop - data$pos_start)/2

chroms <- unique(data$chrom)
out <- data.frame()
sum <- 0
# calculate cumulative center position for plotting
for (i in 1:length(chroms)) {
  print(chroms[i])
  if (i-1 == 0) {
    len <- 0
  } else  {
    len <- max(out[(data$chrom == chroms[i-1]),4])
  }
  print(len)
  sum <- len + sum
  print(sum)
  d <- data[(data$chrom == chroms[i]),]
  d$center_cum <- d$pos_center + sum
  out <- rbind(out, d)
}
data <- out
comps <- unique(data$comp)


###############
# COUNT PEAKS #
###############
unused <- c()
for (i in 1:length(comps)) {
  print(comps[i])
  x <- data[(data$comp == comps[i] & data$nSNPs_avail > maxsnps),11]
  x <- sum(x - maxsnps)
  print(x)
  unused <- append(x, unused)
}
unusedt <- sum(unused)
unusedm1 <- mean(unused[c(1,8)])
unusedm2 <- mean(unused[2:7])

datar <- data[(data$nSNPs >= wgmin),]  # adjust for dataset
keep <- c(expression(data$comp == comps[1] | data$comp == comps[8]),
          expression(data$comp != comps[1] & data$comp != comps[8]))
keepn <- c("NY-PA & PRE-POST", "All Others")
if (ds == "10Kb") {
  data$SNPbins <- cut(data$nSNPs, breaks = append(seq(0,400,25),725), 
                      labels = c("1-25","26-50","51-75","76-100","101-125","126-150",
                                 "151-175","176-200","201-225","226-250","251-275",
                                 "276-300","301-325","326-350","351-375","376-400","401-500"))
} else if (ds == "5Kb") {
  data$SNPbins <- cut(data$nSNPs, breaks = append(seq(0,225,25),420), 
                      labels = c("1-25","26-50","51-75","76-100","101-125","126-150",
                                 "151-175","176-200","201-225","226-250"))
}
bins <- levels(data$SNPbins)

#                                   #
## Create binning diagnostic plots ##
#                                   #
td <- data.frame(Comparisons = keepn, vals = c(1,2))
td$Comparisons <- factor(td$Comparisons, levels = keepn)
l1 <- get_legend(ggplot(td, aes(x = vals, group = Comparisons, fill = Comparisons)) +
                   geom_bar() +
                   scale_fill_manual(values = c("black","grey40")) +
                   theme(legend.position = "top"))
hist <- ggplot() +
  geom_bar(data[(eval(keep[2])),], mapping = aes(x = SNPbins), fill = "grey40") +
  geom_bar(data[(eval(keep[1])),], mapping = aes(x = SNPbins), fill = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dp1 <- plot_grid(l1, hist, ncol = 1, rel_heights = c(0.1, 0.9))
l2 <- get_legend(ggplot(td, aes(x = vals, y = vals, group = Comparisons, linetype = Comparisons)) +
                   geom_line() +
                   scale_linetype_manual(values = c("dashed","solid")) +
                   theme(legend.position = "top", legend.key=element_rect(fill= NA)))
dens <- ggplot() +
  geom_density(data[(eval(keep[1])),], mapping = aes(x = xpclr_norm, group = SNPbins, color = SNPbins)) +
  geom_density(data[(eval(keep[2])),], mapping = aes(x = xpclr_norm, group = SNPbins, color = SNPbins), linetype = "dashed") +
  theme_classic() +
  scale_x_sqrt()
dp2 <- plot_grid(l2, dens, ncol = 1, rel_heights = c(0.1, 0.9))
#pdf(paste("../../../Figures/xpclr_norm_qual_filt_",ds,"_binning_diagnostics.pdf",sep=""), width = 10, height = 5)
#plot_grid(dp1, dp2, rel_widths = c(0.47,0.53))
#dev.off()
cat(paste("Mean # SNPs/window:\n",keepn[1]," ",mean(data[(eval(keep[1])),10]),"\n",keepn[2]," ",mean(data[(eval(keep[2])),10])))
cat(paste("Mean # SNPs unused per comp:\n",keepn[1]," ",unusedm1,"\n",keepn[2]," ",unusedm2))

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
  datr <- data[(eval(keep[i]) & data$nSNPs >= wgmin),]
  t <- nrow(datr)
  v5 <- tail(head(datr[order(datr[,13], decreasing=T),13],nrow(datr)*.05),1)
  nv5 <- nrow(datr[(datr$xpclr_norm >= v5),])
  v1 <- tail(head(datr[order(datr[,13], decreasing=T),13],nrow(datr)*.01),1)
  nv1 <- nrow(datr[(datr$xpclr_norm >= v1),])
  sd5 <- mean(datr$xpclr_norm) + sd(datr$xpclr_norm)*5
  nsd5 <- nrow(datr[(datr$xpclr_norm >= sd5),])
  wg <- data.frame("SNPs/win" = "Whole genome*", "total # win" = t, "top 5% xpclr_norm cutoff" = round(v5,4),
                   "# win in top 5%" = nv5, "top 1% xpclr_norm cutoff" = round(v1,4), "% win in top 1%" = nv1, 
                   "5 SD xpclr_norm cutoff" = round(sd5,4), "# win > 5 SD" = nsd5)
  compsr <- unique(datr$comp)
  for (j in 1:length(compsr)) {
    if (j == 1) {
    }
    v5n <- nrow(datr[datr$comp == compsr[j] & datr$xpclr_norm >= v5,])
    v1n <- nrow(datr[datr$comp == compsr[j] & datr$xpclr_norm >= v1,])
    sd5n <- nrow(datr[datr$comp == compsr[j] & datr$xpclr_norm >= sd5,])
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
    m <- mean(dat[,13])
    sd5n <- m + sd(dat[,13])*5
    nrsd5n <- nrow(dat[(dat$xpclr_norm > sd5n),])
    t5 <- round(nr*.05,0)
    t1 <- round(nr*.01,0)
    v5 <- tail(head(dat[order(dat[,13], decreasing=T),13],t5),1)
    v1 <- tail(head(dat[order(dat[,13], decreasing=T),13],t1),1)
    bin <- data.frame("SNPs/win" = bins[j], "total # win" = nr, "top 5% xpclr_norm cutoff" = round(v5,4),
                      "# win in top 5%" = t5, "top 1% xpclr_norm cutoff" = round(v1,4), "% win in top 1%" = t1, 
                      "5 SD xpclr_norm cutoff" = round(sd5n,4), "# win > 5 SD" = nrsd5n)
    dfbin <- rbind(dfbin, bin)
    for (k in 1:length(compsr)) {
      d <- dat[(dat$comp == compsr[k]),]
      cv5 <- nrow(d[(d$xpclr_norm >= v5),])
      cv1 <- nrow(d[(d$xpclr_norm >= v1),])
      csd5n <- nrow(d[(d$xpclr_norm > sd5n),])
      df <- data.frame("comp" = compsr[k], "bin" = bins[j], "top5" = cv5, "top1" = cv1, "5SD" = csd5n)
      pop_sum <- rbind(df, pop_sum)
    }
    if (bins[j] != "1-25") {
      datp <- dat[(dat$xpclr_norm > sd5n),]
      data5sd <- rbind(datp, data5sd)
    }
  }
  binsum <- data.frame("SNPs/win" = "Sum", "total # win" = sum(dfbin$total...win[2:length(dfbin$total...win)]), "top 5% xpclr_norm cutoff" = "",
                       "# win in top 5%" = sum(dfbin$X..win.in.top.5.[2:length(dfbin$total...win)]), "top 1% xpclr_norm cutoff" = "", 
                       "% win in top 1%" = sum(dfbin$X..win.in.top.1.[2:length(dfbin$total...win)]), "5 SD xpclr_norm cutoff" = "", 
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
  d5sdwg <- rbind(d5sdwg, datr[(datr$xpclr_norm > wg[1,4]),])
}


##############
# MAKE PLOTS #
##############

axis_set <- dat %>% group_by(chrom) %>% summarise(center = mean(center_cum), end = max(center_cum))
axis_set$chr_e <- c(seq(1,3),seq(5,13),"  14",15,16,"",18,"","20 ","","  22")
n <- 1
for (i in comps) {
  cat(i,"\n")
  sub <- gsub('_',"-",i)
  dat <- data[(data$SNPbins != "1-25" & data$comp == i),]
  sig <- d5sdb[(d5sdb$comp == i),]
  p <- ggplot(dat, aes(x = center_cum, y = xpclr_norm, color = as.factor(chrom), alpha = xpclr_norm)) +
    geom_point(size = 0.7, shape = 16) +
    geom_point(sig, mapping = aes(x = center_cum, y = xpclr_norm), 
               size = 0.7, shape = 16, color = "#CC0000", alpha = 0.6) +
    scale_color_manual(values = rep(c("grey30", "grey70"), unique(length(axis_set$chrom)))) +
    scale_alpha_continuous(range = c(0.1, 0.5)) +
    scale_x_continuous(label = axis_set$chr_e, breaks = axis_set$center, expand = c(0.0025,0.0025)) +
    scale_y_continuous(limits = c(0,ceiling(max(d5sdb$xpclr_norm*2))/2)) +
    labs(y = expression('XP-CLR'[norm]), subtitle = sub) +
    theme(legend.position = "none", panel.background = element_blank(), 
          axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.text.x = element_text(size = 8), plot.subtitle = element_text(size = 10),
          axis.line.y = element_line(color = "#545454", linewidth = 0.5))
  #print(p)
  assign(paste("p",n,sep = ""), p)
  n <- n + 1
}

pdf("../../../Figures/XP-CLR_qual_filt_10Kb_5SDbinned_manhattan.pdf", width = 11, height = 8)
plot_grid(p4,p7,p5,p3,p8,p6,p2,p1, ncol = 2)
dev.off()

getwd()
###############
# GET REGIONS #
###############

## Binned cutoffs (> 5 SD from mean)
d5sdb$scaff <- paste("SUPER__", d5sdb$chrom, sep = "")
write.table(d5sdb, file = paste("../qual_filt_",ds,"_binned5SD_peaks.tsv", sep = ""), quote = F, sep = "\t", row.names = F)

## Whole-genome cutoffs (> 5 SD from mean)
d5sdwg$scaff <- paste("SUPER__", d5sdwg$chrom, sep = "")
write.table(d5sdwg, file = paste("../qual_filt_",ds,"_WG5SD_peaks.tsv", sep = ""), quote = F, sep = "\t", row.names = F)

# !! NOTE: !! need to join consecutive windows so as not to get multiple hits/gene
