#install.packages("rehh", dependencies = TRUE)
# Using REHH to look at haplotype homozygosity info following the guide at https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
# Also using the tutorial https://speciationgenomics.github.io/haplotypes/
library(rehh)
library(vcfR)
library(tidyverse)
# The ggman package is used throughout to make Manhattan plots. See https://github.com/drveera/ggman
library(ggman)
library(patchwork)
# Read in the phased/imputed vcf just to get a list of the unique chromosomes in the file
vcf <- read.vcfR("imputed_0.9present.vcf.gz")
chroms <- unique(vcf@fix[,1])
colnames(vcf@fix)
# Create a map file from the vcf for use with rehh
map_file <- tibble(marker_ID = vcf@fix[,3], chromosome = vcf@fix[,1], position = vcf@fix[,2], ref = vcf@fix[,4], alt = vcf@fix[,5])
write_delim(map_file, "rehh_mapfile.inp", delim = "\t", col_names = FALSE)
rm(vcf)

# Loop through the Lake Winnipeg data using the list of chromosomes to scan_hh over each chromosome. This loop has to read over the whole vcf for each chromosome, so it is inefficient. But it works.
for(i in 1:length(chroms)) {
  # create internal representation
  hh <- data2haplohh(hap_file = "imputed_0.9_LkWinnipeg.recode.vcf",
                     map_file = "rehh_mapfile.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    Winnipeg_wgscan <- scan
  } else {
    Winnipeg_wgscan <- rbind(Winnipeg_wgscan, scan)
  }
}

# Loop through the Lake Manitoba data using the list of chromosomes to scan_hh over each chromosome
for(i in 1:length(chroms)) {
  # create internal representation
  hh <- data2haplohh(hap_file = "imputed_0.9_LkManitoba.recode.vcf",
                     map_file = "rehh_mapfile.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    Manitoba_wgscan <- scan
  } else {
    Manitoba_wgscan <- rbind(Manitoba_wgscan, scan)
  }
}

# Loop through Kansas data using the list of chromosomes to scan_hh over each chromosome
for(i in 1:length(chroms)) {
  # create internal representation
  hh <- data2haplohh(hap_file = "imputed_0.9_Kansas.recode.vcf",
                     map_file = "rehh_mapfile.inp", 
                     chr.name = chroms[i],
                     polarize_vcf = FALSE,
                     vcf_reader = "vcfR")
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh, polarized = FALSE)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    Kansas_wgscan <- scan
  } else {
    Kansas_wgscan <- rbind(Kansas_wgscan, scan)
  }
}

# Calculate the XP-EHH statistic between each pair of populations, First, Lakes Winnipeg & Manitoba
xpehh.Win_Man <- ies2xpehh(scan_pop1 = Winnipeg_wgscan,
                        scan_pop2 = Manitoba_wgscan,
                        popname1 = "Lake Winnipeg",
                        popname2 = "Lake Manitoba",
                        p.adjust.method = "fdr")
# Plot the XP-EHH values, filtering for -log10(q) > 1.30103 because -log10(0.05) = 1.30103
head(xpehh.Win_Man)
xpehh.Win_Man <- xpehh.Win_Man %>% 
  mutate(SNP = map_file$marker_ID, `XPEHH_Lake Winnipeg_Lake Manitoba` = tidyr::replace_na(`XPEHH_Lake Winnipeg_Lake Manitoba`, 0), chrom = CHR, CHR = as.numeric(as.factor(CHR)))


# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 3 'extreme' markers of q < 0.05
WinMan_candidates <- calc_candidate_regions(xpehh.Win_Man, threshold = 1.30103, pval = TRUE, window_size = 1E5, overlap = 1E4, min_n_extr_mrk = 2)
WinMan_candidates

# Filter for significant hits (q < 0.05) among the Win-Man XP-EHH SNPs
sig_xpehh_ManWin <- filter(xpehh.Win_Man, LOGPVALUE > 1.30103)
# Check which SNPs appear within a candidate region
sig_xpehh_ManWin$cand_region <- ifelse(sapply(sig_xpehh_ManWin$POSITION, function(p) 
  any(WinMan_candidates$START <= p & WinMan_candidates$END >= p)),"cand", NA)

# Pull positive (selection in Lake Winnipeg) and negative (selection in Lake Manitoba) SNPs for highlighting on the manhattan plot
xpehh_ManWin_pos <- filter(sig_xpehh_ManWin, LOGPVALUE > 1.30103 & `XPEHH_Lake Winnipeg_Lake Manitoba` > 0 & cand_region == "cand")
xpehh_ManWin_neg <- filter(sig_xpehh_ManWin, LOGPVALUE > 1.30103 & `XPEHH_Lake Winnipeg_Lake Manitoba` < 0 & cand_region == "cand")


# Use the ggman package to plot a Manhattan plot with relative chromosome and SNP positions intact
xpehh_win_man_plot <- ggman(gwas = xpehh.Win_Man, snp = "SNP", bp = "POSITION", chrom = "CHR", pvalue = "XPEHH_Lake Winnipeg_Lake Manitoba", logTransform = FALSE, ymin = -9, ymax = 9, relative.positions = TRUE, sigLine = NA, title = "Lake Winnipeg versus Lake Manitoba", ylabel = "XP-EHH", xlabel = "") + theme_bw() 
# The colour #4575b4 is blue, representing Lake Winnipeg here
highlight_xpehh_win_man_plot <- ggmanHighlight(xpehh_win_man_plot, highlight = xpehh_ManWin_pos$SNP, colour = "#4575b4", size = 0.7)
# The colour #d73027 is red, representing Lake Manitoba
twohighlight_highlight_xpehh_win_man_plot <- ggmanHighlight(highlight_xpehh_win_man_plot, highlight = xpehh_ManWin_neg$SNP, colour = "#d73027", size = 0.7)
twohighlight_highlight_xpehh_win_man_plot


# XP-EHH between Lake Winnipeg and Kansas
xpehh.Win_Kan <- ies2xpehh(scan_pop1 = Winnipeg_wgscan,
                           scan_pop2 = Kansas_wgscan,
                           popname1 = "Lake Winnipeg",
                           popname2 = "Kansas",
                           p.adjust.method = "fdr")
head(xpehh.Win_Kan)
xpehh.Win_Kan <- xpehh.Win_Kan %>% 
  mutate(SNP = map_file$marker_ID, `XPEHH_Lake Winnipeg_Kansas` = tidyr::replace_na(`XPEHH_Lake Winnipeg_Kansas`, 0), chrom = CHR, CHR = as.numeric(as.factor(CHR)))

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 3 'extreme' markers of q < 0.05
WinKan_candidates <- calc_candidate_regions(xpehh.Win_Kan, threshold = 1.30103, pval = TRUE, window_size = 1E5, overlap = 1E4, min_n_extr_mrk = 2)
WinKan_candidates

# Filter for significant hits (q < 0.05) among the Win-Kan XP-EHH SNPs
sig_xpehh_WinKan <- filter(xpehh.Win_Kan, LOGPVALUE > 1.30103)
# Check which SNPs appear within a candidate region
sig_xpehh_WinKan$cand_region <- ifelse(sapply(sig_xpehh_WinKan$POSITION, function(p) 
  any(WinKan_candidates$START <= p & WinKan_candidates$END >= p)),"cand", NA)

# Pull positive (selection in Lake Winnipeg) and negative (selection in Cedar Bluff) SNPs for highlighting on the manhattan plot
xpehh_WinKan_pos <- filter(sig_xpehh_WinKan, LOGPVALUE > 1.30103 & `XPEHH_Lake Winnipeg_Kansas` > 0 & cand_region == "cand")
xpehh_WinKan_neg <- filter(sig_xpehh_WinKan, LOGPVALUE > 1.30103 & `XPEHH_Lake Winnipeg_Kansas` < 0 & cand_region == "cand")

# Use the ggman package to plot a Manhattan plot with relative chromosome and SNP positions intact
xpehh_win_kan_plot <- ggman(gwas = xpehh.Win_Kan, snp = "SNP", bp = "POSITION", chrom = "CHR", pvalue = "XPEHH_Lake Winnipeg_Kansas", logTransform = FALSE, ymin = -9, ymax = 9, relative.positions = TRUE, sigLine = NA, title = "Lake Winnipeg versus Cedar Bluff Reservoir", ylabel = "XP-EHH", xlabel = "") + theme_bw()
# The colour #4575b4 is blue, representing Winnipeg here
highlight_xpehh_win_kan_plot <- ggmanHighlight(xpehh_win_kan_plot, highlight = xpehh_WinKan_pos$SNP, colour = "#4575b4", size = 0.7)
# The colour #fdae61 is pastel, representing Cedar Bluff. Commented out because no significant SNPs were found for Cedar Bluff
twohighlight_highlight_xpehh_win_kan_plot <- ggmanHighlight(highlight_xpehh_win_kan_plot, highlight = xpehh_WinKan_neg$SNP, colour = "#fdae61", size = 0.7)
twohighlight_highlight_xpehh_win_kan_plot




# XP-EHH between Lake Manitoba and Kansas
xpehh.Man_Kan <- ies2xpehh(scan_pop1 = Manitoba_wgscan,
                           scan_pop2 = Kansas_wgscan,
                           popname1 = "Lake Manitoba",
                           popname2 = "Kansas",
                           p.adjust.method = "fdr")
head(xpehh.Man_Kan)
xpehh.Man_Kan <- xpehh.Man_Kan %>% 
  mutate(SNP = map_file$marker_ID, `XPEHH_Lake Manitoba_Kansas` = tidyr::replace_na(`XPEHH_Lake Manitoba_Kansas`, 0), chrom = CHR, CHR = as.numeric(as.factor(CHR)))

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 3 'extreme' markers of q < 0.05
ManKan_candidates <- calc_candidate_regions(xpehh.Man_Kan, threshold = 1.30103, pval = TRUE, window_size = 1E5, overlap = 1E4, min_n_extr_mrk = 2)
ManKan_candidates
# Filter for significant hits (q < 0.05) among the Man-Kan XP-EHH SNPs
sig_xpehh_ManKan <- filter(xpehh.Man_Kan, LOGPVALUE > 1.30103)
# Check which SNPs appear within a candidate region
sig_xpehh_ManKan$cand_region <- ifelse(sapply(sig_xpehh_ManKan$POSITION, function(p) 
  any(ManKan_candidates$START <= p & ManKan_candidates$END >= p)),"cand", NA)

xpehh_ManKan_pos <- filter(sig_xpehh_ManKan, LOGPVALUE > 1.30103 & `XPEHH_Lake Manitoba_Kansas` > 0 & cand_region == "cand")
xpehh_ManKan_neg <- filter(sig_xpehh_ManKan, LOGPVALUE > 1.30103 & `XPEHH_Lake Manitoba_Kansas` < 0 & cand_region == "cand")

# Use the ggman package to plot a Manhattan plot with relative chromosome and SNP positions intact
xpehh_man_kan_plot <- ggman(gwas = xpehh.Man_Kan, snp = "SNP", bp = "POSITION", chrom = "CHR", pvalue = "XPEHH_Lake Manitoba_Kansas", logTransform = FALSE, ymin = -9, ymax = 9, relative.positions = TRUE, sigLine = NA, title = "Lake Manitoba versus Cedar Bluff Reservoir", ylabel = "XP-EHH", xlabel = "Chromosome") + theme_bw()
# The colour #d73027 is red, representing Lake Manitoba
highlight_xpehh_man_kan_plot <- ggmanHighlight(xpehh_man_kan_plot, highlight = xpehh_ManKan_pos$SNP, colour = "#d73027", size = 0.7)
# The colour #fdae61 is pastel, representing Cedar Bluff. Commented out because no significant SNPs were found for Cedar Bluff
#twohighlight_xpehh_man_kan_plot <- ggmanHighlight(highlight_xpehh_man_kan_plot, highlight = xpehh_ManKan_neg$SNP, colour = "#fdae61", size = 0.7)
#twohighlight_xpehh_man_kan_plot
highlight_xpehh_man_kan_plot
# Plot all three Manhattan plots for XP-EHH
xpehh_plots_combined <- twohighlight_highlight_xpehh_win_man_plot / twohighlight_highlight_xpehh_win_kan_plot/ highlight_xpehh_man_kan_plot
ggsave("combined_XP-EHH.pdf", plot = xpehh_plots_combined, dpi = 2000)


############### 
# Read in the chromosome IDs used with magma, attach them to the XP-EHH results, and write out significant SNPs in candidate regions for analysis with genes using magma
# First, the magma chromosome IDs
magma_chrom_ids <- read_delim("magma_chrom_IDs.txt", delim = "\t", col_names = c("chrom", "magma_chrom_ID")) %>% 
  distinct(chrom, .keep_all = TRUE)

# Merge the various XP-EHH significant candidate SNP tables with the magma chromosomes. They automatically merge by chromosome because column names match
# Winnipeg & Manitoba
xpehh_ManWin_pos <- left_join(xpehh_ManWin_pos, magma_chrom_ids)
xpehh_ManWin_neg <- left_join(xpehh_ManWin_neg, magma_chrom_ids)
# Winnipeg & Kansas
xpehh_WinKan_pos <- left_join(xpehh_WinKan_pos, magma_chrom_ids)
xpehh_WinKan_neg <- left_join(xpehh_WinKan_neg, magma_chrom_ids)
# Manitoba & Kansas
xpehh_ManKan_pos <- left_join(xpehh_ManKan_pos, magma_chrom_ids)
xpehh_ManKan_neg <- left_join(xpehh_ManKan_neg, magma_chrom_ids)

# Magma annotation input files have SNP ID, magma chromosome ID, and base pair position with no column names
# Name these text files with the comparison and what negative or positive selection implies for a site. Here, pos values are relevant to Lake Winnipeg
write_delim(as.data.frame(cbind(xpehh_ManWin_pos$SNP, xpehh_ManWin_pos$magma_chrom_ID, xpehh_ManWin_pos$POSITION)), "xpehh_WinMan_PosWinnipeg_SNPs_magma_0.9present.txt", col_names = FALSE)
# From the same XP-EHH run, but negative values and hence relevant to Lake Manitoba
write_delim(as.data.frame(cbind(xpehh_ManWin_neg$SNP, xpehh_ManWin_neg$magma_chrom_ID, xpehh_ManWin_neg$POSITION)), "xpehh_WinMan_NegManitoba_SNPs_magma_0.9present.txt", col_names = FALSE)
# Writing Winnipeg & Kansas
write_delim(as.data.frame(cbind(xpehh_WinKan_pos$SNP, xpehh_WinKan_pos$magma_chrom_ID, xpehh_WinKan_pos$POSITION)), "xpehh_WinKan_PosWinnipeg_SNPs_magma_0.9present.txt", col_names = FALSE)
write_delim(as.data.frame(cbind(xpehh_WinKan_neg$SNP, xpehh_WinKan_neg$magma_chrom_ID, xpehh_WinKan_neg$POSITION)), "xpehh_WinKan_NegKansas_SNPs_magma_0.9present.txt", col_names = FALSE)
# Writing Lake Manitoba and Kansas
write_delim(as.data.frame(cbind(xpehh_ManKan_pos$SNP, xpehh_ManKan_pos$magma_chrom_ID, xpehh_ManKan_pos$POSITION)), "xpehh_ManKan_PosManitoba_SNPs_magma_0.9present.txt", col_names = FALSE)
write_delim(as.data.frame(cbind(xpehh_ManKan_neg$SNP, xpehh_ManKan_neg$magma_chrom_ID, xpehh_ManKan_neg$POSITION)), "xpehh_ManKan_NegKansas_SNPs_magma_0.9present.txt", col_names = FALSE)

###########################################################################################################
# Moving to iHS in regions of the genome now
# calculate genome-wide iHS values for Lake Winnipeg. Set freqbin to one here, because our data are unpolarized
Winnipeg_wgscan.ihs <- ihh2ihs(Winnipeg_wgscan, min_maf = 0.05, freqbin = 1, include_freq = TRUE)
head(Winnipeg_wgscan.ihs$ihs)
# Check the top 1% of iHS values. Pull the minimum of the top 1% of absolute IHS values as the minimum filtering threshold for highlighting SNPs on the Manhattan plot
Winnipeg_top_IHS <- Winnipeg_wgscan.ihs$ihs %>% 
  mutate(abs_IHS = abs(IHS)) %>% 
  slice_max(abs_IHS, prop = 0.01)
# The minimum of the top 1% of IHS is 2.902531
min(Winnipeg_top_IHS$abs_IHS)
ggplot(Winnipeg_top_IHS, aes(POSITION, IHS)) + geom_point()
ggplot(Winnipeg_wgscan.ihs$ihs, aes(POSITION, IHS)) + geom_point()
ggplot(Winnipeg_wgscan.ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 2 'extreme' markers of q < 0.05
Winnipeg_candidates <- calc_candidate_regions(Winnipeg_wgscan.ihs, threshold = 2.9, pval = FALSE, window_size = 1E5, overlap = 1E4, min_n_extr_mrk = 2)
Winnipeg_candidates


# Plot iHS for Lake Winnipeg. First, pull the IHS dataframe and stick on SNP IDs plus the full SNP dataset. Replace missing significance values with 0
Lk_Winnipeg_IHS <- Winnipeg_wgscan.ihs$ihs %>% 
  rename(position = POSITION, chromosome = CHR) %>% 
  mutate(position = as.character(position)) %>% 
  left_join(map_file, .) %>% 
  mutate(IHS = tidyr::replace_na(IHS, 0), chromosome = as.numeric(as.factor(chromosome)), abs_IHS = abs(IHS))

# Filter for significant hits (q < 0.05) among the Man-Kan XP-EHH SNPs
sig_IHS_Winnipeg <- filter(Lk_Winnipeg_IHS, IHS >= 2.902531)
# Check which SNPs appear within a candidate region
sig_IHS_Winnipeg$cand_region <- ifelse(sapply(sig_IHS_Winnipeg$position, function(p) 
  any(Winnipeg_candidates$START <= p & Winnipeg_candidates$END >= p)),"cand", NA)

# Filter for the top 1% of IHS values using the threshold found earlier  
Winnipeg_IHS_highlight <- filter(sig_IHS_Winnipeg, cand_region == "cand")

# Plot the base plot
Winnipeg_IHS_plot <- ggman(gwas = Lk_Winnipeg_IHS, snp = "marker_ID", bp = "position", chrom = "chromosome", pvalue = "abs_IHS", logTransform = FALSE, ymin = 0, ymax = 6, relative.positions = TRUE, sigLine = 2.902531, title = "iHS in Lake Winnipeg", ylabel = "iHS", xlabel = "") + theme_bw()
# The colour #4575b4 is blue, representing Winnipeg here
highlight_Winnipeg_IHS_plot <- ggmanHighlight(Winnipeg_IHS_plot, highlight = Winnipeg_IHS_highlight$marker_ID, colour = "#4575b4")
highlight_Winnipeg_IHS_plot



# Follow the same process for Lake Manitoba
Manitoba_wgscan.ihs <- ihh2ihs(Manitoba_wgscan, min_maf = 0.05, freqbin = 1, p.adjust.method = "fdr", include_freq = TRUE)
head(Manitoba_wgscan.ihs$ihs)
# Check the top 1% of iHS values
Manitoba_top_IHS <- Manitoba_wgscan.ihs$ihs %>% 
  mutate(abs_IHS = abs(IHS)) %>% 
  slice_max(abs_IHS, prop = 0.01)
# The minimum of the top 1% of IHS is 2.839257
min(Manitoba_top_IHS$abs_IHS)

# Look at candidate regions by using windows of size 100 kb overlapping by 10 kb with at least 2 'extreme' markers of q < 0.05
Manitoba_candidates <- calc_candidate_regions(Manitoba_wgscan.ihs, threshold = 2.839, pval = FALSE, window_size = 1E5, overlap = 1E4, min_n_extr_mrk = 2)
Manitoba_candidates


# Plot iHS for Lake Manitoba. First, pull the IHS dataframe and stick on SNP IDs plus the full SNP dataset. Replace missing significance values with 0
Lk_Manitoba_IHS <- Manitoba_wgscan.ihs$ihs %>% 
  rename(position = POSITION, chromosome = CHR) %>% 
  mutate(position = as.character(position)) %>% 
  left_join(map_file, .) %>% 
  mutate(IHS = tidyr::replace_na(IHS, 0), chromosome = as.numeric(as.factor(chromosome)), abs_IHS = abs(IHS))



# Filter for significant hits (q < 0.05) among the Man-Kan XP-EHH SNPs
sig_IHS_Manitoba <- filter(Lk_Manitoba_IHS, IHS >= 2.839257)
# Check which SNPs appear within a candidate region
sig_IHS_Manitoba$cand_region <- ifelse(sapply(sig_IHS_Manitoba$position, function(p) 
  any(Manitoba_candidates$START <= p & Manitoba_candidates$END >= p)),"cand", NA)

# Filter for the top 1% of IHS values using the threshold found earlier  
Lk_Manitoba_IHS_highlight <- filter(sig_IHS_Manitoba, cand_region == "cand")


# Plot the base plot
Manitoba_IHS_plot <- ggman(gwas = Lk_Manitoba_IHS, snp = "marker_ID", bp = "position", chrom = "chromosome", pvalue = "abs_IHS", logTransform = FALSE, ymin = 0, ymax = 6, relative.positions = TRUE, sigLine = 2.839257, title = "iHS in Lake Manitoba", ylabel = "iHS", xlabel = "") + theme_bw() + scale_colour_grey()
# The colour #d73027 is red, representing Lake Manitoba. Commented out because no candidate regions were significant
#highlight_Manitoba_IHS_plot <- ggmanHighlight(Manitoba_IHS_plot, highlight = Lk_Manitoba_IHS_highlight$marker_ID, colour = "#d73027")
Manitoba_IHS_plot





distribplot(Manitoba_wgscan.ihs$ihs$IHS, xlab = "iHS", qqplot = TRUE)
manhattanplot(Manitoba_wgscan.ihs,pval = TRUE, threshold = 4, main = "iHS (Lake Manitoba walleye)")
manhattanplot(Manitoba_wgscan.ihs, pval = FALSE, threshold = c(-4,4), main = "iHS (Lake Manitoba walleye)")

# iHS in Kansas walleye
Kansas_wgscan.ihs <- ihh2ihs(Kansas_wgscan, min_maf = 0.05, freqbin = 1, p.adjust.method = "fdr", include_freq = TRUE)
# Check the top 1% of iHS values
Kansas_top_IHS <- Kansas_wgscan.ihs$ihs %>% 
  mutate(abs_IHS = abs(IHS)) %>% 
  slice_max(abs_IHS, prop = 0.01)
# The minimum of the top 1% of IHS is 3.022571
min(Kansas_top_IHS$abs_IHS)


Kansas_candidates <- calc_candidate_regions(Kansas_wgscan.ihs, threshold = 3.02, pval = FALSE, window_size = 1E5, overlap = 1E4, min_n_extr_mrk = 1)
Kansas_candidates

# Plot iHS for Kansas. First, pull the IHS dataframe and stick on SNP IDs plus the full SNP dataset. Replace missing significance values with 0
Kansas_IHS <- Kansas_wgscan.ihs$ihs %>% 
  rename(position = POSITION, chromosome = CHR) %>% 
  mutate(position = as.character(position)) %>% 
  left_join(map_file, .) %>% 
  mutate(IHS = tidyr::replace_na(IHS, 0), chromosome = as.numeric(as.factor(chromosome)), abs_IHS = abs(IHS))

# Filter for significant hits (q < 0.05) among the Man-Kan XP-EHH SNPs
sig_IHS_Kansas <- filter(Kansas_IHS, IHS >= 3.022571)
# Check which SNPs appear within a candidate region
sig_IHS_Kansas$cand_region <- ifelse(sapply(sig_IHS_Kansas$position, function(p) 
  any(Kansas_candidates$START <= p & Kansas_candidates$END >= p)),"cand", NA)

# Filter for the top 1% of IHS values using the threshold found earlier  
Lk_Manitoba_IHS_highlight <- filter(sig_IHS_Manitoba, cand_region == "cand")
# Filter for the top 1% of IHS values using the threshold found earlier  
Kansas_IHS_IHS_highlight <- filter(Kansas_IHS, abs(IHS) >= 2.828413)

# Plot the base plot
Kansas_IHS_plot <- ggman(gwas = Kansas_IHS, snp = "marker_ID", bp = "position", chrom = "chromosome", pvalue = "abs_IHS", logTransform = FALSE, ymin = 0, ymax = 7, relative.positions = TRUE, sigLine = NA, title = "iHS in Cedar Bluff Reservoir", ylabel = "iHS", xlabel = "Chromosome") + theme_bw()
# The colour #fdae61 is pastel, representing Cedar Bluff
highlight_Kansas_IHS_plot <- ggmanHighlight(Kansas_IHS_plot, highlight = Kansas_IHS_IHS_highlight$marker_ID, colour = "#fdae61")
highlight_Kansas_IHS_plot


# put the IHS plots together
highlight_Winnipeg_IHS_plot / highlight_Manitoba_IHS_plot / highlight_Kansas_IHS_plot




distribplot(Kansas_wgscan.ihs$ihs$IHS, xlab = "iHS", qqplot = TRUE)
manhattanplot(Kansas_wgscan.ihs,pval = TRUE, threshold = 4, main = "iHS (Kansas walleye)")
manhattanplot(Kansas_wgscan.ihs, pval = FALSE, threshold = c(-4,4), main = "iHS (Kansas walleye)")



furcation <- calc_furcation(wgscan, mrk = "211582:147:-")
plot(furcation, hap.names = hap.names(hh), xlim = c(1.58E+7, 1.6E+7))

haplen <- calc_haplen(furcation)
haplen$mrk.name
plot(haplen)

ihs <- ihh2ihs(scan, min_maf = 0.01, freqbin = 0)
head(ihs$ihs)
cr.cgu <- calc_candidate_regions(ihs,
                                 threshold = 4,
                                 pval = TRUE,
                                 window_size = 1E6,
                                 overlap = 1E5,
                                 min_n_extr_mrk = 2)
cr.cgu
freqbinplot(ihs)
distribplot(ihs$ihs$IHS, xlab = "iHS")
manhattanplot(ihs,
              main = "iHS (overall walleye)")

###############################################
# Take a look at a region of interest in chromosome 38 for each population
Winnipeg_chr38_hh <- data2haplohh(hap_file = "LkWinnipeg_phased.recode.vcf",
                                  map_file = "rehh_mapfile.inp", 
                                  chr.name = "NC_041338.1",
                                  polarize_vcf = FALSE,
                                  vcf_reader = "vcfR")
Manitoba_chr38_hh <- data2haplohh(hap_file = "LkManitoba_phased.recode.vcf",
                                  map_file = "rehh_mapfile.inp", 
                                  chr.name = "NC_041338.1",
                                  polarize_vcf = FALSE,
                                  vcf_reader = "vcfR")
Kansas_chr38_hh <- data2haplohh(hap_file = "Kansas_phased.recode.vcf",
                                map_file = "rehh_mapfile.inp", 
                                chr.name = "NC_041338.1",
                                polarize_vcf = FALSE,
                                vcf_reader = "vcfR")
overall_chr38_hh <- data2haplohh(hap_file = "beagle_imputed.vcf.gz",
                                 map_file = "rehh_mapfile.inp", 
                                 chr.name = "NC_041338.1",
                                 polarize_vcf = FALSE,
                                 vcf_reader = "vcfR")
admixed_chr38_hh <- data2haplohh(hap_file = "Admixed_phased.recode.vcf",
                                 map_file = "rehh_mapfile.inp", 
                                 chr.name = "NC_041338.1",
                                 polarize_vcf = FALSE,
                                 vcf_reader = "vcfR")
par(mfrow=c(3,1))
plot(calc_ehhs(Winnipeg_chr38_hh,
               mrk = "211567:185:+",
               phased = TRUE))
plot(calc_ehhs(Manitoba_chr38_hh,
               mrk = "211582:147:-",
               phased = TRUE))
plot(calc_ehhs(Kansas_chr38_hh,
               mrk = "211582:147:-",
               phased = TRUE))

# Look at bifurcation plots
furcation_Admixed <- calc_furcation(admixed_chr38_hh,
                                    mrk = "211567:185:+")
furcation_Overall <- calc_furcation(overall_chr38_hh,
                                mrk = "211567:185:+")
furcation_Win <- calc_furcation(Winnipeg_chr38_hh,
                            mrk = "211567:185:+")
furcation_Man <- calc_furcation(Manitoba_chr38_hh,
                            mrk = "211582:147:-")
furcation_Kan <- calc_furcation(Kansas_chr38_hh,
                             mrk = "211582:147:-")
plot(furcation_Admixed,
     hap.names = hap.names(admixed_chr38_hh),
     xlim = c(1.49E+7, 1.69E+7))
plot(furcation_Overall,
     hap.names = hap.names(overall_chr38_hh),
     xlim = c(1.49E+7, 1.69E+7))
plot(furcation_Win,
     hap.names = hap.names(Winnipeg_chr38_hh),
     xlim = c(1.49E+7, 1.69E+7))
plot(furcation_Man,
     hap.names = hap.names(Manitoba_chr38_hh),
     xlim = c(1.58E+7, 1.6E+7))
plot(furcation_Kan,
     hap.names = hap.names(Kansas_chr38_hh),
     xlim = c(1.58E+7, 1.6E+7))
