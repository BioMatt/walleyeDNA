library(tidyverse)
library(vcfR)
library(adegenet)
library(hierfstat)
library(ggman)
library(patchwork)
library(RColorBrewer)
library(cowplot)
###################################################################

# Using Hierfstat to look at Fst 
# Reformatting the vcfR object into a genlight object
vcf_hf <- read.vcfR("HDplot_singletons.recode.vcf")

genlight <- vcfR2genlight(vcf_hf)

# Reformatting this to a matrix for input into hierfstat
matrix <- as.matrix(genlight)
# Re-order SNPs to match the assignment order downstream
matrix <- matrix[order(rownames(matrix)),]


# Prepare 0, 1, 2 matrix in hierfstat format
# Use locations to test for population differentiation
matrix[matrix==0] <- 11
matrix[matrix==1] <- 12
matrix[matrix==2] <- 22

assignments <- read_delim("dadi_assignments_Dec31.txt", delim = "\t", col_names = c("Ind", "Pop"))
# Re-order assignments to the same order as the SNPs
assignments <- assignments[order(assignments$Ind),]
identical(rownames(matrix), assignments[,1]) 

snp_df <- rownames_to_column(as.data.frame(matrix), var = "Ind")
identical(snp_df[,1], assignments[,1]) 
snp_df <- left_join(assignments, snp_df) %>% 
  arrange(desc(Pop)) %>% 
  column_to_rownames(var = "Ind")
snp_df[,1:3]


# Convert all values to numeric because while the Fst and betas functions work on the population names, the boot Fst does not
snp_df_numeric <- snp_df
snp_df_numeric$Pop <- as.numeric(as.factor(snp_df_numeric$Pop))
snp_df[,1:3]

snp_df_numeric[,1:3]
# Looking at pairwise Fst with hierfstat. 4=Lake Winnipeg, 3=Lake Manitoba, 2=Cedar Bluff, 1=Admixed
winman_wc_basic <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==3:4,], diploid = TRUE)
plot(x = 1:length(winman_wc_basic$perloc$Fstp), y = winman_wc_basic$perloc$Fstp)

winman_wc_fst <- wc(snp_df_numeric[snp_df_numeric$Pop==3:4,], diploid = TRUE)
winman_wc_fst$per.loc$FST
winman_wc_fst$FST
plot(x = 1:length(winman_wc_fst$per.loc$FST), y = winman_wc_fst$per.loc$FST)


winkan_wc_fst <- wc(snp_df_numeric[(snp_df_numeric$Pop==4 | snp_df_numeric$Pop==2),], diploid = TRUE)

winkan_wc_basic <- basic.stats(snp_df_numeric[(snp_df_numeric$Pop==4 | snp_df_numeric$Pop==2),], diploid = TRUE)
plot(x = 1:length(winkan_wc_basic$perloc$Fstp), y = winkan_wc_basic$perloc$Fstp)

winkan_wc_fst$per.loc$FST
winkan_wc_fst$FST
plot(x = 1:length(winkan_wc_fst$per.loc$FST), y = winkan_wc_fst$per.loc$FST)


mankan_wc_fst <- wc(snp_df_numeric[(snp_df_numeric$Pop==3 | snp_df_numeric$Pop==2),], diploid = TRUE)

mankan_wc_basic <- basic.stats(snp_df_numeric[(snp_df_numeric$Pop==3 | snp_df_numeric$Pop==2),], diploid = TRUE)
plot(x = 1:length(mankan_wc_basic$perloc$Fstp), y = mankan_wc_basic$perloc$Fstp)

mankan_wc_fst$per.loc$FST
mankan_wc_fst$FST
plot(x = 1:length(mankan_wc_fst$per.loc$FST), y = mankan_wc_fst$per.loc$FST)

# Pulling SNP information for plotting. Re-format NaN into NA
snp_table <- as_tibble(vcf_hf@fix) %>% 
  mutate(winman_fst = winman_wc_basic$perloc$Fstp, winkan_fst = winkan_wc_basic$perloc$Fstp, mankan_fst = mankan_wc_basic$perloc$Fstp, CHROM_ID = as.numeric(as.factor(CHROM))) %>% 
  mutate(winman_fst = na_if(winman_fst, "NaN"), winkan_fst = na_if(winkan_fst, "NaN"), mankan_fst = na_if(mankan_fst, "NaN")) %>% 
  filter(CHROM_ID <= 24)

# Read in pcadapt results for highlighting significant SNPs on different axes
pcadapt <- read_delim("PCAdapt_SNPs.txt", delim = " ") %>% 
  mutate(snp_id = paste0(CHROM, "_", POS), PC = as.factor(PC))
# Break up the pcadapt results into significant SNPs along PC1 and PC2
pcadapt_pc1 <- filter(pcadapt, PC == 1)
pcadapt_pc2 <- filter(pcadapt, PC == 2)


# Plot Fst between each of the three populations, highlighting the top 5% of points
winman_plot <- ggman(gwas = snp_table, snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "winman_fst", logTransform = FALSE, ymin = min(winman_wc_basic$perloc$Fstp, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = winman_wc_basic$overall[8], lineColour = "#fdae61", title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[ST]^"'")) + ggtitle(bquote(italic(F)[ST]^"'"*~between~Lake~Winnipeg~'&'~Lake~Manitoba~(Global~italic(F)[ST]^"'"~"="~0.0325)))
# These commented out rows highlight the top 0.01% of Fst outliers
#winman_highlights <- snp_table %>% slice_max(winman_fst, prop = 0.001)
#nrow(winman_highlights)
#winman_highlight_plot <- ggmanHighlight(winman_plot, highlight = winman_highlights$ID, colour = "#4575b4", size = 0.7)
# Instead, highlight the PC1 and PC2 SNPs
winman_highlight_plot <- ggmanHighlight(winman_plot, highlight = pcadapt_pc1$ID, colour = "#d73027", size = 0.7)
winman_highlight_plot
winman_highlight_plot2 <- ggmanHighlight(winman_highlight_plot, highlight = pcadapt_pc2$ID, colour = "#4575b4", size = 0.7)
winman_highlight_plot2

winkan_plot <- ggman(gwas = snp_table, snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "winkan_fst", logTransform = FALSE, ymin = min(winkan_wc_basic$perloc$Fstp, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = winkan_wc_basic$overall[8], lineColour = "#fdae61", title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[ST]^"'")) + ggtitle(bquote(italic(F)[ST]^"'"*~between~Lake~Winnipeg~'&'~Cedar~Bluff~Reservoir~(Global~italic(F)[ST]^"'"~"="~0.0826)))
# These commented out rows highlight the top 0.01% of Fst outliers
#winkan_highlights <- snp_table %>% slice_max(winkan_fst, prop = 0.001)
#nrow(winkan_highlights)
#winkan_highlight_plot <- ggmanHighlight(winkan_plot, highlight = winkan_highlights$ID, colour = "#4575b4", size = 0.7)
#winkan_highlight_plot
# Instead, highlight the PC1 and PC2 SNPs
winkan_highlight_plot <- ggmanHighlight(winkan_plot, highlight = pcadapt_pc1$ID, colour = "#d73027", size = 0.7)
winkan_highlight_plot
winkan_highlight_plot2 <- ggmanHighlight(winkan_highlight_plot, highlight = pcadapt_pc2$ID, colour = "#4575b4", size = 0.7)
winkan_highlight_plot2

mankan_plot <- ggman(gwas = snp_table, snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "mankan_fst", logTransform = FALSE, ymin = min(mankan_wc_basic$perloc$Fstp, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = mankan_wc_basic$overall[8], lineColour = "#fdae61", title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[ST]^"'")) + ggtitle(bquote(italic(F)[ST]^"'"*~between~Cedar~Bluff~Reservoir~'&'~Lake~Manitoba~(Global~italic(F)[ST]^"'"~"="~0.0859)))
# These commented out rows highlight the top 0.01% of Fst outliers
#mankan_highlights <- snp_table %>% slice_max(mankan_fst, prop = 0.001)
#nrow(mankan_highlights)
#mankan_highlight_plot <- ggmanHighlight(mankan_plot, highlight = mankan_highlights$ID, colour = "#4575b4", size = 0.7)
#mankan_highlight_plot

# Instead, highlight the PC1 and PC2 SNPs
mankan_highlight_plot <- ggmanHighlight(mankan_plot, highlight = pcadapt_pc1$ID, colour = "#d73027", size = 0.7)
mankan_highlight_plot
mankan_highlight_plot2 <- ggmanHighlight(mankan_highlight_plot, highlight = pcadapt_pc2$ID, colour = "#4575b4", size = 0.7)
mankan_highlight_plot2

# Put the plots together
combined_fst_plot <- winman_highlight_plot2 / winkan_highlight_plot2 / mankan_highlight_plot2

ggsave(filename = "combined_fst_plot.pdf", combined_fst_plot, dpi = 2000)

##########################################################################################################################
# Look at absolute allele frequency differences
#winman_freq <- winman_wc_basic$pop.freq

# Read in allele frequency data calculated with vcftools --freq option for each population
# Start with Lake Winnipeg individuals
winnipeg_freq <- read_delim("lkwinnipeg_freq.frq", delim = "\t", col_names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "Allele_1", "Allele_2"), skip = 1) %>% 
  separate(Allele_1, into = c("Allele1", "Allele1_freq"), sep = ":") %>% 
  separate(Allele_2, into = c("Allele2", "Allele2_freq"), sep = ":") %>% 
  mutate(Allele1_freq = as.numeric(Allele1_freq), Allele2_freq = as.numeric(Allele2_freq)) %>% 
  mutate(minor_allele_freq_lkWinnipeg = pmin(Allele1_freq, Allele2_freq)) %>% 
  mutate(minor_allele_lkwinnipeg = ifelse((minor_allele_freq_lkWinnipeg == Allele2_freq), Allele2, Allele1)) # In addition, check the minor allele in the Lake Winnipeg population. This will be used to calculate absolute allele frequency differences among all populations

# Read in allele frequency data for Lake Manitoba indidividuals
manitoba_freq <- read_delim("lkmanitoba_freq.frq", delim = "\t", col_names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "Allele_1", "Allele_2"), skip = 1) %>% 
  separate(Allele_1, into = c("Allele1", "Allele1_freq"), sep = ":") %>% 
  separate(Allele_2, into = c("Allele2", "Allele2_freq"), sep = ":") %>% 
  mutate(Allele1_freq = as.numeric(Allele1_freq), Allele2_freq = as.numeric(Allele2_freq)) %>% 
  mutate(minor_allele = ifelse(Allele2 == winnipeg_freq$minor_allele_lkwinnipeg, Allele2, Allele1)) %>%  # Find the minor allele based on the minor allele in the Lake Winnipeg population, for calculating absolute allele frequency differences 
  mutate(minor_allele_freq = ifelse(minor_allele == Allele2, Allele2_freq, Allele1_freq)) # Look up the minor allele frequency from the Lake Manitoba data

# Read in allele frequency data for Cedar Bluff Reservoir indidividuals
cedar_bluff_freq <- read_delim("cedar_bluff_freq.frq", delim = "\t", col_names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "Allele_1", "Allele_2"), skip = 1) %>% 
  separate(Allele_1, into = c("Allele1", "Allele1_freq"), sep = ":") %>% 
  separate(Allele_2, into = c("Allele2", "Allele2_freq"), sep = ":") %>% 
  mutate(Allele1_freq = as.numeric(Allele1_freq), Allele2_freq = as.numeric(Allele2_freq)) %>% 
  mutate(minor_allele = ifelse(Allele2 == winnipeg_freq$minor_allele_lkwinnipeg, Allele2, Allele1)) %>%  # Find the minor allele based on the minor allele in the Lake Winnipeg population, for calculating absolute allele frequency differences 
  mutate(minor_allele_freq = ifelse(minor_allele == Allele2, Allele2_freq, Allele1_freq)) # Look up the minor allele frequency from the Cedar Bluff data

# Combine the minor allele frequency data into one tibble. Turn NaN values into NA. Calculate absolute allele frequency differences between populations
combined_freq <- tibble(CHROM = winnipeg_freq$CHROM, POS = winnipeg_freq$POS, lkwinnipeg_minorallele = winnipeg_freq$minor_allele_lkwinnipeg, lkwinnipeg_maf = winnipeg_freq$minor_allele_freq_lkWinnipeg, lkmanitoba_maf = manitoba_freq$minor_allele_freq, cedarbluff_maf = cedar_bluff_freq$minor_allele_freq) %>% 
  mutate(cedarbluff_maf = na_if(cedarbluff_maf, "NaN"), lkmanitoba_maf = na_if(lkmanitoba_maf, "NaN"), lkwinnipeg_maf = na_if(lkwinnipeg_maf, "NaN")) %>% 
  mutate(winman_abs = abs((lkwinnipeg_maf - lkmanitoba_maf)), winkan_abs = abs((lkwinnipeg_maf - cedarbluff_maf)), mankan_abs = abs((lkmanitoba_maf - cedarbluff_maf))) %>% 
  mutate(SNP_ID = paste0(CHROM, "_", POS), chrom_id = as.numeric(as.factor(CHROM))) %>%  # Create a SNP identifier column for plotting, along with numeric chromosome IDs
 filter(chrom_id <= 24) # Filter for SNPs only on assembled chromosomes

# Plot all 3 comparisons of absolute allele frequency differences
winman_abs_plot <- ggman(gwas = combined_freq, snp = "SNP_ID", bp = "POS", chrom = "chrom_id", pvalue = "winman_abs", logTransform = FALSE, ymin = 0, ymax = 1, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Chromosome", pointSize = 0.2) + theme_bw() + ylab("Absolute Allele\n Frequency Difference") + ggtitle("Lake Winnpeg vs Lake Manitoba") + scale_colour_brewer(type = "qual", palette = 6)
winman_abs_plot

# Highlight the PC1 and PC2 SNPs
winman_abs_plot_highlight <- ggmanHighlight(winman_abs_plot, highlight = pcadapt_pc1$snp_id, colour = "#d73027", size = 0.7)
winman_abs_plot_highlight
winman_abs_plot_highlight2 <- ggmanHighlight(winman_abs_plot_highlight, highlight = pcadapt_pc2$snp_id, colour = "#4575b4", size = 0.7)
winman_abs_plot_highlight2


winkan_abs_plot <- ggman(gwas = combined_freq, snp = "SNP_ID", bp = "POS", chrom = "chrom_id", pvalue = "winkan_abs", logTransform = FALSE, ymin = 0, ymax = 1, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Chromosome", pointSize = 0.2) + theme_bw() + ylab("Absolute Allele\n Frequency Difference") + ggtitle("Lake Winnpeg vs Cedar Bluff Reservoir") + scale_colour_brewer(type = "qual", palette = 6)

# Highlight the PC1 and PC2 SNPs
winkan_abs_plot_highlight <- ggmanHighlight(winkan_abs_plot, highlight = pcadapt_pc1$snp_id, colour = "#d73027", size = 0.7)
winkan_abs_plot_highlight
winkan_abs_plot_highlight2 <- ggmanHighlight(winkan_abs_plot_highlight, highlight = pcadapt_pc2$snp_id, colour = "#4575b4", size = 0.7)
winkan_abs_plot_highlight2

mankan_abs_plot <- ggman(gwas = combined_freq, snp = "SNP_ID", bp = "POS", chrom = "chrom_id", pvalue = "mankan_abs", logTransform = FALSE, ymin = 0, ymax = 1, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Chromosome", pointSize = 0.2) + theme_bw() + ylab("Absolute Allele\n Frequency Difference") + ggtitle("Lake Manitoba vs Cedar Bluff Reservoir") + scale_colour_brewer(type = "qual", palette = 6)

# Highlight the PC1 and PC2 SNPs
mankan_abs_plot_highlight <- ggmanHighlight(mankan_abs_plot, highlight = pcadapt_pc1$snp_id, colour = "#d73027", size = 0.7)
mankan_abs_plot_highlight
mankan_abs_plot_highlight2 <- ggmanHighlight(mankan_abs_plot_highlight, highlight = pcadapt_pc2$snp_id, colour = "#4575b4", size = 0.7)
mankan_abs_plot_highlight2


# Put the plots together
combined_abs_plot <- winman_abs_plot_highlight2 / winkan_abs_plot_highlight2 / mankan_abs_plot_highlight2

ggsave(filename = "combined_abs_plot.pdf", combined_abs_plot, dpi = 2000)


# Use cowplot to retrieve a legend for the red and blue labels. This is a placeholder plot only meant to create a legend with
legend <- get_legend(ggplot(filter(pcadapt, !is.na(PC)), aes(x = CHROM, y = POS, colour = PC, group = PC)) +
  geom_point() +
  scale_colour_manual(values = c("#d73027", "#4575b4")) +
  theme(legend.key=element_blank()))

ggdraw(plot_grid(legend))
ggsave(filename = "legend.pdf", plot = legend)
# Focus in on just chromosome 8
combined_freq_chr8 <- combined_freq %>% 
  filter(chrom_id == 8) %>% 
  filter(POS > 15100000, POS < 16300000)

ggman(gwas = combined_freq_chr8, snp = "SNP_ID", bp = "POS", chrom = "chrom_id", pvalue = "winman_abs", logTransform = FALSE, ymin = 0, ymax = 1, relative.positions = TRUE, sigLine = 0.9, title = "", ylabel = 'q', xlabel = "Chromosome", pointSize = 0.8) + theme_bw() + ylab("Absolute Allele\n Frequency Difference") + ggtitle("Lake Winnpeg vs Lake Manitoba") + scale_colour_manual(values = "#4575b4") + scale_x_continuous(breaks = round(seq(min(combined_freq_chr8$POS), max(combined_freq_chr8$POS), by = 10000),1))


##################################################################################################
# Plot the candidate regions and the possible inversion from Chromosome 8

# Take a look at a zoomed in plot from the overall plot
#View(dplyr::filter(combined_freq, CHROM == "NC_041338.1", POS >= 15260000, POS <= 15900000))
ggmanZoom(ggmanPlot = winman_abs_plot, chromosome = "8", start.position = 14960000, end.position = 16200000, gene.tracks = FALSE)

# Filter the combined frequency table for a region around the putative inversion
# Add rows for the beginning and end of XP-EHH regions and for genes of interest. These rows represent 'fake' SNPs not used for plotting, with NA values for absolute allele frequency differences so no erroneous new SNPs are added to the figure.
putative_inversion <- dplyr::filter(combined_freq, CHROM == "NC_041338.1", POS >= 14960000, POS <= 16200000) %>% 
  select(CHROM, chrom_id, POS, SNP_ID, winman_abs) %>% 
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15260000, SNP_ID = "start_inv1", winman_abs = NA) %>% 
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15450000, SNP_ID = "end_inv1", winman_abs = NA) %>% 
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15710000, SNP_ID = "start_inv2", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15900000, SNP_ID = "end_inv2", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15281766, SNP_ID = "start_PDHX", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15314272, SNP_ID = "end_PDHX", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15318324, SNP_ID = "start_EHF", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15324968, SNP_ID = "end_EHF", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15722749, SNP_ID = "start_LRRC4C", winman_abs = NA) %>%
  add_row(CHROM = "NC_041338.1", chrom_id = 8, POS = 15724683, SNP_ID = "end_LRRC4C", winman_abs = NA) %>%
  arrange(POS)

# Make an initial plot with no geom rectangles, for getting indices for highlighted regions
winman_chr8 <- ggman(gwas = putative_inversion, snp = "SNP_ID", bp = "POS", chrom = "chrom_id", pvalue = "winman_abs", logTransform = FALSE, ymin = 0, ymax = 1, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Position on Chromosome 8", pointSize = 0.2) + theme_bw() + ylab("Absolute Allele\n Frequency Difference") + ggtitle("Lake Winnpeg vs Lake Manitoba") + scale_colour_brewer(type = "qual", palette = 6)

# Create a dataframe of candidate regions, based on XP-EHH
chr8_cand_regions <- tibble(chr = c("NC_041338.1", "NC_041338.1"), start = c(15260000, 15710000), end = c(15450000, 15900000), start_index = c(NA, NA), end_index = c(NA, NA))

# Using the Candidate region data, find the points in the plot Index that the beginning and end of each plot index corresponds to
# Add an index column to the putative inversion dataframe
putative_inversion <- putative_inversion %>% 
  mutate(SNP_index = rownames(putative_inversion))

# Loop through the candidate region dataframe, and based on the index of point X coordinates on the Manhattan plot, create beginning and end values for each candidate region. Note that even though this comes before the ggman code, that code needs to be run at least once *without* the geom_rect() object to get the dataframe of index values calculated by ggman.
for (region in 1:nrow(chr8_cand_regions)){
  cand_region <- chr8_cand_regions[region,]
  cand_snps <- filter(putative_inversion, CHROM == cand_region$chr, POS >= cand_region$start, POS <= cand_region$end)
  chr8_cand_regions$start_index[region] <- winman_chr8[["data"]][["index"]][as.numeric(cand_snps$SNP_index[1])]
  chr8_cand_regions$end_index[region] <- winman_chr8[["data"]][["index"]][as.numeric(cand_snps$SNP_index[nrow(cand_snps)])]
}
rm(region, cand_region, cand_snps)

# Follow the same process for highlighting a region, but for adding gene locations to the zoomed in plot. Refer back to the winman_chr8 object, the putative_inversion table, and the gene locations to identify the approximate plot coordinates to be highlighted
inversion_genes <- tibble(chr = c("NC_041338.1", "NC_041338.1", "NC_041338.1"), start = c(15281766, 15318324, 15722749), end = c(15314272, 15324968, 15724683), start_index = c(3.203866, 3.422982, 6.161746), end_index = c(3.215444, 3.559671, 6.339457), gene = c("PDHX", "EHF", "LRRC4C")) %>% 
  mutate(length = end - start)

# Loop through the candidate region dataframe, and based on the index of point X coordinates on the Manhattan plot, create beginning and end values for each candidate region. Note that even though this comes before the ggman code, that code needs to be run at least once *without* the geom_rect() object to get the dataframe of index values calculated by ggman.
for (region in 1:nrow(inversion_genes)){
  cand_region <- inversion_genes[region,]
  cand_snps <- filter(putative_inversion, CHROM == cand_region$chr, POS >= cand_region$start, POS <= cand_region$end)
  inversion_genes$start_index[region] <- winman_chr8[["data"]][["index"]][as.numeric(cand_snps$SNP_index[1])]
  inversion_genes$end_index[region] <- winman_chr8[["data"]][["index"]][as.numeric(cand_snps$SNP_index[nrow(cand_snps)])]
}
rm(region, cand_region, cand_snps)

# Remake the plots with the highlighted XP-EHH candidate regions
winman_chr8_xpehh <- ggman(gwas = putative_inversion, snp = "SNP_ID", bp = "POS", chrom = "chrom_id", pvalue = "winman_abs", logTransform = FALSE, ymin = 0, ymax = 1, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Position on Chromosome 8 (Mb)", pointSize = 0.2) + theme_bw() + ylab("Absolute Allele\n Frequency Difference") + ggtitle("Lake Winnipeg versus Lake Manitoba") + scale_colour_brewer(type = "qual", palette = 6) + geom_rect(data = chr8_cand_regions, mapping = aes(xmin = start_index, xmax = end_index, ymin = 0, ymax = 1), fill = "gray", inherit.aes = FALSE, alpha = 0.2) + geom_rect(data = inversion_genes, mapping = aes(xmin = start_index, xmax = end_index, ymin = 0.0, ymax = 0.2), fill = "orange", inherit.aes = FALSE, alpha = 0.4) + scale_x_continuous(breaks = seq(1, 9), labels = format(round((seq(14960000, 16200000, by = 137777.8)/1000000), digits = 2), nsmall = 2))
winman_chr8_xpehh
winman_chr8_highlight <- ggmanHighlight(winman_chr8_xpehh, highlight = pcadapt_pc1$snp_id, colour = "#d73027", size = 0.7)
winman_chr8_highlight
winman_chr8_highlight2 <- ggmanHighlight(winman_chr8_highlight, highlight = pcadapt_pc2$snp_id, colour = "#4575b4", size = 0.7)
winman_chr8_highlight2
ggsave(filename = "winman_chr8.pdf", winman_chr8_highlight2, dpi = 2000, height = 4, width = 11)


###############################################################################################################################
# Look at Fis within the three waterbodies, especially around chromosome 8
# Looking at pairwise Fst with hierfstat. 4=Lake Winnipeg, 3=Lake Manitoba, 2=Cedar Bluff, 1=Admixed
winnipeg_basic <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==4,], diploid = TRUE)
manitoba_basic <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==3,], diploid = TRUE)
kansas_basic <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==2,], diploid = TRUE)

# Add SNP metadata to the per locus stats
winnipeg_perloc <- cbind(winnipeg_basic$perloc, vcf_hf@fix)
manitoba_perloc <- cbind(manitoba_basic$perloc, vcf_hf@fix)
kansas_perloc <- cbind(kansas_basic$perloc, vcf_hf@fix)

# Plot Fis
winnipeg_fis <- ggman(gwas = dplyr::filter(winnipeg_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(winnipeg_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS])) + ggtitle("Lake Winnipeg")

manitoba_fis <- ggman(gwas = dplyr::filter(manitoba_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(manitoba_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS])) + ggtitle("Lake Manitoba")

kansas_fis <- ggman(gwas = dplyr::filter(kansas_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(kansas_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS])) + ggtitle("Cedar Bluff Reservoir")

# Plot Fst between Lake Winnipeg & Lake Manitoba for comparison
winman_chr8_fst <- ggman(gwas = dplyr::filter(snp_table, CHROM_ID == 8), snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "winman_fst", logTransform = FALSE, ymin = min(winman_wc_basic$perloc$Fstp, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = winman_wc_basic$overall[8], lineColour = "#fdae61", title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[ST]^"'")) + ggtitle(bquote(italic(F)[ST]^"'"*~between~Lake~Winnipeg~'&'~Lake~Manitoba~(Global~italic(F)[ST]^"'"~"="~0.0325)))

combined_fis <- winman_chr8_fst / winnipeg_fis / manitoba_fis / kansas_fis
ggsave(filename = "fis_per_waterbody.pdf", combined_fis, dpi = 2000)


##################################################
# Check the distribution of outliers for absolute allele frequency differences between Lake Winnipeg & Lake Manitoba
View(combined_freq %>% arrange(desc(winman_abs)))
