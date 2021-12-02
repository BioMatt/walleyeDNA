library(tidyverse)
library(vcfR)
library(adegenet)
library(hierfstat)
library(ggman)
###################################################################

# Using Hierfstat to look at Fst 
# Reformatting the vcfR object into a genlight object
vcf_hf <- read.vcfR("HDplot_singletons_pruned.recode.vcf")

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
assignments <- assignments[order(assignments$Indiv),]
identical(rownames(matrix), assignments[,1]) 

snp_df <- rownames_to_column(as.data.frame(matrix), var = "Ind")
identical(snp_df[,1], assignments[,1]) 
snp_df <- left_join(assignments, snp_df) %>% 
  filter(Pop != "Admixed") %>% 
  column_to_rownames(var = "Ind")


#hf_pca <- indpca(snp_df)
#plot(hf_pca,col=snp_df[,1],cex=0.7)
# Playing with betas, or coalescent Fst
# 1 is Kansas, 2 is Lake Manitoba, 3 is Lake Winnipeg
betas_reassigned <- betas(snp_df, nboot = 1000)
betas_reassigned$betaiovl
betas_reassigned$ci

# Increase the memory size before running betas individually
# None of this code worked on a PC with 32 gb of ram, commenting out the code.
# memory.size(max = TRUE)
# memory.limit(size = 32679)
# snp_df_ind <- cbind(1:nrow(snp_df),snp_df[,-1])
# gc(reset = TRUE)
# betas_ind <- betas(snp_df_ind, nboot = 0, betaijT = TRUE)

# Looking at pairwise Fst with hierfstat
wc_fst <- pairwise.WCfst(snp_df, diploid = TRUE)
wc_fst

# Convert all values to numeric because while the Fst and betas functions work on the population names, the boot Fst does not
snp_df_numeric <- snp_df
snp_df_numeric$Pop <- as.numeric(as.factor(snp_df_numeric$Pop))

boot_fst <- boot.ppfst(snp_df_numeric, nboot = 1000)
boot_fst$ll
boot_fst$ul
#######################################################################################
# Take a look at fst across the genome
canada_basic_stats <- basic.stats(snp_df[snp_df$V1==2:3,])
canada_perloc <- canada_basic_stats$perloc
canada_perloc$snp_index <- 1:nrow(canada_perloc)
# Trying to make my own Manhattan plot that I can colour
ggplot(data = canada_perloc, aes(x = snp_index, y = Fst)) +
  geom_point()
per_locus_fst <- pp_fst$vc.per.loc


# Basic stats
basic_stats <- basic.stats(snp_df, diploid = TRUE)
basic_stats$overall

# Basic stats for Kansas walleye
basic_kansas <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==1,], diploid = TRUE)
basic_kansas$overall


# Basic stats for Lake Manitoba walleye
basic_manitoba <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==2,], diploid = TRUE)
basic_manitoba$overall

# Basic stats for Lake Winnipeg walleye
basic_winnipeg <- basic.stats(snp_df_numeric[snp_df_numeric$Pop==3,], diploid = TRUE)
basic_winnipeg$overall

# Check fis bootstraps for each population
kansas_boot_fis <- boot.ppfis(snp_df_numeric[snp_df_numeric$Pop==1,], diploid = TRUE, nboot = 1000)
kansas_boot_fis$fis.ci
manitoba_boot_fis <- boot.ppfis(snp_df_numeric[snp_df_numeric$Pop==2,], diploid = TRUE, nboot = 1000)
manitoba_boot_fis$fis.ci
winnipeg_boot_fis <- boot.ppfis(snp_df_numeric[snp_df_numeric$Pop==3,], diploid = TRUE, nboot = 1000)
winnipeg_boot_fis$fis.ci


# Allelic richness
overall_allelicrichness <- allelic.richness(snp_df, diploid = TRUE)
overall_allelicrichness$min.all
overall_allelicrichness$Ar

mean(overall_allelicrichness$Ar$Kansas, na.rm = TRUE)
mean(overall_allelicrichness$Ar$Manitoba, na.rm = TRUE)
mean(overall_allelicrichness$Ar$Winnipeg, na.rm = TRUE)

##############################################################################
# Write out fstat files for each population for Ne Estimator v2
# First, pull individual data frames for each population
Winnipeg_hf <- as_tibble(t(snp_df_numeric[snp_df_numeric$Pop==3,1:1000]))
manitoba_hf <- snp_df_numeric[snp_df_numeric$Pop==2,]
kansas_hf <- snp_df_numeric[snp_df_numeric$Pop==1,]
# Filter for SNPs with complete data within each dataframe. Ignore the population column.
test <- Winnipeg_hf %>% 
  filter_all(any_vars(!is.na(.)))

write.fstat(snp_df_numeric[snp_df_numeric$Pop==1,], "Kansas_fstat.dat")
write.fstat(snp_df_numeric[snp_df_numeric$Pop==2,], "Manitoba_fstat.dat")
write.fstat(snp_df_numeric[snp_df_numeric$Pop==3,], "Winnipeg_fstat.dat")

# Writing a short version of the Kansas data
write.fstat(snp_df_numeric[snp_df_numeric$Pop==1,c(1:500)], "Kansas_short_fstat.dat")


# Write out a fstat file with all groups, including 'admixed' individuals
fstat_df <- rownames_to_column(as.data.frame(matrix), var = "Ind")
identical(fstat_df[,1], assignments[,1]) 
fstat_df <- left_join(assignments, fstat_df) %>% 
  arrange(desc(Pop)) %>% 
  column_to_rownames(var = "Ind")
fstat_df[,1:3]
fstat_df$Pop <- as.numeric(as.factor(fstat_df$Pop))
fstat_df <- arrange(fstat_df, Pop)
write.fstat(fstat_df, "all_fstat.dat")

###################################################################################
# Take a look at admixed individuals
admixed_indivs <- rownames_to_column(as.data.frame(matrix), var = "Ind")
identical(admixed_indivs[,1], assignments[,1]) 
admixed_indivs <- left_join(assignments, admixed_indivs) %>% 
  filter(Pop == "Admixed") %>% 
  column_to_rownames(var = "Ind")

admixed_indivs$Pop <- as.numeric(as.factor(admixed_indivs$Pop))
# Basic stats
admixed_basic_stats <- basic.stats(admixed_indivs, diploid = TRUE)
admixed_basic_stats$overall

admixed_basic_stats$perloc$Fis

# Create a combined dataframe of all SNPs with info for the admixed fish
admixed_perloc <- cbind(admixed_basic_stats$perloc, vcf_hf@fix)
admixed_perloc$index <- 1:nrow(admixed_perloc)
admixed_perloc$POS <- as.numeric(admixed_perloc$POS)

ggplot(data = admixed_perloc, aes(x = index, y = Fis, group = CHROM)) +
  geom_point() +
  facet_wrap(~CHROM)

ggman(gwas = admixed_perloc, snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(admixed_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS]))


ggman(gwas = dplyr::filter(admixed_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(admixed_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS]))

ggman(gwas = dplyr::filter(admixed_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Ho", logTransform = FALSE, ymin = min(admixed_perloc$Ho, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Ho', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(H)[O]))

# Read in HWE deviations from vcftools for these admixed individuals
admixed_hwe <- read_tsv("admixed_hwe.hwe") 
colnames(admixed_hwe) <- c("CHROM", "POS", "OBS(HOM1/HET/HOM2)", "E(HOM1/HET/HOM2)", "ChiSq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS")
# Combine with the overall data
admixed_perloc <- left_join(admixed_perloc, admixed_hwe)

ggman(gwas = dplyr::filter(admixed_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "P_HWE", logTransform = TRUE, ymin = min(admixed_perloc$Ho, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'P_HWE', xlabel = "Chromosome") + theme_bw()

# Check Fis in Lake Winnipeg & Lake Manitoba fish
winnipeg_perloc <- cbind(basic_winnipeg$perloc, vcf_hf@fix)
ggman(gwas = dplyr::filter(winnipeg_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(admixed_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS])) + ggtitle("Lake Winnipeg")

manitoba_perloc <- cbind(basic_manitoba$perloc, vcf_hf@fix)
ggman(gwas = dplyr::filter(manitoba_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fis", logTransform = FALSE, ymin = min(admixed_perloc$Fis, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS])) + ggtitle("Lake Manitoba")


canada_perloc <- cbind(canada_perloc, vcf_hf@fix)
ggman(gwas = dplyr::filter(canada_perloc, CHROM == "NC_041338.1"), snp = "ID", bp = "POS", chrom = "CHROM", pvalue = "Fstp", logTransform = FALSE, ymin = min(canada_perloc$Fstp, na.rm = TRUE), ymax = 1, relative.positions = TRUE, sigLine = NA, lineColour = "#fdae61", title = "", ylabel = 'Fis', xlabel = "Chromosome") + theme_bw() + ylab(bquote(italic(F)[IS])) + ggtitle("Lake Manitoba")
