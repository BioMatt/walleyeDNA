# using PCadapt to detect local adaptation. Based on the vignette at https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

library(pcadapt)
library(tidyverse)
library(vcfR)
library(qvalue)

# Read in the data
vcf <- read.pcadapt("HDplot_singletons.recode.vcf", type = "vcf")
vcf_data <- read.vcfR("HDplot_singletons.recode.vcf")
snp_data <- as_tibble(vcf_data@fix)
individual_order <- tibble(Ind = colnames(vcf_data@gt)[-1])

# Run a PCA
x <- pcadapt(input = vcf, K = 20)

# Save the scree plot with 20 PCs for supplemental figures
plot(x, option = "screeplot")

# Read the metadata for all 345 individuals, then match it up to the individuals included in the PCadapt analysis. Re-order the metadata to match the VCF individual order
metadata <- read_delim("combined_popmap2.txt", col_names = c("Ind", "Site"), delim = "\t")
metadata <- left_join(individual_order, metadata)

# Plotting the first two axes, then axes 3 and 4
x$singular.values
plot(x, option = "scores", pop = metadata$Site)
plot(x, option = "scores", i = 3, j = 4, pop = metadata$Site)

# Re-run PC Adapt with only 2 PCs
x <- pcadapt(vcf, K = 2)
summary(x)

plot(x, option = "scores", pop = metadata$Site)

plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")


# Using FDR instead of p values
qval <- qvalue(x$pvalues)$qvalues
# Choose alpha to be 0.05
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)

summary(x)

# Get PCs associated with different markers
which_pc <- get.pc(x, outliers)

# Attach the PC and outlier info together
snp_data$SNP <- seq(1, nrow(snp_data))
snp_data <- left_join(snp_data, which_pc)
snp_data$log_Q <- -log10(qval)
snp_data$Q <- qval
snp_data$PC <- as.factor(snp_data$PC)

snp_data$log_Q <- replace_na(snp_data$log_Q, 0)
snp_data$Q <- replace_na(snp_data$Q, 1)

# Trying to make my own Manhattan plot that I can colour
ggplot(data = snp_data, aes(x = SNP, y = log_Q, colour = PC, group = CHROM)) +
  geom_point(alpha = 0.8)

facet_wrapped_manhattan <- ggplot(data = snp_data, aes(x = POS, y = log_Q, colour = PC, group = CHROM)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ CHROM)


# Showing one or the other PC at a time
ggplot(data = filter(snp_data, PC == 1), aes(x = SNP, y = log_Q, group = CHROM)) +
  geom_point()
ggplot(data = filter(snp_data, PC == 2), aes(x = SNP, y = log_Q, group = CHROM)) +
  geom_point()

# Taking a closer look at chromosome 38
chr38 <- filter(snp_data, CHROM == "NC_041338.1") %>% 
  mutate(POS = as.numeric(chr38$POS))
# Trying to make my own Manhattan plot that I can colour
ggplot(data = chr38, aes(x = POS, y = log_Q, colour = PC, group = PC)) +
  geom_point() #+
  geom_rect(inherit.aes=TRUE, aes(xmin=15722749, xmax=15845548, ymin=min(log_Q),
                                   ymax=max(log_Q)), color="transparent", fill="orange", alpha=0.1)


# Check out the chromosome 38 SNPs that are significant
sig_chr_38 <- filter(chr38, Q < 0.05)

# Write out the PC Adapt results
write_delim(snp_data, "PCAdapt_SNPs.txt")

######################################
# Looking at plots of individual chromosomes
chr54 <- filter(snp_data, CHROM == "NC_041354.1")
# Trying to make my own Manhattan plot that I can colour
ggplot(data = chr54, aes(x = POS, y = log_Q, colour = PC, group = PC)) +
  geom_point()

chr32 <- filter(snp_data, CHROM == "NC_041332.1") %>% 
  mutate(POS = as.numeric(chr32$POS))
# Trying to make my own Manhattan plot that I can colour
ggplot(data = chr32, aes(x = POS, y = log_Q, colour = PC, group = PC)) +
  geom_point() +
  geom_rect(inherit.aes=TRUE, aes(xmin=15965078, xmax=16067633, ymin=min(log_Q),
                                               ymax=max(log_Q), group=PC), color="transparent", fill="orange", alpha=0.3)


chr36 <- filter(snp_data, CHROM == "NC_041346.1")
chr36 <-  mutate(chr36, POS = as.numeric(chr36$POS))
# Trying to make my own Manhattan plot that I can colour
ggplot(data = chr36, aes(x = POS, y = log_Q, colour = PC, group = PC)) +
  geom_point() +
  geom_rect(inherit.aes=TRUE, aes(xmin=19999287, xmax=20158372, ymin=min(log_Q),
                                  ymax=max(log_Q), group=PC), color="transparent", fill="orange", alpha=0.3)

###############################################
library(ggman)
library(patchwork)
# Create a PCA and Manhattan plot for manuscript figures

# Use the ggman package to plot a Manhattan plot with relative chromosome and SNP positions intact
# Try turning Q values and correlations positive or negative depending on which direction they go. Commented out this line because I will go with just separating the plots, instead
# Also remove scaffolds (NW_) for nicer plotting
snps_plot_data <- snp_data %>% 
  #mutate(log_Q = ifelse(PC == 2, (log_Q * -1), log_Q)) %>% 
  mutate(log_Q = ifelse(is.na(log_Q), 0, log_Q), CHROM_ID = as.numeric(as.factor(CHROM))) %>% 
  filter(CHROM_ID <= 24)

# Filter out PC1 and 2 info from the overall outlier test
sig_PCAdapt_PC1 <- filter(snps_plot_data, PC == 1 & Q < 0.05)
all_PCAdapt_PC1 <- filter(snps_plot_data, PC == 1 | is.na(PC))
sig_PCAdapt_PC2 <- filter(snps_plot_data, PC == 2 & Q < 0.05)
all_PCAdapt_PC2 <- filter(snps_plot_data, PC == 2 | is.na(PC))

# This creates a combined Manhattan plot with negative values for PC2 (using the commented out mutate line above), but it doesn't look great
# pcadapt_plot <- ggman(gwas = snps_plot_data, snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "log_Q", logTransform = FALSE, ymin = -60, ymax = 15, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() 
# # The colour #d73027 is red, representing PC1 SNPs
# pcadapt_plot_PC1 <- ggmanHighlight(pcadapt_plot, highlight = sig_PCAdapt_PC1$ID, colour = "#d73027")
# # The colour #4575b4 is blue, representing PC2 SNPs
# pcadapt_plot_PC1_2 <- ggmanHighlight(pcadapt_plot_PC1, highlight = sig_PCAdapt_PC2$ID, colour = "#4575b4") + ylab(bquote('-log'[10]~italic(q)))
# pcadapt_plot_PC1_2 

# Create separate Manhattan plots for PC1 and PC2. First, PC1
pcadapt_plot_PC1only <- ggman(gwas = all_PCAdapt_PC1, snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "log_Q", logTransform = FALSE, ymin = 0, ymax = 15, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() + ylab(bquote(PC1~'-log'[10]~italic(q)))
pcadapt_plot_PC1only <- ggmanHighlight(pcadapt_plot_PC1only, highlight = sig_PCAdapt_PC1$ID, colour = "#d73027", size = 0.7)
pcadapt_plot_PC1only

pcadapt_plot_PC2only <- ggman(gwas = all_PCAdapt_PC2, snp = "ID", bp = "POS", chrom = "CHROM_ID", pvalue = "log_Q", logTransform = FALSE, ymin = 0, ymax = 60, relative.positions = TRUE, sigLine = NA, title = "", ylabel = 'q', xlabel = "Chromosome") + theme_bw() + ylab(bquote(PC2~'-log'[10]~italic(q)))
pcadapt_plot_PC2only <- ggmanHighlight(pcadapt_plot_PC2only, highlight = sig_PCAdapt_PC2$ID, colour = "#4575b4", size = 0.7)
pcadapt_plot_PC2only

# Read in population assignments to see if assignments as shapes look good on the PCA figure
assignments <- read_delim("dadi_assignments_Dec31.txt", delim = "\t", col_names = c("Ind", "Assignment"))
metadata <- left_join(metadata, assignments)
# Create a plot of the PCDapt PCA results
# Pull the variance explained by each of PC1 and PC2
(x$singular.values^2)*100
pcadapt_pca <- ggplot(as.data.frame(x$scores), aes(x=x$scores[,1], y=x$scores[,2], colour = factor(metadata$Site), shape = factor(metadata$Assignment))) + 
  geom_point(size = 2, alpha= 0.80) + 
  labs(x = "PC1 (5.76% variance)", y = "PC2 (2.28% variance)") +
  labs(color = "Site Collected", shape = "Population Assignments") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 1),
        text = element_text(size=15)
  ) + 
  scale_colour_manual(values = c("#fdae61", "#d73027", "#abd9e9", "#2c7bb6", "#313695"),
                      breaks = c("Kansas", "Swan_Creek", "Dauphin_River", "Matheson", "Red_River"),
                      labels = c("Cedar Bluff Reservoir", "Swan Creek Hatchery", "Dauphin River", "Matheson Island", "Red River")) +
  scale_shape_manual(breaks = c("Admixed", "Kansas", "Manitoba", "Winnipeg"),
                     labels = c("Unassigned", "Kansas", "Lake Manitoba", "Lake Winnipeg"), 
                     values = c(15, 17, 1, 16)) +
  guides(shape = guide_legend(order = 2),col = guide_legend(order = 1))

pcadapt_pca

# Combine the two plots to be next to each other. I never got this to quite look very good
#combined_plots <- (pcadapt_pca + plot_spacer()) / pcadapt_plot_PC1_2 + plot_layout(widths = c(1, 3))
#combined_plots
#ggsave("combined_PCAdapt_plots.pdf", combined_plots, dpi = 2000)

# Try combining the PCA and the two Manhattan plots. This one looks pretty good.
combined_PC_separate <- (pcadapt_pca+ggtitle("Principal Components Analysis")) | ((pcadapt_plot_PC1only+ggtitle("Significant PC1 SNPs")) / (pcadapt_plot_PC2only+ggtitle("Significant PC2 SNPs")))
#combined_PC_separate
ggsave("combined_PCAdapt_plots_separatePCs.pdf", combined_PC_separate, dpi = 2000, height = 6, width = 10)

# Write each plot individually
ggsave("PCAdapt_PCA.pdf", pcadapt_pca, dpi = 2000)
ggsave("PCAdapt_Manhattan.pdf", pcadapt_plot_PC1_2, dpi = 2000)
ggsave("PC1only_Manhattan.pdf", pcadapt_plot_PC1only, dpi = 2000)
ggsave("PC2only_Manhattan.pdf", pcadapt_plot_PC2only, dpi = 2000)

# Write out chromosome codes and IDs for Yue
chroms <- snp_data %>% 
  mutate(log_Q = ifelse(is.na(log_Q), 0, log_Q), CHROM_ID = as.numeric(as.factor(CHROM))) %>% 
  select(CHROM, CHROM_ID) %>% 
  distinct(CHROM_ID, .keep_all = TRUE)
write_tsv(x = chroms, file = "chrom_IDs.txt")

##################################################
# Check the distribution of outliers along PC2
View(sig_PCAdapt_PC2 %>% arrange(desc(log_Q)))
