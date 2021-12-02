library(plinkQC)
library(tidyverse)

# Check for relatedness using the plink tool
relatedness_QC <- evaluate_check_relatedness(
  qcdir=".",
  name="HDplot_pruned",
  highIBDTh = 0.1875,
  imissTh = 0.03,
  interactive = FALSE,
  verbose = TRUE
)


evaluate_check_relatedness(
  qcdir=".",
  name="HDplot_pruned",
  highIBDTh = 0.1875,
  imissTh = 0.03,
  interactive = TRUE,
  verbose = TRUE
)

plink.genome <- read_delim("HDplot_pruned.genome", delim = " ")

# Also check for relatedness using the method of moments in SNPRelate/Plink
# Code to install SNPRelate
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# BiocManager::install("gdsfmt")
# BiocManager::install("SeqArray")
library(gdsfmt)
library(SNPRelate)

vcf.fn <- "HDplot_singletons_pruned.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "wall.gds", method="biallelic.only")

snpgdsSummary("wall.gds")

genofile <- SNPRelate::snpgdsOpen("wall.gds")

# Getting a list of the walleye in the SNP relate dataset
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Using the LD-pruned SNPs for kinship analysis using the method of moments (MoM)
# Based on the guide provided in the vignette here: https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#estimating-ibd-using-plink-method-of-moments-mom
ibd <- snpgdsIBDMoM(genofile, sample.id=NULL, snp.id=NULL,
                    maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE, kinship = TRUE)
# Make a data.frame. This contains pairwise relatedness measures by the method of moments.
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

snpgdsClose("wall.gds")
