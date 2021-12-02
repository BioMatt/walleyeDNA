library(tidyverse)
library(pophelper)
metadata <- read.table("combined_popmap2.txt")

admix_files <- list.files("admixture_outputs//Q_only/", full.names = T)
admix <- readQ(filetype = "auto", files = admix_files)

# Re-order the individuals in each data frame by site, so plotQ can show individuals and sites in the order I want
admix_ordered <- lapply(admix, function(df){
        df[order(metadata$V2),]
})

# Get some summary tables
tabulated <- tabulateQ(admix)
summary <- summariseQ(tabulated)

# Re-order the site metadata in exactly the same way the individuals were reordered
sites <- metadata[order(metadata$V2),2, drop=FALSE]

colnames(sites) <- "Site"
# Plot the barplots
plotQ(admix_ordered[c(2:6)], imgoutput = "join", grplab = sites, sortind = "all", subsetgrp = c("Red_River", "Matheson", "Dauphin_River", "Swan_Creek", "Kansas"), sharedindlab = FALSE, outputfilename = "Admixture_K2-6.png", dpi = 900, showtitle = TRUE, titlelab = "Admixture from K=2-6", splab = c("K=2", "K=3", "K=4", "K=5", "K=6"))


############################################################################################
# Assuming K=3 is the most biologically meaningful value, pulling out ancestry estimates for that group
k3 <- admix$HDplot_pruned.3.Q
# Using case_when to identify the column with the highest ancestry estimates
k3 <- k3 %>% 
        mutate(ancestry = case_when(
                Cluster1 >= 0.85 ~ "Kansas",
                Cluster2 >= 0.85 ~ "Winnipeg",
                Cluster3 >= 0.85 ~ "Manitoba",
                Cluster1 < 0.85 & Cluster1 > 0.5 ~ "Kansas_admixed",
                TRUE ~ "Admixed"
        ))

# By combining the metadata with clusters and ancestry estimates, we see that cluster 1 corresponds to Kansas fish, cluster 2 to Swan Creek/Dauphin River, cluster 3 to Red River/Matheson/Dauphin River (south Lake Winnipeg), and admixed individuals are in the Lake Winnipeg/Lake Manitoba system
admixed_metadata <- cbind(metadata, k3)
admixed_metadata <- admixed_metadata[order(admixed_metadata$ancestry),]



# Write out the four groups for extraction by Plink from the whole dataset, and separation into individual chromosomes/contigs. These will be the input files for Elai
# Because Plink requires family and individual information on each line, we write out each individual twice
# Write out the 49 admixed invididuals for filtering out of the dataset, for use with dadi
admixed_indivs <- dplyr::filter(admixed_metadata, ancestry == "Admixed")
write_delim(as.data.frame(paste0(admixed_indivs$V1, "_", admixed_indivs$V1)), "admixed_indivs.txt", delim = "\t", col_names = FALSE)


kansas_indivs <- dplyr::filter(admixed_metadata, ancestry == "Cluster1")
write_delim(as.data.frame(cbind(kansas_indivs$V1, kansas_indivs$V1)), "Kansas_indivs.txt", delim = "\t", col_names = FALSE)

LkManitoba_indivs <- dplyr::filter(admixed_metadata, ancestry == "Cluster2")
write_delim(as.data.frame(cbind(LkManitoba_indivs$V1, LkManitoba_indivs$V1)), "LkManitoba_indivs.txt", delim = "\t", col_names = FALSE)

LkMWinnipeg_indivs <- dplyr::filter(admixed_metadata, ancestry == "Cluster3")
write_delim(as.data.frame(cbind(LkMWinnipeg_indivs$V1, LkMWinnipeg_indivs$V1)), "LkWinnipeg_indivs.txt", delim = "\t", col_names = FALSE)

# Break down cluster 3 into the various Lake Winnipeg components- Red River, Matheson, Dauphin River
red_river <- dplyr::filter(admixed_metadata, ancestry == "Cluster3", V2 == "Red_River")
write_delim(as.data.frame(cbind(red_river$V1, red_river$V1)), "RedRiver_indivs.txt", delim = "\t", col_names = FALSE)

matheson <- dplyr::filter(admixed_metadata, ancestry == "Cluster3", V2 == "Matheson")
write_delim(as.data.frame(cbind(matheson$V1, matheson$V1)), "Matheson_indivs.txt", delim = "\t", col_names = FALSE)

dauphin <- dplyr::filter(admixed_metadata, ancestry == "Cluster3", V2 == "Dauphin_River")
write_delim(as.data.frame(cbind(dauphin$V1, dauphin$V1)), "Dauphin_indivs.txt", delim = "\t", col_names = FALSE)



# Do the same thing but write out the different individuals only for the different populations, with just one individual per line and no other information
kansas_indivs <- dplyr::filter(admixed_metadata, ancestry == "Kansas")
write_delim(as.data.frame(kansas_indivs$V1), "Kansas_indivs.txt", delim = "\t", col_names = FALSE)

LkManitoba_indivs <- dplyr::filter(admixed_metadata, ancestry == "Manitoba")
write_delim(as.data.frame(LkManitoba_indivs$V1), "LkManitoba_indivs.txt", delim = "\t", col_names = FALSE)

LkMWinnipeg_indivs <- dplyr::filter(admixed_metadata, ancestry == "Winnipeg")
write_delim(as.data.frame(LkMWinnipeg_indivs$V1), "LkWinnipeg_indivs.txt", delim = "\t", col_names = FALSE)

# Take the k=3 table and filter out admixed individuals, for use with population assignments in dadi
dadi_metadata <- filter(admixed_metadata, ancestry != "Admixed")
dadi_metadata <- dadi_metadata[order(dadi_metadata$ancestry),]
# Add in a new individual ID column to match the plink-adjusted names
dadi_metadata$plink_IDs <- paste0(dadi_metadata$V1, "_", dadi_metadata$V1)

write.table(cbind(as.character(dadi_metadata$plink_IDs), dadi_metadata$ancestry), "dadi_assignments.txt", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(cbind(as.character(dadi_metadata$V1), dadi_metadata$ancestry), "dadi_assignments_noplink.txt", quote = F, sep = "\t", col.names = F, row.names = F)

# Just in case, re-order the VCF with the non-admixed indivs using a table of just their names
write.table(cbind(as.character(dadi_metadata$plink_IDs)), "dadi_indivs.txt", quote = F, sep = "\t", col.names = F, row.names = F)
write.table(cbind(as.character(admixed_metadata$V1), admixed_metadata$ancestry), "dadi_assignments_Dec31.txt", quote = F, sep = "\t", col.names = F, row.names = F)
# Also write out a table of assignments with no admixed individuals for diveRsity
write.table(cbind(as.character(dadi_metadata$V1), dadi_metadata$ancestry), "dadi_assignments_Dec31_noAdmixed.txt", quote = F, sep = "\t", col.names = F, row.names = F)

# Write out all the admixed metadata, possibly to include as covariates in movement associations
write.table(admixed_metadata, "k3_ancestry_metadata.txt", quote = F, sep = "\t", row.names = F)
