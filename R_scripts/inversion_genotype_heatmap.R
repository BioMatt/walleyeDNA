# Script based on code written by Yue Shi

library(vcfR)
library(pheatmap)
library(tidyverse)
#####################################
####Get genotype and make heatmap####
#####################################
#load vcf file with only 111 topFstp loci
ersn.topfst.vcf<-read.vcfR("chr8_region.recode.vcf")

#extract genotypes in numeric format
ersn.topfst_genotypes_numeric<-extract.gt(ersn.topfst.vcf, element = "GT", mask = FALSE, as.numeric = FALSE,
                                          return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
                                          convertNA = TRUE)

#function to convert genotypes to 0,1,2 format used for heatmap
convertGenos<-function(genotype){
  alleles<-as.numeric(str_split_fixed(genotype,"/",2))
  genoCode<-alleles[1]+alleles[2]
  return(genoCode)
}

#convert to 0,1,2 format
ersn.topfst_genosCode<-apply(ersn.topfst_genotypes_numeric,1:2,convertGenos)
ersn.topfst_genosCode<-t(ersn.topfst_genosCode)

#format for heatmap2
ersn.topfst_heatmap2<-ersn.topfst_genosCode
ersn.topfst_heatmap2[is.na(ersn.topfst_heatmap2)] <- -1 #replace missing info with -1

# Read in population assignments
assignments <- read_tsv("dadi_assignments_Dec31.txt", col_names = c("Ind", "Assignment"))

# Read in metadata
metadata <- read_tsv("combined_metadata.txt", col_names = c("Ind", "Site", "Waterbody")) %>% 
  left_join(., assignments) %>% 
  arrange(desc(Waterbody), desc(Site), desc(Assignment))


ersn.topfst_heatmap2 <- dplyr::as_tibble(ersn.topfst_heatmap2, rownames = NA)
ersn.topfst_heatmap2 <- ersn.topfst_heatmap2 %>% 
  rownames_to_column(var = "Ind")
ersn.topfst_heatmap2 <- left_join(metadata, ersn.topfst_heatmap2) %>% 
  select(-Site, -Waterbody, -Assignment) %>% 
  column_to_rownames(var = "Ind")
ersn.topfst_heatmap2<-apply(ersn.topfst_heatmap2,1:2,function(x) as.numeric(x))




metadata <- column_to_rownames(metadata, var = "Ind")

genotype_heatmap <- pheatmap(ersn.topfst_heatmap2,
         annotation_row = metadata,
         #annotation_colors  = class_colors,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         show_rownames = FALSE,
         legend_breaks = c(-1,0,1,2),
         legend_labels =c("NA","homo1","het","homo2"), cluster_rows = FALSE, cluster_cols = FALSE)
genotype_heatmap
ggsave(filename = "chr8_heatmap.pdf", plot = genotype_heatmap, dpi = 2000)
