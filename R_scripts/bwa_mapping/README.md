This files were used to analyze mapping percentages of walleye RAPTURE data to the yellow perch genome by bwa. 

Samtools flagstat (flagstat.sh) was used in an array to get mapping percentages.

The outputs from samtools flagstat are saved in slurm_out. 

The slurm outputs were analyzed with the R script bwa_mapping.R, which looped through all the files and calculated a mean and standard deviation for mapping percentages.

A mean 97.8% of walleye reads were mapped to the yellow perch genome succesfully, with 0.11% standard deviation.
