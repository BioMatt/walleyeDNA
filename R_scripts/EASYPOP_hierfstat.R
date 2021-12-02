library(tidyverse)
library(vcfR)
library(adegenet)
library(hierfstat)

# A for loop in a function to read in the different replicate dataframes and calculate Fst for each
fst_loop <- function(iterations, base_name) {
  fst_list <- c() # this is an empty vector of Fst values
  for (i in 1:iterations) {
    # Read the fstat, taking in a base name, the replicate for this run, and the .dat file extension. The sprintf function makes the replicate 3-digit to match the easypop replicate naming scheme
    temp_fstat <- hierfstat::read.fstat(paste0(base_name, sprintf('%0.3d', 1:iterations)[i], ".dat"))
    temp_fst <- pairwise.WCfst(temp_fstat, diploid = TRUE) # Calculate Fst
    print(paste0("Replicate #", i," Fst: ", temp_fst[1,2])) # Print Fst to keep track of function progress
    fst_list[i] <- temp_fst[1,2] # Add the Fst value to the list of Fst values for this run
  }
  assign(base_name, fst_list, envir = .GlobalEnv) # Assign the list of Fst values to the base name in the global R environment. This will come out as a vector
}

# Run the function on the 0.01 migration model
fst_loop(100, "1kSNPs_0.01mig_Ne")
# Check mean and standard deviation
mean(`1kSNPs_0.01mig_Ne`)
sd(`1kSNPs_0.01mig_Ne`)

# Run the function with 0.001 migration model
fst_loop(100, "1kSNPs_0.001mig_Ne")
mean(`1kSNPs_0.001mig_Ne`)
sd(`1kSNPs_0.001mig_Ne`)

# Run the function for the 0.0001 migration model
fst_loop(100, "1kSNPs_0.0001migNe")

# Run the function with 0 migration model
fst_loop(100, "1kSNPs_0mig_Ne")
mean(`1kSNPs_0mig_Ne`)

# Run the function with 0.03 migration model
fst_loop(100, "1kSNPs_0.03mig_Ne")
mean(`1kSNPs_0.03mig_Ne`)

# Run the function with 0.02 migration model
fst_loop(100, "1kSNPs_0.02mig_Ne")
mean(`1kSNPs_0.02mig_Ne`)

# Run the function with 0.02 migration model
fst_loop(100, "1kSNPs_0.005mig_Ne")
mean(`1kSNPs_0.005mig_Ne`)

fst_loop(100, "1kSNPs_0.002mig_Ne")
mean(`1kSNPs_0.002mig_Ne`)

fst_loop(100, "1kSNPs_0.003mig_Ne")
mean(`1kSNPs_0.003mig_Ne`)

fst_loop(100, "1kSNPs_0.004mig_Ne")
mean(`1kSNPs_0.004mig_Ne`)

# Create a tibble of combined Fst values
#mig_0 = `1kSNPs_0mig_Ne`,  # Removed the 0-migration simulated data from plotting because Fst was much too high
combined_fst <- tibble(mig_0 = `1kSNPs_0mig_Ne`, mig_0.0001 = `1kSNPs_0.0001migNe`, mig_0.001 = `1kSNPs_0.001mig_Ne`, mig_0.002 = `1kSNPs_0.002mig_Ne`, mig_0.003 = `1kSNPs_0.003mig_Ne`, mig_0.004 = `1kSNPs_0.004mig_Ne`, mig_0.005 = `1kSNPs_0.005mig_Ne`, mig_0.01 = `1kSNPs_0.01mig_Ne`, mig_0.02 = `1kSNPs_0.02mig_Ne`, mig_0.03 = `1kSNPs_0.03mig_Ne`) %>% 
  gather()

easypop_fst <- ggplot(combined_fst, aes(x = key, y = value, group = key)) + 
  geom_boxplot(notch = FALSE) +
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.4) +
  geom_hline(yintercept = 0.027) +
  geom_hline(yintercept = 0.026, linetype = 2) +
  geom_hline(yintercept = 0.029, linetype = 2) +
  scale_x_discrete(labels = c("0", "0.0001", "0.001","0.002", "0.003", "0.004", "0.005", "0.01", "0.02", "0.03")) +
  ylab(bquote(italic(F)[ST])) +
  xlab("Migration Rate") +
  theme_bw() +
  theme(text = element_text(size=18))
easypop_fst
ggsave(filename = "easypop_fst.pdf", plot = easypop_fst, dpi = 2000)
