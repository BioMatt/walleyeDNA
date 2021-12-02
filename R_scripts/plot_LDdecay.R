# An R script to plot LD decay results from PopLDdecay
# Not used in the final manuscript
library(tidyverse)
library(patchwork)


# Read in Lake Winnipeg LD
winnipeg_ld <- read_tsv("winnipeg_LDdecay.stat") %>% 
  dplyr::rename(dist_kb = "#Dist", mean_rsq = "Mean_r^2", sum_rsq = "Sum_r^2", num_pairs = "NumberPairs") %>% 
  select(-`Mean_D'`, -`Sum_D'`)
ggplot(data = winnipeg_ld, aes(x = dist_kb, y = mean_rsq)) +
  geom_point() +
  theme_bw()




# Read in Lake Manitoba LD
manitoba_ld <- read_tsv("manitoba_LDdecay.stat") %>% 
  dplyr::rename(dist_kb = "#Dist", mean_rsq = "Mean_r^2", sum_rsq = "Sum_r^2", num_pairs = "NumberPairs") %>% 
  select(-`Mean_D'`, -`Sum_D'`)
ggplot(data = manitoba_ld, aes(x = dist_kb, y = mean_rsq)) +
  geom_point() +
  theme_bw()


# Read in Cedar Bluff LD
kansas_ld <- read_tsv("kansas_LDdecay.stat") %>% 
  dplyr::rename(dist_kb = "#Dist", mean_rsq = "Mean_r^2", sum_rsq = "Sum_r^2", num_pairs = "NumberPairs") %>% 
  select(-`Mean_D'`, -`Sum_D'`)

ggplot(data = kansas_ld, aes(x = dist_kb, y = mean_rsq)) +
  geom_point() +
  theme_bw()



# Trying pairwise LD plots in Cedar Bluff
# For using non-linear regression to look at estimated linkage decay. Use the pikeperch genome and linkage map (2725.54 cM linkage map over 896,481,186 bp of placed scaffolds in genome) from Rios-Perez et al. 2020. There is thus 896481196 bp per 2725.54 centimorgan, thus 328918.7 bp/centimorgan on average. Thus the distance between SNPs is converted to centimorgans by distance * 1 cm/328918.7 bp, or distance/328918.7
# For recombination fraction r, use the equation r = 0.5 * (1 - e ^ (-2 * d)), where d is distance in centimorgans, from https://jvanderw.une.edu.au/05_Basics_of_linkagemapping.PDF
# C = 4 * Ne * recombination fraction, where Ne was calculated with Ne Estimator v2
# *** NOTE *** These steps are unnecessary, and possibly even misleading with all the assumptions made. nls() can be used to estimate C in addition to expected r^2 from the data
 # pairwise_kansas <- read_tsv("kansas_LDdecay.LD") %>% 
 #   dplyr::rename(rsq = `r^2`, chr = `#chr`) %>% 
 #   mutate(recomb_fraction = (0.5 * (1 - (2.71828 ^ (-2 * Dist/328918.7))))) %>% 
 #   mutate(C = 4 * 762 * recomb_fraction, uncorrected_exp_rsq = (1/(1 + C)), expected_rsq = (((10 + C)/((2 + C)*(11 + C))) * (1 + ((3 + C)*(12 + (12*C) + (C^2)))/(60*(2 + C)*(11 + C)))))


# A much simpler funciton for looking at pairwise r^2 in Cedar Bluff
pairwise_kansas <- read_tsv("kansas_LDdecay.LD") %>% 
  dplyr::rename(rsq = `r^2`, chr = `#chr`) %>% 
  mutate(pop = rep("kansas", nrow(.)))




# Using nonlinear least squares to plot LD decay across increasing distances between SNPs. Code borrowed from https://www.r-bloggers.com/2011/08/estimate-decay-of-linkage-disequilibrium-with-distance/
# Run the NLS model, following the formula for expected r^2
kansas_nls <- nls(formula = (rsq ~ ((10+C*Dist)/((2+C*Dist)*(11+C*Dist)))*(1+((3+C*Dist)*(12+12*C*Dist+(C*Dist)^2))/(60*(2+C*Dist)*(11+C*Dist)))), start=c(C = 0.1), control=nls.control(maxiter=100), data = pairwise_kansas)
# Estimate the population recombination parameter, which is C/distance
kansas_rho <- summary(kansas_nls)$parameters[1]
# Pull the distance values for each row
kansas_dist <- c(pairwise_kansas$Dist)
# Use the population recombination parameter rho with the distance values to estimate expected r^2 for each distance, controlling for sample size (60 here). Stick that info onto the main LD table for kansas
n_kansas <- 60
pairwise_kansas$fpoints<-((10+kansas_rho*kansas_dist)/((2+kansas_rho*kansas_dist)*(11+kansas_rho*kansas_dist)))*(1+((3+kansas_rho*kansas_dist)*(12+12*kansas_rho*kansas_dist+(kansas_rho*kansas_dist)^2))/(n_kansas*(2+kansas_rho*kansas_dist)*(11+kansas_rho*kansas_dist)))

# Check where LD decays below 0.5 (strong LD), where it decays to 0.2 (moderate LD), and where it decays to half its maximum value
# Half decay and strong LD terms end up not being informative, likely because of the low coverage data here. Moderate LD seems informative though
kansas_half_decay <- pairwise_kansas$Dist[which.min(abs(pairwise_kansas$fpoints - max(pairwise_kansas$fpoints/2)))]
kansas_strongLD <- min(pairwise_kansas$Dist[pairwise_kansas$fpoints <= 0.50])
kansas_moderateLD <- min(pairwise_kansas$Dist[pairwise_kansas$fpoints <= 0.20])

# Take a look at the plot
ggplot(pairwise_kansas, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x="Distance (Base Pairs)", y=expression(Pairwise~italic(r)^{2}), title = "Cedar Bluff Reservoir") +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  xlim(0, 50000) +
  ylim(0, 0.5) +
  theme_bw()

# Plot pairwise LD between points along with the model.
kansas_pairwise_plot <- ggplot(pairwise_kansas, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x=element_blank(), y=expression(Pairwise~italic(r)^{2}), title = "Cedar Bluff Reservoir") +
  geom_hline(yintercept = 0.5) + 
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  theme_bw()
kansas_pairwise_plot


# Plotting script adapted from user rmf's response in the forum thread https://www.biostars.org/p/300381/
# Group data into 20 kb intervals
pairwise_kansas$distc <- cut(pairwise_kansas$Dist,breaks=seq(from=min(pairwise_kansas$Dist)-1,to=max(pairwise_kansas$Dist)+1,by=1000))


# Compute mean r2 within each 20 kb interval
pairwise_kansas_grouped <- pairwise_kansas %>% 
  group_by(distc) %>% 
  summarise(mean_20kb = mean(rsq),median_20kb = median(rsq))

# A helper step to get mid points of our distance intervals for plotting.
pairwise_kansas_grouped <- pairwise_kansas_grouped %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                                  end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                                  mid=start+((end-start)/2))
# Add in a column that just describes the population 
pairwise_kansas_grouped$pop <- rep("kansas", nrow(pairwise_kansas_grouped))


ggplot()+
  geom_point(data = pairwise_kansas_grouped, aes(x = start, y = mean_20kb), size = 0.4,colour="grey20") +
  geom_line(data = pairwise_kansas_grouped, aes(x = start, y = mean_20kb), size=0.3, alpha=0.5, colour="grey40") +
  labs(x="Distance (Kilobases)", y=expression(LD~(r^{2})))+
  theme_bw()



# Trying pairwise LD plots in Lake Manitoba
# Add in a column that just describes the population 
pairwise_manitoba <- read_tsv("manitoba_LDdecay.LD") %>% 
  dplyr::rename(rsq = `r^2`, chr = `#chr`) %>% 
  mutate(pop = rep("manitoba", nrow(.)))


#ggsave("manitoba_pairwiseLD.png", manitoba_pairwise_plot, dpi = 300)


# Using nonlinear least squares to plot LD decay across increasing distances between SNPs. Code borrowed from https://www.r-bloggers.com/2011/08/estimate-decay-of-linkage-disequilibrium-with-distance/
# Run the NLS model, following the formula for expected r^2
manitoba_nls <- nls(formula = (rsq ~ ((10+C*Dist)/((2+C*Dist)*(11+C*Dist)))*(1+((3+C*Dist)*(12+12*C*Dist+(C*Dist)^2))/(60*(2+C*Dist)*(11+C*Dist)))), start=c(C = 0.1), control=nls.control(maxiter=100), data = pairwise_manitoba)
# Estimate the population recombination parameter, which is C/distance
manitoba_rho <- summary(manitoba_nls)$parameters[1]
# Pull the distance values for each row
manitoba_dist <- c(pairwise_manitoba$Dist)
# Use the population recombination parameter rho with the distance values to estimate expected r^2 for each distance, controlling for sample size (60 here). Stick that info onto the main LD table for kansas
n_manitoba <- 50
pairwise_manitoba$fpoints<-((10+manitoba_rho*manitoba_dist)/((2+manitoba_rho*manitoba_dist)*(11+manitoba_rho*manitoba_dist)))*(1+((3+manitoba_rho*manitoba_dist)*(12+12*manitoba_rho*manitoba_dist+(manitoba_rho*manitoba_dist)^2))/(n_manitoba*(2+manitoba_rho*manitoba_dist)*(11+manitoba_rho*manitoba_dist)))

# Take a look at the plot
ggplot(pairwise_manitoba, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x="Distance (Base Pairs)", y=expression(Pairwise~italic(r)^{2}), title = "Lake Manitoba") +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  xlim(0, 50000) +
  ylim(0, 0.5) +
  theme_bw()

# Plot pairwise LD between points along with the model.
manitoba_pairwise_plot <- ggplot(pairwise_manitoba, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x="Distance (Base Pairs)", y=expression(Pairwise~italic(r)^{2}), title = "Lake Manitoba") +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  geom_hline(yintercept = 0.5) + 
  theme_bw()

# Check where LD decays below 0.5 (strong LD), where it decays to 0.2 (moderate LD), and where it decays to half its maximum value
# Half decay and strong LD terms end up not being informative, likely because of the low coverage data here. Moderate LD seems informative though
manitoba_half_decay <- pairwise_manitoba$Dist[which.min(abs(pairwise_manitoba$fpoints - max(pairwise_manitoba$fpoints/2)))]
manitoba_strongLD <- min(pairwise_manitoba$Dist[pairwise_manitoba$fpoints <= 0.50])
manitoba_moderateLD <- min(pairwise_manitoba$Dist[pairwise_manitoba$fpoints <= 0.20])

# Plotting script adapted from user rmf's response in the forum thread https://www.biostars.org/p/300381/
# Group data into 20 kb intervals
pairwise_manitoba$distc <- cut(pairwise_manitoba$Dist,breaks=seq(from=min(pairwise_manitoba$Dist)-1,to=max(pairwise_manitoba$Dist)+1,by=1000))

# Compute mean r2 within each 20 kb interval
pairwise_manitoba_grouped <- pairwise_manitoba %>% 
  group_by(distc) %>% 
  summarise(mean_20kb = mean(rsq),median_20kb = median(rsq))

# A helper step to get mid points of our distance intervals for plotting.
pairwise_manitoba_grouped <- pairwise_manitoba_grouped %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                                  end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                                  mid=start+((end-start)/2))
# Add in a column that just describes the population 
pairwise_manitoba_grouped$pop <- rep("manitoba", nrow(pairwise_manitoba_grouped))


ggplot()+
  geom_point(data = pairwise_manitoba_grouped, aes(x = start, y = mean_20kb), size = 0.4,colour="grey20") +
  geom_line(data = pairwise_manitoba_grouped, aes(x = start, y = mean_20kb), size=0.3, alpha=0.5, colour="grey40") +
  labs(x="Distance (Kilobases)", y=expression(LD~(r^{2})))+
  theme_bw()

# Trying pairwise LD plots in Lake Winnipeg
# Add in a column that just describes the population 
pairwise_winnipeg <- read_tsv("winnipeg_LDdecay.LD") %>% 
  dplyr::rename(rsq = `r^2`, chr = `#chr`) %>% 
  mutate(pop = rep("winnipeg", nrow(.)))



# Using nonlinear least squares to plot LD decay across increasing distances between SNPs. Code borrowed from https://www.r-bloggers.com/2011/08/estimate-decay-of-linkage-disequilibrium-with-distance/
# Run the NLS model, following the formula for expected r^2
winnipeg_nls <- nls(formula = (rsq ~ ((10+C*Dist)/((2+C*Dist)*(11+C*Dist)))*(1+((3+C*Dist)*(12+12*C*Dist+(C*Dist)^2))/(60*(2+C*Dist)*(11+C*Dist)))), start=c(C = 0.1), control=nls.control(maxiter=100), data = pairwise_winnipeg)
# Estimate the population recombination parameter, which is C/distance
winnipeg_rho <- summary(winnipeg_nls)$parameters[1]
# Pull the distance values for each row
winnipeg_dist <- c(pairwise_winnipeg$Dist)
# Use the population recombination parameter rho with the distance values to estimate expected r^2 for each distance, controlling for sample size (60 here). Stick that info onto the main LD table for kansas
n_winnipeg <- 235
pairwise_winnipeg$fpoints<-((10+winnipeg_rho*winnipeg_dist)/((2+winnipeg_rho*winnipeg_dist)*(11+winnipeg_rho*winnipeg_dist)))*(1+((3+winnipeg_rho*winnipeg_dist)*(12+12*winnipeg_rho*winnipeg_dist+(winnipeg_rho*winnipeg_dist)^2))/(n_winnipeg*(2+winnipeg_rho*winnipeg_dist)*(11+winnipeg_rho*winnipeg_dist)))

# Take a look at the plot
ggplot(pairwise_winnipeg, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x="Distance (Base Pairs)", y=expression(Pairwise~italic(r)^{2}), title = "Lake Winnipeg") +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  xlim(0, 50000) +
  ylim(0, 0.5) +
  theme_bw()

# Plot pairwise LD with the model
winnipeg_pairwise_plot <- ggplot(pairwise_winnipeg, aes(x = Dist, y = rsq)) +
  labs(x=element_blank(), y=expression(Pairwise~italic(r)^{2}), title = "Lake Winnipeg") +
  geom_point(alpha = 0.4) +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  geom_hline(yintercept = 0.5) + 
  theme_bw()

# Check where LD decays below 0.5 (strong LD), where it decays to 0.2 (moderate LD), and where it decays to half its maximum value
# Half decay and strong LD terms end up not being informative, likely because of the low coverage data here. Moderate LD seems informative though
winnipeg_half_decay <- pairwise_winnipeg$Dist[which.min(abs(pairwise_winnipeg$fpoints - max(pairwise_winnipeg$fpoints/2)))]
winnipeg_strongLD <- min(pairwise_winnipeg$Dist[pairwise_winnipeg$fpoints <= 0.50])
winnipeg_moderateLD <- min(pairwise_winnipeg$Dist[pairwise_winnipeg$fpoints <= 0.20])

# Plotting script adapted from user rmf's response in the forum thread https://www.biostars.org/p/300381/
# Group data into 1 kb intervals
pairwise_winnipeg$distc <- cut(pairwise_winnipeg$Dist,breaks=seq(from=min(pairwise_winnipeg$Dist)-1,to=max(pairwise_winnipeg$Dist)+1,by=1000))

# Compute mean r2 within each 1 kb interval
pairwise_winnipeg_grouped <- pairwise_winnipeg %>% 
  group_by(distc) %>% 
  summarise(mean_20kb = mean(rsq),median_20kb = median(rsq))

# A helper step to get mid points of our distance intervals for plotting.
pairwise_winnipeg_grouped <- pairwise_winnipeg_grouped %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                      end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                      mid=start+((end-start)/2))

# Add in a column that just describes the population 
pairwise_winnipeg_grouped$pop <- rep("winnipeg", nrow(pairwise_winnipeg_grouped))

ggplot()+
  geom_point(data = pairwise_winnipeg_grouped, aes(x = start, y = mean_20kb), size = 0.4,colour="grey20") +
  geom_line(data = pairwise_winnipeg_grouped, aes(x = start, y = mean_20kb), size=0.3, alpha=0.5, colour="grey40") +
  labs(x="Distance (Kilobases)", y=expression(LD~(r^{2})))+
  theme_bw()
#ggsave("winnipeg_pairwiseLD.png", winnipeg_pairwise_plot, dpi = 300)

# Trying pairwise LD plots in Unassigned individuals
# Add in a column that just describes the population 
pairwise_unassigned <- read_tsv("unassigned_LDdecay.LD") %>% 
  dplyr::rename(rsq = `r^2`, chr = `#chr`) %>% 
  mutate(pop = rep("unassigned", nrow(.)))



# Using nonlinear least squares to plot LD decay across increasing distances between SNPs. Code borrowed from https://www.r-bloggers.com/2011/08/estimate-decay-of-linkage-disequilibrium-with-distance/
# Run the NLS model, following the formula for expected r^2
unassigned_nls <- nls(formula = (rsq ~ ((10+C*Dist)/((2+C*Dist)*(11+C*Dist)))*(1+((3+C*Dist)*(12+12*C*Dist+(C*Dist)^2))/(60*(2+C*Dist)*(11+C*Dist)))), start=c(C = 0.1), control=nls.control(maxiter=100), data = pairwise_unassigned)
# Estimate the population recombination parameter, which is C/distance
unassigned_rho <- summary(unassigned_nls)$parameters[1]
# Pull the distance values for each row
unassigned_dist <- c(pairwise_unassigned$Dist)
# Use the population recombination parameter rho with the distance values to estimate expected r^2 for each distance, controlling for sample size (60 here). Stick that info onto the main LD table for kansas
n_unassigned <- 29
pairwise_unassigned$fpoints<-((10+unassigned_rho*unassigned_dist)/((2+unassigned_rho*unassigned_dist)*(11+unassigned_rho*unassigned_dist)))*(1+((3+unassigned_rho*unassigned_dist)*(12+12*unassigned_rho*unassigned_dist+(unassigned_rho*unassigned_dist)^2))/(n_unassigned*(2+unassigned_rho*unassigned_dist)*(11+unassigned_rho*unassigned_dist)))

# Take a look at the plot
ggplot(pairwise_unassigned, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x="Distance (Base Pairs)", y=expression(Pairwise~italic(r)^{2}), title = "Unassigned Individuals") +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  xlim(0, 50000) +
  ylim(0, 0.5) +
  theme_bw()

# Plot pairwise LD between points along with the model.
unassigned_pairwise_plot <- ggplot(pairwise_unassigned, aes(x = Dist, y = rsq)) +
  geom_point(alpha = 0.4) +
  labs(x="Distance (Base Pairs)", y=expression(Pairwise~italic(r)^{2}), title = "Unassigned Individuals") +
  geom_line(aes(x = Dist, y = fpoints), colour = "red") +
  geom_hline(yintercept = 0.5) + 
  theme_bw()
# Check where LD decays below 0.5 (strong LD), where it decays to 0.2 (moderate LD), and where it decays to half its maximum value
# Half decay and strong LD terms end up not being informative, likely because of the low coverage data here. Moderate LD seems informative though
unassigned_half_decay <- pairwise_unassigned$Dist[which.min(abs(pairwise_unassigned$fpoints - max(pairwise_unassigned$fpoints/2)))]
unassigned_strongLD <- min(pairwise_unassigned$Dist[pairwise_unassigned$fpoints <= 0.50])
unassigned_moderateLD <- min(pairwise_unassigned$Dist[pairwise_unassigned$fpoints <= 0.20])

# Plotting script adapted from user rmf's response in the forum thread https://www.biostars.org/p/300381/
# Group data into 1 kb intervals
pairwise_unassigned$distc <- cut(pairwise_unassigned$Dist,breaks=seq(from=min(pairwise_unassigned$Dist)-1,to=max(pairwise_unassigned$Dist)+1,by=1000))

# Compute mean r2 within each 1 kb interval
pairwise_unassigned_grouped <- pairwise_unassigned %>% 
  group_by(distc) %>% 
  summarise(mean_20kb = mean(rsq),median_20kb = median(rsq))

# A helper step to get mid points of our distance intervals for plotting.
pairwise_unassigned_grouped <- pairwise_unassigned_grouped %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                                                      end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                                                      mid=start+((end-start)/2))

# Add in a column that just describes the population 
pairwise_unassigned_grouped$pop <- rep("unassigned", nrow(pairwise_unassigned_grouped))

ggplot()+
  geom_point(data = pairwise_unassigned_grouped, aes(x = start, y = mean_20kb), size = 0.4,colour="grey20") +
  geom_line(data = pairwise_unassigned_grouped, aes(x = start, y = mean_20kb), size=0.3, alpha=0.5, colour="grey40") +
  labs(x="Distance (Kilobases)", y=expression(LD~(r^{2})))+
  theme_bw()


#
# combined_pairwise_plot <- ggplot() +
#   geom_point(data = pairwise_winnipeg, aes(x = Dist, y = rsq), alpha = 0.1) +
#   geom_point(data = pairwise_manitoba, aes(x = Dist, y = rsq), alpha = 0.1, colour = "green") +
#   geom_point(data = pairwise_kansas, aes(x = Dist, y = rsq), alpha = 0.1, colour = "red") +
#     theme_bw()

combined_pairwise_plot <- winnipeg_pairwise_plot / manitoba_pairwise_plot | kansas_pairwise_plot / unassigned_pairwise_plot

ggsave("combined_pairwiseLD.pdf", combined_pairwise_plot, dpi = 600)


# Combine the dataframes of the three populations into one
combined_ld_grouped <- rbind(pairwise_winnipeg_grouped, pairwise_manitoba_grouped, pairwise_kansas_grouped, pairwise_unassigned_grouped)


combined_ld <- rbind(pairwise_winnipeg, pairwise_manitoba, pairwise_kansas, pairwise_unassigned)

# Plot the combined figure
combined_ld_plot <- ggplot(data = combined_ld_grouped, aes(x = start, y = mean_20kb, colour = pop)) +
  scale_colour_manual(values = c("#fdae61", "#d73027", "#4575b4", "black"),
                      breaks = c("kansas", "manitoba", "winnipeg", "unassigned"),
                      labels = c("Cedar Bluff Reservoir", "Lake Manitoba", "Lake Winnipeg", "Unassigned Individuals")) +
  scale_fill_manual(values = c("#fdae61", "#d73027", "#4575b4", "black"),
                      breaks = c("kansas", "manitoba", "winnipeg", "unassigned"),
                      labels = c("Cedar Bluff Reservoir", "Lake Manitoba", "Lake Winnipeg", "Unassigned Individuals")) +
  geom_point(alpha = 0.22) +
  labs(x="Distance (Kilobases)", y=expression(Mean~italic(r)^{2}~Kilobase^-1), colour = "Population", fill = "Population") +
  geom_smooth(se = TRUE, method = "glm", formula = (y ~ log(x)), linetype = 1, size = 0.8, alpha = 0.50, aes(fill = pop)) +
  scale_x_continuous(labels = c("0" = "0", "250000" = "250", "500000" = "500", "750000" = "750", "1000000" = "1000")) +
  theme_bw() +
  theme(text = element_text(size=15), legend.justification = c("right", "top"), legend.position = c(0.99, 0.99), legend.background = element_blank())
combined_ld_plot



#ggsave("combined_ld_plot.pdf", combined_ld_plot, dpi = 2000)  


# Plot curves for LD Decay with nonlinear regression models for each site
combined_ld_curve_plot <- ggplot(data = dplyr::filter(combined_ld, Dist <= 250000), aes(x = Dist, y = fpoints, colour = pop)) +
  scale_colour_manual(values = c("#fdae61", "#d73027", "#4575b4", "black"),
                      breaks = c("kansas", "manitoba", "winnipeg", "unassigned"),
                      labels = c("Cedar Bluff Reservoir", "Lake Manitoba", "Lake Winnipeg", "Unassigned Individuals")) +
  scale_fill_manual(values = c("#fdae61", "#d73027", "#4575b4", "black"),
                      breaks = c("kansas", "manitoba", "winnipeg", "unassigned"),
                      labels = c("Cedar Bluff Reservoir", "Lake Manitoba", "Lake Winnipeg", "Unassigned Individuals")) +
  #geom_point(alpha = 0.10) +
  geom_line() +
  labs(x="Distance between SNPs (Kilobases)", y=expression(Expected~Pairwise~italic(r)^{2}), colour = "Population", fill = "Population") +
  scale_x_continuous(labels = c("0" = "0", "50000" = "50", "100000" = "100", "150000" = "150", "200000" = "200", "250000" = "250")) +
  theme_bw() +
  theme(text = element_text(size=15), legend.justification = c("right", "top"), legend.position = c(0.99, 0.99), legend.background = element_blank())
combined_ld_curve_plot
ggsave("combined_ld_curve_plot.pdf", combined_ld_curve_plot, dpi = 2000)  



ggplot(data = combined_ld, aes(x = Dist, y = fpoints, colour = pop)) +
  scale_colour_manual(values = c("#fdae61", "#d73027", "#4575b4", "black"),
                      breaks = c("kansas", "manitoba", "winnipeg", "unassigned"),
                      labels = c("Cedar Bluff Reservoir", "Lake Manitoba", "Lake Winnipeg", "Unassigned Individuals")) +
  scale_fill_manual(values = c("#fdae61", "#d73027", "#4575b4", "black"),
                    breaks = c("kansas", "manitoba", "winnipeg", "unassigned"),
                    labels = c("Cedar Bluff Reservoir", "Lake Manitoba", "Lake Winnipeg", "Unassigned Individuals")) +
  #geom_point(alpha = 0.10) +
  geom_line() +
  labs(x="Distance between SNPs (Kilobases)", y=expression(Expected~Pairwise~italic(r)^{2}), colour = "Population", fill = "Population") +
  theme_bw() +
  theme(text = element_text(size=15), legend.justification = c("right", "top"), legend.position = c(0.99, 0.99), legend.background = element_blank())
######################################################################################
# Playing with looking at r^2 across each region of the genome
ggplot(data = pairwise_winnipeg, aes(x = Site2 - ((Site2 - Site1)/2), y = rsq), group = chr) +
  geom_point() +
  facet_grid(~chr) +
  theme_bw()
  
chr8_winnipeg <- ggplot(data = filter(pairwise_winnipeg, chr == "NC_041338.1"), aes(x = Site2 - ((Site2 - Site1)/2), y = rsq), group = chr) +
  geom_point() +
  ggtitle("Pairwise R^2 in Lake Winnipeg Walleye") +
  theme_bw()


chr8_unassigned <- ggplot(data = filter(pairwise_unassigned, chr == "NC_041338.1"), aes(x = Site2 - ((Site2 - Site1)/2), y = rsq), group = chr) +
  geom_point() +
  ggtitle("Pairwise R^2 in Unassigned Walleye", subtitle = "Admixed Between Lakes Manitoba & Winnipeg") +
  theme_bw()

chr8_manitoba <- ggplot(data = filter(pairwise_manitoba, chr == "NC_041338.1"), aes(x = Site2 - ((Site2 - Site1)/2), y = rsq), group = chr) +
  geom_point() +
  ggtitle("Pairwise R^2 in Lake Manitoba Walleye") +
  theme_bw()

chr8_kansas <- ggplot(data = filter(pairwise_kansas, chr == "NC_041338.1"), aes(x = Site2 - ((Site2 - Site1)/2), y = rsq), group = chr) +
  geom_point() +
  ggtitle("Pairwise R^2 in Cedar Bluff Reservoir Walleye") +
  theme_bw()

combined_chr8 <- (chr8_winnipeg | chr8_manitoba) / (chr8_unassigned | chr8_kansas)
ggsave(filename = "pairwise_rsq.pdf", plot = combined_chr8, dpi = 2000)
