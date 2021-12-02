# This script is for plotting CPM for PDHX, EHF, and LRRC4C for walleye sampled in 2017 and 2018

library(tidyverse)
library(patchwork)

# Read in count data, then metadata from Thorstensen et al. 2020
cpm <- read_tsv("inversion_genes.txt")

metadata <- read_tsv("wall_samples.txt")

# Combine the datasets. Make Year collected into a factor instead of a numeric variable. Re-arrange site collected so plots go from south to north, left to right.
combined_data <- left_join(cpm, metadata, by = c("walleye_ID" = "Sample_ID")) %>% 
  mutate(Year = as.factor (Year), Location = fct_reorder(Location, desc(Location)))


# Plot all three genes
pdhx <- ggplot(combined_data, aes(x = Location, y = pdhx, colour = Year)) +
  geom_jitter(width = 0.15, size = 1.5) +
  scale_colour_manual(values = c("#d73027", "#4575b4")) +
  ggtitle(bquote(italic(PDHX)~counts~per~million)) +
  ylab("CPM") +
  xlab(element_blank()) +
  scale_x_discrete(labels = element_blank()) +
  theme_bw() +
  theme(
    legend.position = c(.98, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_blank(),
    legend.key = element_blank(),
    text = element_text(size = 16),
    axis.ticks.x = element_blank()
  )


ehf <- ggplot(combined_data, aes(x = Location, y = ehf, colour = Year)) +
  geom_jitter(width = 0.15, size = 1.5) +
  scale_colour_manual(values = c("#d73027", "#4575b4")) +
  ggtitle(bquote(italic(EHF)~counts~per~million)) +
  ylab("CPM") +
  xlab(element_blank()) +
  scale_x_discrete(labels = element_blank()) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    axis.ticks.x = element_blank()
  )

lrrc4c <- ggplot(combined_data, aes(x = Location, y = lrrc4c, colour = Year)) +
  geom_jitter(width = 0.15, size = 1.5) +
  scale_colour_manual(values = c("#d73027", "#4575b4")) +
  ggtitle(bquote(italic(LRRC4C)~counts~per~million)) +
  ylab("CPM") +
  scale_x_discrete(labels = c("Red River", "Matheson Island", "Dauphin River")) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 16)
  )

combined_plot <- pdhx / ehf/ lrrc4c
ggsave(filename = "inversion_genes.pdf", plot = combined_plot, dpi = 2000)
