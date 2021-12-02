library(tidyverse)

admixture <- read_tsv("k3_ancestry_metadata.txt")

admixture %>% group_by(V2, ancestry) %>% tally()

inversion <- read_delim("chr8_inv_pca_cluster.txt", delim = " ")

combined_data <- left_join(admixture, inversion, by = c("V1" = "name")) %>% 
  mutate(cluster = as.factor(cluster))

#View(filter(combined_data, V2 == "Dauphin_River"))

dauphin_river <- filter(combined_data, V2 == "Dauphin_River")
dauphin_tally <- dauphin_river %>% group_by(ancestry, cluster) %>%  tally()



ggplot(dauphin_tally, aes(x=ancestry, y=n, fill=cluster, group = ancestry)) +
  geom_bar(stat="identity", width=0.5) +
  ggtitle("Ancestry and Inversion Cluster for Dauphin River Walleye") +
  theme_bw()


combined_tally <- combined_data %>% group_by(V2, ancestry, cluster) %>%  tally()

freq_plot <- ggplot(combined_tally, aes(x=ancestry, y=n, fill=cluster, group = V2)) +
  geom_bar(stat="identity", width=0.5) +
  facet_wrap(facets = vars(V2)) +
  xlab("Assignment") +
  theme_bw()

ggsave(filename = "inversion_assignments.pdf", plot = freq_plot, dpi = 2000, height = 10, width = 10)


ggplot(combined_tally, aes(x=V2, y=n, fill=cluster, group = V2)) +
  geom_bar(stat="identity", width=0.5) +
  xlab("Assignment") +
  theme_bw()

tally_inversion <- combined_data %>% group_by(ancestry, pop, cluster) %>%  tally()
ancestry_inversion <- combined_data %>% group_by(ancestry, cluster) %>%  tally()
