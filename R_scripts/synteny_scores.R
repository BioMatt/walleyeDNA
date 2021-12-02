library(tidyverse)

# A script to read in weighted and relaxed synteny scores, plot histograms, and calculate confidence intervals

weighted_scores <- read_tsv("weighted_scores.txt", col_names = c("index", "score"))
relaxed_scores <- read_tsv("relaxed_scores.txt", col_names = c("index", "score"))

weighted_plot <- ggplot(data = weighted_scores, aes(score)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = "Weighted Synteny Scores", y = "Count", title = "Histogram of Weighted Synteny Scores") +
  theme_bw()

relaxed_plot <- ggplot(data = relaxed_scores, aes(score)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = "Relaxed Synteny Scores", y = "Count", title = "Histogram of Relaxed Synteny Scores") +
  theme_bw()

# Get quick confidence intervals
# Using bootstraps because the data are obviously not normal
# Notes for this process taken from http://www.stat.ucla.edu/~rgould/110as02/bsci
weighted_bstrap <- c()
for (i in 1:10000){
  b_sample <- sample(weighted_scores$score, 100, replace=T)
  b_estimate <-mean(b_sample)
  weighted_bstrap <- c(weighted_bstrap, b_estimate)
}
rm(b_sample, b_estimate, i)
#Mean
mean(weighted_scores$score)
#lower bound
quantile(weighted_bstrap,.05)
#upper bound
quantile(weighted_bstrap,.95)

# The same process for relaxed synteny scores
relaxed_bstrap <- c()
for (i in 1:10000){
  b_sample <- sample(relaxed_scores$score, 100, replace=T)
  b_estimate <-mean(b_sample)
  relaxed_bstrap <- c(relaxed_bstrap, b_estimate)
}
rm(b_sample, b_estimate, i)
#Mean
mean(relaxed_scores$score)
#lower bound
quantile(relaxed_bstrap,.05)
#upper bound
quantile(relaxed_bstrap,.95)

# Write the plots out
ggsave("relaxed_syntenyplot.pdf", relaxed_plot, dpi = 900)
ggsave("weighted_syntenyplot.pdf", weighted_plot, dpi = 900)
