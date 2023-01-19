################################################### -
## Title: Script for analyzing simulation results
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-10
################################################### -

library(tidyverse)
library(here)

# load previously-saved results
results <- readRDS(here("DataRaw/simulation_results_2023-01-18.rds"))

# create "error" matrix
errors <- results[, 10:16] - results[, "Rho"]

# calculate average error, or "bias", and mean squared error
df_perf_by_method <- data.frame(errors) %>%
  bind_cols(select(data.frame(results), Rho, M)) %>%
  pivot_longer(cols=Max.conserv:Bayes.unif, names_to="Method") %>%
  group_by(M, Rho, Method) %>%
  summarise(Bias=mean(value), Mean.sqd.error=mean(value^2))
    
# plot bias as a function of true correlation
df_perf_by_method %>%
  filter(! Method %in% c("Max.conserv", "Bayes.unif")) %>%
  ggplot(aes(col=Method, x=Rho, y=Bias)) +
  facet_wrap(~ M) +
  geom_line()

# plot MSE as a function of true correlation
df_perf_by_method %>%
  filter(! Method %in% c("Max.conserv", "Bayes.unif")) %>%
  ggplot(aes(col=Method, x=Rho, y=Mean.sqd.error)) +
  facet_wrap(~ M) +
  geom_line()


