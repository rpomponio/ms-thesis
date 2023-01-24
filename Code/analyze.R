################################################### -
## Title: Script for analyzing simulation results
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-10
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(tidyr)
library(ggplot2)
theme_set(ggpubr::theme_classic2(base_size=12))
library(here)

# load previously-saved results
results <- readRDS(here("DataRaw/simulation_results_2023-01-21.rds"))

# create "error" matrix
errors <- results[, 15:26] - results[, "Rho"]

# calculate "failure" rate, by frequency of NAs
df_failures <- data.frame(results) %>%
  filter(Delta==0) %>%
  pivot_longer(cols=Max.conserv:Bayes.arcsine, names_to="Method") %>%
  group_by(Method, Distribution, Rho, Delta, N, M) %>%
  summarise(Failures=sum(is.na(value)), Failure.rate=mean(is.na(value))) %>%
  filter(Failures > 0)

# calculate average error, or "bias", and mean squared error
df_performance <- data.frame(errors) %>%
  bind_cols(select(data.frame(results), Distribution:M)) %>%
  pivot_longer(cols=Max.conserv:Bayes.arcsine, names_to="Method") %>%
  group_by(Method, Distribution, Rho, Delta, N, M) %>%
  summarise(Bias=mean(value, na.rm=T), Mean.sqd.error=mean(value^2, na.rm=T))

saveRDS(df_performance, here("DataProcessed/performance_2023-01-21.rds"))
    
# plot bias as a function of true correlation
df_performance %>%
  filter(Distribution==1, Delta==0, N==50,
         Method %in% c("EM.alg", "Bayes.unif", "Bayes.Jeffreys",
                       "Bayes.arcsine", "Pearson"),
          M %in% 0:15) %>%
  ggplot(aes(col=Method, x=Rho, y=Bias)) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_performance$Rho)) +
  labs(title="Bias of estimators by number of matched samples (0 to 15)",
       subtitle="Measured by mean difference (in all cases total sample size was 50)",
       x="True correlation",
       y="Bias (Estimated minus true value)",
       caption="Results averaged over 1,000 datasets at each point.")

ggsave(filename="~/Downloads/Sim_Results_Bias.png", width=10.5, height=7.5)

# plot MSE as a function of true correlation
df_performance %>%
  filter(Distribution==1, Delta==0, N==50,
         Method %in% c("EM.alg", "Bayes.unif", "Bayes.Jeffreys",
                       "Bayes.arcsine", "Pearson"),
         M %in% 0:15) %>%
  ggplot(aes(col=Method, x=Rho, y=Mean.sqd.error)) +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_performance$Rho)) +
  labs(title="Variance of estimators by number of matched samples (0 to 15)",
       subtitle="Measured by Mean Sqd. Error (in all cases total sample size was 50)",
       x="True correlation",
       caption="Results averaged over 1,000 datasets at each point.")

ggsave(filename="~/Downloads/Sim_Results_Variance.png", width=10.5, height=7.5)








