################################################### -
## Title: Script for analyzing simulation results
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-01-10
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(ggplot2)
theme_set(ggpubr::theme_classic2(base_size=12))
library(here)

# load processed data
df_results_long <- readRDS(here("DataProcessed/all_2023-02-03.rds"))
df_failures <- readRDS(here("DataProcessed/failures_2023-02-03.rds"))
df_inference <- readRDS(here("DataProcessed/inference_2023-02-03.rds"))
df_performance <- readRDS(here("DataProcessed/performance_2023-02-03.rds"))
    
# plot bias as a function of true correlation
df_performance %>%
  filter(!(Method=="Max.conserv" & Bias < -1) & !(Method=="EM.alg" & M==0)) %>%
  filter(Distribution==1, Delta==0, N==10,
         Method %in% c("Rho.hat", "EM.alg", "Bayes.Jeffreys", "Pearson",
                       "Freq.20th.quantile", "Max.conserv"),
         Prop.matched %in% c(0, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  ggplot(aes(col=Method, x=Rho, y=Bias)) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_performance$Rho)) +
  labs(title="Bias of estimators by number of matched samples",
       subtitle="Varied from 0 to 10 matched samples (n=10)",
       x="True correlation",
       y="Bias (Estimated minus true value)",
       caption="Results averaged over 1,000 datasets at each point.")

# plot MSE as a function of true correlation
df_performance %>%
  filter(!(Method=="Max.conserv" & Mean.sqd.error > 1) & !(Method=="EM.alg" & M==0)) %>%
  filter(Distribution==1, Delta==0, N==10,
         Method %in% c("Rho.hat", "EM.alg", "Bayes.Jeffreys", "Pearson",
                       "Freq.20th.quantile", "Max.conserv"),
         Prop.matched %in% c(0, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  ggplot(aes(col=Method, x=Rho, y=Mean.sqd.error)) +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_performance$Rho)) +
  labs(title="MSE of estimators by number of matched samples",
       subtitle="Varied from 0 to 10 matched samples (n=10)",
       x="True correlation",
       caption="Results averaged over 1,000 datasets at each point.")

# plot association with between any two estimators
df_results_long %>%
  filter(!(Method=="EM.alg" & M==0)) %>%
  filter(Distribution==1, Delta==0, N==10,
         Method %in% c("Rho.hat", "EM.alg"),
         Prop.matched %in% c(0, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  tidyr::pivot_wider(names_from=Method, values_from=Estimate) %>%
  sample_frac(0.1) %>%
  ggplot(aes(x=Rho.hat, y=EM.alg)) +
  facet_wrap(~ M) +
  geom_point(alpha=0.25) +
  geom_smooth(method="lm", formula="y ~ x", se=F) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title="Scatterplot of association between estimators",
       subtitle="Varied from 0 to 10 matched samples (n=10)",
       caption="Results sampled from 1,000 datasets for efficiency.")

# plot Type-I error rate as function of correlation
df_inference %>%
  filter(!(Method=="EM.alg" & M==0)) %>%
  filter(Distribution==1, Delta==0, N==10,
         Method %in% c("Oracle", "Independent", "Pearson",
                       "EM.alg", "Bayes.Jeffreys"),
         Prop.matched %in% c(0, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  ggplot(aes(col=Method, x=Rho, y=Rejection.rate)) +
  geom_hline(yintercept=0.05, linetype="dashed") +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_inference$Rho)) +
  labs(title="Type-I error rates by number of matched samples",
       subtitle="Varied from 0 to 10 matched samples (n=10)",
       x="True correlation",
       y="Type-I error (rejected null hypotheses)",
       caption="Results averaged over 1,000 datasets at each point.")

