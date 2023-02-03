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
results <- readRDS(here("DataRaw/simulation_results_2023-02-03.rds"))

# temporary fix for NaNs
results[is.nan(results)] <- NA

# pivot results to longer, flag all "Failures"
df_results_long <- data.frame(results) %>%
  pivot_longer(cols=Max.conserv:Bayes.arcsine, names_to="Method") %>%
  rename(Estimate=value) %>%
  mutate(Failure=ifelse(is.na(Estimate) | Estimate >= 1 | Estimate <= -1, 1, 0))

# calculate "failure" rate, by frequency of NAs
df_failures <- df_results_long %>%
  group_by(Method, Distribution, M) %>%
  summarise(
    Num.iterations = n(),
    Num.failures=sum(Failure),
    Failure.rate=mean(Failure))

df_failures %>%
  filter(Failure.rate > 0.01, M > 3) %>%
  View()

# oracle method knows exact correlation
df_oracle <- df_results_long %>%
  filter(Method=="Max.conserv") %>%
  mutate(Estimate=Rho, Method="Oracle", Falure=0)

# naive method assumes zero correlation
df_naive <- df_results_long %>%
  filter(Method=="Max.conserv") %>%
  mutate(Estimate=0, Method="Independence", Failure=0)

# compute (modified) t-statistics and corresponding p-values
df_inference <- df_results_long %>%
  bind_rows(df_oracle, df_naive) %>%
  mutate(
    d.bar = Mu.hat.X - Mu.hat.Y,
    sum.sqs = (Sigma.sqd.hat.X * N) + (Sigma.sqd.hat.Y * N),
    t.ind = d.bar / sqrt( sum.sqs / (N * (N - 1)) ),
    t.mod = ifelse( Failure==1, NA, t.ind / sqrt(1 - Estimate) ),
    p.value = 2 * pt(abs(t.mod), df=(2 * N - 2), lower.tail=F),
    d.lwr = 0, ###how to compute Conf Int of mean?
    d.upr = 0) %>%
  group_by(Method, Distribution, Rho, Delta, N, M) %>%
  summarise(Rejection.rate=mean(p.value < 0.05, na.rm=T))

# create "error" matrix
errors <- results[, 15:26] - results[, "Rho"]

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

# plot correlation with pseudo-"oracle" estimator
df_results_long %>%
  filter(Distribution==1, Delta==0, N==20,
         Method %in% c("EM.alg", "Bayes.unif", "Bayes.Jeffreys",
                       "Bayes.arcsine", "Pearson", "Max.conserv"),
         M == 10) %>%
  sample_frac(0.1) %>%
  ggplot(aes(x=Rho.hat, y=Estimate)) +
  facet_wrap(~ Method) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype="dashed", col="blue") +
  labs(title="Correlation between estimators and 'oracle' estimator",
       x="Pearson correlation of matched samples",
       caption="Results sampled from 1,000 datasets for efficiency.")

# plot Type-I error rate as function of correlation
df_inference %>%
  filter(Distribution==1, Delta==0, N==50,
         Method %in% c("EM.alg", "Max.conserv", "Freq.20th.quantile",
                       "Bayes.arcsine", "Pearson", "Oracle"),
         M %in% 0:15) %>%
  ggplot(aes(col=Method, x=Rho, y=Rejection.rate)) +
  geom_hline(yintercept=0.05, linetype="dashed") +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_inference$Rho)) +
  labs(title="Type-I error of methods by number of matched samples (0 to 15)",
       subtitle="In all cases total sample size was 50",
       x="True correlation",
       y="Type-I error (rejected null hypotheses)",
       caption="Results averaged over 1,000 datasets at each point.")

