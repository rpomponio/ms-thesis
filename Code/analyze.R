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
library(kableExtra, warn.conflicts=FALSE)
library(here)

# load processed data
df_failures <- readRDS(here("DataProcessed/failures_2023-02-20.rds"))
df_inference <- readRDS(here("DataProcessed/inference_2023-02-20.rds"))
df_performance <- readRDS(here("DataProcessed/performance_2023-02-20.rds"))

##### Bias plots #####

# plot bias versus true correlation
df_performance %>%
  filter(Distribution==1,
         Delta==0,
         N==20,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile", "Pearson"),
         Prop.matched==0.2) %>%
  mutate(lower=Bias-Bias.mce, upper=Bias+Bias.mce) %>%
  ggplot(aes(col=Method, x=Rho, y=Bias)) +
  geom_hline(yintercept=0, linetype="dashed") +
  # facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05, alpha=0.25) +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Bias of estimators by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Bias (Estimated minus true value)",
       caption="Results averaged over 10,000 datasets at each point.")

# plot bias versus number of matched samples
df_performance %>%
  filter(Distribution==1,
         Delta==0,
         N==20,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile", "Pearson")) %>%
  mutate(Rho = round(Rho, 2)) %>%
  ggplot(aes(col=Method, x=M, y=Bias)) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~ Rho) +
  geom_line() +
  geom_point(size=2) +
  scale_x_continuous(breaks=seq(1, 19, 2)) +
  labs(title="Bias of estimators by number of matched samples",
       subtitle="Simulated with n=20, m=4",
       x="Number of matched samples",
       y="Bias (Estimated minus true value)",
       caption="Results averaged over 10,000 datasets at each point.")

##### MSE plots #####

df_performance %>%
  filter(Distribution==1,
         Delta==0,
         N==20,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile", "Pearson"),
         Prop.matched==0.2) %>%
  mutate(
    lower=Mean.sqd.error-Mean.sqd.error.mce,
    upper=Mean.sqd.error+Mean.sqd.error.mce) %>%
  ggplot(aes(col=Method, x=Rho, y=Mean.sqd.error)) +
  geom_hline(yintercept=0, linetype="dashed") +
  # facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05, alpha=0.25) +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="MSE of estimators by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Mean squared error",
       caption="Results averaged over 10,000 datasets at each point.")

##### SE Ratio plots #####

df_oracle <- df_inference %>%
  filter(Distribution==1,
         Delta==0,
         N==20,
         Method=="Oracle",
         Prop.matched==0.2) %>%
  rename(SE.mean.oracle=SE.mean, Rejection.rate.oracle=Rejection.rate) %>%
  select(Distribution:M, SE.mean.oracle, Rejection.rate.oracle)

# plot the oracle SE by true correlation
df_oracle %>%
  ggplot(aes(x=Rho, y=SE.mean.oracle)) +
  geom_hline(yintercept=0, linetype="dashed") +
  # facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Oracle SE by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Standard Error",
       caption="Results averaged over 10,000 datasets at each point.")

# plot the ratio of SE:Oracle SE by true correlation
df_inference %>%
  filter(Distribution==1,
         Delta==0,
         N==20,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile", "Pearson"),
         Prop.matched==0.2) %>%
  left_join(df_oracle) %>%
  mutate(SE.mean.ratio=SE.mean / SE.mean.oracle) %>%
  ggplot(aes(col=Method, x=Rho, y=SE.mean.ratio)) +
  geom_hline(yintercept=1, linetype="dashed") +
  # facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Ratio of SE:Oracle SE by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Standard Error Ratio",
       caption="Results averaged over 10,000 datasets at each point.")

##### Rejection rate plots #####

# plot the oracle Type I error rate by true correlation
df_oracle %>%
  ggplot(aes(x=Rho, y=Rejection.rate.oracle)) +
  geom_hline(yintercept=0.05, linetype="dashed") +
  # facet_wrap(~ Delta, scales="free_y") +
  ylim(c(0, 0.2)) +
  geom_line() +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Rejection rate of oracle test by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Prop. rejected null hypotheses",
       caption="Results averaged over 10,000 datasets at each point.")

# plot the Type I error rate by true correlation
df_inference %>%
  filter(Distribution==1,
         Delta %in% c(0, 0.25),
         N==20,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile", "Pearson"),
         Prop.matched==0.2) %>%
  ggplot(aes(col=Method, x=Rho, y=Rejection.rate)) +
  geom_hline(yintercept=0.05, linetype="dashed") +
  facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Rejection rate of hypoth. tests by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Prop. rejected null hypotheses",
       caption="Results averaged over 10,000 datasets at each point.")

# compare rejection rate to oracle rate
df_inference %>%
  filter(Distribution==1,
         Delta==0,
         N==20,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile", "Pearson"),
         Prop.matched==0.2) %>%
  left_join(df_oracle) %>%
  mutate(Rejection.rate.ratio=Rejection.rate / Rejection.rate.oracle) %>%
  ggplot(aes(col=Method, x=Rho, y=Rejection.rate.ratio)) +
  geom_hline(yintercept=1, linetype="dashed") +
  # facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Ratio of Rej. rate:Oracle rej. rate by true correlation",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Ratio of rejected null hypotheses",
       caption="Results averaged over 10,000 datasets at each point.")


##### Summarize failure rates by simulation setting #####

tab_failures <- df_failures %>%
  filter(Distribution==1,
         Delta==0,
         N <= 50,
         Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile",
                       "Pearson")) %>%
  select(Method, Rho, N, Prop.matched, M, Failure.rate) %>%
  tidyr::pivot_wider(id_cols=c(Method, N, Prop.matched, M),
                     names_from=Rho,
                     values_from=Failure.rate)

tab_failures %>%
  kable(digits=3) %>%
  kable_styling(fixed_thead=T, full_width=F, font_size=16) %>%
  column_spec(5, background=spec_color(tab_failures$`-0.9`, alpha=0.5)) %>%
  column_spec(6, background=spec_color(tab_failures$`-0.5`, alpha=0.5)) %>%
  column_spec(7, background=spec_color(tab_failures$`-0.25`, alpha=0.5)) %>%
  column_spec(8, background=spec_color(tab_failures$`0`, alpha=0.5)) %>%
  column_spec(9, background=spec_color(tab_failures$`0.25`, alpha=0.5)) %>%
  column_spec(10, background=spec_color(tab_failures$`0.5`, alpha=0.5)) %>%
  column_spec(11, background=spec_color(tab_failures$`0.9`, alpha=0.5))

##### Plot bias as function of matched samples #####

# requires loading full results set (may exceed memory capacity)
df_results_long <- readRDS(here("DataProcessed/all_2023-02-20.rds"))

# Calculate simulated P(r > rho) for various values of m for each method
# If r > rho, we'd expect SEs to be too low 
df_results_long %>%
  filter(Method %in% c("Bayes.arcsine", "EM.alg", "Freq.20th.quantile",
                       "Pearson")) %>%
  filter(Distribution==1, Delta==0, N==20,
         Prop.matched %in% c(0, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  group_by(Method, Rho, Prop.matched, M) %>%
  summarise(Prop.overestimate=mean(Estimate > Rho, na.rm=T)) %>%
  filter(!is.na(Prop.overestimate)) %>%
  ggplot(aes(col=Method, x=M, y=Prop.overestimate)) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  facet_wrap(~ Rho) +
  geom_line() +
  geom_point(size=2) +
  ylim(c(0, 1)) +
  scale_x_continuous(breaks=seq(1, 20, 2)) +
  labs(subtitle="Simulated with n=20, number of matched samples varying.",
       x="Number of matched samples",
       y="mBias",
       caption="mBias is the probability P(r > Rho)")


##### Adjusted (re-scaled) T statistic #####

# compare rates of un-scaled T versus scaled T (using @Ryan's correction)
df_inference %>%
  filter(Distribution==1,
         Delta%in%c(0, 0.25),
         N==20,
         Method %in% c("Pearson", "Pearson.scaled"),
         Prop.matched==0.2) %>%
  ggplot(aes(shape=Method, x=Rho, y=Rejection.rate)) +
  geom_hline(yintercept=0.05, linetype="dashed") +
  facet_wrap(~ Delta, scales="free_y") +
  geom_line() +
  geom_point(size=5) +
  scale_x_continuous(breaks=unique(df_performance$Rho[df_performance$Distribution==1])) +
  labs(title="Comparison of Type I error rate before/after scaling T statistic",
       subtitle="Simulated with n=20, m=4",
       x="True correlation",
       y="Prop. rejected null hypotheses",
       caption="Results averaged over 10,000 datasets at each point.")

# requires loading full results set (may exceed memory capacity)
df_stats <- readRDS(here("DataProcessed/stats_2023-02-20.rds"))



