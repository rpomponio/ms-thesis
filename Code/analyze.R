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
# df_results_long <- readRDS(here("DataProcessed/all_2023-02-20.rds"))
# df_failures <- readRDS(here("DataProcessed/failures_2023-02-20.rds"))
# df_inference <- readRDS(here("DataProcessed/inference_2023-02-20.rds"))
# df_performance <- readRDS(here("DataProcessed/performance_2023-02-20.rds"))

##### Summarize failure rates by simulation setting #####

tab_failures <- df_failures %>%
  filter(Distribution==2,
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

# ggsave("~/Downloads/mBias.png", width=10, height=8)





##### Reproduce Shiny App plots #####

# plot bias as a function of true correlation
df_performance %>%
  filter(!(Method=="Max.conserv" & Bias < -1)) %>%
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
       caption="Results averaged over 10,000 datasets at each point.")

# plot MSE as a function of true correlation
df_performance %>%
  filter(!(Method=="Max.conserv" & Mean.sqd.error > 1)) %>%
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
       caption="Results averaged over 10,000 datasets at each point.")

# plot association with between any two estimators
df_results_long %>%
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
       caption="Results sampled from 10,000 datasets for efficiency.")

# plot Type-I error rate as function of correlation
df_inference %>%
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
       caption="Results averaged over 10,000 datasets at each point.")


# plot SE average as function of correlation
df_inference %>%
  filter(Distribution==1, Delta==0, N==10,
         Method %in% c("Oracle", "Independent", "Pearson",
                       "EM.alg", "Bayes.Jeffreys"),
         Prop.matched %in% c(0, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  ggplot(aes(col=Method, x=Rho, y=SE.mean)) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~ M, scales="free_y") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks=unique(df_inference$Rho)) +
  labs(title="Standard error averages by number of matched samples",
       subtitle="Varied from 0 to 10 matched samples (n=10)",
       x="True correlation",
       y="Standard Error (mean)",
       caption="Results averaged over 10,000 datasets at each point.")

# how to present this data in a table
df_inference %>%
  filter(Distribution==1, Delta==0, N==10,
         Method %in% c("Oracle", "Independent", "Pearson",
                       "EM.alg", "Bayes.Jeffreys"),
         Prop.matched == 0.5) %>%
  select(Delta, Method, Rho, N, M, SE.mean) %>%
  tidyr::pivot_wider(names_prefix="Rho=", names_from=Rho, values_from=SE.mean) %>%
  DT::datatable(options=list(pageLength=15, dom='tip'), caption="SE") %>%
  DT::formatRound(columns=5:11)
