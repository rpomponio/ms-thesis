################################################### -
## Title: Process simulation results for down-stream analysis
## Author: Ray Pomponio
## Email: raymond.pomponio@cuanschutz.edu
## Date Created: 2023-02-10
################################################### -

options(stringsAsFactors=FALSE)
library(dplyr, warn.conflicts=FALSE)
options(dplyr.summarise.inform=FALSE)
library(tidyr)
library(here)

# load previously-saved simulation results
results <- readRDS(here("DataRaw/simulation_results_2023-02-15.rds"))

# pivot to long data frame
df_results_long <- data.frame(results) %>%
  pivot_longer(cols=Rho.hat:Bayes.arcsine, names_to="Method") %>%
  rename(Estimate=value)

# warn of invalid estimates
n_ests_invalid <- sum(
  df_results_long$Estimate > 1 | df_results_long$Estimate < -1, na.rm=T)
warning(paste0("Number of invalid correlation estimates: ", n_ests_invalid))

# flag all "Failures" as cases where valid correlation was not estimated
df_results_long <- df_results_long %>%
  mutate(Failure.type=case_when(
    !is.finite(Estimate) ~ "Undefined",
    Estimate < -1 | Estimate > 1 ~ "Invalid",
    abs(Estimate - 1) < 0.000001 ~ "Boundary",
    abs(Estimate + 1) < 0.000001 ~ "Boundary",
    TRUE ~ "Nonfailure")) %>%
  mutate(Failure=ifelse(Failure.type=="Nonfailure", 0, 1))

# calculate failure rate per method
df_failures <- df_results_long %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  summarise(
    Num.iterations=n(),
    Num.failures=sum(Failure),
    Failure.rate=mean(Failure)) %>%
  ungroup()

# for 'oracle' inference method: use true correlation (or effective correlation)
df_oracle <- df_results_long %>%
  filter(Method=="Rho.hat") %>%
  mutate(Estimate=case_when(
    Distribution==1 ~ Rho,
    Rho==-0.9  ~ -0.579453,
    Rho==-0.5  ~ -0.3591659,
    Rho==-0.25 ~ -0.1890325,
    Rho==0     ~ 0.0009432011,
    Rho==0.25  ~ 0.203865,
    Rho==0.5   ~ 0.425501,
    Rho==0.9   ~ 0.8407945)) %>%
  mutate(Method="Oracle")

# for naive inference method: assume zero correlation
df_naive <- df_results_long %>%
  filter(Method=="Rho.hat") %>%
  mutate(Estimate=0) %>%
  mutate(Method="Independent")

# compute (modified) t-statistics and corresponding p-values
df_inference <- df_results_long %>%
  bind_rows(df_oracle, df_naive) %>%
  mutate(Estimate=ifelse(Failure.type=="Boundary", Estimate * 0.9999, Estimate)) %>%
  mutate(Failure=ifelse(Failure.type=="Boundary", 0, Failure)) %>%
  mutate(
    d.bar = Mu.hat.X - Mu.hat.Y,
    sum.sqs = (Sigma.sqd.hat.X * N) + (Sigma.sqd.hat.Y * N),
    t.ind = d.bar / sqrt( sum.sqs / (N * (N - 1)) ),
    t.mod = suppressWarnings(ifelse(Failure==1, NA, t.ind / sqrt(1 - Estimate))),
    p.value = 2 * pt(abs(t.mod), df=(2 * N - 2), lower.tail=F),
    # To-do: calculate SEs using denominator of modified t stat
    d.lwr = NA, # @Ryan: how to compute Conf Int of mean for modified test?
    d.upr = NA) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  summarise(Rejection.rate=mean(p.value < 0.05, na.rm=T)) %>%
  ungroup()

# calculate bias and MSE (accounting for effective correlation in ordinal data)
df_performance <- df_results_long %>%
  mutate(Error=case_when(
    Distribution==1 ~ Estimate - Rho,
    Rho==-0.9  ~ Estimate - (-0.579453),
    Rho==-0.5  ~ Estimate - (-0.3591659),
    Rho==-0.25 ~ Estimate - (-0.1890325),
    Rho==0     ~ Estimate - (0.0009432011),
    Rho==0.25  ~ Estimate - (0.203865),
    Rho==0.5   ~ Estimate - (0.425501),
    Rho==0.9   ~ Estimate - (0.8407945))) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  summarise(Bias=mean(Error, na.rm=T), Mean.sqd.error=mean(Error^2, na.rm=T)) %>%
  ungroup()

# save processed data
saveRDS(df_results_long, here("DataProcessed/all_2023-02-15.rds"))
saveRDS(df_failures, here("DataProcessed/failures_2023-02-15.rds"))
saveRDS(df_inference, here("DataProcessed/inference_2023-02-15.rds"))
saveRDS(df_performance, here("DataProcessed/performance_2023-02-15.rds"))

# save 'lite' version of long data for improved speed
set.seed(2023)
df_results_lite <- df_results_long %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  slice_sample(prop=0.05) %>%
  ungroup()
saveRDS(df_results_lite, here("DataProcessed/lite_2023-02-15.rds"))
