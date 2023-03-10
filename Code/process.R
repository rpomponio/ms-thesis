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
results <- rbind(
  readRDS(here("DataRaw/simulation_results_2023-02-19_seed2020.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-19_seed2021.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-19_seed2022.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-15_seed2023.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-15_seed2024.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-16_seed2025.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-17_seed2026.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-17_seed2027.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-18_seed2028.rds")),
  readRDS(here("DataRaw/simulation_results_2023-02-18_seed2029.rds"))
  )

# temp fix : assign NA wherever EM algorithm had zero matched samples
results[(results[, "M"] == 0), "EM.alg"] <- NA

# temp fix : correct overestimates in EM algorithm (ordinal datasets)
results[(results[, "EM.alg"] > 1), "EM.alg"] <- 1

# pivot to long data frame
df_results_long <- data.frame(results) %>%
  pivot_longer(cols=Rho.hat:Bayes.arcsine, names_to="Method") %>%
  rename(Estimate=value)

# for ordinal datasets: change true correlation to "effective" correlation
df_results_long <- df_results_long %>%
  mutate(Rho=case_when(
    Distribution==1 ~ Rho,
    Rho==-0.9  ~ -0.5797045,
    Rho==-0.5  ~ -0.3609232,
    Rho==-0.25 ~ -0.1890938,
    Rho==0     ~ -8.124838e-05,
    Rho==0.25  ~ 0.2047809,
    Rho==0.5   ~ 0.4266981,
    Rho==0.9   ~ 0.840682)) %>%
  mutate(Sigma.X=case_when(
    Distribution==1 ~ Sigma.X,
    Distribution==2 ~ 2.051202)) %>%
  mutate(Sigma.Y=case_when(
    Distribution==1 ~ Sigma.Y,
    Distribution==2 ~ 1.955840))

# warn of invalid estimates, arising from Bayesian estimators with ordinal data
n_ests_invalid <- sum(
  df_results_long$Estimate > 1 | df_results_long$Estimate < -1, na.rm=T)
warning(paste0("Number of invalid correlation estimates: ", n_ests_invalid))

# flag all estimation failures (cases where valid correlation was not estimated)
df_results_long <- df_results_long %>%
  mutate(Failure.type=case_when(
    !is.finite(Estimate) ~ "Undefined",
    Estimate < -1 | Estimate > 1 ~ "Invalid",
    abs(Estimate - 1) < 0.000001 ~ "Boundary",
    abs(Estimate + 1) < 0.000001 ~ "Boundary",
    TRUE ~ "Nonfailure"))

# in ordinal, correct undefined failures due to zero variance in matched samples
df_pearson_mod <- df_results_long %>%
  filter(Method %in% c("Pearson", "Max.conserv")) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  mutate(Index=1:n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(Index, Repetition:Sigma.sqd.hat.Y), names_from=Method,
              values_from=c(Estimate, Failure.type)) %>%
  mutate(Estimate_Pearson=case_when(
    Distribution==2 & M>=2 & Failure.type_Pearson=="Undefined" ~ Estimate_Max.conserv,
    TRUE ~ Estimate_Pearson)) %>%
  mutate(Failure.type_Pearson=case_when(
    Distribution==2 & M>=2 & Failure.type_Pearson=="Undefined" ~ "Nonfailure",
    TRUE ~ Failure.type_Pearson)) %>%
  select(-contains("Max.conserv"), -Index) %>%
  rename(Estimate=Estimate_Pearson, Failure.type=Failure.type_Pearson) %>%
  mutate(Method="Pearson.mod")

df_quantile_mod <- df_results_long %>%
  filter(Method %in% c("Freq.20th.quantile", "Max.conserv")) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  mutate(Index=1:n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(Index, Repetition:Sigma.sqd.hat.Y), names_from=Method,
              values_from=c(Estimate, Failure.type)) %>%
  mutate(Estimate_Freq.20th.quantile=case_when(
    Distribution==2 & M>=4 & Failure.type_Freq.20th.quantile=="Undefined" ~
      Estimate_Max.conserv,
    TRUE ~ Estimate_Freq.20th.quantile)) %>%
  mutate(Failure.type_Freq.20th.quantile=case_when(
    Distribution==2 & M>=4 & Failure.type_Freq.20th.quantile=="Undefined" ~
      "Nonfailure",
    TRUE ~ Failure.type_Freq.20th.quantile)) %>%
  select(-contains("Max.conserv"), -Index) %>%
  rename(Estimate=Estimate_Freq.20th.quantile,
         Failure.type=Failure.type_Freq.20th.quantile) %>%
  mutate(Method="Freq.20th.quantile.mod")

# add "modified" estimators
df_results_long <- df_results_long %>%
  bind_rows(df_pearson_mod, df_quantile_mod)

# calculate failure rate per method and scenario
df_failures <- df_results_long %>%
  mutate(Failure=ifelse(Failure.type=="Nonfailure", 0, 1)) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  summarise(
    Num.iterations=n(),
    Num.failures=sum(Failure),
    Failure.rate=mean(Failure)) %>%
  ungroup()

# save processed data and clear from memory
saveRDS(df_failures, here("DataProcessed/failures_2023-02-20.rds"))
rm(df_failures)

# for 'oracle' inference method: use true correlation
df_oracle <- df_results_long %>%
  filter(Method=="Rho.hat") %>%
  mutate(Method="Oracle")

# for naive inference method: assume zero correlation
df_naive <- df_results_long %>%
  filter(Method=="Rho.hat") %>%
  mutate(Estimate=0) %>%
  mutate(Method="Independent")

# compute (modified) t-statistics and corresponding p-values
df_stats <- df_results_long %>%
  bind_rows(df_oracle, df_naive) %>%
  mutate(Estimate=ifelse(Failure.type=="Invalid", NA, Estimate)) %>%
  mutate(Estimate=ifelse(Failure.type=="Boundary", Estimate * 0.9999, Estimate)) %>%
  mutate(Failure=ifelse(Failure.type %in% c("Boundary", "Nonfailure"), 0, 1)) %>%
  mutate(
    d.bar = Mu.hat.X - Mu.hat.Y,
    sum.sqs = (Sigma.sqd.hat.X * N) + (Sigma.sqd.hat.Y * N),
    se.ind = sqrt( sum.sqs / (N * (N - 1)) ),
    se.mod = suppressWarnings(ifelse(Failure==1, NA, se.ind * sqrt(1 - Estimate))),
    t.ind = d.bar / se.ind,
    t.mod = suppressWarnings(ifelse(Failure==1, NA, t.ind / sqrt(1 - Estimate))),
    p.value = 2 * pt(abs(t.mod), df=(2 * N - 2), lower.tail=F),
    d.lwr = d.bar - 1.96 * se.mod,
    d.upr = d.bar + 1.96 * se.mod)

# apply @Ryan's correction for scaling the t statistic for Pearson method
df_pearson_scaled_stats <- df_stats %>%
  filter(Method=="Pearson") %>%
  mutate(t.mod=case_when(
    !is.na(t.mod) & M>=3 ~
      t.mod / sqrt( 0.5 * (1 + ((1-Estimate)/(1+Estimate)) * exp(2/(M-3))) ),
    TRUE ~ t.mod)) %>%
  mutate(p.value = 2 * pt(abs(t.mod), df=(2 * N - 2), lower.tail=F)) %>%
  mutate(Method="Pearson.scaled")

# compute rejection rates, average standard errors and coverage
df_inference <- df_stats %>%
  bind_rows(df_pearson_scaled_stats) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  summarise(
    Rejection.rate=mean(p.value < 0.05, na.rm=T),
    SE.mean=mean(se.mod, na.rm=T),
    Coverage.rate=mean(Delta >= d.lwr & Delta <= d.upr, na.rm=T)) %>%
  ungroup()

# save processed data and clear from memory
saveRDS(df_inference, here("DataProcessed/inference_2023-02-20.rds"))
rm(df_inference)

# calculate bias and MSE, including "monte carlo error" (estimates)
df_performance <- df_results_long %>%
  mutate(Error=Estimate - Rho) %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  summarise(
    Bias=mean(Error, na.rm=T),
    Bias.mce=sqrt(var(Error, na.rm=T)),
    Mean.sqd.error=mean(Error^2, na.rm=T),
    Mean.sqd.error.mce=sqrt(var(Error^2, na.rm=T))) %>%
  ungroup()

# save processed data and clear from memory
saveRDS(df_performance, here("DataProcessed/performance_2023-02-20.rds"))
rm(df_performance)

# save processed data (commented out because too big)
# saveRDS(df_results_long, here("DataProcessed/all_2023-02-20.rds"))

# save 'lite' version of long data for improved speed
set.seed(2023)
df_results_lite <- df_results_long %>%
  group_by(Method, Distribution, Rho, Delta, N, Prop.matched, M) %>%
  slice_sample(prop=0.05) %>%
  ungroup()
saveRDS(df_results_lite, here("DataProcessed/lite_2023-02-20.rds"))
