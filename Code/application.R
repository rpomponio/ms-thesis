library(mvtnorm)
library(here)

source(here("Code/estimators.R"))

load(here("DataRaw/score.RData"))

# how many repeated IDs in the control group?
control <- subset(score, cluster=="Control")
table(control$time)

# how many repeated IDs in the intervention group?
intervention <- subset(score, cluster=="Intervention")
table(intervention$time)

# define helper function to run all tests on application dataset
run.all.tests <- function(df) {
  df_t1 <- subset(df, time==1)
  df_t2 <- subset(df, time==2)
  # re order
  df_t1 <- df_t1[order(df_t1$anonymous_identifier), ]
  df_t2 <- df_t2[order(df_t2$anonymous_identifier), ]
  # estimate correlation using all candidate estimators
  matched_ids <- intersect(df_t1$anonymous_identifier,
                           df_t2$anonymous_identifier)
  n_matched <- length(matched_ids)
  X.matched <- subset(df_t1, anonymous_identifier %in% matched_ids)$sum_score
  Y.matched <- subset(df_t2, anonymous_identifier %in% matched_ids)$sum_score
  X.unmatched <- subset(df_t1, !anonymous_identifier %in% matched_ids)$sum_score
  Y.unmatched <- subset(df_t2, !anonymous_identifier %in% matched_ids)$sum_score
  X <- c(X.matched, X.unmatched)
  Y <- c(Y.matched, Y.unmatched)
  res_df <- data.frame(
    "Method"=c("Max.conserv", "Pearson", "20th.quantile", "Bayes", "EM",
               "Independent"),
    "Rho.estimate"=c(
      est.cor.conserv(X[1:length(Y)], Y),
      est.cor.pearson(X, Y, n.matched=n_matched),
      est.cor.quantile(X, Y, n.matched=n_matched),
      est.cor.bayesian.asin(X, Y, n.matched=n_matched),
      NA, # est.cor.emalg(X, Y, n.matched=n_matched)
      0),
    check.names=F)
  res_df$Delta.estimate <- mean(X) - mean(Y)
  res_df$T.student <- t.test(X, Y, var.equal=T)$stat
  res_df$T.df <- t.test(X, Y, var.equal=T)$parameter
  res_df$T.modified <- res_df$T.student / sqrt(1 - res_df$Rho.estimate)
  res_df$SE <- t.test(X, Y, var.equal=T)$stderr * sqrt(1 - res_df$Rho.estimate)
  res_df$CI.lwr <- res_df$Delta.estimate - 1.96 * res_df$SE
  res_df$CI.upr <- res_df$Delta.estimate + 1.96 * res_df$SE
  res_df$p.value <- 2 * pt(abs(res_df$T.modified), df=res_df$T.df, lower.tail=F)
  
  return(res_df)
}

# sanity check with traditional t-tests
t.test(sum_score ~ time, data=control, var.equal=T)
t.test(sum_score ~ time, data=intervention, var.equal=T)


run.all.tests(control)
run.all.tests(intervention)







