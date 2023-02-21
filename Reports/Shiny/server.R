#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(here)
library(mvtnorm)
library(dplyr)
library(ggpubr)
library(ggplot2)
theme_set(theme_classic2(base_size=18))

# load processed results
df_results_long <- readRDS(here("DataProcessed/lite_2023-02-20.rds"))
df_failures <- readRDS(here("DataProcessed/failures_2023-02-20.rds"))
df_inference <- readRDS(here("DataProcessed/inference_2023-02-20.rds"))
df_performance <- readRDS(here("DataProcessed/performance_2023-02-20.rds"))

# Define server logic required to draw a histogram
function(input, output, session) {
  
  # reactive values continually update
  vals <- reactiveValues()
  observe({
    vals$prop.matched <- as.numeric(gsub("%", "", input$perc.matched)) / 100
    vals$n.matched <- floor(
      (as.numeric(gsub("%", "", input$perc.matched)) / 100) * as.numeric(input$n))
  })
  
  output$n.matched <- renderText({
    vals$n.matched
  })

  output$scatterplot <- renderPlot({
    
    set.seed(2023)
    
    Sigma.X <- 1
    Sigma.Y <- 1
    
    # generate a dataset
    S <- diag(2)
    diag(S) <- c(Sigma.X, Sigma.Y)
    S[upper.tri(S) | lower.tri(S)] <- rep(as.numeric(input$rho)) * Sigma.X * Sigma.Y
    sim_mat <- rmvnorm(
      n=as.numeric(input$n),
      mean=c(0, as.numeric(input$delta)[1]),
      sigma=S,
      checkSymmetry=FALSE)
    X <- sim_mat[, 1]
    Y <- sim_mat[, 2]
    if (input$distribution=="Ordinal") {
      cut_seq <- c(-999, c(-2, 1, 0, 0.7, 1.3, 2.0), 999)
      X <- as.numeric(cut(X, cut_seq))
      Y <- as.numeric(cut(Y, cut_seq))
    }
    
    data.frame(X=X, Y=Y, matched=ifelse(1:input$n<=vals$n.matched, "Y", "N")) %>%
      ggplot(aes(X, Y)) +
      geom_point(aes(col=matched, shape=matched), size=5) +
      geom_smooth(formula="y ~ x", method="lm", se=F, linetype="dashed") +
      labs(title="Scatterplot of dummy dataset",
           x="Pre-intervention",
           y="Post-intervention",
           caption="Regression line added for illustration only.")
  })
  
  output$failures <- DT::renderDataTable({
    
    df_failures %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             Rho == as.numeric(input$rho),
             N == as.numeric(input$n),
             Method %in% input$method,
             Prop.matched == vals$prop.matched) %>%
      select(Method, Rho, Delta, N, M, Num.iterations, Num.failures, Failure.rate)
  })
  
  output$biasplot <- renderPlot({
    
    sim_note <- paste0("Simulated with n=", as.numeric(input$n),
                       ", m=", vals$n.matched)
    
    df_performance %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% input$method,
             Prop.matched == vals$prop.matched) %>%
      ggplot(aes(col=Method, x=Rho, y=Bias)) +
      geom_hline(yintercept=0, linetype="dashed") +
      facet_wrap(~ Delta, scales="free_y") +
      geom_line() +
      geom_point(size=5) +
      scale_x_continuous(breaks=unique(df_performance$Rho)) +
      labs(title="Bias of estimators by true correlation",
           subtitle=sim_note,
           x="True correlation",
           y="Bias (Estimated minus true value)",
           caption="Results averaged over 10,000 datasets at each point.")
  })
  
  output$biastable <- DT::renderDataTable({
    
    df_performance %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% input$method,
             Prop.matched == vals$prop.matched) %>%
      select(Delta, Method, Rho, N, M, Bias) %>%
      tidyr::pivot_wider(names_from=Rho, values_from=Bias) %>%
      DT::datatable(options=list(pageLength=15, dom='tip'),
                    caption="Bias for values of true correlation") %>%
      DT::formatRound(columns=5:11)
  })
  
  output$varianceplot <- renderPlot({
    
    sim_note <- paste0("Simulated with n=", as.numeric(input$n),
                       ", m=", vals$n.matched)
    
    df_performance %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% input$method,
             Prop.matched == vals$prop.matched) %>%
      ggplot(aes(col=Method, x=Rho, y=Mean.sqd.error)) +
      geom_hline(yintercept=0, linetype="dashed") +
      facet_wrap(~ Delta, scales="free_y") +
      geom_line() +
      geom_point(size=5) +
      scale_x_continuous(breaks=unique(df_performance$Rho)) +
      labs(title="MSE of estimators by true correlation",
           subtitle=sim_note,
           x="True correlation",
           y="Mean squared error",
           caption="Results averaged over 10,000 datasets at each point.")
  })
  
  output$variancetable <- DT::renderDataTable({
    
    df_performance %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% input$method,
             Prop.matched == vals$prop.matched) %>%
      select(Delta, Method, Rho, N, M, Mean.sqd.error) %>%
      tidyr::pivot_wider(names_from=Rho, values_from=Mean.sqd.error) %>%
      DT::datatable(options=list(pageLength=15, dom='tip'),
                    caption="MSE for values of true correlation") %>%
      DT::formatRound(columns=5:11)
  })
  
  output$seplot <- renderPlot({
    
    addl.methods <- c()
    if (input$plot.oracle) {
      addl.methods <- c(addl.methods, "Oracle")
    }
    if (input$plot.indep) {
      addl.methods <- c(addl.methods, "Independent")
    }
    
    sim_note <- paste0("Simulated with n=", as.numeric(input$n),
                       ", m=", vals$n.matched)
    
    df_inference %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% c(addl.methods, input$method),
             Prop.matched == vals$prop.matched) %>%
      ggplot(aes(col=Method, x=Rho, y=SE.mean)) +
      geom_hline(yintercept=0, linetype="dashed") +
      facet_wrap(~ Delta, scales="free_y") +
      geom_line() +
      geom_point(size=5) +
      scale_x_continuous(breaks=unique(df_inference$Rho)) +
      labs(title="Avgerage standard errors by true correlation",
           subtitle=sim_note,
           x="True correlation",
           y="Standard Error",
           caption="Results averaged over 10,000 datasets at each point.")
  })
  
  output$setable <- DT::renderDataTable({
    
    addl.methods <- c()
    if (input$plot.oracle) {
      addl.methods <- c(addl.methods, "Oracle")
    }
    if (input$plot.indep) {
      addl.methods <- c(addl.methods, "Independent")
    }
    
    df_inference %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% c(addl.methods, input$method),
             Prop.matched == vals$prop.matched) %>%
      select(Delta, Method, Rho, N, M, SE.mean) %>%
      tidyr::pivot_wider(names_from=Rho, values_from=SE.mean) %>%
      DT::datatable(options=list(pageLength=15, dom='tip'),
                    caption="Mean SEs for values of true correlation") %>%
      DT::formatRound(columns=5:11)
  })

  output$powerplot <- renderPlot({
    
    addl.methods <- c()
    if (input$plot.oracle) {
      addl.methods <- c(addl.methods, "Oracle")
    }
    if (input$plot.indep) {
      addl.methods <- c(addl.methods, "Independent")
    }
    
    sim_note <- paste0("Simulated with n=", as.numeric(input$n),
                       ", m=", vals$n.matched)
    
    df_inference %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% c(addl.methods, input$method),
             Prop.matched == vals$prop.matched) %>%
      ggplot(aes(col=Method, x=Rho, y=Rejection.rate)) +
      geom_hline(yintercept=0.05, linetype="dashed") +
      facet_wrap(~ Delta, scales="free_y") +
      geom_line() +
      geom_point(size=5) +
      scale_x_continuous(breaks=unique(df_inference$Rho)) +
      labs(title="Rejection rate of hypoth. tests by true correlation",
           subtitle=sim_note,
           x="True correlation",
           y="Prop. rejected null hypotheses",
           caption="Results averaged over 10,000 datasets at each point.")
  })
  
  output$powertable <- DT::renderDataTable({
    
    addl.methods <- c()
    if (input$plot.oracle) {
      addl.methods <- c(addl.methods, "Oracle")
    }
    if (input$plot.indep) {
      addl.methods <- c(addl.methods, "Independent")
    }
    
    df_inference %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% c(addl.methods, input$method),
             Prop.matched == vals$prop.matched) %>%
      select(Delta, Method, Rho, N, M, Rejection.rate) %>%
      tidyr::pivot_wider(names_from=Rho, values_from=Rejection.rate) %>%
      DT::datatable(options=list(pageLength=15, dom='tip'),
                    caption="Rejection rate for values of true correlation") %>%
      DT::formatRound(columns=5:11)
  })
  
  output$comparisonplot <- renderPlot({
    
    df_results_long %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% c(input$method[1:2]),
             Prop.matched == vals$prop.matched) %>%
      tidyr::pivot_wider(names_from=Method, values_from=Estimate) %>%
      ggplot(aes_string(x=input$method[1], y=input$method[2])) +
      facet_wrap(~ Delta) +
      geom_point(size=5, alpha=0.25) +
      geom_smooth(method="lm", formula="y ~ x", se=F) +
      geom_abline(intercept=0, slope=1, linetype="dashed") +
      labs(title="Scatterplot of association between estimators",
           caption="Results sampled from 10,000 datasets for efficiency.")
    
  })
  
  output$methods <- DT::renderDataTable({
    
    data.frame(
      Method.Name=c("Rho.hat", "Pearson", "Max.conserv", "EM.alg",
                    "Shrunken", "Unbiased", "Freq.20th.quantile",
                    "Bayes.arcsine", "Bayes.Jeffreys", "Bayes.unif"),
      Description=c(
        "Maximum likelihood estimate under fully matched data (for comparison only)",
        "Pearson correlation coefficient of matched samples",
        "Minimum possible correlation given all samples (including unmatched)",
        "EM-based maximum likelihood estimate of correlation given unmatched data",
        "The root of adjusted R-squared, computed on matched samples",
        "Approximately unbiased estimate of correlation, computed on matched samples",
        "Lower bound of an 80% frequentist interval for the correlation coefficient",
        "Approximate Bayesian posterior mean correlation, assuming arcsine prior",
        "Approximate Bayesian posterior mean correlation, assuming Jeffreys prior",
        "Approximate Bayesian posterior mean correlation, assuming Uniform prior")
    )
  }, options=list(pageLength=15, dom='tip'))
  
  
}
