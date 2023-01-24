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

# eventually: load processed results for better performance
df_performance <- readRDS(here("DataProcessed/performance_2023-01-21.rds"))

# Define server logic required to draw a histogram
function(input, output, session) {
  
  output$scatterplot <- renderPlot({
    
    set.seed(2023)
    
    prop.matched <- as.numeric(gsub("%", "", input$perc.matched)) / 100
    
    Sigma.X <- 1
    Sigma.Y <- 1
    
    # generate a dataset
    S <- diag(2)
    diag(S) <- c(Sigma.X, Sigma.Y)
    S[upper.tri(S) | lower.tri(S)] <- rep(input$rho) * Sigma.X * Sigma.Y
    sim_mat <- rmvnorm(
      n=as.numeric(input$n),
      mean=c(0, as.numeric(input$delta)),
      sigma=S,
      checkSymmetry=FALSE)
    X <- sim_mat[, 1]
    Y <- sim_mat[, 2]
    if (input$distribution=="Ordinal") {
      cut_seq <- c(-999, c(-2, 1, 0, 0.7, 1.3, 2.0), 999)
      X <- as.numeric(cut(X, cut_seq))
      Y <- as.numeric(cut(Y, cut_seq))
    }
    n_matched <- floor(prop.matched * as.numeric(input$n))

    data.frame(X=X, Y=Y, matched=ifelse(1:input$n<=n_matched, "Y", "N")) %>%
      ggplot(aes(X, Y)) +
      geom_point(aes(col=matched, shape=matched), size=5) +
      geom_smooth(formula="y ~ x", method="lm", se=F, linetype="dashed") +
      labs(title="Scatterplot of dummy dataset",
           x="Pre-intervention",
           y="Post-intervention",
           caption="Regression line added for illustration only.")
  })
  
  output$biasplot <- renderPlot({
    
    prop.matched <- as.numeric(gsub("%", "", input$perc.matched)) / 100
    n_matched <- floor(prop.matched * as.numeric(input$n))
    
    df_performance %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% input$method,
             M == n_matched) %>%
      ggplot(aes(col=Method, x=Rho, y=Bias)) +
      geom_hline(yintercept=0, linetype="dashed") +
      # facet_wrap(~ M, scales="free_y") +
      geom_line() +
      geom_point(size=5) +
      scale_x_continuous(breaks=unique(df_performance$Rho)) +
      labs(title="Bias of estimators by true correlation",
           x="True correlation",
           y="Bias (Estimated minus true value)",
           caption="Results averaged over 1,000 datasets at each point.")
  })
  
  output$varianceplot <- renderPlot({
    
    prop.matched <- as.numeric(gsub("%", "", input$perc.matched)) / 100
    n_matched <- floor(prop.matched * as.numeric(input$n))
    
    df_performance %>%
      filter(Distribution==ifelse(input$distribution=="Normal", 1, 2),
             Delta %in% as.numeric(input$delta),
             N == as.numeric(input$n),
             Method %in% input$method,
             M == n_matched) %>%
      ggplot(aes(col=Method, x=Rho, y=Mean.sqd.error)) +
      # facet_wrap(~ M, scales="free_y") +
      geom_line() +
      geom_point(size=5) +
      scale_x_continuous(breaks=unique(df_performance$Rho)) +
      labs(title="Variance of estimators by true correlation",
           x="True correlation",
           caption="Results averaged over 1,000 datasets at each point.")
  })
  

}
