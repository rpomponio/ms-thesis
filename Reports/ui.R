#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
fluidPage(
  titlePanel("Simulation Results"),
  sidebarLayout(
    sidebarPanel(
      selectInput("distribution", "Distribution type", c("Normal", "Ordinal")),
      checkboxGroupInput("delta", "Mean difference", c(0, 0.25, 0.5), c(0)),
      selectInput("n", "(Group) sample size", c(10, 20, 50, 100, 200)),
      selectInput("perc.matched", "Prop. matched", c("0%", "5%", "10%", "15%",
                                                     "20%", "30%", "50%", "100%"),
                  c("30%")),
      selectizeInput("method", "Estimators",
                     c("Pearson", "Max.conserv", "EM.alg",
                       "Shrunken", "Unbiased", "Boot.mean",
                       "Boot.5th.quantile", "Boot.20th.quantile",
                       "Freq.20th.quantile", "Bayes.arcsine", "Bayes.Jeffreys",
                       "Bayes.unif"),
                     c("Pearson", "Max.conserv", "EM.alg"), multiple=T),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Dummy dataset",
                 verticalLayout(
                   sliderInput("rho", "Correlation", -0.95, 0.95, 0.25, 0.05),
                   plotOutput("scatterplot", width="80%"))),
        tabPanel("Bias", plotOutput("biasplot")),
        tabPanel("Variance", plotOutput("varianceplot"))
      )
    )
  )
)
