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
      selectInput("n", "Sample size (per group)", c(10, 20, 50, 100, 200), c(20)),
      selectInput("perc.matched", "Prop. matched (b/w groups)",
                  c("0%", "5%", "10%", "15%", "20%", "30%", "50%", "100%"),
                  c("30%")),
      wellPanel("Number of matched samples:", textOutput("n.matched")),
      selectInput("rho", "Correlation",
                  c(-0.9, -0.5, -0.25, 0, 0.25, 0.5, 0.9),
                  c(0.5)),
      selectizeInput("method", "Estimators",
                     c("Rho.hat", "Pearson", "Max.conserv", "EM.alg",
                       "Shrunken", "Unbiased", "Boot.mean",
                       "Boot.5th.quantile", "Boot.20th.quantile",
                       "Freq.20th.quantile", "Bayes.arcsine", "Bayes.Jeffreys",
                       "Bayes.unif"),
                     c("Pearson", "EM.alg", "Bayes.Jeffreys"), multiple=T),
      wellPanel(
        "Power plot settings",
        checkboxInput("plot.oracle", "Plot 'oracle' test", value=T),
        checkboxInput("plot.indep", "Plot independent t-test", value=F)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Dummy dataset", plotOutput("scatterplot", width="90%", height="500px")),
        tabPanel("Failures", verticalLayout(
          p("Failures were recorded when an invalid correlation estimate was obtained."),
          DT::dataTableOutput("failures"))),
        tabPanel("Bias", plotOutput("biasplot", height="600px")),
        tabPanel("Variance", plotOutput("varianceplot", height="600px")),
        tabPanel("Comparison", verticalLayout(
          p("Compares the first two estimators selected."),
          plotOutput("comparisonplot", height="600px"))),
        tabPanel("Power", plotOutput("powerplot", height="600px")),
        tabPanel("Estimators", DT::dataTableOutput("methods"))
      )
    )
  )
)