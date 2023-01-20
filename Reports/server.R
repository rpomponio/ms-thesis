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

# eventually: load processed results for better performance
results <- readRDS(here("DataRaw/simulation_results_2023-01-18.rds"))

# Define server logic required to draw a histogram
function(input, output, session) {

    output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        x    <- results[, "Rho.hat"]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')

    })
    
    output$moreControls <- renderUI({
      sliderInput("N",
                  "Sample size",
                  min = min(results[, "N"]),
                  max = max(results[, "N"]),
                  value = 10)
    })

}
