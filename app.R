#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
load("SMARTdat.rda")
source("funcs.R")
true.dat <- t


metric_list <- c("MSE", "Coverage Probability")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("SMART Missingness Simulator"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            # sliderInput("pct_missing",
            #             "% Missing Data:",
            #             min = 0,
            #             max = 99,
            #             value = 10),
            radioButtons("pct_missing",
                         "% Missing Data:",
                         c("30", "60"),
                         selected = "30"),
            radioButtons("mechanism",
                         "Missingness Mechanism",
                         c("MCAR", "MAR", "MNAR"),
                         selected = "MAR"),
            numericInput("n",
                        "Sample Size:",
                        min = 100,
                        max = 1000,
                        value = 200),
            numericInput("n_sims",
                         "Number of Simulations:",
                         value = 50,
                         min = 2,
                         max = 100),
            actionButton("simulate", "Simulate!"),
            checkboxGroupInput("metrics",
                               "Which metrics?",
                               metric_list,
                               selected = metric_list)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           textOutput("test_message")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    # Run simulations whenever simulate button is pressed
    re <- reactive({
        input$simulate
        data.frame(t(sapply(1:input$n_sims, main, 
                            pct_mis = input$pct_missing, 
                            mis_mec = input$mechanism, 
                            samp_size = input$n)))
    })

    output$test_message <- renderText(
        paste0("You have selected ",
               input$pct_missing, "% missingness, ",
               input$mechanism, " as the missingness mechanism, and ",
               input$n, " for the sample size. ",
               "To be carried out over ",
               input$n_sims, " simulations with the following metrics displayed: ",
               paste0(input$metrics, collapse = ", "), ".")
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
