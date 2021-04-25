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
true.dat <- t
source("funcs.R")


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
            # numericInput("n",
            #             "Sample Size:",
            #             min = 100,
            #             max = 1000,
            #             value = 200),
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
        # Run simulation ----
        sim.res <- data.frame(t(sapply(1:input$n_sims, main, 
                            pct_mis = input$pct_missing, 
                            mis_mec = input$mechanism, 
                            samp_size = 1000)))
        # Compute summary measures ----
        res <- sapply(sim.res, mean)
        names(res)=c("O1","O2","R1","R2","Q1","Q2","Q3","C.O1","C.O2","C.R1","C.R2","C.Q1","C.Q2","C.Q3")
        
        #Metrics send to user
        bias <- res[1:7]-true.means
        se <- sapply(sim.res[1:7],sd)
        cov.probability <- res[8:14]
        mse <- c()
        for (i in 1:7){
            m <- mean( (sim.res[,i] - true.means[i])^2 )
            mse <- c(m,mse)}
        
        #Plot
        Regime <- c(rep("DTR1",n_sims),rep("DTR2",n_sims),rep("DTR3",n_sims),
                    rep("DTR4",n_sims),rep("DTR5",n_sims),rep("DTR6",n_sims),
                    rep("DTR7",n_sims))
        
        reg.means <- c(sim.res[,1],sim.res[,2],sim.res[,3],sim.res[,4],sim.res[,5],
                       sim.res[,6],sim.res[,7])
        
        plot.df <- data.frame(Regime, reg.means)
        
        plot = ggplot(plot.df, aes(Regime, reg.means, fill=Regime)) +
            geom_boxplot() +
            labs(y = "Expected PANSS score", x="Embedded DTR",
                 title="Comparing embedded DTRs")
        
        #Add the true means to the plot (red dots)
        plot +
            annotate("point", x = "DTR1", y = true.means[1], colour = "red",size=3)+
            annotate("point", x = "DTR2", y = true.means[2], colour = "red",size=3)+
            annotate("point", x = "DTR3", y = true.means[3], colour = "red",size=3)+
            annotate("point", x = "DTR4", y = true.means[4], colour = "red",size=3)+
            annotate("point", x = "DTR5", y = true.means[5], colour = "red",size=3)+
            annotate("point", x = "DTR6", y = true.means[6], colour = "red",size=3)+
            annotate("point", x = "DTR7", y = true.means[7], colour = "red",size=3)+
            theme(legend.position = "none")
        
        # Return all results in a list
        list(
            bias = bias,
            se = se,
            cov.probability = cov.probability,
            mse = mse,
            plot = plot
        )
        
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
