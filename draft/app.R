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
library(boot)
library(ggplot2)

source("funcs.R")


metric_list <- c("Bias", "SE", "MSE")


# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("SMART Missingness Simulator"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      radioButtons("pct_missing",
        "% Missing Data:",
        c("30", "60"),
        selected = "30"
      ),
      radioButtons("mechanism",
        "Missingness Mechanism",
        c("MCAR", "MAR", "MNAR"),
        selected = "MAR"
      ),
      numericInput("n_sims",
        "Number of Simulations:",
        value = 2,
        min = 2,
        max = 100
      ),
      actionButton("simulate", "Simulate!"),
      h3("How to use this app"),
      p(
        "Use the input panel (above) to select
         your desired missingness mechanism (missing completely at random,
         missing at random, or missing not at random), percentage of
         missing data, and the number of simulations to generate. Then,
         press 'Simulate!' to begin simulations. Once completed (it takes
         a while to run, be patient!), results will display on the 'Plot'
         and 'Table' tabs."
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel(
          "About",
          h3("Background"),
          p(
            "Missing data can compromise the validity of inference,
             particularly when the missingness is not random (MNAR).
             The complex structure of sequential multiple assignment
             randomized trials (SMART) presents unique challenges
             when handling missing data."
          ),
          p(
            "In other longitudinal contexts, multiple imputation
             (MI) has been widely accepted as a flexible framework
             to handle missingness when the missingness mechanism is
             missing completely at random (MCAR) or missing at
             random (MAR)."
          ),
          p(
            "While multiple imputation may provide some
             protection when the MAR assumption is violated in some
             contexts, it is unclear how these violations will
             affect an analysis in a SMART setting."
          ),
          p(
            "This Shiny application aims to evaluate the performance of MI under
             various randomness mechanisms
             and levels of missingness."
          ),
          h3("Methods"),
          p(
            "We consider a two-stage SMART study with three
             study visits with data that mimics the Clinical
             Antipsychotic Trials of Intervention Effectiveness
             (CATIE) SMART study for patients with Schizophrenia."
          ),
          p(
            "In this study Subjects were randomized to three possible first line
             treatments: Olanzapine, Risperidone, and Quetiapine (second
             generation or atypical antipyshcotics). For the second treatment
             stage, responders will remain on their first line treatment and
             non-responders will be assigned treatment according to a
             randomization scheme depending on the reason for non-response
             (efficacy or tolerability). Both atypical antipsychotics and
             typical (first generation) antipsychotics are possible second line
             treatments for non-responders. This gives seven total treatment
             regimes."
          ),
          p(
            "The study's response of interest is the subject's symptom severity
             score  (PANSS). A higher score indicates a worsening of symptoms.
             We will look at the average PANSS score over the three study visits
             as our response."
          ),
          p(
            "In all scenarios, we assumed monotone missingness. To generate
             MAR data we let the probability of a patient dropping out at time
             t be a function of their BMI, years on prescription medication
             prior to the study, proportion of study pills taken (adherence) in
             the past treatment stage, and their most recently observed PANSS
             score. For MNAR data, the probability of dropping out at time t is
             a function of the complete data. The primary factors that dictate
             this probability are response status, reason for non-response,
             BMI,  adherence, and their final (unobserved) PANSS score."
          ),
          p(
            "To control the percentage of missing data, we varied the
              coefficients in our missing data models such that there was
              10% or 20% dropout during stage 1, and then an additional
              10% or 20% dropout at the end of stage 1, during stage 2,
              and at the end of stage 2, for a total dropout rate of 30% or 60%."
          ),
          p(
            "The multiple imputation procedure was based on a time-ordered
           nested multiple imputation stategy proposed for SMART studies in
           Shortreed et al. (2014)."
          ),
          h3("Future Work"),
          p(
            "In later versions of this app, we hope to allow the user to...",
            tags$div( # HTML bulleted list
              tags$ul(
                tags$li("vary the sample size"),
                tags$li("control which metrics are reported and which are hidden"),
                tags$li("select missingness percentages other than 30/60%"),
                tags$li("see the results of a complete case analysis for comparison"),
                tags$li("input multiple lists of paramaters and compare between settings"),
                tags$li("utilize parallel computing to simulate faster")
              )
            )
          )
        ),
        tabPanel("Plot", plotOutput("plot")),
        tabPanel("Table", tableOutput("table1"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  # Run simulations whenever simulate button is pressed
  re <- eventReactive(
    input$simulate,
    {
      # Run simulation ----
      load("SMARTdat.rda")
      true.dat <- t
      true.dat$SWITCH <- ifelse(true.dat$SWITCH == "SWITCHED", 1, 0)

      sim.res <- data.frame(t(sapply(1:input$n_sims, main,
        pct_mis = input$pct_missing,
        mis_mec = input$mechanism,
        samp_size = 1000,
        true.dat = true.dat
      )))

      # Compute summary measures ----
      res <- sapply(sim.res, mean)
      names(res) <- c(
        "O1", "O2", "R1", "R2", "Q1", "Q2", "Q3", "C.O1", "C.O2",
        "C.R1", "C.R2", "C.Q1", "C.Q2", "C.Q3"
      )

      # Metrics send to user
      bias <- res[1:7] - true.means
      se <- sapply(sim.res[1:7], sd)
      cov.probability <- res[8:14]
      mse <- c()
      for (i in 1:7) {
        m <- mean((sim.res[, i] - true.means[i])^2)
        mse <- c(m, mse)
      }

      # Plot
      Regime <- c(
        rep("DTR1", input$n_sims), rep("DTR2", input$n_sims), rep("DTR3", input$n_sims),
        rep("DTR4", input$n_sims), rep("DTR5", input$n_sims), rep("DTR6", input$n_sims),
        rep("DTR7", input$n_sims)
      )

      reg.means <- c(
        sim.res[, 1], sim.res[, 2], sim.res[, 3], sim.res[, 4], sim.res[, 5],
        sim.res[, 6], sim.res[, 7]
      )

      plot.df <- data.frame(Regime, reg.means)

      plot <- ggplot(plot.df, aes(Regime, reg.means, fill = Regime)) +
        geom_boxplot() +
        labs(
          y = "Expected PANSS score", x = "Embedded DTR",
          title = "Comparing embedded DTRs"
        )

      # Add the true means to the plot (red dots)
      plot <-
        plot +
        annotate("point", x = "DTR1", y = true.means[1], colour = "red", size = 3) +
        annotate("point", x = "DTR2", y = true.means[2], colour = "red", size = 3) +
        annotate("point", x = "DTR3", y = true.means[3], colour = "red", size = 3) +
        annotate("point", x = "DTR4", y = true.means[4], colour = "red", size = 3) +
        annotate("point", x = "DTR5", y = true.means[5], colour = "red", size = 3) +
        annotate("point", x = "DTR6", y = true.means[6], colour = "red", size = 3) +
        annotate("point", x = "DTR7", y = true.means[7], colour = "red", size = 3) +
        theme(legend.position = "none")

      # Return all results in a list
      list(
        bias = bias,
        se = se,
        mse = mse,
        cov.probability = cov.probability,
        plot = plot
      )
    }
  )

  output$plot <- renderPlot(
    re()$plot
  )

  output$table1 <- renderTable(
    expr = {
      rbind(re()$bias, re()$se, re()$mse) %>%
        `rownames<-`(c("Bias", "SE", "MSE")) %>%
        `colnames<-`(paste0("DTR", 1:7))
    },
    striped = TRUE,
    hover = TRUE,
    rownames = TRUE
  )
}

# Run the application
shinyApp(ui = ui, server = server)
