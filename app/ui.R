#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)

#options(shiny.maxRequestSize=1000*1024^2)
shinyUI(
    
    
    # Sidebar with a slider input for number of bins
    navbarPage(
        # Application title
        title="Surveillance test accuracy in an epidemic",
        tabPanel("Population-level statistics",
               
                 sidebarPanel(
                     helpText("This page plots the incidence, prevalence,
                          percentage positive from different testing strategies,
                          the positive predictive values and negative predictive 
                          values. Please refer to the accompanying manuscript,
                          Figure 1 legend, for further detail."),
                     sliderInput("growth_rate",
                                 "Exponential growth rate:",
                                 min = -0.1,
                                 max = 0.1,
                                 value = -0.05,
                                 step=0.01),
                 fluidRow(
                     column(4,numericInput("cumu_inc",
                                 "Cumulative incidence:",
                                 min = 0,
                                 max = 0.1,
                                 value = 0.01,
                                 step=0.001)),
                     column(4,numericInput("N",
                                 "Population size:",
                                 min = 10000,
                                 max = 10000000,
                                 value = 1500000,
                                 step=10000)),
                     column(4,numericInput("sim_duration",
                                 "Duration of simulation:",
                                 min = 10,
                                 max = 100,
                                 value = 55,
                                 step=5))
                     ),
                 
                 sliderInput("test_frequency",
                             "Test frequency (days):",
                             min = 1,
                             max = 10,
                             value = 3,
                             step=1),
                 fluidRow(
                     column(6,numericInput("strategy_start_time",
                                  "Start of routine screening:",
                                  min = 0,
                                  max = 100,
                                  value = 9,
                                  step=1)),
                     column(6,numericInput("isolation_period",
                                  "Isolation period:",
                                  min = 1,
                                  max = 30,
                                  value = 10,
                                  step=1))
                 ),
                     hr(),
                     width=3
                 ),
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("all_plots",width="100%")
        )
    ),
    tabPanel("Simulated cohorts",
             sidebarPanel(
                 helpText("This page simulates an agent-based model (with a smaller population size), 
                          equivalent to the previous page but tracking
                          all individual infection statues. The top row shows
                          individual-level infection statues over time
                          as observed by different surveillance and screening 
                          strategies. The middle row shows the number of individuals 
                          in each infection states on the specified observation day 
                          (red vertical line). The bottom row shows the number
                          of false and true positives on each day of the routine screening
                          strategy, as well as the number of individuals in isolation over
                          time."),
                 hr(),
                 helpText("NOTE: this page runs an agent-based model, so will take a while to load if 
                          the subsample fraction is set too high (>0.01). Increasing this value will make the summary statistics more accurate. Please
                          be patient!"),
                 hr(),
             fluidRow(
                 column(6,numericInput("subsamp_frac",
                          "Subsample population for cohort simulation:",
                          min = 0.01,
                          max = 1,
                          value = 0.1,
                          step=0.01)),
                 column(6,  numericInput("t_choose_index",
                                         "Cross section index (testing rounds):",
                                         min = 1,
                                         max = 100,
                                         value = 4,
                                         step=1))
             ), 
             hr(),
             width=3
    ),
             # Show a plot of the generated distribution
             mainPanel(
                 #helpText("Due to memory constraints, this part of the simulation cannot be run online. Please download the github repository containing this code and run the app locally, following the instructions on the github README.")
                 plotOutput("linelist_plot",width="100%")
             )
    ),  
    tabPanel("Test characteristics",
             sidebarPanel(
                 helpText("Vary assumptions regarding qPCR and LFT test characteristics
                          underpinning all simulations."),
                 hr(),
                 sliderInput("duration_pos",
                             "Duration of \"true\" infection:",
                             min = 5,
                             max = 30,
                             value = 21,
                             step=1),
                 sliderInput("pcr_gamma_mode",
                             "Day of peak sensitivity for both tests:",
                             min = 2,
                             max = 15,
                             value = 5,
                             step=1),
                 fluidRow(
                     column(4,numericInput("pcr_gamma_sd",
                                          "PCR width:",
                                          min = 1,
                                          max = 25,
                                          value = 10,
                                          step=1)),
                     bsTooltip("pcr_gamma_sd", "The standard deviation of the underlying gamma distribution. Increase to increase the persistence of PCR sensitivity", placement = "right", trigger = "hover",
                               options = list(container="body")),
                     column(4,numericInput("pcr_gamma_max_sens",
                                           "Peak PCR sensitivity:",
                                           min = 0,
                                           max = 1,
                                           value = 0.95,
                                           step=0.01)),
                     tags$head(tags$style(HTML('.irs-from, .irs-to, .irs-min, .irs-max, .irs-single {
            visibility: hidden !important;
    }'))),
                     column(4,sliderInput("pcr_gamma_height",
                                          "Scale PCR:",
                                          min = 0,
                                          max = 2,
                                          value = 0.95,
                                          ticks=FALSE,
                                          step=0.01)),
                     bsTooltip("pcr_gamma_height", "This slider arbitrarily scales the PCR sensitivity curve up, maintaining the maximum specified sensitivity", placement = "right", trigger = "hover",
                               options = list(container="body"))
                     ),
                 hr(),
                 fluidRow(
                     column(4,numericInput("lft_gamma_sd",
                                          "LFT width:",
                                          min = 1,
                                          max = 25,
                                          value = 5,
                                          step=1)),
                     bsTooltip("lft_gamma_sd", "The standard deviation of the underlying gamma distribution. Increase to increase the persistence of LFT sensitivity", placement = "right", trigger = "hover",
                               options = list(container="body")),
                     column(4, numericInput("lft_gamma_max_sens",
                                           "Peak LFT sensitivity:",
                                           min = 0,
                                           max = 1,
                                           value = 0.95,
                                           step=0.01)),
                     tags$head(tags$style(HTML('.irs-from, .irs-to, .irs-min, .irs-max, .irs-single {
            visibility: hidden !important;
    }'))),
                     column(4,sliderInput("lft_gamma_height",
                                           "Scale LFT:",
                                           min = 0,
                                           max = 2,
                                           value = 0.95,
                                           ticks=FALSE,
                                           step=0.01)),
                     bsTooltip("lft_gamma_height", "This slider arbitrarily scales the LFT sensitivity curve up, maintaining the maximum specified sensitivity", placement = "right", trigger = "hover",
                               options = list(container="body"))
                     ),
                 hr(),
                 fluidRow(
                     column(6,numericInput("pcr_spec","PCR specificity:",
                                           min=0.5,max=1,value=0.9999,step=0.0001)),
                     column(6,numericInput("lft_spec","LFT specificity:",
                                           min=0.5,max=1,value=0.9997,step=0.0001))
                     
                 ),
                 width=3),
             mainPanel(
                 textOutput('test_characteristics1'),
                 textOutput('test_characteristics2'),
                plotOutput("assumptions_plot",width="100%")
             )
    )
    )
)
