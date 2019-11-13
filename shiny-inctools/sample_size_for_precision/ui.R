# Created by and Copyright (C) 2017-2018 Lamin Juwara (McGill).
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.


library(shiny)
# Define UI

shinyUI(fluidPage(

  # Application title
  titlePanel("Sample Size for Precision"),
  br(),
  sidebarLayout(
    sidebarPanel(

      wellPanel(fluidPage(
        fluidRow(
          column(4, downloadButton('downloadPlot', 'Save Plot')),
          column(4, downloadButton('downloadData', 'Download Table'))
        )
        # fluidRow(
        #   column(6, downloadButton('downloadPlot', 'Save Plot')))
      )),
      #hr(),
      wellPanel(fluidPage(
        h3("Sample Size"),
        fluidRow(column(12,
                        sliderInput("n_range", 
                                    label = "Sample size range",
                                    min = 0, max = 50000,
                                    value = c(0, 50000), 
                                    step = 100, animate = FALSE)
        ))
      )),
      
      # fluid page for the assay parameters
      wellPanel(fluidPage(
        
        h3("Assay Parameters"),
        
        fluidRow(
          column(6,
                 numericInput("MDRI",
                              label = h5("MDRI estimate (days)"),
                              min = 0,
                              max = 720,
                              step = 1,
                              value = 211),
                 numericInput("FRR",
                              label = h5("FRR estimate (%)"),
                              min = 0,
                              max = 100,
                              step = 0.1,
                              value = 0.9)),
          column(6,
                 numericInput("RSE_MDRI",
                              label = h5("RSE of MDRI estimate (%)"),
                              min = 0,
                              max = 100,
                              value = 5,
                              step = 0.1),
                 numericInput("RSE_FRR",
                              label = h5("RSE of FRR estimate (%)"),
                              min = 0,
                              max = 100,
                              value = 20,
                              step = 0.1)),
          column(8,
                 numericInput("BigT",
                              label = h5("Cut-off time T (days)"),
                              value = 720,
                              step = 1)
          )
          
        )
        
      )),
      wellPanel(# starting a new fluid page
        fluidPage(
          h3("Survey Parameters"),
          fluidRow(
            column(6,
                   numericInput("I", 
                                label = h5("Reference Incidence (% p.a)"),
                                min = 0,
                                max = 100,
                                value = 1.7,
                                step = 0.1)),
            column(6,
                   numericInput("PrevH",
                                label = h5("Reference prevalence (%)"),
                                min = 0,
                                max = 100,
                                value = 10,
                                step = 0.1))),
          fluidRow(
            column(9,
                   sliderInput("CR", # the coverage rate probability
                               label = h5("Coverage rate (%)"),
                               min = 0,
                               max = 100, 
                               value = 50,
                               step = 0.1)
                   )
          ))
      ),
      wellPanel(#Design Effect parameters
        fluidPage(
          h3("Design Effect Parameters"),
          fluidRow(
            column(6,
                   numericInput("DE_H",
                                label = h5("Design effect for HIV infection prevalence"),
                                value = 0.9, 
                                step = 0.1)),
            column(6,
                   numericInput("DE_R",
                                label = h5("Design effect for Recent infection prevalence among positives"),
                                value = 0.9,
                                step = 0.1)
            )
          )
        )),

      wellPanel( fluidPage(
        fluidRow(
          column(12,
                 sliderInput("RSE_req_Inc",
                             label = h4("RSE required for incidence estimate"),
                             value = .25, min = 0, max = 0.8,
                             step = 0.01))
        ),
        fluidRow(
          column(12,
                 numericInput("RSE_infSS",
                              label = h4("RSE of the incidence estimate at infinite sample size"),
                              value = 0.0761, min = 0, max = 1,
                              step = 0.0001))
        )
      ))
      
     
      ),
    mainPanel(
      fluidRow(
        column(12,
               img(src='SACEMA_logo.jpg', align = "right", height = "75px")
               #img(src='mcgill.png', align = "right", height = "40px"),
        )),

      #plotOutput("plot")
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot", width = "100%", height = "600px")),
                  tabPanel("Table", tableOutput("tab"),
                           p(""),
                           p(strong('Definition of Parameters')),
                           br("I: Expected Incidence."),
                           br("PrevH: Prevalence of HIV."),
                           br("CR: Coverage rate."),
                           br("MDRI: mean duration of recent infection in days "),
                           br("RSE_MDRI: Relative standard error of MDRI"),
                           br("FRR: False recent rate "),
                           br("n:  The Sample Size")),
                  tabPanel("User Guide", value='tab3_val', id = 'tab3',
                           wellPanel(includeHTML("sample_size_calculator.html")),
                           p("")),
                  tabPanel("About", value='tab4_val', id = 'tab4',
                           wellPanel( p(""),
                                      p(HTML("Calculates the minimum sample size required for a desired relative 
                                               standard error (RSE) of the incidence estimat given assay characteristics,
                                               reference epidemic state, design effects and recency test coverage.")),
                                      p("Contributors:"),
                                      tags$ul(
                                        tags$li("Lamin Juwara"),
                                        tags$li("Eduard Grebe"),
                                        tags$li("Stefano Ongarello"),
                                        tags$li("Cari van Schalkwyk"),
                                        tags$li("Alex Welte")
                                      ),
                                      p(em("Built using", a(strong("inctools"), href = "https://cran.r-project.org/web/packages/inctools/index.html", target = "_blank")))
                           )
                  )
                  
      )
    )
  )
))
