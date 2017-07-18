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
  titlePanel("Incidence/Prevalence Calculator"),
  br(),
  sidebarLayout(
    sidebarPanel(
      fluidPage(
        fluidRow(column(12, downloadButton('downloadData', 'Download Estimates'))),

      hr(),
      fluidPage(
        h3("Sample Counts"),
        fluidRow(column(12,
                        numericInput("N", 
                                    label = "Total Sample Size",
                                    value = 5000, 
                                    step = 1)
        )),
        fluidRow(
        column(12,
               numericInput("N_H", 
                           label = "Number of HIV positive among total",
                           # min = 0,
                           # max = "input.N",
                           value = 1000, 
                           step = 1)
        )
        ),
        fluidRow(column(12,
                        numericInput("N_testR", 
                                    label = "HIV positives tested for recency",
                                    # min = 0,
                                    # max = "input.N_H",
                                    value = 1000, 
                                    step = 1)
        )),
        fluidRow(
        column(12,
               numericInput("N_R", 
                           label = "Number of recent cases",
                           # min = 0,
                           # max = "input.N_testR",
                           value = 50, 
                           step = 1 )
        )
        )
      ),

      # fluid page for the assay parameters
      #hr(),
      fluidPage(
        h3("Assay Parameters"),
        
        fluidRow(
          column(6,
                 numericInput("MDRI",
                             label = h5("MDRI estimate (days)"),
                             step = 1,
                             value = 210),
                 numericInput("FRR",
                             label = h5("FRR estimate (%)"),
                             min = 0,
                             max = 100,
                             step = 0.1,
                             value = 0.5)),
          column(6,
                 numericInput("RSE_MDRI",
                                label = h5("RSE of MDRI estimate (%)"),
                              min = 0,
                              max = 100,
                                value = 5,
                                step = 0.1),
                 numericInput("RSE_FRR",
                              label = h5("RSE of FRR estimate (%) "),
                              min = 0,
                              max = 100,
                              value = 19,
                              step = 0.1)),
          column(10,
                 numericInput("BigT",
                              label = h5("Cut-off time T (days)"),
                              value = 700,
                              step = 1))
                 )
        
      ),

      #Design Effect parameters
      #hr(),
      fluidPage(
        h3("Design Effect Parameters"),
        fluidRow(
          column(6,
                 numericInput("DE_H",
                              label = h5("Design effect for HIV infection prevalence"),
                              value = 1, 
                              step = 0.1)),
          column(6,
                 numericInput("DE_R",
                              label = h5("Design effect for Recent infection prevalence among positives "),
                              value = 1,
                              step = 0.1)
                 )
                 )
     
         
        ))
    )
      ,
    mainPanel(
      img(src='SACEMA_logo.jpg', align = "right", height = "75px"),
      #img(src='mcgill.png', align = "right", height = "40px"),
      br(),
      #plotOutput("plot")
      tabsetPanel(type = "tabs",
                  tabPanel("Estimated Prevalence", tableOutput("tab1"), plotOutput("plot1"),
                           br(),
                           p(""),
                           p(strong('Definition of Parameters')),
                           br("PrevH: Prevalence of HIV."),
                           br("RSE_PrevH: Relative standard error of PrevH"),
                           br("PrevR: Prevalence of recency"),
                           br("RSE_PrevR: Relative standard error of PrevR"),
                           br("x: Sample Count")),
                  tabPanel("Estimated Incidence", tableOutput("tab2"),
                           br(),
                           p(""),
                           p(strong('Definition of Parameters')),
                           br("Incidence: Estimated incidence"),
                           br("CI.low: Confidence interval(lower limit)"),
                           br("CI.up: Confidence interval(upper limit)"),
                           br("RSE: Relative standard error of incidence estimate")),
                  tabPanel("Risk of Infection", tableOutput("tab3"),
                           br(),
                           p(""),
                           p(strong('Definition of Parameters')),
                           br("ARI: Annual Risk of Infection"),
                           br("ARI.CI.low: Lower confidence limit of Annual Risk of Infection"),
                           br("ARI.CI.up: Upper confidence limit of Annual Risk of Infection")),
                  tabPanel("About", value='tab4_val', id = 'tab4',
                           wellPanel( p("Calculates the point estimates and confidence intervals for prevalence,
                                        incidence and annual risk of infection."),
                                      p(HTML("")),
                                      p("Contributors:"),
                                      tags$ul(
                                        tags$li("Eduard Grebe"),
                                        tags$li("Stefano Ongarello"),
                                        tags$li("Cari van Schalkwyk"),
                                        tags$li("Alex Welte"),
                                        tags$li("Lamin Juwara")
                                      ),
                                      p(em("Built using", a(strong("inctools"), href = "https://cran.r-project.org/web/packages/inctools/index.html", target = "_blank")))
                           )
                  )
      )
    )
  )
))
