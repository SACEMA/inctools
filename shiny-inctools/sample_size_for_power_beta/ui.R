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
  titlePanel("Sample Size for Power"),
  br(),
  sidebarLayout(
    sidebarPanel(
      
      wellPanel(
        fluidPage(
        fluidRow(
          column(5, downloadButton('downloadData', 'Download Table')),
          column(5, downloadButton('downloadPlot', 'Save Plot')))
      )),

        wellPanel( # fluid page for the assay parameters
          fluidPage(
            h3("Assay Parameters"),
            ###
            fluidRow(
              column(6,
                       numericInput("MDRI_1",
                                    label = h5("MDRI estimate for survey 1 (days)"),
                                    min = 0,
                                    max = 720,
                                    step = 1,
                                    value = 240),
                       numericInput("MDRI_2",
                                    label = h5("MDRI estimate for survey 2 (days)"),
                                    min = 0,
                                    max = 720,
                                    step = 1,
                                    value = 240),
                       numericInput("FRR_1",
                                    label = h5("FRR estimate for survey 1 (%)"),
                                    min = 0,
                                    max = 100,
                                    step = 0.1,
                                    value = 1),
                       numericInput("FRR_2",
                                    label = h5("FRR estimate for survey 2 (%)"),
                                    min = 0,
                                    max = 100,
                                    step = 0.1,
                                    value = 1),
                     selectInput("COV_MDRI",
                                 label = h5("Covariance of MDRI values"),
                                 choices = c(0,1),
                                 selected = 0
                                  )
                     ),
              
              column(6,
                       numericInput("RSE_MDRI_1", label = h5("RSE of MDRI estimate for survey 1 (%)"), value = 5, step = 0.1),
                       numericInput("RSE_MDRI_2", label = h5("RSE of MDRI estimate for survey 2(%)"), value = 5, step = 0.1),
                       numericInput("RSE_FRR_1", label = h5("RSE of FRR estimate for survey 1 (%)"), value = 20, step = 0.1),
                       numericInput("RSE_FRR_2", label = h5("RSE of FRR estimate for survey 2 (%)"), value = 20, step = 0.1),
                       selectInput("COV_FRR",
                                  label = h5("Covariance of FRR values"),
                                  choices = c(0,1),
                                  selected = 0
                       )
                     
                     
                     )
                     
              ),
            fluidRow(column(9,
                            numericInput("BigT", label = h5("Cut-off BigT T (days)"), value = 730, step = 10)
            ))
            )
          ),

      wellPanel(fluidPage(
        h3("Survey Parameters"),
        fluidRow(
          column(6,
                 numericInput("I1", label = h5("Incidence in survey 1 (%)"), value = 5, step = 0.1),
                 numericInput("PrevH1", label = h5("Prevalence in survey 1 (%)"), value = 15, step = 0.1),
                 numericInput("CR_1", label = h5("Percentage of HIV positives tested for recency in survey 1 (%)"), value = 100, step = 1)
                 
          ),
          
          column(6,
                 numericInput("I2", label = h5("Incidence in survey 2 (%)"), value = 2.5, step = 0.1),
                 numericInput("PrevH2", label = h5("Prevalence in survey 2 (%)"), value = 12, step = 0.1),
                 numericInput("CR_2", label = h5("Percentage of HIV positives tested for recency in survey 2 (%)"), value = 100, step = 1)
          )
        )
      )
      ),

      wellPanel(
        #Design Effect parameters
        fluidPage(
          h3("Design Effect Parameters"),
          fluidRow(
            column(6,
                   numericInput("DE_H1", label = h5("Infection prevalence in survey 1"), value = 1, step = 0.1),
                   numericInput("DE_R1", label = h5("Recent infection prevalence among positives in survey 1"), value = 1, step = 0.1)),
            column(6,
                   numericInput("DE_H2", label = h5("Infection prevalence in survey 2"), value = 1, step = 0.1),
                   numericInput("DE_R2", label = h5("Recent infection prevalence among positives in survey 2"), value = 1, step = 0.1))
          )
        )
      ),

               
      wellPanel(
        fluidPage(
        fluidRow(numericInput("alpha",
                              label = "Significance level (alpha)",
                              value = 0.05,
                              step = 0.1)),
        fluidRow(
          numericInput("statPower",
                       label = "Statistical Power (%)",
                       min = 0,
                       max = 100,
                       step = 0.01,
                       value = 80)
        ),
        fluidRow(sliderInput("Power_range",
                             label = "Plot Power Range",
                             min = 0.01, max = 0.99,
                             value = c(0.2, 0.9),
                             step = 0.01))
      )
      )
#)
      ),
    mainPanel(
      img(src='SACEMA_logo.jpg', align = "right", height = "75px"),
      br(),
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot", width = "100%", height = "600px")),
                  tabPanel("Table", tableOutput("tab"),
                           p(""),
                           p(strong('Definition of Parameters')),
                           br("MDRI: Mean Duration of Recent Infection"),
                           br("FRR: False Recent Rate"),
                           br("BigT: Cut-off Time"),
                           br("RSE_FRR: Covariance of FRR estimate"), br("RSE_MDRI: Covariance of MDRI estimate"),
                           br("DE_H1: Design effect for prevalence in survey 1"),
                           br("DE_H2: Design effect for prevalence in survey 2"),
                           br("DE_R1: Design effect for recent infection prevalence among positives in survey 1"),
                           br("DE_R2: Design effect for recent infection prevalence among positives in survey 2"),
                           br("I1 (%): Incidence (%) in survey 1"),
                           br("I2 (%): Incidence (%) in survey 2"),
                           br("PrevH1 (%): Prevalence (%) in survey 1"),
                           br("PrevH2 (%): Prevalence (%) in survey 2"),
                           br("ss: Minimal number of subjects required")),
                  tabPanel("User Guide", value='tab3_val', id = 'tab3',
                           wellPanel( p(""),
                                      wellPanel(includeHTML("sample_size_calculator.html"))
                                      )
                  ),
                  tabPanel("About", value='tab4_val', id = 'tab4',
                           wellPanel( p(""),
                                      p(HTML("Calculates the minimum sample size (common to the two surveys) required for a desired probability
                                      of detecting difference in incidence (with the correct sign).")),
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
