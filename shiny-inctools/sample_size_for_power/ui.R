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
      fluidPage(
        fluidRow(column(6, downloadButton('downloadData', 'Download table'))),
        fluidRow(
                 column(6, downloadButton('downloadPlot', 'Save plot'))),

        fluidRow(column(9,
                        radioButtons("scenario_case", label = h3("Scenario Type:"),
                                     c("Same MDRI, same FRR estimates in the two surveys" = 1,
                                       "Same MDRI, but different FRR estimates in the two surveys" = 2,
                                       "Different MDRI and FRR estimates in the two surveys" = 3)),
                        selected = 1)
        )
      ),

      # fluid page for the assay parameters
      #hr(),
      fluidPage(
        h3("Assay Parameters"),
        fluidRow(
          column(6,
                 conditionalPanel(
                   condition = "input.scenario_case != 3",
                   sliderInput("MDRI",
                               label = h5("MDRI estimate (days)"),
                               min = 0,
                               max = 720,
                               step = 1,
                               value = 240)
                 ),
                 conditionalPanel(
                   condition = "input.scenario_case == 3",
                   sliderInput("MDRI_1",
                               label = h5("MDRI estimate for survey 1 (days)"),
                               min = 0,
                               max = 720,
                               step = 1,
                               value = 240),
                   sliderInput("MDRI_2",
                               label = h5("MDRI estimate for survey 2 (days)"),
                               min = 0,
                               max = 720,
                               step = 1,
                               value = 240)
                 ),
                 conditionalPanel(
                   condition = "input.scenario_case == 1",
                   sliderInput("frrhat",
                               label = h5("FRR estimate (%)"),
                               min = 0,
                               max = 10,
                               step = 0.1,
                               value = 1)
                 ),
                 conditionalPanel(
                   condition = "input.scenario_case != 1",
                   sliderInput("frrhat_1",
                               label = h5("FRR estimate for survey 1 (%)"),
                               min = 0,
                               max = 10,
                               step = 0.1,
                               value = 1),
                   sliderInput("frrhat_2",
                               label = h5("FRR estimate for survey 2 (%)"),
                               min = 0,
                               max = 10,
                               step = 0.1,
                               value = 1)
                 ),
        
                 numericInput("TIME", label = h5("Cut-off time T (days)"), value = 730, step = 10)
          ),

          column(6,
                 conditionalPanel(
                   condition = "input.scenario_case != 3",
                   numericInput("mdrihatcov", label = h5("RSE of MDRI estimate (%)"), value = 5, step = 0.1)),
                 conditionalPanel(
                   condition = "input.scenario_case == 3",
                   numericInput("mdrihatcov_1", label = h5("RSE of MDRI estimate for survey 1 (%)"), value = 5, step = 0.1),
                   numericInput("mdrihatcov_2", label = h5("RSE of MDRI estimate for survey 2(%)"), value = 5, step = 0.1)),
                 conditionalPanel(
                   condition = "input.scenario_case == 1",
                   numericInput("frrhatcov", label = h5("RSE of FRR estimate (%)"), value = 20, step = 0.1)),
                 conditionalPanel(
                   condition = "input.scenario_case != 1",
                   numericInput("frrhatcov_1", label = h5("RSE of FRR estimate for survey 1 (%)"), value = 20, step = 0.1),
                   numericInput("frrhatcov_2", label = h5("RSE of FRR estimate for survey 2 (%)"), value = 20, step = 0.1))
                
          )
        )),
      #hr(),
      fluidPage(
        h3("Survey Parameters"),
        fluidRow(
          column(6,
                 numericInput("inc_1", label = h5("Incidence in survey 1 (%)"), value = 5, step = 0.1),
                 numericInput("p_pos_1", label = h5("Prevalence in survey 1 (%)"), value = 15, step = 0.1),
                 numericInput("rec_test_coverage_1", label = h5("Percentage of HIV positives tested for recency in survey 1 (%)"), value = 100, step = 1),
                 #numericInput("power", label = h5("Power required"), value = 0.8, step = 0.01)
                 sliderInput("power_range", label = "Power Range", min = 0.01, max = 0.99, value = c(0.2, 0.9), step = 0.01)
          ),

          column(6,
                 numericInput("inc_2", label = h5("Incidence in survey 2 (%)"), value = 2.5, step = 0.1),
                 numericInput("p_pos_2", label = h5("Prevalence in survey 2 (%)"), value = 12, step = 0.1),
                 numericInput("rec_test_coverage_2", label = h5("Percentage of HIV positives tested for recency in survey 2 (%)"), value = 100, step = 1),
                 numericInput("alpha", label = h5("Significance level (alpha)"), value = 0.05, step = 0.1)
          )
        )
      ),
      #Design Effect parameters
      #hr(),
      fluidPage(
        h3("Design Effect Parameters"),
        fluidRow(
          column(6,
                 numericInput("DE_prev_1", label = h5("Infection prevalence in survey 1"), value = 1, step = 0.1),
                 numericInput("DE_RgivenTested_1", label = h5("Recent infection prevalence among positives in survey 1"), value = 1, step = 0.1)),
          column(6,
                 numericInput("DE_prev_2", label = h5("Infection prevalence in survey 2"), value = 1, step = 0.1),
                 numericInput("DE_RgivenTested_2", label = h5("Recent infection prevalence among positives in survey 2"), value = 1, step = 0.1))
        )
      ),
      fluidPage(
        #h3(""),
        fluidRow(
          sliderInput("statPower",
                      label = h3("Statistical Power"),
                      min = 0,
                      max = 1,
                      step = 0.01,
                      value = 0.80)
        )
         
      )
      ),
    mainPanel(
      img(src='SACEMA_logo.jpg', align = "right", height = "75px"),
      #img(src='mcgill.png', align = "right", height = "40px"),
      br(),
      #plotOutput("plot")
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot", width = "100%", height = "600px")),
                  tabPanel("Table", tableOutput("tab"),
                           p(""),
                           p(strong('Definition of Parameters')),
                           br("MDRI: Mean Duration of Recent Infection"),
                           br("FRR: False Recent Rate"),
                           br("BigT: Cut-off time"),
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
                  tabPanel("About", value='tab4_val', id = 'tab4',
                           wellPanel( p("Calculates the minimum sample size (common in two surveys) requited to achieve a desired 
                                        probability of detecting a difference in incidence (with the correct sign)."),
                                      p(HTML("")),
                                      # p("Contributors:"),
                                      # tags$ul(
                                      #   tags$li("Eduard Grebe"),
                                      #   tags$li("Stefano Ongarello"),
                                      #   tags$li("Cari van Schalkwyk"),
                                      #   tags$li("Alex Welte"),
                                      #   tags$li("Lamin Juwara")
                                      # ),
                             p(em("Built using", a(strong("inctools"), href = "https://cran.r-project.org/web/packages/inctools/index.html", target = "_blank")))
                           )
                  )
      )
    )
  )
))
