# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND).
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

#ui.R
library(shiny)
#library(shinyURL)
# Define UI

shinyUI(fluidPage(

  # Application title
  titlePanel("Sample size calculator - Based on SACEMA 3 ABIE"),
  br(),
  sidebarLayout(
    sidebarPanel(
      fluidPage(
        fluidRow(column(6, downloadButton('downloadData', 'Download table'))),
        fluidRow(column(6, downloadButton('downloadPlot', 'Save plot'))),
        fluidRow(
          h3("Analysis parameters"),
          column(9, selectInput("x_variable", label = "Choose analysis type",
                                choices = list("Explore x-axis MDRI" = "MDRI",
                                               "Exolore x-axis FRR" = "FRR",
                                               "Simulation x: MDRI (Case 1 only)" = "simulate",
                                               "Simulation x: FRR (Case 1 only)" = "simulate_FRR"
                                ), selected = "MDRI"))),
        fluidRow(column(9,
                        radioButtons("scenario_case", label = h5("Scenario case type:"),
                                     c("Case 1 (Same MDRI, same FRR estimates in the two surveys)" = 1,
                                       "Case 2 (Same MDRI, but different FRR estimates in the two surveys)" = 2,
                                       "Case 3 (Different MDRI and FRR estimates in the two surveys)" = 3)),
                        selected = 1)
        ),
        fluidRow(column(9,
                        h4("Plot options"),
                        checkboxInput("checkbox_plotparams", label = "Display parameters values on the plot", value = FALSE))),
        conditionalPanel(
          condition = "input.scenario_case == 1 & input.x_variable == 'FRR'",
          sliderInput("FRR_range", label = "FRR range", min = 0, max = 15, value = c(0, 10), step = 0.25, animate = TRUE)
        ),
        conditionalPanel(
          condition = "input.scenario_case != 3 & input.x_variable == 'MDRI'",
          sliderInput("MDRI_range", label = "MDRI range", min = 60, max = 720, value = c(120, 360), step = 10, animate = FALSE)
        )
      ),

      hr(),
      fluidPage(
        h3("Assay parameters"),
        fluidRow(
          column(6,
                 conditionalPanel(
                   condition = "input.scenario_case != 3 & input.x_variable == 'FRR'",
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
                   condition = "input.scenario_case == 1 & input.x_variable == 'MDRI'",
                   sliderInput("FRR",
                               label = h5("FRR estimate (%)"),
                               min = 0,
                               max = 10,
                               step = 0.1,
                               value = 1)
                 ),
                 conditionalPanel(
                   condition = "input.scenario_case != 1",
                   sliderInput("FRR_1",
                               label = h5("FRR estimate for survey 1 (%)"),
                               min = 0,
                               max = 10,
                               step = 0.1,
                               value = 1),
                   sliderInput("FRR_2",
                               label = h5("FRR estimate for survey 2 (%)"),
                               min = 0,
                               max = 10,
                               step = 0.1,
                               value = 1)
                 ),
                 conditionalPanel(
                   condition = "input.x_variable == 'simulate'",
                   sliderInput("MDRI_range_sim1", label = "MDRI range", min = 60, max = 720, value = c(120, 360), step = 10),
                   sliderInput("FRR_range_simul", label = "FRR range", min = 0, max = 15, value = c(0, 3), step = 0.5)
                 ),
                 conditionalPanel(
                   condition = "input.x_variable == 'simulate_FRR'",
                   sliderInput("MDRI_range_simul", label = "MDRI range", min = 60, max = 720, value = c(120, 360), step = 10),
                   sliderInput("FRR_range_sim2", label = "FRR range", min = 0, max = 15, value = c(0, 10), step = 0.25)
                 ),
                 numericInput("BigT", label = h5("Cut-off BigT T (days)"), value = 730, step = 10)
          ),

          column(6,
                 conditionalPanel(
                   condition = "input.scenario_case != 3",
                   numericInput("RSE_MDRI", label = h5("MDRI estimate standard error (%)"), value = 5, step = 0.1)),
                 conditionalPanel(
                   condition = "input.scenario_case == 3",
                   numericInput("RSE_MDRI_1", label = h5("MDRI estimate standard error for survey 1 (%)"), value = 5, step = 0.1),
                   numericInput("RSE_MDRI_2", label = h5("MDRI estimate standard error for survey 2(%)"), value = 5, step = 0.1)),
                 conditionalPanel(
                   condition = "input.scenario_case == 1",
                   numericInput("RSE_FRR", label = h5("FRR estimate standard error (%)"), value = 20, step = 0.1)),
                 conditionalPanel(
                   condition = "input.scenario_case != 1",
                   numericInput("RSE_FRR_1", label = h5("FRR estimate standard error for survey 1 (%)"), value = 20, step = 0.1),
                   numericInput("RSE_FRR_2", label = h5("FRR estimate standard error for survey 2 (%)"), value = 20, step = 0.1))
          )
        )),
      hr(),
      fluidPage(
        h3("Survey parameters"),
        fluidRow(
          column(6,
                 numericInput("I_1", label = h5("Incidence in survey 1 (%)"), value = 5, step = 0.1),
                 numericInput("PrevH_1", label = h5("Prevalence in survey 1 (%)"), value = 15, step = 0.1),
                 numericInput("CR_1", label = h5("Percentage of HIV positives tested for recency in survey 1"), value = 100, step = 1),
                 numericInput("Power", label = h5("Power required"), value = 0.8, step = 0.01)
          ),

          column(6,
                 numericInput("I_2", label = h5("Incidence in survey 2 (%)"), value = 2.5, step = 0.1),
                 numericInput("PrevH_2", label = h5("Prevalence in survey 2 (%)"), value = 12, step = 0.1),
                 numericInput("CR_2", label = h5("Percentage of HIV positives tested for recency in survey 2"), value = 100, step = 1),
                 numericInput("alpha", label = h5("Significance level (alpha)"), value = 0.05, step = 0.1)
          )
        )
      ),
      #Design Effect parameters
      hr(),
      fluidPage(
        h3("Design Effect parameters"),
        fluidRow(
          column(6,
                 numericInput("DE_H_1", label = h5("Infection prevalence in survey 1"), value = 1, step = 0.1),
                 numericInput("DE_R_1", label = h5("Recent infection prevalence among positives in survey 1"), value = 1, step = 0.1)),
          column(6,
                 numericInput("DE_H_2", label = h5("Infection prevalence in survey 2"), value = 1, step = 0.1),
                 numericInput("DE_R_2", label = h5("Recent infection prevalence among positives in survey 2"), value = 1, step = 0.1))
        )
      )
      ),
    mainPanel(
      #img(src='FIND_logo.jpg', align = "right", height = "75px"),
      img(src='SACEMA_logo.jpg', align = "right", height = "75px"),
      br(),
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("plot", width = "100%", height = "600px")),
                  tabPanel("Table", tableOutput("tab"),
                           p(""),
                           p("Table legend"),
                           br("MDRI: Mean Duration of Recent Infection"),
                           br("FRR: False Recent Rate"),
                           br("BigT: Cut-off BigT"),
                           br("RSE_FRR: Covariance of FRR estimate"), br("RSE_MDRI: Covariance of MDRI estimate"),
                           br("DE_H_1: Design effect for prevalence in survey 1"),
                           br("DE_H_2: Design effect for prevalence in survey 2"),
                           br("DE_R_1: Design effect for recent infection prevalence among positives in survey 1"),
                           br("DE_R_2: Design effect for recent infection prevalence among positives in survey 2"),
                           br("I_1 (%): Incidence (%) in survey 1"),
                           br("I_2 (%): Incidence (%) in survey 2"),
                           br("PrevH_1 (%): Prevalence (%) in survey 1"),
                           br("PrevH_2 (%): Prevalence (%) in survey 2"),
                           br("N: Minimal number of subjects required"))
      )
    )
  )
))
