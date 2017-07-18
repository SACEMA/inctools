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
  titlePanel("Incidence Difference Calculator"),
  br(),
  sidebarLayout(
    sidebarPanel(
      fluidPage(
        fluidRow(column(6, downloadButton('downloadData', 'Download table'))),
        fluidRow(column(9,
                        radioButtons("case", label = h3("Scenario Type:"),
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
                   condition = "input.case != 3",
                   sliderInput("MDRI",
                               label = h5("MDRI estimate (days)"),
                               min = 0,
                               max = 720,
                               step = 1,
                               value = 240)
                 ),
                 conditionalPanel(
                   condition = "input.case == 3",
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
                               value = 200)
                 ),
                 conditionalPanel(
                   condition = "input.case == 1",
                   sliderInput("FRR",
                               label = h5("FRR estimate (%)"),
                               min = 0,
                               max = 100,
                               step = 0.1,
                               value = 1)
                 ),
                 conditionalPanel(
                   condition = "input.case != 1",
                   sliderInput("FRR_1",
                               label = h5("FRR estimate for survey 1 (%)"),
                               min = 0,
                               max = 100,
                               step = 0.1,
                               value = 1),
                   sliderInput("FRR_2",
                               label = h5("FRR estimate for survey 2 (%)"),
                               min = 0,
                               max = 100,
                               step = 0.1,
                               value = 2)
                 ),
                 
                 numericInput("BigT", label = h5("Cut-off time T (days)"), value = 730, step = 10)
          ),
          
          column(6,
                 conditionalPanel(
                   condition = "input.case != 3",
                   numericInput("RSE_MDRI", label = h5("RSE of MDRI estimate (%)"), value = 5, step = 0.1)),
                 conditionalPanel(
                   condition = "input.case == 3",
                   numericInput("RSE_MDRI_1", label = h5("RSE of MDRI estimate for survey 1 (%)"), value = 5, step = 0.1),
                   numericInput("RSE_MDRI_2", label = h5("RSE of MDRI estimate for survey 2 (%)"), value = 5, step = 0.1)),
                 conditionalPanel(
                   condition = "input.case == 1",
                   numericInput("RSE_FRR", label = h5("RSE of FRR estimate (%)"), value = 20, step = 0.1)),
                 conditionalPanel(
                   condition = "input.case != 1",
                   numericInput("RSE_FRR_1", label = h5("RSE of FRR estimate for survey 1 (%)"), value = 20, step = 0.1),
                   numericInput("RSE_FRR_2", label = h5("RSE of FRR estimate for survey 2 (%)"), value = 30, step = 0.1))
                 
          )
        )),
      #hr(),
      fluidPage(
        h3("Survey Parameters"),
        fluidRow(
          column(6,
                 numericInput("I1", label = h5("Incidence in survey 1 (%)"), value = 5, step = 0.1),
                 numericInput("PrevH1", label = h5("Prevalence in survey 1 (%)"), value = 20, step = 0.1),
                 numericInput("CR_1", label = h5("Percentage of HIV positives tested for recency in survey 1 (%)"), value = 100, step = 0.1),
                 sliderInput("Power", label = "Power", min = 0, max = 1, value = 0.8, step = 0.01)
          ),
          
          column(6,
                 numericInput("I2", label = h5("Incidence in survey 2 (%)"), value = 2.5, step = 0.1),
                 numericInput("PrevH2", label = h5("Prevalence in survey 2 (%)"), value = 15, step = 0.1),
                 numericInput("CR_2", label = h5("Percentage of HIV positives tested for recency in survey 2 (%)"), value = 100, step = 0.1),
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
                 numericInput("DE_H_1", label = h5("Infection prevalence in survey 1"), value = 1, step = 0.1),
                 numericInput("DE_R_1", label = h5("Recent infection prevalence among positives in survey 1"), value = 1, step = 0.1)),
          column(6,
                 numericInput("DE_H_2", label = h5("Infection prevalence in survey 2"), value = 1, step = 0.1),
                 numericInput("DE_R_2", label = h5("Recent infection prevalence among positives in survey 2"), value = 1, step = 0.1))
        )
      )

    ),
    mainPanel(
      img(src='SACEMA_logo.jpg', align = "right", height = "75px"),
      #img(src='mcgill.png', align = "right", height = "40px"),
      br(),
      tabsetPanel(type = "tabs",
                  tabPanel("Incidence Difference", tableOutput("tab"),
                           br(),
                           p(""),
                           p(style = "color:black", strong('Parameter Definitions')),
                           #p(strong('Parameter Definitions')),
                           br(style = "color:grey","deltaI_Est: point estimate of the estimated difference (p.a)"),
                           br(style = "color:grey","RSE_deltaI: RSE of difference estimate"),
                           br(style = "color:grey","RSE_deltaI.infSS: RSE of the difference estimate at the infinite sample size"),
                           br(style = "color:grey","Power: P-value for incidence estimate "),
                           br(style = "color:grey","Power.infSS: P-value at infinite sample size"),
                           br(style = "color:grey","CI.low: lower limit of confidence interval "),
                           br(style = "color:grey","CI.up: Upper limit of confidence interval")),
                 
                  tabPanel("About", value='tab4_val', id = 'tab4',
                           wellPanel(
                             p("This tool calculates the point estimate of incidence difference, and the 
                               95% CI and p-value for incidence differenc from two surveys."),
                             # p("Contributors:"),
                             # tags$ul(
                             #   tags$li("Eduard Grebe"),
                             #   tags$li("Stefano Ongarello"),
                             #   tags$li("Cari van Schalkwyk"),
                             #   tags$li("Alex Welte"),
                             #   tags$li("Lamin Juwara")
                             # ),
                             p(em("Built using ", a(strong("inctools"), href = "https://cran.r-project.org/web/packages/inctools/index.html", target = "_blank")))
                           )
                  )
      )
    )
  )
))
