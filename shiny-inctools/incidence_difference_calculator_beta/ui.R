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
         wellPanel(
         fluidRow(column(9, downloadButton('downloadData', 'Download Results')))
         ),
        # wellPanel(fluidRow(
        #     column(12,
        #            radioButtons("case", label = h3("Scenario Type:"),
        #                         c(" Same MDRI, same FRR estimates in the two surveys" = 1,
        #                           " Same MDRI, but different FRR estimates in the two surveys" = 2,
        #                           " Different MDRI and FRR estimates in the two surveys" = 3)),
        #            selected = 1)
        #   )),
       
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
        
       # wellPanel(
       #    # fluid page for the assay parameters
       #    fluidPage(
       #      h3("Assay Parameters"),
       #        fluidRow(
       #          column(6,
       #                 conditionalPanel(
       #                   condition = "input.case != 3",
       #                   numericInput("MDRI",
       #                                label = h5("MDRI estimate (days)"),
       #                                min = 0,
       #                                max = 720,
       #                                step = 1,
       #                                value = 240)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case == 3",
       #                   numericInput("MDRI_1",
       #                                label = h5("MDRI estimate for survey 1 (days)"),
       #                                min = 0,
       #                                max = 720,
       #                                step = 1,
       #                                value = 240)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case == 1",
       #                   numericInput("FRR",
       #                                label = h5("FRR estimate (%)"),
       #                                min = 0,
       #                                max = 100,
       #                                step = 0.1,
       #                                value = 1)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case != 1",
       #                   numericInput("FRR_1",
       #                                label = h5("FRR estimate for survey 1 (%)"),
       #                                min = 0,
       #                                max = 100,
       #                                step = 0.1,
       #                                value = 1)
       #                 )
       #          ),
       #          column(6,
       #                 conditionalPanel(
       #                   condition = "input.case != 3",
       #                   numericInput("RSE_MDRI", label = h5("RSE of MDRI estimate (%)"), value = 5, step = 0.1)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case == 3",
       #                   numericInput("RSE_MDRI_1", label = h5("RSE of MDRI estimate for survey 1 (%)"), value = 5, step = 0.1)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case == 1",
       #                   numericInput("RSE_FRR", label = h5("RSE of FRR estimate (%)"), value = 20, step = 0.1)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case != 1",
       #                   numericInput("RSE_FRR_1", label = h5("RSE of FRR estimate for survey 1 (%)"), value = 20, step = 0.1)
       #                 )
       #          )
       #        ),
       #        fluidRow(
       #          column(6,
       #                 conditionalPanel(
       #                   condition = "input.case == 3",
       #                   numericInput("MDRI_2",
       #                                label = h5("MDRI estimate for survey 2 (days)"),
       #                                min = 0,
       #                                max = 720,
       #                                step = 1,
       #                                value = 240)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case != 1",
       #                   numericInput("FRR_2",
       #                                label = h5("FRR estimate for survey 2 (%)"),
       #                                min = 0,
       #                                max = 100,
       #                                step = 0.1,
       #                                value = 1)
       #                 )
       #          ),
       #          column(6,
       #                 conditionalPanel(
       #                   condition = "input.case == 3",
       #                   numericInput("RSE_MDRI_2", label = h5("RSE of MDRI estimate for survey 2(%)"), value = 5, step = 0.1)
       #                 ),
       #                 conditionalPanel(
       #                   condition = "input.case != 1",
       #                   numericInput("RSE_FRR_2", label = h5("RSE of FRR estimate for survey 2 (%)"), value = 20, step = 0.1)
       #                 )
       #          )
       # 
       #        ),
       #     
       #        fluidRow(
       #        column(10,numericInput("BigT", label = h5("Cut-off time T (days)"), value = 730, step = 10)
       #        )))
       # ),
       
       wellPanel(
         fluidPage(
           h3("Survey Parameters"),
           fluidRow(
             
             column(6, 
                    numericInput("PrevH_1", label = h5("Prevalence of HIV infection in survey 1 (%)"), value = 20, step = 0.1, min=0, max = 100),
                    numericInput("PrevR_1", label = h5("Prevalence of recent infections among positives in survey 1 (%)"), value = 10, step = 0.1, min=0, max = 100)
             ),
             column(6,
                    numericInput("RSE_PrevH_1", label = h5("RSE of Prevalence HIV infection in survey 1 (%)"), value = 2.8, step = 0.1, min=0, max = 100),
                    numericInput("RSE_PrevR_1", label = h5("RSE of Prevalence of recent infections among positives 1 (%)"), value = 9.8, step = 0.1, min=0, max = 100))
             # )
           ),
           # wellPanel(
           fluidRow(
             
             column(6, 
                    numericInput("PrevH_2", label = h5("Prevalence of HIV infection in survey 2 (%)"), value = 21, step = 0.1, min=0, max = 100),
                    numericInput("PrevR_2", label = h5("Prevalence of recent infections among positives in survey 2 (%)"), value = 13, step = 0.1, min=0, max = 100)
             ),
             column(6, 
                    numericInput("RSE_PrevH_2", label = h5("RSE of Prevalence HIV infection in survey 2 (%)"), value = 3, step = 0.1, min=0, max = 100),
                    numericInput("RSE_PrevR_2", label = h5("RSE of Prevalence of recent infections among positives 2 (%)"), value = 9.5, step = 0.1, min=0, max = 100))
             # )
           )
          
         )
       )
          )
      #  )
    ),
    mainPanel(
      fluidRow(
      column(12,
             img(src='SACEMA_logo.jpg', align = "right", height = "75px")
             #img(src='mcgill.png', align = "right", height = "40px"),
      )),
      fluidRow(
        tabsetPanel(type = "tabs",
                    tabPanel("Incidence Difference", tableOutput("tab"),
                             p(""),
                             p(style = "color:black", strong('Parameter Definitions')),
                             br(style = "color:grey","compare: surveys compared"),
                             br(style = "color:grey","Diff: point estimate of the estimated difference (p.a)"),
                             br(style = "color:grey","CI.Diff.low: lower limit of confidence interval "),
                             br(style = "color:grey","CI.Diff.up: Upper limit of confidence interval"),
                             br(style = "color:grey","RSE.Diff: RSE of difference estimate"),
                             br(style = "color:grey","RSE.Diff.inf.SS: RSE of the difference estimate at the infinite sample size"),
                             br(style = "color:grey","p.value: P-value for incidence estimate "),
                             br(style = "color:grey","p.value.inf.SS: P-value at infinite sample size"),
                             
                             br(),
                             br()),
                    tabPanel("User Guide", value='tab3_val', id = 'tab3',
                             wellPanel(includeHTML("incidence_difference_calculator.html"))
                             #includeHTML("incidence_difference_calculator.html")
                             ),
                    tabPanel("About", value='tab4_val', id = 'tab4',
                             wellPanel( p(""),
                                        p(HTML("calculates the point estimate of incidence difference, the 
                                               95% and p-value for incidence difference, from two surveys.")),
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
  )
))
