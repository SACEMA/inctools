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

  # Application title  # formerly called the incidence/Prevalence calculator
  titlePanel("Incidence Calculator"),
  br(),
  sidebarLayout(
    sidebarPanel(


      wellPanel(
        fluidRow(column(9,
                        radioButtons("data_type", label = h3("Data Type:"),
                                     c("Sample Counts" = 2,
                                       "Proportions" = 1
                                       )),
                        selected = 2)
        )
      ),
      conditionalPanel(
        condition = "input.data_type == 1",
        wellPanel(
          fluidPage(
            h3("Sample Proportions"),
            fluidRow(
              column(6, 
                     numericInput("PrevH",
                                  label = h5("Prevalence of HIV infection (%)"),
                                  value = 20, step = 0.1, min=0, max = 100),
                     numericInput("PrevR",
                                  label = h5("Prevalence of recent infections among positives (%)"), value = 10, step = 0.1, min=0, max = 100)
              ),
              column(6,
                     numericInput("RSE_PrevH",
                                  label = h5("RSE of Prevalence HIV infection (%)"),
                                  value = 2.8, 
                                  step = 0.1, min=0, max = 100),
                     numericInput("RSE_PrevR",
                                  label = h5("RSE of Prevalence of recent infections among positives (%)"),
                                  value = 9.8,
                                  step = 0.1, min=0, max = 100))
            )
          ))
      ),
      
      conditionalPanel(
        condition = "input.data_type == 2",
        wellPanel(
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
                                  value = 1000, 
                                  step = 1)
              )
            ),
            fluidRow(column(12,
                            numericInput("N_testR", 
                                         label = "HIV positives tested for recency",
                                         value = 1000, 
                                         step = 1)
            )),
            fluidRow(
              column(12,
                     numericInput("N_R", 
                                  label = "Number of recent cases",
                                  value = 50, 
                                  step = 1 )
              )
            )
          ))
        
      ),
      conditionalPanel(
        condition = "input.data_type == 1 | input.data_type == 2 ",
      wellPanel(
        # fluid page for the assay parameters
        fluidPage(
          h3("Assay Parameters"),
          fluidRow(
            column(6,
                   numericInput("MDRI",
                                label = h5("MDRI estimate (days)"),
                                step = 1,
                                value = 210)),
            column(6,
                   numericInput("FRR",
                                label = h5("FRR estimate (%)"),
                                min = 0,
                                max = 100,
                                step = 0.1,
                                value = 0.5))),
          fluidRow(
                   column(6,
                          numericInput("RSE_MDRI",
                                       label = h5("RSE of MDRI estimate (%)"),
                                       min = 0,
                                       max = 100,
                                       value = 5,
                                       step = 0.1)),
                   column(6,
                          numericInput("RSE_FRR",
                                       label = h5("RSE of FRR estimate (%) "),
                                       min = 0,
                                       max = 100,
                                       value = 19,
                                       step = 0.1))),
          fluidRow(
            column(10,
                   numericInput("BigT",
                                label = h5("Cut-off time T (days)"),
                                value = 700,
                                step = 1))
          )
          
        )),
      wellPanel(
        fluidPage(
        #Design Effect parameters
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
          ))
      ))
      )
      
      ),
    mainPanel(
      fluidRow(
        column(12,
               img(src='SACEMA_logo.jpg', align = "right", height = "75px")
        )),
      conditionalPanel(
        condition = "input.data_type == 2",
        tabsetPanel(type = "tabs",
                    tabPanel("Incidence Estimates",
                             br(""),
                             fluidRow(column(12,
                                               downloadButton('downloadData1', 'Download Estimates'))
                             ),
                             br(""),
                             tableOutput("tab4"),
                             fluidRow(
                               plotOutput("plot1")
                             )),
                    tabPanel("User Guide", value='tab3_val', id = 'tab3',
                             wellPanel( p(""),
                                        wellPanel(includeHTML("Incidence_Prevalence_Calculator.html"))
                             )
                    ),
                    tabPanel("About", 
                             #value='tab4_val', id = 'tab4',
                             wellPanel( p(""),
                                        p(HTML("Calculates the point estimate and confidence interval for incidence, prevalence and
                                             annual risk of infection.")),
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
        ),
      conditionalPanel(
        condition = "input.data_type == 1",
        tabsetPanel(type = "tabs",
                    tabPanel("Incidence Estimates",
                             br(""),
                             fluidRow(column(12,
                                               downloadButton('downloadData2', 'Download Estimates'))
                             ),
                             br(""),
                             tableOutput("tab5b")
                           
                             ),

                    tabPanel("User Guide", value='tab3_val', id = 'tab3',
                             wellPanel( p(""),
                                        wellPanel(includeHTML("Incidence_Prevalence_Calculator.html"))
                             )
                    ),
                    tabPanel("About", value='tab4_val', id = 'tab10',
                             wellPanel( p(""),
                                        p(HTML("Calculates the point estimate and confidence interval for incidence, prevalence and
                                             annual risk of infection.")),
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
                    ))
      )

    )
  )
))
