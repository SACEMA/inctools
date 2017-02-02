library(shiny)

fluidPage(
  titlePanel("Required sample size for power (baseline survey and cohort)"),
  fluidRow(
    column(3,
           wellPanel(
             h4("Context"),
             numericInput("Inc","Baseline incidence (cases/100PY):", value = 1, step = 0.1, min = 0),
             numericInput("Prev","Baseline prevalence (%):", value = 10, step = 0.5, min = 0, max = 100),
             conditionalPanel("input.tabs == 'tab2_val'",
                              sliderInput("FracIncRed","Incidence reduction (%):",min = 5, max = 95, value = 50, step = 5)
             )
             #sliderInput("FracIncRed","Incidence reduction (%):",min = 5, max = 95, value = 50, step = 5)
           ),
           wellPanel(
             h4("Recency test"),
             numericInput("MDRI","Mean Duration of Recent Infection (days):", value = 180, step = 1, min = 0, max = 730),
             numericInput("RSE_MDRI","Relative Standard Error on MDRI (%):", value = 10, step = 0.5, min = 0, max = 100),
             numericInput("FRR","False-Recent Rate (%):", value = 0.5, step = 0.1, min = 0, max = 100),
             numericInput("RSE_FRR","Relative Standard Error on FRR (%):", value = 25, step = 0.5, min = 0, max = 100),
             sliderInput("BigT", "Time cutoff T (days):", min = 180, max = 1096, value = 730, step = 1)
           ),
           wellPanel(
             h4("Significance"),
             numericInput("Power","Required power:", value = 0.8, step = 0.05, max = 1, min = 0),
             #sliderInput("Power","Required power:", min = 0, max = 1, value = 0.8, step = 0.1),
             numericInput("alpha","Level of significance (α):", value = 0.05, step = 0.01, max = 0.2, min = 0)
             # sliderInput("alpha","Level of significance (α):", min = 0.01, max = 0.1,
             #             value = 0.05, step = 0.01)
           )
    ),
    column(3,
           wellPanel(
             h4("Baseline survey"),
             sliderInput("CR", "Recency test coverage rate (%):", min = 0, max = 100, value = 100, step = 1),
             numericInput("DE_H", "Design effect on HIV prevalence:", min = 0, max = 50, value = 1, step = 0.1),
             numericInput("DE_R", "Design effect on prevalence of recency:", min = 0, max = 50, value = 1, step = 0.1)
           ),
           wellPanel(
             h4("Cohort"),
             sliderInput("CohortCR","Proportion of negatives recruited (%):", min = 0, max = 100, value = 100, step = 1),
             numericInput("FUT","Follow-up time (years):", value = 1, step = 0.5, min = 0, max = 5),
             numericInput("DE_C","Design effect on cohort incidence:", min = 0, max = 50, value = 1, step = 0.1)
           ),
           wellPanel(
             div(img(src='SACEMA_logo.png', width = "100%"), style = "text-align: center;"),
             br(),
             div(img(src='FIND_logo.jpg', width = "100%"), style = "text-align: center;")
           )
           # wellPanel(
           #   conditionalPanel("input.tabs == 'tab1_val'",
           #                        sliderInput('tab1_slider', 'tab1 slider', min=2,max=7,value=2)
           #   )
           #   # conditionalPanel("input.tabs == 'tab2_val'",
           #   #                  sliderInput('tab2_slider', 'tab2 slider', min=2,max=7,value=2)
           #   # )
           # )

    ),
    column(6,
           tabsetPanel(id = "tabs", type = "tabs",
                       tabPanel("Incidence reduction", value='tab1_val', id = 'tab1',
                                wellPanel(
                                  p("Select an incidence reduction range to see the minimum sample sizes in a plot"),
                                  sliderInput("IncRedRange", "Incidence reduction range (%):",
                                              min = 5, max = 95, step = 5, value = c(50, 75))
                                  #actionButton("plot_incred","Plot")
                                ),
                                wellPanel(
                                  plotOutput("incredplot", height = "500px", width = "100%"),
                                  tableOutput("incred_result_table")
                                )
                       ),
                       tabPanel("Power vs sample size", value='tab2_val', id = 'tab2',
                                wellPanel(
                                  span(strong(textOutput("text"))),
                                  br(),
                                  plotOutput("powerplot", height = "500px", width = "100%"),
                                  tableOutput("table_result_power_ss")
                                  #tableOutput("incred_result_table")
                                )
                       ),
                       #new
                       tabPanel("Power vs sample size and incidence % reduction", value='tab3_val', id = 'tab3',
                                wellPanel(
                                  p("Select an incidence reduction range and compare power and sample size required"),
                                  sliderInput("IncRedRange_grid", "Incidence Reduction Range (%):",
                                              min = 5, max = 95, step = 5, value = c(50, 75))
                                  #actionButton("plot_incred_many", "Plot")
                                ),
                                wellPanel(
                                  plotOutput("incred_result_many", height = "500px", width = "100%"),
                                  #renderDataTable(incred_table_power_ss)
                                  tableOutput("incred_table_power_ss")
                                )
                       ),
                       #endnew
                       tabPanel("About", value='tab4_val', id = 'tab4',
                                wellPanel(
                                  p("This tool calculates the minimum required sample size to achieve a specified
                                    probability of correctly inferring an incidence decline (power) in the special
                                    case where a baseline cross-sectional HIV incidence survey is conducted and HIV-
                                    survey respondents are recruited into a cohort to observe an expected decline
                                    in incidence. It also allows to evaluate the impact on power and required sample size
                                    for different values of % incidence reduction."),
                                  p("Contributors:"),
                                  tags$ul(
                                    tags$li("Eduard Grebe"),
                                    tags$li("Stefano Ongarello"),
                                    tags$li("Cari van Schalkwyk"),
                                    tags$li("Alex Welte")
                                  ),
                                  p(em("Built using ", a(strong("inctools"), href = "https://cran.r-project.org/web/packages/inctools/index.html", target = "_blank")))
                                  )
                                )
           )
    )
  )

  )
