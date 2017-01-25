library(shiny)

fluidPage(
  titlePanel("Required sample size for power (baseline survey and cohort)"),
  fluidRow(
    column(3,
           wellPanel(
             h4("Significance"),
             sliderInput("Power","Required power:",min=0, max=1, value=0.8, step=0.1),
             sliderInput("alpha","Level of significance (α):",min=0.01, max=0.1, value=0.05, step=0.01)
           ),
           wellPanel(
             h4("Context"),
             numericInput("Inc","Baseline incidence (cases/100PY):", value=1, step=0.1, min=0),
             numericInput("Prev","Baseline prevalence (%):", value=10, step=0.5, min=0, max = 100),
             sliderInput("FracIncRed","Incidence reduction (%):",min=5, max=95, value=50, step=5)
           ),
           wellPanel(
             h4("Recency test"),
             numericInput("MDRI","Mean Duration of Recent Infection (days):", value=180, step=1, min=0, max = 730),
             numericInput("RSE_MDRI","Relative Standard Error on MDRI (%):", value=10, step=0.5, min=0, max = 100),
             numericInput("FRR","False-Recent Rate (%):", value=0.5, step=0.1, min=0, max = 100),
             numericInput("RSE_FRR","Relative Standard Error on FRR (%):", value=25, step=0.5, min=0, max = 100)
           )
    ),
    column(3,
           wellPanel(
             h4("Baseline survey"),
             sliderInput("CR","Recency test coverage rate (%):", min=0, max=100, value=100, step=1),
             numericInput("DE_H","Design effect on HIV prevalence:", min=0, max=50, value=1, step=0.1),
             numericInput("DE_R","Design effect on prevalence of recency:", min=0, max=50, value=1, step=0.1),
             sliderInput("BigT","Time cutoff T (days):", min=180, max=1096, value=730, step=1)
           ),
           wellPanel(
             h4("Cohort"),
             sliderInput("CohortCR","Proportion of negatives recruited (%):", min=0, max=100, value=100, step=1),
             numericInput("FUT","Follow-up time (years):", value=1, step=0.5, min=0, max = 5),
             numericInput("DE_C","Design effect on cohort incidence:", min=0, max=50, value=1, step=0.1)
           ),
           wellPanel(
             div(img(src='SACEMA_logo.png', width = "100%"), style="text-align: center;"),
             br(),
             div(img(src='FIND_logo.jpg', width = "100%"), style="text-align: center;")
           )
           
    ),
    column(6,
           tabsetPanel(id = "tabs", type = "tabs", 
                       tabPanel("Results",
                                wellPanel(
                                  span(strong(textOutput("text"))),
                                  br(),
                                  plotOutput("powerplot", height = "500px", width = "100%")
                                )
                       ),
                       tabPanel("Incidence Reduction Plot",
                                wellPanel(
                                  p("Select an incidence reduction range to see the minimum sample sizes in a plot"),
                                  sliderInput("IncRedRange", "Incidence Reduction Range (%):",
                                              min = 5, max = 95, step = 5, value = c(50,75)),
                                  actionButton("plot_indred","Plot")
                                ),
                                wellPanel(
                                  plotOutput("incredplot", height = "500px", width = "100%")
                                )
                       ),
                       tabPanel("About",
                                wellPanel(
                                  p("This tool calculates the minimum required sample size to achieve a specified 
                                       probability of correctly inferring an incidence decline (power) in the special 
                                       case where a baseline cross-sectional HIV incidence survey is conducted and HIV- 
                                       survey respondents are recruited into a cohort to observe an expected decline 
                                       in incidence."),
                                  p("Contributors:"),
                                  tags$ul(
                                    tags$li("Eduard Grebe"), 
                                    tags$li("Cari van Schalkwyk"), 
                                    tags$li("Stefano Ongarello"),
                                    tags$li("Alex Welte")
                                  ),
                                  p(em("Built using ", a(strong("inctools"), href = "https://cran.r-project.org/web/packages/inctools/index.html", target = "_blank")))
                                  )
                       )
           )
    )
  )
  
)

# # Define UI for application that draws a histogram
# shinyUI(fluidPage(
#   
#    , 
#   sidebarLayout(
#     sidebarPanel(
#       h1("Inputs")
#       
#       
#       
#       
#       #img(src = "SACEMA_logo.png", height = 100, width = 200), img(src = "US_logo.png", height = 100, width = 200)
#     ),
#     mainPanel(
#       
#       #p(em("Incidence difference from baseline cross-sectional survey and cohort recruited from those screened HIV-")),
#       #br(),br(),br(),
#       
#     )
#   )
# ))