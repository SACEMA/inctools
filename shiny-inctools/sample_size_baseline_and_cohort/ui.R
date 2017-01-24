library(shiny)

fluidPage(
  titlePanel("Required sample size for power (baseline survey and cohort)"),
  fluidRow(
    column(12,
           wellPanel(
             p(em("This tool calculates the minimum required sample size to achieve a specified 
                  probability of correctly inferring an incidence decline (power) in the special 
                  case where a baseline cross-sectional HIV incidence survey is conducted and HIV- 
                  survey respondents are recruited into a cohort to observe an expected decline 
                  in incidence."))
           )
           )
  ),
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
             sliderInput("FracIncRed","Incidence reduction (%):",min=0, max=100, value=50, step=5)
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
           )
    ),
    column(6,
           wellPanel(
             h4("Minimum sample size required:"),
             span(strong(textOutput("text_ss"))),
             h4("Power vs sample size:"),
             plotOutput("plot1", height = 400, width = 400)
           ),
           wellPanel(
             div(img(src='SACEMA_logo.png', height = "100px"), style="text-align: center;"),
             br(),
             div(em("Created by Eduard Grebe, Cari Van Schalkwyk and Alex Welte at the ", a("South African Centre for Epidemiological Modelling and Analysis (SACEMA)", href = "http://www.sacema.org", target="_blank"), "."), style="text-align: center;")
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