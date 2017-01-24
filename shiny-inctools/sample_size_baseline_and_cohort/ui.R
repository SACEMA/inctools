library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  titlePanel("Required sample size for power (baseline survey and cohort)") , 
  sidebarLayout(
    sidebarPanel(
      h2("Significance"),
      sliderInput("Power","Required power (probability of correctly diagnosing incidence difference):",min=0, max=1, value=0.8, step=0.1),
      sliderInput("alpha","Level of significance (α):",min=0.01, max=0.1, value=0.05, step=0.01),
      h2("Context"),
      numericInput("Inc","Baseline incidence (cases/100PY):", value=1, step=0.1, min=0),
      numericInput("Prev","Baseline prevalence (%):", value=10, step=0.5, min=0, max = 100),
      sliderInput("FracIncRed","Incidence reduction (%):",min=0, max=100, value=50, step=5),
      h2("Baseline survey"),
      numericInput("MDRI","Mean Duration of Recent Infection (days):", value=180, step=1, min=0, max = 730),
      numericInput("RSE_MDRI","Relative Standard Error on MDRI (%):", value=10, step=0.5, min=0, max = 100),
      numericInput("FRR","False-Recent Rate (%):", value=0.5, step=0.1, min=0, max = 100),
      numericInput("RSE_FRR","Relative Standard Error on FRR (%):", value=25, step=0.5, min=0, max = 100),
      sliderInput("CR","Recency test coverage rate (%):", min=0, max=100, value=100, step=1),
      numericInput("DE_H","Design effect on HIV prevalence:", min=0, max=50, value=1, step=0.1),
      numericInput("DE_R","Design effect on prevalence of recency:", min=0, max=50, value=1, step=0.1),
      sliderInput("BigT","Time cutoff T (days):", min=180, max=1096, value=730, step=1),
      h2("Cohort"),
      sliderInput("CohortCR","Proportion of negatives recruited (%):", min=0, max=100, value=100, step=1),
      numericInput("FUT","Follow-up time (years):", value=1, step=0.5, min=0, max = 5),
      numericInput("DE_C","Design effect on cohort incidence:", min=0, max=50, value=1, step=0.1),
      actionButton("button","Calculate")
      #img(src = "SACEMA_logo.png", height = 100, width = 200), img(src = "US_logo.png", height = 100, width = 200)
    ),
    mainPanel(
      img(src='SACEMA_logo.png', align = "right", height = "75px"),
      #p(em("Incidence difference from baseline cross-sectional survey and cohort recruited from those screened HIV-")),
      #br(),br(),br(),
      p(strong(em(textOutput("text_desc")))),
      p(strong(textOutput("text_ss"))),
      br(),
      plotOutput("plot1", height = 600, width = 600),
      br(),br(),br(),
      p(em("Created by Alex Welte, Eduard Grebe and Cari Van Schalkwyk."))
    )
  )
))