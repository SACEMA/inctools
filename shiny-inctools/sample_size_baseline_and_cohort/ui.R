library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  titlePanel("Required sample size for power") , 
  headerPanel("Power to detect incidence difference from baseline cross-sectional survey and cohort"),
  
  sidebarLayout(
    sidebarPanel(
      h1("Input"),
      numericInput("Inc","Baseline incidence (cases/100PY):", value=1, step=0.1, min=0),
      numericInput("Prev","Baseline prevalence (%):", value=10, step=0.5, min=0, max = 100),
      sliderInput("FracIncRed","Incidence reduction (%):",min=0, max=100, value=50, step=5),
      sliderInput("Power","Required power (probability of correctly diagnosing incidence difference):",min=0, max=1, value=0.8, step=0.1),
      sliderInput("alpha","Level of significance (α):",min=0.01, max=0.1, value=0.05, step=0.01),
      numericInput("MDRI","Mean Duration of Recent Infection (days):", value=180, step=1, min=0, max = 730),
      numericInput("RSE_MDRI","Relative Standard Error on MDRI (%):", value=10, step=0.5, min=0, max = 100),
      numericInput("FRR","False-Recent Rate (%):", value=0.5, step=0.1, min=0, max = 100),
      numericInput("RSE_FRR","Relative Standard Error on FRR (%):", value=25, step=0.5, min=0, max = 100),
      sliderInput("CR","Recency test coverage rate (%):", min=0, max=100, value=100, step=1),
      numericInput("DE_H","Design effect on HIV prevalence:", min=0, max=50, value=1, step=0.1),
      numericInput("DE_R","Design effect on prevalence of recency:", min=0, max=50, value=1, step=0.1),
      sliderInput("BigT","Time cutoff T (days):", min=180, max=1096, value=730, step=1),
      sliderInput("CohortCR","Proportion of negatives recruited (%):", min=0, max=100, value=100, step=1),
      numericInput("FUT","Follow-up time (years):", value=1, step=0.5, min=0, max = 5),
      numericInput("DE_C","Design effect on cohort incidence:", min=0, max=50, value=1, step=0.1),
      br(),
      actionButton("button","Calculate")
      #img(src = "SACEMA_logo.png", height = 100, width = 200), img(src = "US_logo.png", height = 100, width = 200)
    ),
    mainPanel(
      h1("Output"),
      #verbatimTextOutput("text1"),
      p(strong(textOutput("text1"))),
      br(),
      plotOutput("plot1", height = 600, width = 600),
      tags$head(tags$style("#text1{color: blue;
                                 font-size: 20px;
                                
                                 }"))
      
    )
  )
))