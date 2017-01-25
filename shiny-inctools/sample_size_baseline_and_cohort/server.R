library(shiny)

source("ss_baseline_cohort.R")

shinyServer(function(input, output) {
  result <- reactive({
    ss_baseline_cohort(Inc = input$Inc/100, 
                                         Prev = input$Prev/100, 
                                         FracIncRed = input$FracIncRed/100,
                                         Power = input$Power, 
                                         alpha = input$alpha, 
                                         MDRI = input$MDRI, 
                                         RSE_MDRI = input$RSE_MDRI/100, 
                                         FRR = input$FRR/100, 
                                         RSE_FRR = input$RSE_FRR/100, 
                                         CR = input$CR/100, 
                                         DE_H = input$DE_H,
                                         DE_R = input$DE_R,
                                         BigT = input$BigT, 
                                         CohortCR = input$CohortCR/100, 
                                         FUT = input$FUT,
                                         DE_C = input$DE_C,
                                         produce_plot = TRUE)
    })
  
  incred_result <- eventReactive(
    input$plot_indred,{
      incred_plot(IncRedRange=input$IncRedRange/100,
                  steps = 50,
                  Inc = input$Inc/100, 
                  Prev = input$Prev/100,
                  Power = input$Power, 
                  alpha = input$alpha, 
                  MDRI = input$MDRI, 
                  RSE_MDRI = input$RSE_MDRI/100, 
                  FRR = input$FRR/100, 
                  RSE_FRR = input$RSE_FRR/100, 
                  CR = input$CR/100, 
                  DE_H = input$DE_H,
                  DE_R = input$DE_R,
                  BigT = input$BigT, 
                  CohortCR = input$CohortCR/100, 
                  FUT = input$FUT,
                  DE_C = input$DE_C
      )
    }
  )
  
  output$powerplot <- renderPlot({
    print(result()$plot)
  })
  output$text <- renderText({
    if (!is.na(result()$RequiredN)) {
      return(paste("Minimum sample size required:", result()$RequiredN))
    } else {
      return("Desired power level cannot be achieved!")
    }
  })
  output$incredplot <- renderPlot({
    print(incred_result())
  })
  
  
})