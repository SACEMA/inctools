library(shiny)

source("ss_baseline_cohort.R")

shinyServer(function(input, output) {
  do_it <- reactive({ss_baseline_cohort(Inc = input$Inc/100, 
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
                                                           DE_C = input$DE_C)})
  
  output$text_ss <- renderText({
    ss_do <- do_it()
    
    if (!is.na(ss_do$RequiredN)) {
      return(ss_do$RequiredN)
    } else {
      return("Desired power level cannot be achieved!")
    }
  })
  
  output$plot1 <- renderPlot({
    plot_do <- do_it()
    plot_do$plot
  })
  
})