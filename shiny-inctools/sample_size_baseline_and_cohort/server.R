library(shiny)
library(ggplot2)
library(scales)
library(dplyr)
library(plyr)
library(inctools)
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
  
  table_result <- reactive({
    do_table_data_ss_bl_cohort(Inc = input$Inc/100, 
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
                               DE_C = input$DE_C)
  })
  
  #incred_result <- eventReactive(
  incred_result <- reactive(  
    incred_plot(incred_data(Inc = input$Inc/100, 
                            Prev = input$Prev/100, 
                            IncRedRange = c(input$IncRedRange[1]/100, input$IncRedRange[2]/100),
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
                            DE_C = input$DE_C))
  )
  
  incred_result_datatable <- reactive(
    incred_data(Inc = input$Inc/100, 
                Prev = input$Prev/100, 
                IncRedRange = c(input$IncRedRange[1]/100, input$IncRedRange[2]/100),
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
                DE_C = input$DE_C)
  ) 
  
  incred_result_many <- reactive(
    do_inc_red_many(
      FracIncRed = seq(input$IncRedRange_grid[1]/100, input$IncRedRange_grid[2]/100, 0.05),
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
  )
  incred_result_table_ss_power <- reactive(
    generate_table_power_ss_incred(
      FracIncRed = seq(input$IncRedRange_grid[1]/100, input$IncRedRange_grid[2]/100, 0.05),         
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
      DE_C = input$DE_C,
      power_threshold = input$Power
    )
  )
  
  output$powerplot <- renderPlot({
    print(result()$plot)
  })
  output$table_result_power_ss <- renderTable({
    tab <- table_result()
    tab
    tab$Ns <- comma(as.integer(tab$Ns))
    tab$powers <- 100*tab$powers
    colnames(tab) <- c("Sample size", "Power")
    tab[seq(1, dim(tab)[1], by = 5), c("Power", "Sample size")]
  })
  output$text <- renderText({
    if (!is.na(result()$RequiredN)) {
      return(paste("Minimum sample size required:", result()$RequiredN))
    } else {
      return("Desired power level cannot be achieved!")
    }
  })
  output$incred_result_many <- renderPlot({
    print(incred_result_many())
  })
  output$incredplot <- renderPlot({
    print(incred_result())
  })
  output$incred_table_power_ss <- renderTable({
    tab <- incred_result_table_ss_power()
    tab$N <- comma(as.integer(tab$N))
    tab$FracIncRed <- 100*tab$FracIncRed
    tab$Power <- 100*tab$Power
    colnames(tab) <- c("Incidence reduction (%)", "Sample size", "Power (%)")
    tab[c("Incidence reduction (%)", "Power (%)", "Sample size")]
  }
  )
  output$incred_result_table <- renderTable({
    tab <- incred_result_datatable()
    tab$Ns <- comma(as.integer(tab$Ns))
    colnames(tab) <- c("Incidence reduction (%)", "Sample size")
    tab
  })
})