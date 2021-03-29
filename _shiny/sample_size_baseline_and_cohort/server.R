# Incidence Estimation Tools (Shiny webapps).
# Copyright (C) 2017-2019, Eduard Grebe, Stefano Ongarello, individual 
# contributors and Stellenbosch University.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
  
  ss_for_power_and_precision <- reactive({
    ss_baseline_cohort_precision(Inc = input$Inc/100, 
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
                                 min_precision = "delta", 
                                 RSE_required = input$RSEdelta_required_perc/100
                                 )
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
  
  output$ss_power_text <- renderText({
    result <- ss_for_power_and_precision()
    if (is.na(result$RequiredN_power)) {
      ss_for_power_string <- "Required power level cannot be achieved."
    } else {
      ss_for_power_string <- paste0("Minimum sample size to detect incidence reduction: ", prettyNum(result$RequiredN_power, big.mark = ","))
    }
    return(ss_for_power_string)
  })
  
  output$ss_precision_text <- renderText({
    result <- ss_for_power_and_precision()
    if (is.na(result$RequiredN_precision)) {
      ss_for_precision_string <- "Required RSE on incidence difference cannot be achieved."
    } else {
      ss_for_precision_string <- paste0("Minimum sample size to achieve desired precision: ", prettyNum(result$RequiredN_precision, big.mark = ","))
    }
    return(ss_for_precision_string)
  })
  
  output$deltaI_ci_header <- renderText({
    result <- ss_for_power_and_precision()
    if (is.na(result$RequiredN_precision)) {
      header_string <- ""
    } else {
      header_string <- paste0("Precision with sample size of ", prettyNum(result$RequiredN, big.mark = ","),":")
    }
    return(header_string)
  })
  
  
  
  output$deltaI_ci_text <- renderText({
    result <- ss_for_power_and_precision()
    if (is.na(result$RequiredN_precision)) {
      ci_string <- ""
    } else {
      #browser()
      deltaI <- input$Inc*(1 - input$FracIncRed/100) - input$Inc
      SE <- result$achieved_RSE * abs(deltaI)
      LB <- round(stats::qnorm(input$alpha/2, mean = deltaI, sd = SE), 2)
      UB <- round(stats::qnorm(1-input$alpha/2, mean = deltaI, sd = SE), 2)
      ci_string <- paste0("Incidence difference: ", round(deltaI, 2), " (", round((1-input$alpha)*100), "% CI: ", LB, ",", UB, ") cases/100PY")
    }
    return(ci_string)
  })
  
})