  # Created by Lamin Juwara (McGill) 2017/18 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

#server.R
library(shiny)
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)
library(grid)
library(inctools)
source('prev_inc_calc.R')


shinyServer(function(input, output, session){

  #shinyURL.server(session)
  data_prevalence <- reactive({ # for prevalence calculation
    validate(
      need(input$N>0,"Please enter a valid total population sample size"),
      need(input$N>=input$N_H,"HIV-positive samples should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive samples tested for recency should be less than HIV-positive samples among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive samples tested for recency"),
      need(input$RSE_FRR >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR <= 100, 'Please provide a valid RSE for FRR'),
      need(!(input$RSE_FRR == "" ), 'Please provide a value for RSE_FRR'),
      need(input$RSE_MDRI >= 0, 'Please provide a valid RSE for MDRI'),
      need(input$RSE_MDRI <= 100, 'Please provide a valid RSE for MDRI'),
      need(!(input$RSE_MDRI == "" ), 'Please provide a value for RSE_MDRI'),
      need(input$MDRI >= 0, 'Please provide a valid value for MDRI'),
      need(input$FRR >= 0, 'Please provide a valid value for FRR'),
      need(input$FRR <= 100, 'Please provide a valid value for FRR'),
      need(input$BigT, 'Please provide a value for the cut-off time'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off time (>120)')
    )
    temp <- prevalence_calc(N = input$N, N_H = input$N_H, 
                            N_testR = input$N_testR, N_R = input$N_R, 
                            DE_H = input$DE_H, DE_R = input$DE_R) 
    return(temp)
  })
  
  data_incidence <- reactive({ # for incidence calculiation
    validate(
      need(input$N>0,"Please enter a valid total population sample size"),
      need(input$N>=input$N_H,"HIV-positive samples should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive samples tested for recency should be less than HIV-positive samples among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive samples tested for recency"),
      need(input$N_R>=0,"Please enter a valid number for recent HIV cases"),
      need(input$RSE_FRR >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR <= 100, 'Please provide a valid RSE for FRR'),
      need(!(input$RSE_FRR == "" ), 'Please provide a value for RSE_FRR'),
      need(input$RSE_MDRI >= 0, 'Please provide a valid RSE for MDRI'),
      need(input$RSE_MDRI <= 100, 'Please provide a valid RSE for MDRI'),
      need(!(input$RSE_MDRI == "" ), 'Please provide a value for RSE_MDRI'),
      need(input$MDRI >= 0, 'Please provide a valid value for MDRI'),
      need(input$FRR >= 0, 'Please provide a valid value for FRR'),
      need(input$FRR <= 100, 'Please provide a valid value for FRR'),
      need(input$BigT, 'Please provide a value for the cut-off time'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off time (>120)')
    )
    temp <- incidence_calc(N = input$N, N_H = input$N_H, 
                           N_testR = input$N_testR, N_R = input$N_R, 
                           DE_H = input$DE_H, DE_R = input$DE_R,
                           MDRI = input$MDRI, RSE_MDRI = input$RSE_MDRI/100,
                           FRR = input$FRR/100, RSE_FRR = input$RSE_FRR/100,
                           BigT = input$BigT) 
    return(temp)
  })
  
  data_risk <- reactive({ # for risk of infection calculiation
    validate(
      need(input$N>0,"Please enter a valid total population sample size"),
      need(input$N>=input$N_H,"HIV-positive samples should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive samples tested for recency should be less than HIV-positive samples among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive samples tested for recency"),
      need(input$RSE_FRR >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR <= 100, 'Please provide a valid RSE for FRR'),
      need(!(input$RSE_FRR == "" ), 'Please provide a value for RSE_FRR'),
      need(input$RSE_MDRI >= 0, 'Please provide a valid RSE for MDRI'),
      need(input$RSE_MDRI <= 100, 'Please provide a valid RSE for MDRI'),
      need(!(input$RSE_MDRI == "" ), 'Please provide a value for RSE_MDRI'),
      need(input$MDRI >= 0, 'Please provide a valid value for MDRI'),
      need(input$FRR >= 0, 'Please provide a valid value for FRR'),
      need(input$FRR <= 100, 'Please provide a valid value for FRR'),
      need(input$BigT, 'Please provide a value for the cut-off time'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off time (>120)')
    )
    temp <- risk_of_infection_calc(N = input$N, N_H = input$N_H, 
                           N_testR = input$N_testR, N_R = input$N_R, 
                           DE_H = input$DE_H, DE_R = input$DE_R,
                           MDRI = input$MDRI, RSE_MDRI = input$RSE_MDRI/100,
                           FRR = input$FRR/100, RSE_FRR = input$RSE_FRR/100,
                           BigT = input$BigT) 
    return(temp)
  })
  
  data_pie<-reactive({
    Sample.Count<-c((input$N-input$N_H),input$N_R,
         (input$N_H-input$N_R),(input$N_H-input$N_testR))
    piepercent<-paste(round(100*Sample.Count/sum(Sample.Count),1),sep = "","%") 
    label<-c("HIV-negative","HIV-positive and 'recent' ",
             "HIV-positve and 'not recent' ", 
             "HIV-positive and not tested for recency")
    color<-c("blue","orange","yellow","violet")
    pie(x = Sample.Count,labels = piepercent, col = color,
                         main = "Sample Counts")
    options(error = NULL)
    legend("bottomleft",legend = label,cex = 1.0,fill = color)
  })
  
# Output value  
   

  # Produce an output table value.
  output$tab1 <- renderTable({
    data_prevalence()
    
  })
  output$plot1<-renderPlot({
    data_pie()
  })
  
  
  output$tab2 <- renderTable({
    data_incidence()
    
  })
  output$tab3 <- renderTable({
    data_risk()
    
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("resultTable-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(c(data_prevalence(),data_incidence(),data_risk()) , file)
    }
  )

 

})
