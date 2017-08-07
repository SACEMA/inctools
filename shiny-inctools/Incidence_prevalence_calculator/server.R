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
      need(input$N>=input$N_H,"HIV-positive subjects should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive subjects tested for recency should be less than HIV-positive subjects among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive subjects tested for recency"),
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
      need(input$N>=input$N_H,"HIV-positive subjects should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive subjects tested for recency should be less than HIV-positive subjects among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive subjects tested for recency"),
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
      need(input$N>=input$N_H,"HIV-positive subjects should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive subjects tested for recency should be less than HIV-positive subjects among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive subjects tested for recency"),
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
  
  
  data_incidence_count <- reactive({ # for the general wrapper from above
    validate(
      need(input$N>0,"Please enter a valid total population sample size"),
      need(input$N>=input$N_H,"HIV-positive subjects should be less than total sample size"),
      need(input$N_H>=input$N_testR,"HIV-positive subjects tested for recency should be less than HIV-positive subjects among total sample size"),
      need(input$N_testR>=input$N_R,"The number of recent HIV cases should be less than HIV-positive subjects tested for recency"),
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
    temp <- prev_inc_calc_counts(N = input$N, N_H = input$N_H, 
                                   N_testR = input$N_testR, N_R = input$N_R, 
                                   DE_H = input$DE_H, DE_R = input$DE_R,
                                   MDRI = input$MDRI, RSE_MDRI = input$RSE_MDRI/100,
                                   FRR = input$FRR/100, RSE_FRR = input$RSE_FRR/100,
                                   BigT = input$BigT) 
    return(temp)
  })
  
  data_pie<-reactive({
    validate(
      need(input$N>0,""),
      need(input$N>=input$N_H,""),
      need(input$N_H>=input$N_testR,""),
      need(input$N_testR>=input$N_R,""))
    Sample.Count<-c((input$N-input$N_H),input$N_R,
         (input$N_testR-input$N_R),(input$N_H-input$N_testR))
    piepercent<-paste(round(100*Sample.Count/sum(Sample.Count),1),sep = "","%") 
    label<-c("HIV-negative",
             "HIV-positive and 'recent' ",
             "HIV-positve and 'not recent' ", 
             "HIV-positive and not tested for recency")
    color<-c("blue","orange","yellow","violet")
    pie(x = Sample.Count,labels = piepercent, col = color,
                         main = "Sample Counts")
    #options(error = NULL)
    legend("bottomleft",legend = label,cex = 1.0,fill = color)
  })
  # create a function for incidence from the incprops calculator
  data_prev_inc_calc_incprop<- reactive({
    validate(
      need(input$PrevH >= 0, 'Please provide a valid value for HIV prevalence'),
      need(input$PrevH <= 100, 'Please provide a valid HIV prevalence'),
      need(!(input$PrevH == "" ), 'Please provide a value HIV prevalence'),
      need(input$RSE_PrevH >= 0, 'Please provide a valid RSE for HIV prevalence'),
      need(input$RSE_PrevH <= 100, 'Please provide a valid RSE for Hiv prevalence'),
      need(!(input$RSE_PrevH == "" ), 'Please provide a value for RSE for HIV prevalence'),
      need(input$PrevR >= 0, 'Please provide a valid value for recency among HIV prevalence'),
      need(input$PrevR <= 100, 'Please provide a valid recency among HIV prevalence'),
      need(!(input$PrevR == "" ), 'Please provide a value recency among HIV prevalence'),
      need(input$RSE_PrevR >= 0, 'Please provide a valid RSE for recency of HIV positives'),
      need(input$RSE_PrevR <= 100, 'Please provide a valid RSE for recency of Hiv prevalence'),
      need(!(input$RSE_PrevR == "" ), 'Please provide a value RSE for recency of HIV prevalence'),
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
    temp<-prev_inc_calc_incprop(PrevH = input$PrevH/100, RSE_PrevH = input$RSE_PrevH/100,
                                PrevR = input$PrevR/100, RSE_PrevR = input$RSE_PrevR/100,
                                MDRI = input$MDRI, RSE_MDRI = input$RSE_MDRI/100,
                                FRR = input$FRR/100, RSE_FRR = input$RSE_FRR/100,
                                BigT = input$BigT)
    return(temp)
  }) 
  
  # the plot output rendered
  output$plot1<-renderPlot({
    data_pie()
  })
  # Produce an output table value. redundant on this version
  output$tab1 <- renderTable({
    data_prevalence()
    
  })
  output$tab2 <- renderTable({
    data_incidence()
  })
  output$tab3 <- renderTable({
    data_risk()
  })
  
  # need a tab that displays the results in a dataframe form.
  output$tab5<-renderTable(
    data_prev_inc_calc_incprop()
  )
  
 # A tab to display the all the incidence estimates from count data in the form of a dataframe
  output$tab4<-renderTable({
    temp<-data_incidence_count()
    data.frame("Parameter" = c("Prevalence of HIV (PrevH)","Prevalence of recency (PrevR)",
                      "Relative standard error of PrevH (RSE_PrevH)","Relative standard error of PrevR (RSE_PrevR)",
                      "Estimated incidence (Incidence)","Lower limit of confidence interval (CI.low)",
                      "Upper limit of confidence interval (CI.up)","Relative standard error of incidence estimate (RSE)",
                      "Relative standard error at infinite sample size (RSE.Inf.SS)",
                      "Annual Risk of Infection (ARI)",
                      "Lower confidence limit of Annual Risk of Infection (ARI.CI.low)",
                      "Upper confidence limit of Annual Risk of Infection (ARI.CI.up)"),
               "Value"=c(temp$PrevH,temp$PrevR,temp$RSE_PrevH,temp$RSE_PrevR,
                          temp$Incidence,temp$CI.low,temp$CI.up,temp$RSE,
                          temp$RSE.Inf.SS,temp$ARI,temp$ARI.CI.low,temp$ARI.CI.up))

  })

  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste("resultTable-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      temp<-data_incidence_count()
      temp<-data.frame("Parameter" = names(temp), 
                       "Value"=c(temp$PrevH,temp$PrevR,temp$RSE_PrevH,temp$RSE_PrevR,
                            temp$Incidence,temp$CI.low,temp$CI.up,temp$RSE,
                            temp$RSE.Inf.SS,temp$ARI,temp$ARI.CI.low,temp$ARI.CI.up))
      tt=xtabs(Value~Parameter,data = temp) 
      write.csv(temp, file)
    }
  )
  
  # A tab to display the incidence estimates from proportions in a dataframe form
  output$tab5b<-renderTable({
    temP<-data_prev_inc_calc_incprop()
    data.frame("Parameter" = c("Estimated incidence (Incidence)","Lower limit of confidence interval (CI.low)",
                               "Upper limit of confidence interval (CI.up)","Relative standard error of incidence estimate (RSE)",
                               "Annual Risk of Infection (ARI)",
                               "Lower confidence limit of Annual Risk of Infection (ARI.CI.low)",
                               "Upper confidence limit of Annual Risk of Infection (ARI.CI.up)"),
               "Value"=c(temP$Incidence,temP$CI.low,temP$CI.up,temP$RSE,temP$ARI,temP$ARI.CI.low,temP$ARI.CI.up))
    
  })
  
   	 	 	 	
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste("resultTable-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      temp<-data_prev_inc_calc_incprop()
      temp<-data.frame("Parameter" = names(temp), 
                       "Value"=c(temp$Incidence,temp$CI.low,temp$CI.up,temp$RSE,
                                 temp$ARI,temp$ARI.CI.low,temp$ARI.CI.up))
      tt=xtabs(Value~Parameter,data = temp) 
      write.csv(temp, file)
    }
  )

 

})
