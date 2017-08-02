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
source('rse_calc_ss_v297.R')


shinyServer(function(input, output, session){

  #shinyURL.server(session)
  df <- reactive({
    temp <- mdply(expand.grid(n = seq(input$n_range[1],input$n_range[2], by = 100),
                              I = input$I/100,
                              PrevH = input$PrevH/100,
                              CR = input$CR/100, 
                              MDRI = input$MDRI, RSE_MDRI = input$RSE_MDRI/100,
                              FRR = input$FRR/100, RSE_FRR = input$RSE_FRR/100,
                              BigT = input$BigT),  
                  ss_calc_precision)
    return(temp)

  })
  
  ## Define the function renametable() to renames the variable names to use inctools 
  #  variable names
  
  renameTable<-function(data){ # the structure of the output value changes with the same value.
    #data <- data[,c(1:9,11)]
    names(data)<- c("n","I", "PrevH","CR", "MDRI", "RSE_MDRI","FRR",
                    "RSE_FRR","BigT","RSE_I")
    data.frame(data)
  }
  
# Output value  
  output$plot <- renderPlot({
    validate(
      need(input$RSE_FRR >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR <= 100, 'Please provide a valid RSE for FRR'),
      need(!(input$RSE_FRR == "" ), 'Please provide a value for RSE_FRR'),
      need(input$RSE_MDRI >= 0, 'Please provide a valid RSE for MDRI'),
      need(input$RSE_MDRI <= 100, 'Please provide a valid RSE for MDRI'),
      need(!(input$RSE_MDRI == "" ), 'Please provide a value for RSE_MDRI'),
      need(input$MDRI >= 0, 'Please provide a valid value for MDRI'),
      need(input$FRR >= 0, 'Please provide a valid value for FRR'),
      need(input$FRR <= 100, 'Please provide a valid value for FRR'),
      need(input$I, 'Please provide an incidence value'),
      need(input$I > 0, 'Please provide a valid incidence value for survey (>0)'),
      need(input$BigT, 'Please provide a value for the cut-off time'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off time (>120)'),
      need(input$PrevH, 'Please provide an prevalence value for survey '),
      need(input$DE_H, 'Please provide a D.E. value for survey '),
      need(input$DE_R, 'Please provide a D.E. value for survey '),
      need(input$CR, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$RSE_infSS>=0 & input$RSE_infSS<=1, "Please provide a RSE of incidence estimate at infinite sample size between 0 and 1.")
      
    )

    data <- renameTable(df())
    
    #plot<- plot(data[,1],data[,length(names(data))],
    plot<- plot(data[,"n"],data[,"RSE_I"],
                main = "Relative standard error of incidence estimate as \n a function of sample size",
                xlab = "Sample Size", 
                ylab = "Relative standard error",type = "l",col='blue',ylim = c(0,0.8))
    abline(h = c(input$RSE_infSS,input$RSE_req_Inc), lty=c(2,2), col=c("red","grey"),lwd=c(1,2))
    abline(v = data[which(round(data[,"RSE_I"],2)==input$RSE_req_Inc),"n"],
           lty=2, col="grey",lwd=2)
    print(plot)
  })

  # Produce an output table value.
  output$tab <- renderTable({
    data <- df()
    renameTable(data)  # We call the function rename table here
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("downloadData-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(renameTable(renameTable(df()) ) , file)
    }
  )

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("downloadPlot-", Sys.Date(), ".jpeg", sep="")
    },
    content = function(file) {
      ##
      data <- renameTable(df())
      jpeg(filename = file)
      plot<- plot(data[,"n"],data[,"RSE_I"],
                  main = "Relative standard error of incidence estimate as \n a function of sample size",
                  xlab = "Sample Size", 
                  ylab = "Relative standard error",type = "l",col='blue',ylim = c(0,0.8))
      abline(h = c(input$RSE_infSS,input$RSE_req_Inc), lty=c(2,2), col=c("red","grey"),lwd=c(1,2))
      abline(v = data[which(round(data[,"RSE_I"],2)==round(input$RSE_req_Inc,2)),"n"],
             lty=2, col="grey",lwd=2)
      print(plot)
      dev.off()
 
    }
  )

})
