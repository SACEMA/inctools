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
source('Incdiff_calc.R')


shinyServer(function(input, output, session){

  #shinyURL.server(session)
  incdiff.data <- reactive({  # reactive function for incidence difference
    temp <- incdiff_calc(case = input$case, I1 = input$I1/100, I2 = input$I2/100, 
                         PrevH1 = input$PrevH1/100, PrevH2 = input$PrevH2/100,
                         alpha = input$alpha, 
                         Power = input$Power,
                         DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                         DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
                         FRR = input$FRR/100, RSE_FRR = input$RSE_FRR/100,
                         FRR_1 = input$FRR_2/100, RSE_FRR_1 = input$RSE_FRR_1/100, 
                         FRR_2 = input$FRR_2/100, RSE_FRR_2 = input$RSE_FRR_2/100,
                         MDRI = input$MDRI, RSE_MDRI = input$RSE_MDRI/100,
                         MDRI_1 = input$MDRI_1, MDRI_2 = input$MDRI_2,
                         RSE_MDRI_1 = input$RSE_MDRI_1/100, RSE_MDRI_2 = input$RSE_MDRI_2/100,
                         CR_1 = input$CR_1/100, CR_2 = input$CR_2/100, BigT = input$BigT)
    return(temp)

  })
  


  # Produce an output table value.
  output$tab <- renderTable({
    validate(
      need(!(input$RSE_FRR < 0 ), 'Please provide a value for RSE_FRR'),
      need(!(input$RSE_FRR > 100), 'Please provide a value for RSE_FRR'),
      need(input$RSE_FRR >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR <= 100, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR_1 >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR_1 <= 100, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR_2 >= 0, 'Please provide a valid RSE for FRR'),
      need(input$RSE_FRR_2 <= 100, 'Please provide a valid RSE for FRR'),
      need(input$RSE_MDRI >= 0, 'Please provide a valid RSE for MDRI'),
      need(input$RSE_MDRI <= 100, 'Please provide a valid RSE for MDRI'),
      need(input$BigT, 'Please provide a value for the cut-off time'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off time (>120)'),
      need(input$PrevH1, 'Please provide an prevalence value for survey '),
      need(input$PrevH1 >= 0, 'Please provide a valid prevalence value for survey'),
      need(input$PrevH1 <= 100, 'Please provide a valid prevalence value for survey'),
      need(input$PrevH2, 'Please provide an prevalence value for survey '),
      need(input$PrevH2 >= 0, 'Please provide a valid prevalence value for survey'),
      need(input$PrevH2 <= 100, 'Please provide a valid prevalence value for survey'),
      need(input$I1 >= 0, 'Please provide a valid incidence value for survey'),
      need(input$I1 <= 100, 'Please provide a valid incidence value for survey'),
      need(input$I2 >= 0, 'Please provide a valid incidence value for survey'),
      need(input$I2 <= 100, 'Please provide a valid incidence value for survey'),
      need(input$CR_1 >= 0, 'Please provide a valid value for HIV positives tested for recency in survey 1 (%)'),
      need(input$CR_1 <= 100, 'Please provide a valid value for HIV positives tested for recency in survey 1 (%)'),
      need(input$CR_2 >= 0, 'Please provide a valid value for HIV positives tested for recency in survey 2 (%)'),
      need(input$CR_2 <= 100, 'Please provide a valid value for HIV positives tested for recency in survey 2 (%)'),
      need(input$DE_H_1, 'Please provide a D.E. value for survey '),
      need(input$DE_R_1, 'Please provide a D.E. value for survey '),
      need(input$DE_H_2, 'Please provide a D.E. value for survey '),
      need(input$DE_R_2, 'Please provide a D.E. value for survey ')
    )
    data <- incdiff.data()

  })

  output$downloadData <- downloadHandler(
    filename = function(){
      paste("downloadData-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(incdiff.data() , file)
    }
  )

 

})
