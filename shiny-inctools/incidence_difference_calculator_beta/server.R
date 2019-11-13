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
source('Incdiff_incprop_calc.R')


shinyServer(function(input, output, session){

  #shinyURL.server(session)
  incdiff.data <- reactive({  # reactive function for incidence difference
    temp <- incdiff_calc(PrevH_1 = input$PrevH_1/100, RSE_PrevH_1 = input$RSE_PrevH_1/100,
                         PrevR_1 = input$PrevR_1/100, RSE_PrevR_1 = input$RSE_PrevR_1/100,
                         MDRI_1 = input$MDRI_1, RSE_MDRI_1 = input$RSE_MDRI_1/100,
                         FRR_1 = input$FRR_1/100, RSE_FRR_1 = input$RSE_FRR_1/100,
                         PrevH_2 = input$PrevH_2/100, RSE_PrevH_2 = input$RSE_PrevH_2/100,
                         PrevR_2 = input$PrevR_2/100, RSE_PrevR_2 = input$RSE_PrevR_2/100,
                         MDRI_2 = input$MDRI_2, RSE_MDRI_2 = input$RSE_MDRI_2/100,
                         FRR_2 = input$FRR_2/100, RSE_FRR_2 = input$RSE_FRR_2/100,
                         COV_FRR=input$COV_FRR, COV_MDRI=input$COV_MDRI,
                         BigT = input$BigT )
    return(temp)

  })
  
  # Produce an output table value.
  output$tab <- renderTable({
    validate(
      need(input$BigT, 'Please provide a value for the cut-off time'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off time (>120)'),
      need(!(input$COV_FRR==0 & input$COV_MDRI==1), 'Please provide a valid combination for the covariance between FRR and MDRI'),
      need(input$RSE_FRR_1 >= 0, 'Please provide a valid RSE for FRR of survey 1'),
      need(input$RSE_FRR_1 <= 100, 'Please provide a valid RSE for FRR of survey 1'),
      need(input$RSE_MDRI_1 >= 0, 'Please provide a valid RSE for MDRI of survey 1'),
      need(input$RSE_MDRI_1, 'Please provide a valid prevalence value for survey 1'),
      need(input$PrevH_1 <= 100, 'Please provide a valid prevalence value for survey 1'),
      need(input$PrevR_1, 'Please provide an prevalence of recent infection value for survey 1 '),
      need(input$PrevR_1 >= 0, 'Please provide a valid prevalence of recent infection value for survey 1'),
      need(input$PrevR_1 <= 100, 'Please provide a valid prevalence of recent infection value for survey 1'),
      need(input$RSE_PrevH_1, 'Please provide an RSE of prevalence value for survey 1 '),
      need(input$RSE_PrevH_1 >= 0, 'Please provide a valid RSE of prevalence value for survey 1'),
      need(input$RSE_PrevH_1 <= 100, 'Please provide a valid RSE of prevalence value for survey 1'),
      need(input$RSE_PrevR_1, 'Please provide an RSE of prevalence of recent infection value for survey 1 '),
      need(input$RSE_PrevR_1 >= 0, 'Please provide a valid RSE of prevalence of recent infection value for survey 1'),
      need(input$RSE_PrevR_1 <= 100, 'Please provide a valid RSE of prevalence of recent infection value for survey 1'),
      need(input$RSE_FRR_2 >= 0, 'Please provide a valid RSE for FRR of survey 2'),
      need(input$RSE_FRR_2 <= 100, 'Please provide a valid RSE for FRR of survey 2'),
      need(input$RSE_MDRI_2 >= 0, 'Please provide a valid RSE for MDRI of survey 2'),
      need(input$RSE_MDRI_2, 'Please provide a valid prevalence value for survey 2'),
      need(input$PrevH_2 <= 100, 'Please provide a valid prevalence value for survey 2'),
      need(input$PrevR_2, 'Please provide an prevalence of recent infection value for survey 2 '),
      need(input$PrevR_2 >= 0, 'Please provide a valid prevalence of recent infection value for survey 2'),
      need(input$PrevR_2 <= 100, 'Please provide a valid prevalence of recent infection value for survey 2'),
      need(input$RSE_PrevH_2, 'Please provide an RSE of prevalence value for survey 2'),
      need(input$RSE_PrevH_2 >= 0, 'Please provide a valid RSE of prevalence value for survey 2'),
      need(input$RSE_PrevH_2 <= 100, 'Please provide a valid RSE of prevalence value for survey 2'),
      need(input$RSE_PrevR_2, 'Please provide an RSE of prevalence of recent infection value for survey 2 '),
      need(input$RSE_PrevR_2 >= 0, 'Please provide a valid RSE of prevalence of recent infection value for survey 2'),
      need(input$RSE_PrevR_2 <= 100, 'Please provide a valid RSE of prevalence of recent infection value for survey 2'),
      need(input$MDRI_1, 'Please provide a  value for MDRI for survey 1'),
      need(input$MDRI_1 >= 0, 'Please provide a valid  value for MDRI for survey 1'),
      #need(input$MDRI_1 <= 720, 'Please provide a valid  value for MDRI for survey 1'),
      need(input$MDRI_2, 'Please provide a  value for MDRI for survey 2'),
      need(input$MDRI_2 >= 0, 'Please provide a valid  value for MDRI for survey 2')
      #need(input$MDRI_2 <= 720, 'Please provide a valid  value for MDRI for survey 2'),
    )
    incdiff.data()

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
