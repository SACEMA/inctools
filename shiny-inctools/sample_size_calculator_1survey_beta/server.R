# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND).
# Recoded by Lamin Juwara (McGill)(2017/18)
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

library(shiny)
#library(shinyURL)
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)
library(grid)
library(inctools)
source('incidence_difference_calculator.R')
source('plot_fcn.R')

shinyServer(function(input, output, session) {

  myData <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath, header = TRUE)
    data
  })

  output$contents <- renderTable({
    if (is.null(input$file1))
      return(NULL)
    myData()
  })

  output$tab <- renderTable({
    if (is.null(input$file1))
      return(NULL)

    cd <- myData()
    temp <- do_table(cd = myData(), DE_prev = input$DE_prev, DE_RgivenTested = input$DE_RgivenTested,
             alpha = input$alpha, inc_cov = input$inc_cov/100, rec_cov = input$rec_cov/100, TIME = input$TIME,
             mdrihatcov = input$mdrihatcov/100, frrhatcov = input$frrhatcov/100)
    temp$scenario <- cd$scenario
    temp
  }
  )
output$downloadData <- downloadHandler(
  filename = function() {
    paste0("ss_1survey_", "alpha_", input$alpha,".csv")},
    #paste(input$file1, '.csv', sep='')},
  content = function(file) {

    cd <- myData()
    temp <- do_table(cd = cd, DE_prev = input$DE_prev, DE_RgivenTested = input$DE_RgivenTested,
                     alpha = input$alpha, inc_cov = input$inc_cov/100, rec_cov = input$rec_cov/100, TIME = input$TIME,
                     mdrihatcov = input$mdrihatcov/100, frrhatcov = input$frrhatcov/100)
    temp
    
    temp$scenario <- cd$scenario
    write.csv(temp, file, row.names = FALSE)
  }
)
}
)


