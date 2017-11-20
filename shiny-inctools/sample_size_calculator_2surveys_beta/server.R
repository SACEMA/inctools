# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND).
# Recoded by Lamin Juwara (McGill) to use functions from R Package 'inctools'
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
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)
library(grid)
library(inctools)
source('sample_size_calc_v297.R')
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
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.

  output$tab <- renderTable({
    if (is.null(input$file1))
      return(NULL)
    cd <- myData()

    temp <- do_table(cd, input$percent_reduction,
                 input$DE_p_1, input$DE_p_2,
                 input$DE_R_1, input$DE_R_2,
                 input$power, input$alpha)
    #temp <- temp[, c("Scenario", "prev_1", "inc_1", "FRR", "MDRI", "Sample size")]
    colnames(temp)[1] <- "scenario"
    temp$scenario <- cd$scenario
    temp
  }
  )

output$downloadData <- downloadHandler(
  filename = function() {
    paste0("ss_2surveys_", "alpha_", input$alpha, "_power_", input$power,".csv")},
    #paste(input$file1, '.csv', sep='')},
  content = function(file) {
    cd <- myData()
    res <- do_table(myData(), input$percent_reduction,
                    input$DE_p_1, input$DE_p_2,
                    input$DE_R_1, input$DE_R_2,
                    input$power, input$alpha)
    #res <- res[, c("Scenario", "prev_1", "inc_1", "FRR", "MDRI", "Sample size")]
    colnames(res)[1] <- "scenario"
    res$scenario <- cd$scenario
    write.csv(res, file, row.names = FALSE)
  }
)
}
)


