# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND).
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
#sample size table
library(shiny)
#library(shinyURL)
#main_path <- 'D:/PROJECTS/AIDS/INCIDENCE_ESTIMATION/R/'
#setwd(paste0(main_path, 'abie-r/Sample.Size.Calculator'))
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)
library(grid)
#library(googleVis)
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



#
#   output$contents <- renderTable({
#
#     # input$file1 will be NULL initially. After the user selects
#     # and uploads a file, it will be a data frame with 'name',
#     # 'size', 'type', and 'datapath' columns. The 'datapath'
#     # column will contain the local filenames where the data can
#     # be found.
#
#     inFile <- input$file1
#
#     if (is.null(inFile))
#       return(NULL)
#
#     read.csv(inFile$datapath, header=input$header, sep=input$sep,
#              quote=input$quote)
#
#
#   })

  output$tab <- renderTable({
    #     input_df <- read.csv(inFile$datapath, header=input$header, sep=input$sep,
    #                          quote=input$quote)
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
#     cd$inc_1 <- cd$inc_1 * 1/100
#     cd$p_pos_1 <- cd$p_pos_1 * 1/100
#     cd$frrhat <- cd$frrhat * 1/100
#     df <- assemble_data(cd, percent_reduction = (100 - input$percent_reduction)/100,
#                         DE_p_1 = input$DE_p_1, DE_p_2 = input$DE_p_2,
#                         DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_1,
#                         power = input$power, alpha = input$alpha)
#
#     res_I001 <- mdply(df[, -5], ss_calc, case = 1)
#     res <- cbind(res_I001, country = df[,5])
#     id <- which(colnames(res) == "V1")
#     colnames(res)[id] <- "N"
#     res_clean <- res[, c("country", "inc_1", "p_pos_1", "MDRI", "frrhat", "N")]
#     res_clean[, c(2, 3, 5)] <- 100*res_clean[, c(2, 3, 5)]
#     res_clean$N <- round(res_clean$N, digits = 0)
#     res_clean[, 2:6] <- sapply(res_clean[, 2:6], FUN=function(x) prettyNum(x, big.mark=","))
#     res_clean
  }
  )
#
#   output$downloadData <- downloadHandler(
#
#     filename = function() { paste0(input$file1, "_results", '.csv', sep='') },
#     content = function(file) {
#       cd <- myData()
#       cd$inc_1 <- cd$inc_1 * 1/100
#       cd$p_pos_1 <- cd$p_pos_1 * 1/100
#       cd$frrhat <- cd$frrhat * 1/100
#       #     df <- assemble_data(cd, percent_reduction = 0.5, DE_p_1 = 1.3, DE_p_2 = 1.3, DE_R_1 = 1.3, DE_R_2 = 1.3,
#       #                         power = 0.8, alpha = 0.05)
#       df <- assemble_data(cd, percent_reduction = (100 - input$percent_reduction)/100,
#                           DE_p_1 = input$DE_p_1, DE_p_2 = input$DE_p_2,
#                           DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_1,
#                           power = input$power, alpha = input$alpha)
#
#       res_I001 <- mdply(df[, -5], ss_calc, case = 1)
#       res <- cbind(res_I001, country = df[,5])
#       id <- which(colnames(res) == "V1")
#       colnames(res)[id] <- "N"
#       #res
#       res_clean <- res[, c("country", "inc_1", "p_pos_1", "MDRI", "frrhat", "N")]
#       res_clean[, c(2, 3, 5)] <- 100*res_clean[, c(2, 3, 5)]
#       res_clean$N <- round(res_clean$N, digits = 0)
#       res_clean[, 2:6] <- sapply(res_clean[, 2:6], FUN=function(x) prettyNum(x, big.mark=","))
#       write.csv(res_clean, file)
#     }
#   )
# })

output$downloadData <- downloadHandler(
  filename = function() {
    paste0("ss_2surveys_", "alpha_", input$alpha, "_power_", input$power,".csv")},
    #paste(input$file1, '.csv', sep='')},
  content = function(file) {
#     temp <- do_table(myData(), input$percent_reduction,
#                  input$DE_p_1, input$DE_p_2,
#                  input$DE_R_1, input$DE_R_2,
#                  input$power, input$alpha)
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


# output$downloadData <- downloadHandler(
#   filename = function() {
#     paste(input$x_variable, '.csv', sep='')
#   },
#   content = function(file) {
#     write.csv(format_table(df(), input$scenario_case, input$x_variable), file)
#   }
# )
