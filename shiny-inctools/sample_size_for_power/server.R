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
source('sample_size_calc_v297.R')


shinyServer(function(input, output, session) {

  #shinyURL.server(session)
  df <- reactive({

    if(1 == input$scenario_case) {
      temp <- mdply(expand.grid(power = seq(input$power_range[1],input$power_range[2], by = 0.01), 
                                MDRI = input$MDRI, 
                                frrhat = (1/100)*input$frrhat,
                                TIME = input$TIME,
                                frrhatcov = (1/100)*input$frrhatcov,
                                mdrihatcov = (1/100)*input$mdrihatcov,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
                        
      return(temp)
    }

    if(2 == input$scenario_case) {
      temp <- mdply(expand.grid(power = seq(input$power_range[1],input$power_range[2], by = 0.01), 
                                MDRI = input$MDRI,
                                frrhat_1 = (1/100)*input$frrhat_1, frrhat_2 = (1/100)*input$frrhat_2,
                                TIME = input$TIME,
                                frrhatcov_1 = (1/100)*input$frrhatcov_1, frrhatcov_2 = (1/100)*input$frrhatcov_2,
                                mdrihatcov = (1/100)*input$mdrihatcov,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }
    if(3 == input$scenario_case) {
      temp <- mdply(expand.grid(power = seq(input$power_range[1],input$power_range[2], by = 0.01),
                                MDRI_1 = input$MDRI_1, MDRI_2 = input$MDRI_2,
                                TIME = input$TIME,
                                frrhat_1 = (1/100)*input$frrhat_1, frrhat_2 = (1/100)*input$frrhat_2,
                                frrhatcov_1 = (1/100)*input$frrhatcov_1, frrhatcov_2 = (1/100)*input$frrhatcov_2,
                                mdrihatcov_1 = (1/100)*input$mdrihatcov_1, mdrihatcov_2 = (1/100)*input$mdrihatcov_2,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

  })
  

  
  output$plot <- renderPlot({
    validate(
      need(!(input$frrhatcov_1 < 0  & input$scenario_case >1), 'Please provide a value for FRR_1 covariance'),
      need(!(input$frrhatcov_1 > 100  & input$scenario_case >1), 'Please provide a value for FRR_1 covariance'),
      need(!(input$frrhatcov_2 < 0  & input$scenario_case >1), 'Please provide a value for FRR_2 covariance'),
      need(!(input$frrhatcov_2 > 100  & input$scenario_case >1), 'Please provide a value for FRR_2 covariance'),
      need(!(input$mdrihatcov_1 < 0  & input$scenario_case == 3), 'Please provide a value for MDRI_1 covariance'),
      need(!(input$mdrihatcov_1 > 100  & input$scenario_case == 3), 'Please provide a value for MDRI_1 covariance'),
      need(!(input$mdrihatcov_2 < 0  & input$scenario_case == 3), 'Please provide a value for MDRI_2 covariance'),
      need(!(input$mdrihatcov_2 > 100  & input$scenario_case == 3), 'Please provide a value for MDRI_2 covariance'),
      need(input$frrhatcov >= 0, 'Please provide a valid covariance value for FRR'),
      need(input$frrhatcov <= 100, 'Please provide a valid covariance value for FRR'),
      need(input$mdrihatcov, 'Please provide a covariance value for MDRI'),
      need(input$mdrihatcov >= 0, 'Please provide a valid covariance value for MDRI'),
      need(input$mdrihatcov <= 100, 'Please provide a valid covariance value for MDRI'),
      need(input$inc_1, 'Please provide an incidence value for survey 1'),
      need(input$inc_2, 'Please provide an incidence value for survey 2'),
      need(input$inc_1 > 0, 'Please provide a valid incidence value for survey 1 (>0)'),
      need(input$inc_2 > 0, 'Please provide a valid incidence value for survey 2 (>0)'),
      need(input$inc_2 < input$inc_1, 'Incidence in survey 1 must be >=  Incidence in survey 2'),
      need(input$TIME, 'Please provide a value for the cut-off time'),
      need(input$TIME > 120, 'Please provide a valid value for the cut-off time (>120)'),
      need(input$p_pos_1, 'Please provide an prevalence value for survey 1'),
      need(input$p_pos_2, 'Please provide an prevalence value for survey 2'),
      need(input$p_pos_1 > 0, 'Please provide a valid prevalence value for survey 1 (>0)'),
      need(input$p_pos_2 > 0, 'Please provide a valid prevalence value for survey 2 (>0)'),
      need(input$alpha, 'Please provide a valid value for alpha (0:1)'),
      need(input$alpha > 0, 'Please provide a valid value for alpha (0:1)'),
      need(input$alpha <= 1, 'Please provide a valid value for alpha (0:1)'),
      need(input$DE_prev_1, 'Please provide a D.E. value for survey 1'),
      need(input$DE_RgivenTested_1, 'Please provide a D.E. value for survey 1'),
      need(input$DE_prev_2, 'Please provide a D.E. value for survey 2'),
      need(input$DE_RgivenTested_2, 'Please provide a D.E. value for survey 2'),
      need(input$rec_test_coverage_1, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$rec_test_coverage_2, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$rec_test_coverage_1 > 0, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$rec_test_coverage_1 <= 100, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$rec_test_coverage_2 > 0, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$rec_test_coverage_2 <= 100, 'Please provide a valid value for the % of recently tested for positives')
    )

    data <- format_table(df(), input$scenario_case)
  
    plot<- plot(data[,length(names(data))],data[,1],
                main = "Probability of correcting inferring incidence 1 > incidence 2 \n as a function of sample size",
                xlab = "Sample Size common to survey 1 and survey 2", 
                ylab = "Probability ",type = "l", col='blue')
    #abline(h = 0.8, v = data[which(data[,1]==0.8),"N"])
    abline(h = input$statPower,  lty=2, col="grey", lwd=2)
    abline(v = data[which(round(data[,1],2)==input$statPower),"N"], lty=2, col="grey", lwd=2)
    print(plot)
  })

  ## Define the function renametable() to renames the variable names to use inctools 
  #  variable names
  
  renameTable<-function(data){
    colnames(data)[which(colnames(data)=="power")] <- "Power"
    if("MDRI" %in% colnames(data)){
      data$MDRI <- format(data$MDRI, digits = 0) } 
    if("MDRI_1" %in% colnames(data)) {
      data$MDRI_1 <- format(data$MDRI_1, digits = 0) }
    if("MDRI_2" %in% colnames(data)) {
      data$MDRI_2 <- format(data$MDRI_2, digits = 0) }
    data$N <- format(data$N, digits = 10)
    colnames(data)[which(colnames(data)=="N")] <- "ss"
    if("FRR" %in% colnames(data)) data$FRR <- format(data$FRR, digits = 2)
    if("FRR_1" %in% colnames(data)) data$FRR_1 <- format(data$FRR_1, digits = 2)
    if("FRR_2" %in% colnames(data)) data$FRR_2 <- format(data$FRR_2, digits = 2)
    if("mdrihatcov" %in% colnames(data)) {
      data$mdrihatcov <- format(data$mdrihatcov, digits = 4)
      colnames(data)[which(colnames(data)=="mdrihatcov")] <- "RSE_MDRI" }
    if("mdrihatcov_1" %in% colnames(data)) {
      data$mdrihatcov_1 <- format(data$mdrihatcov_1, digits = 4)
      colnames(data)[which(colnames(data)=="mdrihatcov_1")] <- "RSE_MDRI_1" }
    if("mdrihatcov_2" %in% colnames(data)) {
      data$mdrihatcov_2 <- format(data$mdrihatcov_2, digits = 4)
      colnames(data)[which(colnames(data)=="mdrihatcov_2")] <- "RSE_MDRI_2" }
    if("frrhatcov" %in% colnames(data)){
      data$frrhatcov <- format(data$frrhatcov, digits = 4)
      colnames(data)[which(colnames(data)=="frrhatcov")] <- "RSE_FRR" } 
    if("frrhatcov_1" %in% colnames(data)){
      data$frrhatcov_1 <- format(data$frrhatcov_1, digits = 4)
      colnames(data)[which(colnames(data)=="frrhatcov_1")] <- "RSE_FRR_1" } 
    if("frrhatcov_2" %in% colnames(data)){
      data$frrhatcov_2 <- format(data$frrhatcov_2, digits = 4)
      colnames(data)[which(colnames(data)=="frrhatcov_2")] <- "RSE_FRR_2" } 
    data$DE_prev_1 <- format(data$DE_prev_1, digits = 2)
    colnames(data)[which(colnames(data)=="DE_prev_1")] <- "DE_H1"
    data$DE_prev_2 <- format(data$DE_prev_2, digits = 2)
    colnames(data)[which(colnames(data)=="DE_prev_2")] <- "DE_H2"
    data$DE_RgivenTested_1 <- format(data$DE_RgivenTested_1, digits = 2)
    colnames(data)[which(colnames(data)=="DE_RgivenTested_1")] <- "DE_R1"
    data$DE_RgivenTested_2 <- format(data$DE_RgivenTested_2, digits = 2)
    colnames(data)[which(colnames(data)=="DE_RgivenTested_2")] <- "DE_R2"
    data[, "Incidence_1 (%)"] <- format(data[, "Incidence_1 (%)"], digits = 2)
    colnames(data)[which(colnames(data)=="Incidence_1 (%)")] <- "I1"
    data[, "Incidence_2 (%)"] <- format(data[, "Incidence_2 (%)"], digits = 2)
    colnames(data)[which(colnames(data)=="Incidence_2 (%)")] <- "I2"
    data[, "Prevalence_1 (%)"] <- format(data[, "Prevalence_1 (%)"], digits = 3)
    colnames(data)[which(colnames(data)=="Prevalence_1 (%)")] <- "PrevH1"
    data[, "Prevalence_2 (%)"] <- format(data[, "Prevalence_2 (%)"], digits = 3)
    colnames(data)[which(colnames(data)=="Prevalence_2 (%)")] <- "PrevH2"
    data$TIME <- format(data$TIME, digits = 4)
    colnames(data)[which(colnames(data)=="TIME")] <- "BigT"
    
    data.frame(data)
  }
  
  ### Define the function format_table()  (reference: DS(FIND))
  format_table <- function(df, case_scenario) {
    df$inc_1 <- 100*df$inc_1
    df$inc_2 <- 100*df$inc_2
    df$p_pos_1 <- 100*df$p_pos_1
    df$p_pos_2 <- 100*df$p_pos_2
    if("frrhat" %in% colnames(df)){
      df$frrhat <- 100*df$frrhat}
    if("frrhat_1" %in% colnames(df)){
      df$frrhat_1 <- 100*df$frrhat_1 }
    if("frrhat_2" %in% colnames(df)){
      df$frrhat_2 <- 100*df$frrhat_2 }
    colnames(df)[which(colnames(df) == "frrhat")] <- "FRR"
    colnames(df)[which(colnames(df) == "frrhat_1")] <- "FRR_1"
    colnames(df)[which(colnames(df) == "frrhat_2")] <- "FRR_2"
    colnames(df)[which(colnames(df) == "mdrihat")] <- "MDRI"
    colnames(df)[which(colnames(df) == "mdrihat_1")] <- "MDRI_1"
    colnames(df)[which(colnames(df) == "mdrihat_2")] <- "MDRI_2"
    colnames(df)[which(colnames(df) == "V1")] <- "N"
    colnames(df)[which(colnames(df) == "inc_1")] <- "Incidence_1 (%)"
    colnames(df)[which(colnames(df) == "inc_2")] <- "Incidence_2 (%)"
    colnames(df)[which(colnames(df) == "p_pos_1")] <- "Prevalence_1 (%)"
    colnames(df)[which(colnames(df) == "p_pos_2")] <- "Prevalence_2 (%)"
    df$N <- as.numeric(df$N)
    df$N[which(df$N<0)] <- NA
    return(df)
  }
  
  # Produce an output table value.
  output$tab <- renderTable({
    data <- format_table(df(), input$scenario_case)
    renameTable(data)  # We call the function rename table here
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(renameTable(format_table(df(), input$scenario_case)) , file)
    }
  )

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".jpeg", sep="")
    },
    content = function(file) {
      data <- format_table(df(), input$scenario_case)
      jpeg(filename = file)
      plot<- plot(data[,length(names(data))],data[,1],
                  main = "Probability of correcting inferring incidence 1 > incidence 2 \n as a function of sample size",
                  xlab = "Sample Size common to survey 1 and survey 2", 
                  ylab = "Probability",type = "l",col='blue')
      abline(h = input$statPower,  lty=2, col="grey")
      abline(v = data[which(round(data[,1],2)==input$statPower),"N"])
      print(plot)
      dev.off()
      #write.csv(renameTable(format_table(df(), input$scenario_case)) , file)
      # ggsave(file = file, plot = plot, width = 13, height = 9)
      
    }
  )
  
 

})
