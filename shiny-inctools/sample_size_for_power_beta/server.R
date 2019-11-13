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

    if(0 == input$COV_FRR & 0==input$COV_MDRI) {
      temp <- mdply(expand.grid(Power = seq(input$Power_range[1],input$Power_range[2], by = 0.01), 
                                MDRI_1 = input$MDRI_1, 
                                BigT = input$BigT,
                                FRR_1 = (1/100)*input$FRR_1,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1, 
                                RSE_MDRI_1 = (1/100)*input$RSE_MDRI_1, 
                                DE_H1 = input$DE_H1, DE_H2 = input$DE_H2,
                                DE_R1 = input$DE_R1, DE_R2 = input$DE_R2,
                                I1 = (1/100)*input$I1, I2 = (1/100)*input$I2,
                                PrevH1 = (1/100)*input$PrevH1, PrevH2 = (1/100)*input$PrevH2),
                    alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, COV_FRR=input$COV_FRR, COV_MDRI=input$COV_MDRI)
                        
      return(temp)
    }

    if(0 == input$COV_FRR & 1==input$COV_MDRI) {
      temp <- mdply(expand.grid(Power = seq(input$Power_range[1],input$Power_range[2], by = 0.01), 
                                MDRI_1 = input$MDRI_1,
                                MDRI_2 = input$MDRI_2,
                                RSE_MDRI_1 = (1/100)*input$RSE_MDRI_1,
                                RSE_MDRI_2 = (1/100)*input$RSE_MDRI_2,
                                BigT = input$BigT,
                                FRR_1 = (1/100)*input$FRR_1,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1,
                                DE_H1 = input$DE_H1, DE_H2 = input$DE_H2,
                                DE_R1 = input$DE_R1, DE_R2 = input$DE_R2,
                                I1 = (1/100)*input$I1, I2 = (1/100)*input$I2,
                                PrevH1 = (1/100)*input$PrevH1, PrevH2 = (1/100)*input$PrevH2),
                    alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, COV_FRR=input$COV_FRR, COV_MDRI=input$COV_MDRI)
      return(temp)
    }
    if(1 == input$COV_FRR & 1==input$COV_MDRI) {
      temp <- mdply(expand.grid(Power = seq(input$Power_range[1],input$Power_range[2], by = 0.01),
                                MDRI_1 = input$MDRI_1, MDRI_2 = input$MDRI_2,
                                BigT = input$BigT,
                                FRR_1 = (1/100)*input$FRR_1, FRR_2 = (1/100)*input$FRR_2,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1, RSE_FRR_2 = (1/100)*input$RSE_FRR_2,
                                RSE_MDRI_1 = (1/100)*input$RSE_MDRI_1, RSE_MDRI_2 = (1/100)*input$RSE_MDRI_2,
                                DE_H1 = input$DE_H1, DE_H2 = input$DE_H2,
                                DE_R1 = input$DE_R1, DE_R2 = input$DE_R2,
                                I1 = (1/100)*input$I1, I2 = (1/100)*input$I2,
                                PrevH1 = (1/100)*input$PrevH1, PrevH2 = (1/100)*input$PrevH2),
                    alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, COV_FRR=input$COV_FRR, COV_MDRI=input$COV_MDRI)
      return(temp)
    }
    if(1 == input$COV_FRR & 0==input$COV_MDRI) {
      temp <- mdply(expand.grid(Power = seq(input$Power_range[1],input$Power_range[2], by = 0.01),
                                MDRI_1 = input$MDRI_1, 
                                RSE_MDRI_1 = (1/100)*input$RSE_MDRI_1,
                                BigT = input$BigT,
                                FRR_1 = (1/100)*input$FRR_1,
                                FRR_2 = (1/100)*input$FRR_2,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1,
                                RSE_FRR_2 = (1/100)*input$RSE_FRR_2,
                                DE_H1 = input$DE_H1, DE_H2 = input$DE_H2,
                                DE_R1 = input$DE_R1, DE_R2 = input$DE_R2,
                                I1 = (1/100)*input$I1, I2 = (1/100)*input$I2,
                                PrevH1 = (1/100)*input$PrevH1, PrevH2 = (1/100)*input$PrevH2),
                    alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, COV_FRR=input$COV_FRR, COV_MDRI=input$COV_MDRI)
      return(temp)
    }

  })
  

  
  output$plot <- renderPlot({
    validate(
      need(input$RSE_FRR_1 >= 0, 'Please provide a valid RSE for FRR of survey 1'),
      need(input$RSE_FRR_1 <= 100, 'Please provide a valid RSE for FRR of survey 1'),
      need(input$RSE_MDRI_1 >= 0, 'Please provide a valid RSE for MDRI of survey 1'),
      need(input$RSE_MDRI_1, 'Please provide a valid prevalence value for survey 1'),
      need(input$RSE_FRR_2 >= 0, 'Please provide a valid RSE for FRR of survey 2'),
      need(input$RSE_FRR_2 <= 100, 'Please provide a valid RSE for FRR of survey 2'),
      need(input$RSE_MDRI_2 >= 0, 'Please provide a valid RSE for MDRI of survey 2'),
      need(input$RSE_MDRI_2, 'Please provide a valid prevalence value for survey 2'),
      need(!(input$COV_FRR==0 & input$COV_MDRI==1), 'Please provide a valid combination for the covariance between FRR and MDRI'),
      need(input$MDRI_1, 'Please provide a  value for MDRI for survey 1'),
      need(input$MDRI_1 >= 0, 'Please provide a valid  value for MDRI for survey 1'),
      #need(input$MDRI_1 <= 720, 'Please provide a valid  value for MDRI for survey 1'),
      need(input$MDRI_2, 'Please provide a  value for MDRI for survey 2'),
      need(input$MDRI_2 >= 0, 'Please provide a valid  value for MDRI for survey 2'),
      #need(input$MDRI_2 <= 720, 'Please provide a valid  value for MDRI for survey 2'),
      need(input$I1, 'Please provide an incidence value for survey 1'),
      need(input$I2, 'Please provide an incidence value for survey 2'),
      need(input$I1 > 0, 'Please provide a valid incidence value for survey 1 (>0)'),
      need(input$I2 > 0, 'Please provide a valid incidence value for survey 2 (>0)'),
      need(input$I2 < input$I1, 'Incidence in survey 1 must be >=  Incidence in survey 2'),
      need(input$BigT, 'Please provide a value for the cut-off BigT'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off BigT (>120)'),
      need(input$PrevH1, 'Please provide an prevalence value for survey 1'),
      need(input$PrevH2, 'Please provide an prevalence value for survey 2'),
      need(input$PrevH1 > 0, 'Please provide a valid prevalence value for survey 1 (>0)'),
      need(input$PrevH2 > 0, 'Please provide a valid prevalence value for survey 2 (>0)'),
      need(input$alpha, 'Please provide a valid value for alpha (0:1)'),
      need(input$alpha > 0, 'Please provide a valid value for alpha (0:1)'),
      need(input$alpha <= 1, 'Please provide a valid value for alpha (0:1)'),
      need(input$DE_H1, 'Please provide a D.E. value for survey 1'),
      need(input$DE_R1, 'Please provide a D.E. value for survey 1'),
      need(input$DE_H2, 'Please provide a D.E. value for survey 2'),
      need(input$DE_R2, 'Please provide a D.E. value for survey 2'),
      need(input$CR_1, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_2, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_1 > 0, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_1 <= 100, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_2 > 0, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_2 <= 100, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$statPower>=0 & input$statPower<=100,"Please provide a valid value for statistical Power (%)")
    )

    data <- format_table(df())
    #data<-df()
    plot<- plot(data[,length(names(data))],data[,1],
                main = "Probability of correcting inferring incidence 1 > incidence 2 \n as a function of sample size",
                xlab = "Sample Size common to survey 1 and survey 2", 
                ylab = "Probability ",type = "l", col='blue')
    abline(h = input$statPower/100,  lty=2, col="grey", lwd=2)
    abline(v = data[which(round(data[,1],2)==input$statPower/100),"N"], lty=2, col="grey", lwd=2)
    print(plot)
  })

  ## Define the function renametable() to renames the variable names to use inctools 
  #  variable names
  
  renameTable<-function(data){
    colnames(data)[which(colnames(data)=="Power")] <- "Power"
    if("MDRI_1" %in% colnames(data)) {
      data$MDRI_1 <- format(data$MDRI_1, digits = 0) }
    if("MDRI_2" %in% colnames(data)) {
      data$MDRI_2 <- format(data$MDRI_2, digits = 0) }
    data$N <- format(data$N, digits = 10)
    colnames(data)[which(colnames(data)=="N")] <- "ss"
    if("FRR_1" %in% colnames(data)) data$FRR_1 <- format(data$FRR_1, digits = 2)
    if("FRR_2" %in% colnames(data)) data$FRR_2 <- format(data$FRR_2, digits = 2)
    if("RSE_MDRI_1" %in% colnames(data)) {
      data$RSE_MDRI_1 <- format(data$RSE_MDRI_1, digits = 4)
      colnames(data)[which(colnames(data)=="RSE_MDRI_1")] <- "RSE_MDRI_1" }
    if("RSE_MDRI_2" %in% colnames(data)) {
      data$RSE_MDRI_2 <- format(data$RSE_MDRI_2, digits = 4)
      colnames(data)[which(colnames(data)=="RSE_MDRI_2")] <- "RSE_MDRI_2" }
    if("RSE_FRR_1" %in% colnames(data)){
      data$RSE_FRR_1 <- format(data$RSE_FRR_1, digits = 4)
      colnames(data)[which(colnames(data)=="RSE_FRR_1")] <- "RSE_FRR_1" } 
    if("RSE_FRR_2" %in% colnames(data)){
      data$RSE_FRR_2 <- format(data$RSE_FRR_2, digits = 4)
      colnames(data)[which(colnames(data)=="RSE_FRR_2")] <- "RSE_FRR_2" } 
    data$DE_H1 <- format(data$DE_H1, digits = 2)
    colnames(data)[which(colnames(data)=="DE_H1")] <- "DE_H1"
    data$DE_H2 <- format(data$DE_H2, digits = 2)
    colnames(data)[which(colnames(data)=="DE_H2")] <- "DE_H2"
    data$DE_R1 <- format(data$DE_R1, digits = 2)
    colnames(data)[which(colnames(data)=="DE_R1")] <- "DE_R1"
    data$DE_R2 <- format(data$DE_R2, digits = 2)
    colnames(data)[which(colnames(data)=="DE_R2")] <- "DE_R2"
    data[, "Incidence_1 (%)"] <- format(data[, "Incidence_1 (%)"], digits = 2)
    colnames(data)[which(colnames(data)=="Incidence_1 (%)")] <- "I1"
    data[, "Incidence_2 (%)"] <- format(data[, "Incidence_2 (%)"], digits = 2)
    colnames(data)[which(colnames(data)=="Incidence_2 (%)")] <- "I2"
    data[, "Prevalence_1 (%)"] <- format(data[, "Prevalence_1 (%)"], digits = 3)
    colnames(data)[which(colnames(data)=="Prevalence_1 (%)")] <- "PrevH1"
    data[, "Prevalence_2 (%)"] <- format(data[, "Prevalence_2 (%)"], digits = 3)
    colnames(data)[which(colnames(data)=="Prevalence_2 (%)")] <- "PrevH2"
    data$BigT <- format(data$BigT, digits = 4)
    colnames(data)[which(colnames(data)=="BigT")] <- "BigT"
    
    data.frame(data)
  }
  
  ### Define the function format_table()  (reference: DS(FIND))
  format_table <- function(df) {
    df$I1 <- 100*df$I1
    df$I2 <- 100*df$I2
    df$PrevH1 <- 100*df$PrevH1
    df$PrevH2 <- 100*df$PrevH2
    if("FRR" %in% colnames(df)){
      df$FRR <- 100*df$FRR}
    if("FRR_1" %in% colnames(df)){
      df$FRR_1 <- 100*df$FRR_1 }
    if("FRR_2" %in% colnames(df)){
      df$FRR_2 <- 100*df$FRR_2 }
    colnames(df)[which(colnames(df) == "FRR")] <- "FRR"
    colnames(df)[which(colnames(df) == "FRR_1")] <- "FRR_1"
    colnames(df)[which(colnames(df) == "FRR_2")] <- "FRR_2"
    colnames(df)[which(colnames(df) == "MDRI")] <- "MDRI"
    colnames(df)[which(colnames(df) == "MDRI_1")] <- "MDRI_1"
    colnames(df)[which(colnames(df) == "MDRI_2")] <- "MDRI_2"
    colnames(df)[which(colnames(df) == "V1")] <- "N"
    colnames(df)[which(colnames(df) == "I1")] <- "Incidence_1 (%)"
    colnames(df)[which(colnames(df) == "I2")] <- "Incidence_2 (%)"
    colnames(df)[which(colnames(df) == "PrevH1")] <- "Prevalence_1 (%)"
    colnames(df)[which(colnames(df) == "PrevH2")] <- "Prevalence_2 (%)"
    df$N <- as.numeric(df$N)
    df$N[which(df$N<0)] <- NA
    return(df)
  }
  
  # Produce an output table value.
  output$tab <- renderTable({
    data <- format_table(df())
    renameTable(data)  # We call the function rename table here
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(renameTable(format_table(df())) , file)
    }
  )

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".jpeg", sep="")
    },
    content = function(file) {
      data <- format_table(df())
      jpeg(filename = file)
      plot<- plot(data[,length(names(data))],data[,1],
                  main = "Probability of correcting inferring incidence 1 > incidence 2 \n as a function of sample size",
                  xlab = "Sample Size common to survey 1 and survey 2", 
                  ylab = "Probability",type = "l",col='blue')
      abline(h = input$statPower/100,  lty=2, col="grey")
      abline(v = data[which(round(data[,1],2)==input$statPower/100),"N"])
      print(plot)
      dev.off()
  
    }
  )
  
 

})
