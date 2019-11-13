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

#server.R
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

  #shinyURL.server(session)
  df <- reactive({

    if("MDRI" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = seq(input$MDRI_range[1], input$MDRI_range[2], by = 30), FRR = (1/100)*input$FRR,
                                BigT = input$BigT,
                                RSE_FRR = (1/100)*input$RSE_FRR,
                                RSE_MDRI = (1/100)*input$RSE_MDRI,
                                DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
                                I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
                                PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
                    Power = input$Power, alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, case = input$scenario_case)
                        
      return(temp)
    }

    if("MDRI" == input$x_variable & 2 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = seq(input$MDRI_range[1], input$MDRI_range[2], by = 30),
                                FRR_1 = (1/100)*input$FRR_1, FRR_2 = (1/100)*input$FRR_2,
                                BigT = input$BigT,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1, RSE_FRR_2 = (1/100)*input$RSE_FRR_2,
                                RSE_MDRI = (1/100)*input$RSE_MDRI,
                                DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
                                I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
                                PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
                    Power = input$Power, alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("FRR" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = input$MDRI,
                                FRR = seq((1/100)*input$FRR_range[1], (1/100)*input$FRR_range[2], by = 0.005),
                                BigT = input$BigT,
                                RSE_FRR = (1/100)*input$RSE_FRR, RSE_MDRI = (1/100)*input$RSE_MDRI,
                                DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
                                I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
                                PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
                    Power = input$Power, alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("FRR" == input$x_variable & 2 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = input$MDRI,
                                BigT = input$BigT,
                                FRR_1 = (1/100)*input$FRR_1, FRR_2 = (1/100)*input$FRR_2,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1, RSE_FRR_2 = (1/100)*input$RSE_FRR_2,
                                RSE_MDRI = (1/100)*input$RSE_MDRI,
                                DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
                                I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
                                PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
                    Power = input$Power, alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, case = input$scenario_case)
      #temp$FRR <- (100)*temp$FRR
      return(temp)
    }

    if(3 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI_1 = input$MDRI_1, MDRI_2 = input$MDRI_2,
                                BigT = input$BigT,
                                FRR_1 = (1/100)*input$FRR_1, FRR_2 = (1/100)*input$FRR_2,
                                RSE_FRR_1 = (1/100)*input$RSE_FRR_1, RSE_FRR_2 = (1/100)*input$RSE_FRR_2,
                                RSE_MDRI_1 = (1/100)*input$RSE_MDRI_1, RSE_MDRI_2 = (1/100)*input$RSE_MDRI_2,
                                DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
                                I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
                                PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
                    Power = input$Power, alpha = input$alpha,
                    CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("simulate" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(
        MDRI = seq(input$MDRI_range_sim1[1], input$MDRI_range_sim1[2], by = 30),
        FRR = seq((1/100)*input$FRR_range_simul[1], (1/100)*input$FRR_range_simul[2], by = 0.005),

        BigT = input$BigT,
        RSE_FRR = (1/100)*input$RSE_FRR,
        RSE_MDRI = (1/100)*input$RSE_MDRI,
        DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
        DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
        I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
        PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
        Power = input$Power, alpha = input$alpha,
        CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
        ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("simulate_FRR" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(
        MDRI = seq(input$MDRI_range_simul[1], input$MDRI_range_simul[2], by = 30),
        FRR = seq((1/100)*input$FRR_range_sim2[1], (1/100)*input$FRR_range_sim2[2], by = 0.005),
        BigT = input$BigT,
        RSE_FRR = (1/100)*input$RSE_FRR,
        RSE_MDRI = (1/100)*input$RSE_MDRI,
        DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
        DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2,
        I_1 = (1/100)*input$I_1, I_2 = (1/100)*input$I_2,
        PrevH_1 = (1/100)*input$PrevH_1, PrevH_2 = (1/100)*input$PrevH_2),
        Power = input$Power, alpha = input$alpha,
        CR_1 = (1/100)*input$CR_1, CR_2 = (1/100)*input$CR_2,
        ss_calc, case = input$scenario_case)
      return(temp)
    }

  })
  output$plot <- renderPlot({
    validate(
      need(input$scenario_case <3, 'No plot available for this scenario (only one value available -see table)'),
      need(!(input$scenario_case == 2 & input$x_variable == "FRR"), 'This scenario cannot be represented by a plot (2 FRR values cannot live on one x-axis -see table)'),
      need(!(input$scenario_case != 1 & input$x_variable == "simulate"), 'Simulation is only implemented for case 1'),
      need(!(input$scenario_case != 1 & input$x_variable == "simulate_FRR"), 'Simulation is only implemented for case 1'),
      #need(!(input$RSE_FRR < 0 ), 'Please provide a covariance value for FRR'),
      need(!(input$RSE_FRR_1 < 0  & input$scenario_case >1), 'Please provide a value for FRR_1 covariance'),
      need(!(input$RSE_FRR_1 > 100  & input$scenario_case >1), 'Please provide a value for FRR_1 covariance'),
      need(!(input$RSE_FRR_2 < 0  & input$scenario_case >1), 'Please provide a value for FRR_2 covariance'),
      need(!(input$RSE_FRR_2 > 100  & input$scenario_case >1), 'Please provide a value for FRR_2 covariance'),
      need(!(input$RSE_MDRI_1 < 0  & input$scenario_case == 3), 'Please provide a value for MDRI_1 covariance'),
      need(!(input$RSE_MDRI_1 > 100  & input$scenario_case == 3), 'Please provide a value for MDRI_1 covariance'),
      need(!(input$RSE_MDRI_2 < 0  & input$scenario_case == 3), 'Please provide a value for MDRI_2 covariance'),
      need(!(input$RSE_MDRI_2 > 100  & input$scenario_case == 3), 'Please provide a value for MDRI_2 covariance'),
      need(input$RSE_FRR >= 0, 'Please provide a valid covariance value for FRR'),
      need(input$RSE_FRR <= 100, 'Please provide a valid covariance value for FRR'),
      need(input$RSE_MDRI, 'Please provide a covariance value for MDRI'),
      need(input$RSE_MDRI >= 0, 'Please provide a valid covariance value for MDRI'),
      need(input$RSE_MDRI <= 100, 'Please provide a valid covariance value for MDRI'),
      need(input$I_1, 'Please provide an incidence value for survey 1'),
      need(input$I_2, 'Please provide an incidence value for survey 2'),
      need(input$I_1 > 0, 'Please provide a valid incidence value for survey 1 (>0)'),
      need(input$I_2 > 0, 'Please provide a valid incidence value for survey 2 (>0)'),
      need(input$I_2 < input$I_1, 'Incidence in survey 1 must be >=  Incidence in survey 2'),
      need(input$BigT, 'Please provide a value for the cut-off BigT'),
      need(input$BigT > 120, 'Please provide a valid value for the cut-off BigT (>120)'),
      need(input$PrevH_1, 'Please provide an prevalence value for survey 1'),
      need(input$PrevH_2, 'Please provide an prevalence value for survey 2'),
      need(input$PrevH_1 > 0, 'Please provide a valid prevalence value for survey 1 (>0)'),
      need(input$PrevH_2 > 0, 'Please provide a valid prevalence value for survey 2 (>0)'),
      need(input$Power, 'Please provide a valid value for the Power required (0:1)'),
      need(input$Power > 0, 'Please provide a valid value for the Power required (0:1)'),
      need(input$Power <= 1, 'Please provide a valid value for the Power required (0:1)'),
      need(input$alpha, 'Please provide a valid value for alpha (0:1)'),
      need(input$alpha > 0, 'Please provide a valid value for alpha (0:1)'),
      need(input$alpha <= 1, 'Please provide a valid value for alpha (0:1)'),
      need(input$DE_H_1, 'Please provide a D.E. value for survey 1'),
      need(input$DE_R_1, 'Please provide a D.E. value for survey 1'),
      need(input$DE_H_2, 'Please provide a D.E. value for survey 2'),
      need(input$DE_R_2, 'Please provide a D.E. value for survey 2'),
      need(input$CR_1, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_2, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_1 > 0, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_1 <= 100, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_2 > 0, 'Please provide a valid value for the % of recently tested for positives'),
      need(input$CR_2 <= 100, 'Please provide a valid value for the % of recently tested for positives')
    )

    data <- format_table(df(), input$scenario_case, input$x_variable)
    plot_title <- plot_title(x_var = input$x_variable,
                             Power = input$Power, alpha = input$alpha,
                             I_1 = input$I_1, I_2 = input$I_2, BigT = input$BigT,
                             PrevH_1 = input$PrevH_1, PrevH_2 = input$PrevH_2,
                             FRR = input$FRR, FRR_1 = input$FRR_1, FRR_2 = input$FRR_2,
                             MDRI = input$MDRI, scenario_case = input$scenario_case)
    legend_text <- legendbox_text(scenario_case = input$scenario_case, x_var = input$x_variable,
                                  RSE_FRR = input$RSE_FRR , RSE_FRR_1 = input$RSE_FRR_1, RSE_FRR_2 = input$RSE_FRR_2,
                                  RSE_MDRI = input$RSE_MDRI, RSE_MDRI_1 = input$RSE_MDRI_1, RSE_MDRI_2 = input$RSE_MDRI_2,
                                  CR_1 = input$CR_1, CR_2 = input$CR_2,
                                  DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                  DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2)
    if (input$x_variable != "simulate" & input$x_variable != "simulate_FRR") {
      plot <- myplot(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)}
    if (input$x_variable == "simulate") {
      plot <- myplot_simulation(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)
    }
    if (input$x_variable == "simulate_FRR") {
      plot <- myplot_simulation_FRR(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)
    }
    print(plot)
    #return(plot)
  })

  output$tab <- renderTable({
    validate(
      need(!(input$scenario_case != 1 & input$x_variable == "simulate"), 'Simulation is only implemented for case 1'),
      need(!(input$scenario_case != 1 & input$x_variable == "simulate_FRR"), 'Simulation is only implemented for case 1'))
    data <- format_table(df(), input$scenario_case, input$x_variable)

    if("MDRI" %in% colnames(data)) data$MDRI <- format(data$MDRI, digits = 0)
    if("MDRI_1" %in% colnames(data)) data$MDRI_1 <- format(data$MDRI_1, digits = 0)
    if("MDRI_2" %in% colnames(data)) data$MDRI_2 <- format(data$MDRI_2, digits = 0)
    data$N <- format(data$N, digits = 10)
    if("FRR" %in% colnames(data)) data$FRR <- format(data$FRR, digits = 2)
    if("FRR_1" %in% colnames(data)) data$FRR_1 <- format(data$FRR_1, digits = 2)
    if("FRR_2" %in% colnames(data)) data$FRR_2 <- format(data$FRR_2, digits = 2)
    if("RSE_MDRI" %in% colnames(data)) data$RSE_MDRI <- format(data$RSE_MDRI, digits = 4)
    if("RSE_MDRI_1" %in% colnames(data)) data$RSE_MDRI_1 <- format(data$RSE_MDRI_1, digits = 4)
    if("RSE_MDRI_2" %in% colnames(data)) data$RSE_MDRI_2 <- format(data$RSE_MDRI_2, digits = 4)
    if("RSE_FRR" %in% colnames(data)) data$RSE_FRR <- format(data$RSE_FRR, digits = 4)
    if("RSE_FRR_1" %in% colnames(data)) data$RSE_FRR_1 <- format(data$RSE_FRR_1, digits = 4)
    if("RSE_FRR_2" %in% colnames(data)) data$RSE_FRR_2 <- format(data$RSE_FRR_2, digits = 4)
    data$DE_H_1 <- format(data$DE_H_1, digits = 2)
    data$DE_H_2 <- format(data$DE_H_2, digits = 2)
    data$DE_R_1 <- format(data$DE_R_1, digits = 2)
    data$DE_R_2 <- format(data$DE_R_2, digits = 2)
    data[, "I_1"] <- format(data[, "I_1"], digits = 2)
    data[, "I_2"] <- format(data[, "I_2"], digits = 2)
    data[, "PrevH_1"] <- format(data[, "PrevH_1"], digits = 3)
    data[, "PrevH_2"] <- format(data[, "PrevH_2"], digits = 3)
    data$BigT <- format(data$BigT, digits = 4)
    data.frame(data)
  #}, digits = 1)
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$x_variable, '.csv', sep='')
    },
    content = function(file) {
      write.csv(format_table(df(), input$scenario_case, input$x_variable), file)
    }
  )

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste(input$x_variable, '.pdf', sep='')
    },
    content = function(file) {
      ##begin
      data <- format_table(df(), input$scenario_case, input$x_variable)
      plot_title <- plot_title(x_var = input$x_variable,
                               Power = input$Power, alpha = input$alpha,
                               I_1 = input$I_1, I_2 = input$I_2, BigT = input$BigT,
                               PrevH_1 = input$PrevH_1, PrevH_2 = input$PrevH_2,
                               FRR = input$FRR, FRR_1 = input$FRR_1, FRR_2 = input$FRR_2,
                               MDRI = input$MDRI, scenario_case = input$scenario_case)
      legend_text <- legendbox_text(scenario_case = input$scenario_case, x_var = input$x_variable,
                                    RSE_FRR = input$RSE_FRR , RSE_FRR_1 = input$RSE_FRR_1, RSE_FRR_2 = input$RSE_FRR_2,
                                    RSE_MDRI = input$RSE_MDRI, RSE_MDRI_1 = input$RSE_MDRI_1, RSE_MDRI_2 = input$RSE_MDRI_2,
                                    CR_1 = input$CR_1, CR_2 = input$CR_2,
                                    DE_H_1 = input$DE_H_1, DE_H_2 = input$DE_H_2,
                                    DE_R_1 = input$DE_R_1, DE_R_2 = input$DE_R_2)
      if (input$x_variable != "simulate" & input$x_variable != "simulate_FRR") {
        plot <- myplot(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)}
      if (input$x_variable == "simulate") {
        plot <- myplot_simulation(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)
      }
      if (input$x_variable == "simulate_FRR") {
        plot <- myplot_simulation_FRR(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)
      }
      ggsave(file = file, plot = plot, width = 13, height = 9)
    }
  )

})
