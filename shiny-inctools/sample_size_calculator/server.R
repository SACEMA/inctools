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
      temp <- mdply(expand.grid(MDRI = seq(input$MDRI_range[1], input$MDRI_range[2], by = 30), frrhat = (1/100)*input$frrhat,
                                TIME = input$TIME,
                                frrhatcov = (1/100)*input$frrhatcov,
                                mdrihatcov = (1/100)*input$mdrihatcov,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    power = input$power, alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
                        
      return(temp)
    }

    if("MDRI" == input$x_variable & 2 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = seq(input$MDRI_range[1], input$MDRI_range[2], by = 30),
                                frrhat_1 = (1/100)*input$frrhat_1, frrhat_2 = (1/100)*input$frrhat_2,
                                TIME = input$TIME,
                                frrhatcov_1 = (1/100)*input$frrhatcov_1, frrhatcov_2 = (1/100)*input$frrhatcov_2,
                                mdrihatcov = (1/100)*input$mdrihatcov,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    power = input$power, alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("frrhat" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = input$MDRI,
                                frrhat = seq((1/100)*input$FRR_range[1], (1/100)*input$FRR_range[2], by = 0.005),
                                TIME = input$TIME,
                                frrhatcov = (1/100)*input$frrhatcov, mdrihatcov = (1/100)*input$mdrihatcov,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    power = input$power, alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("frrhat" == input$x_variable & 2 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI = input$MDRI,
                                TIME = input$TIME,
                                frrhat_1 = (1/100)*input$frrhat_1, frrhat_2 = (1/100)*input$frrhat_2,
                                frrhatcov_1 = (1/100)*input$frrhatcov_1, frrhatcov_2 = (1/100)*input$frrhatcov_2,
                                mdrihatcov = (1/100)*input$mdrihatcov,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    power = input$power, alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
      #temp$frrhat <- (100)*temp$frrhat
      return(temp)
    }

    if(3 == input$scenario_case) {
      temp <- mdply(expand.grid(MDRI_1 = input$MDRI_1, MDRI_2 = input$MDRI_2,
                                TIME = input$TIME,
                                frrhat_1 = (1/100)*input$frrhat_1, frrhat_2 = (1/100)*input$frrhat_2,
                                frrhatcov_1 = (1/100)*input$frrhatcov_1, frrhatcov_2 = (1/100)*input$frrhatcov_2,
                                mdrihatcov_1 = (1/100)*input$mdrihatcov_1, mdrihatcov_2 = (1/100)*input$mdrihatcov_2,
                                DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
                                inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
                                p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
                    power = input$power, alpha = input$alpha,
                    rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
                    ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("simulate" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(
        #MDRI = seq(120, 720, by = 120),S
        MDRI = seq(input$MDRI_range_sim1[1], input$MDRI_range_sim1[2], by = 30),
        frrhat = seq((1/100)*input$FRR_range_simul[1], (1/100)*input$FRR_range_simul[2], by = 0.005),

        TIME = input$TIME,
        frrhatcov = (1/100)*input$frrhatcov,
        # frrhatcov_1 = (1/100)*input$frrhatcov_1, frrhatcov_2 = (1/100)*input$frrhatcov_2,
        mdrihatcov = (1/100)*input$mdrihatcov,
        DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
        DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
        inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
        p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
        power = input$power, alpha = input$alpha,
        rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
        ss_calc, case = input$scenario_case)
      return(temp)
    }

    if("simulate_FRR" == input$x_variable & 1 == input$scenario_case) {
      temp <- mdply(expand.grid(
        MDRI = seq(input$MDRI_range_simul[1], input$MDRI_range_simul[2], by = 30),
        frrhat = seq((1/100)*input$FRR_range_sim2[1], (1/100)*input$FRR_range_sim2[2], by = 0.005),
        TIME = input$TIME,
        frrhatcov = (1/100)*input$frrhatcov,
        mdrihatcov = (1/100)*input$mdrihatcov,
        DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
        DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2,
        inc_1 = (1/100)*input$inc_1, inc_2 = (1/100)*input$inc_2,
        p_pos_1 = (1/100)*input$p_pos_1, p_pos_2 = (1/100)*input$p_pos_2),
        power = input$power, alpha = input$alpha,
        rec_test_coverage_1 = (1/100)*input$rec_test_coverage_1, rec_test_coverage_2 = (1/100)*input$rec_test_coverage_2,
        ss_calc, case = input$scenario_case)
      return(temp)
    }

  })
  output$plot <- renderPlot({
    #Input checking to avoid the application crashes due to lack of proper input
    #logic_case2 <- input$scenario_case == 2 & input$x_variable == "frrhat"
    validate(
      need(input$scenario_case <3, 'No plot available for this scenario (only one value available -see table)'),
      need(!(input$scenario_case == 2 & input$x_variable == "frrhat"), 'This scenario cannot be represented by a plot (2 FRR values cannot live on one x-axis -see table)'),
      need(!(input$scenario_case != 1 & input$x_variable == "simulate"), 'Simulation is only implemented for case 1'),
      need(!(input$scenario_case != 1 & input$x_variable == "simulate_FRR"), 'Simulation is only implemented for case 1'),
      #need(!(input$frrhatcov < 0 ), 'Please provide a covariance value for FRR'),
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
      need(input$power, 'Please provide a valid value for the power required (0:1)'),
      need(input$power > 0, 'Please provide a valid value for the power required (0:1)'),
      need(input$power <= 1, 'Please provide a valid value for the power required (0:1)'),
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

    data <- format_table(df(), input$scenario_case, input$x_variable)
    plot_title <- plot_title(x_var = input$x_variable,
                             power = input$power, alpha = input$alpha,
                             inc_1 = input$inc_1, inc_2 = input$inc_2, TIME = input$TIME,
                             p_pos_1 = input$p_pos_1, p_pos_2 = input$p_pos_2,
                             FRR = input$frrhat, FRR_1 = input$frrhat_1, FRR_2 = input$frrhat_2,
                             MDRI = input$MDRI, scenario_case = input$scenario_case)
    legend_text <- legendbox_text(scenario_case = input$scenario_case, x_var = input$x_variable,
                                  frrhatcov = input$frrhatcov , frrhatcov_1 = input$frrhatcov_1, frrhatcov_2 = input$frrhatcov_2,
                                  mdrihatcov = input$mdrihatcov, mdrihatcov_1 = input$mdrihatcov_1, mdrihatcov_2 = input$mdrihatcov_2,
                                  rec_test_coverage_1 = input$rec_test_coverage_1, rec_test_coverage_2 = input$rec_test_coverage_2,
                                  DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                  DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2)
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
    if("mdrihatcov" %in% colnames(data)) data$mdrihatcov <- format(data$mdrihatcov, digits = 4)
    if("mdrihatcov_1" %in% colnames(data)) data$mdrihatcov_1 <- format(data$mdrihatcov_1, digits = 4)
    if("mdrihatcov_2" %in% colnames(data)) data$mdrihatcov_2 <- format(data$mdrihatcov_2, digits = 4)
    if("frrhatcov" %in% colnames(data)) data$frrhatcov <- format(data$frrhatcov, digits = 4)
    if("frrhatcov_1" %in% colnames(data)) data$frrhatcov_1 <- format(data$frrhatcov_1, digits = 4)
    if("frrhatcov_2" %in% colnames(data)) data$frrhatcov_2 <- format(data$frrhatcov_2, digits = 4)
    data$DE_prev_1 <- format(data$DE_prev_1, digits = 2)
    data$DE_prev_2 <- format(data$DE_prev_2, digits = 2)
    data$DE_RgivenTested_1 <- format(data$DE_RgivenTested_1, digits = 2)
    data$DE_RgivenTested_2 <- format(data$DE_RgivenTested_2, digits = 2)
    data[, "Incidence_1 (%)"] <- format(data[, "Incidence_1 (%)"], digits = 2)
    data[, "Incidence_2 (%)"] <- format(data[, "Incidence_2 (%)"], digits = 2)
    data[, "Prevalence_1 (%)"] <- format(data[, "Prevalence_1 (%)"], digits = 3)
    data[, "Prevalence_2 (%)"] <- format(data[, "Prevalence_2 (%)"], digits = 3)
    data$TIME <- format(data$TIME, digits = 4)
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
                               power = input$power, alpha = input$alpha,
                               inc_1 = input$inc_1, inc_2 = input$inc_2, TIME = input$TIME,
                               p_pos_1 = input$p_pos_1, p_pos_2 = input$p_pos_2,
                               FRR = input$frrhat, FRR_1 = input$frrhat_1, FRR_2 = input$frrhat_2,
                               MDRI = input$MDRI, scenario_case = input$scenario_case)
      legend_text <- legendbox_text(scenario_case = input$scenario_case, x_var = input$x_variable,
                                    frrhatcov = input$frrhatcov , frrhatcov_1 = input$frrhatcov_1, frrhatcov_2 = input$frrhatcov_2,
                                    mdrihatcov = input$mdrihatcov, mdrihatcov_1 = input$mdrihatcov_1, mdrihatcov_2 = input$mdrihatcov_2,
                                    rec_test_coverage_1 = input$rec_test_coverage_1, rec_test_coverage_2 = input$rec_test_coverage_2,
                                    DE_prev_1 = input$DE_prev_1, DE_prev_2 = input$DE_prev_2,
                                    DE_RgivenTested_1 = input$DE_RgivenTested_1, DE_RgivenTested_2 = input$DE_RgivenTested_2)
      if (input$x_variable != "simulate" & input$x_variable != "simulate_FRR") {
        plot <- myplot(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)}
      if (input$x_variable == "simulate") {
        plot <- myplot_simulation(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)
      }
      if (input$x_variable == "simulate_FRR") {
        plot <- myplot_simulation_FRR(data, input$x_variable, input$checkbox_plotparams, title = plot_title, legend_text = legend_text)
      }
      ###end
      #write.csv(format_table(df(), input$scenario_case, input$x_variable), file)
      ggsave(file = file, plot = plot, width = 13, height = 9)
    }
  )

})
