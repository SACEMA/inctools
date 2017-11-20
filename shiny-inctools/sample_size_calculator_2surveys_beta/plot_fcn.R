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

#' Simple plotting function to produce the N vs x plot
#' At the moment only takes into account the cases with MDRI or FRR on x
#'
#' @param df input dataframe
#' @param x_var string with name of variable to put on x-axis
#' @return a ggplot object

myplot <- function(df, x_var = "MDRI", checkbox_plotparams = FALSE, title = "XXX", legend_text = "XXX") {
  if("MDRI" == x_var) {
    x_breaks <- c(120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840)
    x_lab <- "MDRI (days)"
    x_variable <- "MDRI"}
  if("frrhat" == x_var) {
    x_breaks <- seq(0, 15, 1)
    x_lab <- "FRR (%)"
    x_variable <- "FRR"}
  #   gp <- ggplot(df, aes_string(y = "V1", x = "MDRI", colour = "factor(frrhat)")) +
  gp <- ggplot(df, aes_string(y = "N", x = x_variable)) +
    geom_point(size = 5) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), alpha = 0, size = 0.6, colour = "darkblue") +
    #stat_smooth(method = "auto", alpha = 0.1, size = 0.6, colour = "darkblue") +
    scale_y_continuous(limits = c(0, max(df$V1))) +
    # geom_hline(yintercept = seq(10000, 50000, by = 20000)) +
    xlab(x_lab) + ylab("Number of subjects required") +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle(paste0("Required sample size")) +
    theme_linedraw() + theme(
      panel.grid.major = element_line(colour = "darkgrey"),
      panel.grid.minor = element_line(colour = "darkgrey")
    ) +
    scale_y_continuous(labels = comma) +
    theme(legend.text = element_text(size = 20, family = "URWHelvetica")) +
    theme(axis.title = element_text(size = 18, family = "URWHelvetica")) +
    theme(text = element_text(size = 20, family = "URWHelvetica")) +
    ggtitle(title)

  if(checkbox_plotparams) {
    if("MDRI" == x_var) x_grob_pos <- 0.65
    if("frrhat" == x_var) x_grob_pos <- 0.05
    my_grob <- grobTree(textGrob(legend_text, x = x_grob_pos,  y = 0.7, hjust = 0,
                                 gp = gpar(col = "darkblue", fontsize = 12)))
    gp <- gp + annotation_custom(my_grob)
  }
  gp
}


myplot_simulation <- function(df, x_var = "MDRI", checkbox_plotparams = FALSE, title = "XXX", legend_text = "XXX") {
  # if("MDRI" == x_var) {
  x_breaks <- c(120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840)
  x_lab <- "MDRI (days)"
  x_variable <- "MDRI"
  var_variable <- "FRR"
  df$FRR <- as.factor(df$FRR)
  #}
  gp <- ggplot(df, aes_string(y = "N", x = x_variable, colour = var_variable)) +
    geom_point(size = 5)  +
    #stat_smooth(method = "auto", alpha = 0.1, size = 1, aes_string(fill = var_variable)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), alpha = 0, size = 1) +
   # stat_smooth(method = "auto", alpha = 0.1, size = 1) +
    scale_y_continuous(limits = c(0, max(df$N))) +
    # geom_hline(yintercept = seq(10000, 50000, by = 20000)) +
    xlab(x_lab) + ylab("Number of subjects required") +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle(title) +
    #  scale_colour_brewer(palette="Paired") +
    #ggtitle(paste0("Required sample size")) +
    theme_linedraw() + theme(
      panel.grid.major = element_line(colour = "darkgrey"),
      panel.grid.minor = element_line(colour = "darkgrey")
    ) +
    scale_y_continuous(labels = comma) +
    theme(legend.text = element_text(size = 20, family = "URWHelvetica")) +
    theme(axis.title = element_text(size = 18, family = "URWHelvetica")) +
    theme(text = element_text(size = 20, family = "URWHelvetica"))
  if(checkbox_plotparams) {
    #  if("MDRI" == x_var) x_grob_pos <- 0.65
    #  if("frrhat" == x_var) x_grob_pos <- 0.05
    x_grob_pos <- 0.65
    my_grob <- grobTree(textGrob(legend_text, x = x_grob_pos,  y = 0.7, hjust = 0,
                                 gp = gpar(col = "darkblue", fontsize = 12)))
    gp <- gp + annotation_custom(my_grob)
  }
  gp
}


myplot_simulation_FRR <- function(df, x_var = "FRR", checkbox_plotparams = FALSE, title = "XXX", legend_text = "XXX") {
  # if("MDRI" == x_var) {
  x_breaks <- seq(0, 15, 1)
  x_lab <- "FRR (%)"
  x_variable <- "FRR"
  var_variable <- "MDRI"
  df$MDRI <- as.factor(df$MDRI)
  df$N <- as.numeric(df$N)
  #}
  gp <- ggplot(df, aes_string(y = "N", x = x_variable, colour = var_variable)) +
    geom_point(size = 5) +
    #stat_smooth(method = "auto", alpha = 0.1, size = 1, aes_string(fill = var_variable)) +
    stat_smooth(method = "auto", size = 1, alpha = 0) +
    scale_y_continuous(limits = c(0, max(df$N))) +
    # geom_hline(yintercept = seq(10000, 50000, by = 20000)) +
    xlab(x_lab) + ylab("Number of subjects required") +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle(title) +
    #  scale_colour_brewer(palette="Paired") +
    #ggtitle(paste0("Required sample size")) +
    theme_linedraw() + theme(
      panel.grid.major = element_line(colour = "darkgrey"),
      panel.grid.minor = element_line(colour = "darkgrey")
    ) +
    scale_y_continuous(labels = comma) +
    theme(legend.text = element_text(size = 20, family = "URWHelvetica")) +
    theme(axis.title = element_text(size = 18, family = "URWHelvetica")) +
    theme(text = element_text(size = 20, family = "URWHelvetica"))
  if(checkbox_plotparams) {
    #  if("MDRI" == x_var) x_grob_pos <- 0.65
    #  if("frrhat" == x_var) x_grob_pos <- 0.05
    x_grob_pos <- 0.05
    my_grob <- grobTree(textGrob(legend_text, x = x_grob_pos,  y = 0.7, hjust = 0,
                                 gp = gpar(col = "darkblue", fontsize = 12)))
    gp <- gp + annotation_custom(my_grob)
  }
  gp
}



format_table <- function(df, case_scenario, x_var) {
  df$inc_1 <- 100*df$inc_1
  df$inc_2 <- 100*df$inc_2
  df$p_pos_1 <- 100*df$p_pos_1
  df$p_pos_2 <- 100*df$p_pos_2
  if("frrhat" %in% colnames(df)){
    df$frrhat <- 100*df$frrhat}
  if("frrhat_1" %in% colnames(df)){
    df$frrhat_1 <- 100*df$frrhat_1
  }
  if("frrhat_2" %in% colnames(df)){
    df$frrhat_2 <- 100*df$frrhat_2
  }
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
  colnames(df)[which(colnames(df) == "p_pos_2")] <- "Prevalence_1 (%)"
  df$N <- as.numeric(df$N)
  df$N[which(df$N<0)] <- NA

  return(df)
}

#' Creates the title plot based on different cases
plot_title <- function(x_var = "MDRI", power = 0.8, alpha = 0.1, inc_1 = 5, inc_2 = 2.5, TIME = 720,
                       p_pos_1 = 15, p_pos_2 = 15, FRR = 1, FRR_1 = 1.5, FRR_2 = 2, MDRI = 360, scenario_case = 1) {
  if("MDRI" == x_var & 1 == scenario_case) {
    title <- paste0("Required sample size \n",
                    "Power: ", power, " Alpha: ", alpha, " FRR: ", FRR, "%, T: ", TIME, " days\n",
                    "Inc_1: ", inc_1, "% Inc_2: ", inc_2, "% Prev_1: ", p_pos_1, "% Prev_2: ", p_pos_2, "%\n")}

  if("MDRI" == x_var & 2 == scenario_case) {
    title <- paste0("Required sample size \n",
                    "Power: ", power, " Alpha: ", alpha, " FRR_1: ", FRR_1, "% FRR_2: ", FRR_2,  "%, T: ", TIME, " days\n",
                    "Inc_1: ", inc_1, "% Inc_2: ", inc_2, "% Prev_1: ", p_pos_1, "% Prev_2: ", p_pos_2, "%\n")}

  if("frrhat" == x_var & 1 == scenario_case) {
    title <-  paste0("Required sample size \n",
                     "Power: ", power, " Alpha: ", alpha, " MDRI: ", MDRI, " days, T: ", TIME, " days\n",
                     "Inc_1: ", inc_1, "% Inc_2: ", inc_2, "% Prev_1: ", p_pos_1, "% Prev_2: ", p_pos_2, "%\n")}

  if(any("simulate" == x_var, "simulate_FRR" == x_var) & 1 == scenario_case) {
    title <- paste0("Required sample size \n",
                    "Power: ", power, " Alpha: ", alpha, "%, T: ", TIME, " days\n",
                    "Inc_1: ", inc_1, "% Inc_2: ", inc_2, "% Prev_1: ", p_pos_1, "% Prev_2: ", p_pos_2, "%\n")}
#
#   if("simulate_FRR" == x_var & 1 == scenario_case) {
#     title <- paste0("Required sample size \n",
#                     "Power: ", power, " Alpha: ", alpha, "%, T: ", TIME, " days\n",
#                     "Inc_1: ", inc_1, "% Inc_2: ", inc_2, "% Prev_1: ", p_pos_1, "% Prev_2: ", p_pos_2, "%\n")}

  return(title)
}

#' Creates additional text to display on the plot
legendbox_text <- function(scenario_case = 1, x_var = "MDRI",
                           frrhatcov = 0.2 , frrhatcov_1 = 0.3, frrhatcov_2 = 0.4,
                           mdrihatcov = 0.05, mdrihatcov_1 = 0.06, mdrihatcov_2 = 0.07,
                           rec_test_coverage_1 = 100, rec_test_coverage_2 = 100,
                           DE_prev_1 = 1, DE_prev_2 = 1, DE_RgivenTested_1 = 1, DE_RgivenTested_2 = 1) {
  if(1 == scenario_case) {
    legend <- paste0("FRR covariance: ", frrhatcov, "%",
                     "\nMDRI covariance:", mdrihatcov, "%",
                     "\nRecency test coverage in survey 1: ", rec_test_coverage_1, "%\n",
                     "Recency test coverage in survey 2: ", rec_test_coverage_2, "%\n",
                     "DE_1 prevalence: ", DE_prev_1,
                     "\nDE_2 prevalence: ", DE_prev_2,
                     "\nDE_1 recent prevalence among positives: ", DE_RgivenTested_1,
                     "\nDE_2 recent prevalence among positives: ", DE_RgivenTested_2)
  }
  if(2 == scenario_case) {
    legend <- paste0("FRR_1 covariance: ", frrhatcov_1, "% ",
                     "FRR_2 covariance: ", frrhatcov_2, "%",
                     "\nMDRI covariance:", mdrihatcov, "%",
                     "\nRecency test coverage in survey 1: ", rec_test_coverage_1, "%\n",
                     "Recency test coverage in survey 2: ", rec_test_coverage_2, "%\n",
                     "DE_1 prevalence: ", DE_prev_1,
                     "\nDE_2 prevalence: ", DE_prev_2,
                     "\nDE_1 recent prevalence among positives: ", DE_RgivenTested_1,
                     "\nDE_2 recent prevalence among positives: ", DE_RgivenTested_2)
  }
  return(legend)
}

# plot_title(x_var = x_var,
#            power = 0.8, alpha = 0.1, inc_1 = 5, inc_2 = 2.5, TIME = 720,
#            p_pos_1 = 15, p_pos_2 = 15, FRR = 1, FRR_1 = 1.5, FRR_2 = 2, MDRI = 360, case_scenario = 1)
