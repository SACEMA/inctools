# Created by and Copyright (C) 2015-2016 Stefano Ongarello (FIND).
# recoded by lamin Juwara (mcgill) 2017
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
  if("FRR" == x_var) {
    x_breaks <- seq(0, 15, 1)
    x_lab <- "FRR (%)"
    x_variable <- "FRR"}
  gp <- ggplot(df, aes_string(y = "N", x = x_variable)) +
    geom_point(size = 5) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), alpha = 0, size = 0.6, colour = "darkblue") +
    scale_y_continuous(limits = c(0, max(df$V1))) +
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
    if("FRR" == x_var) x_grob_pos <- 0.05
    my_grob <- grobTree(textGrob(legend_text, x = x_grob_pos,  y = 0.7, hjust = 0,
                                 gp = gpar(col = "darkblue", fontsize = 12)))
    gp <- gp + annotation_custom(my_grob)
  }
  gp
}


myplot_simulation <- function(df, x_var = "MDRI", checkbox_plotparams = FALSE, title = "XXX", legend_text = "XXX") {
  x_breaks <- c(120, 180, 240, 300, 360, 420, 480, 540, 600, 660, 720, 780, 840)
  x_lab <- "MDRI (days)"
  x_variable <- "MDRI"
  var_variable <- "FRR"
  df$FRR <- as.factor(df$FRR)
  #}
  gp <- ggplot(df, aes_string(y = "N", x = x_variable, colour = var_variable)) +
    geom_point(size = 5)  +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), alpha = 0, size = 1) +
    scale_y_continuous(limits = c(0, max(df$N))) +
    xlab(x_lab) + ylab("Number of subjects required") +
    scale_x_continuous(breaks = x_breaks) +
    ggtitle(title) +
    theme_linedraw() + theme(
      panel.grid.major = element_line(colour = "darkgrey"),
      panel.grid.minor = element_line(colour = "darkgrey")
    ) +
    scale_y_continuous(labels = comma) +
    theme(legend.text = element_text(size = 20, family = "URWHelvetica")) +
    theme(axis.title = element_text(size = 18, family = "URWHelvetica")) +
    theme(text = element_text(size = 20, family = "URWHelvetica"))
  if(checkbox_plotparams) {
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
    x_grob_pos <- 0.05
    my_grob <- grobTree(textGrob(legend_text, x = x_grob_pos,  y = 0.7, hjust = 0,
                                 gp = gpar(col = "darkblue", fontsize = 12)))
    gp <- gp + annotation_custom(my_grob)
  }
  gp
}



format_table <- function(df, case_scenario, x_var) {
  df$I_1 <- 100*df$I_1
  df$I_2 <- 100*df$I_2
  df$PrevH_1 <- 100*df$PrevH_1
  df$PrevH_2 <- 100*df$PrevH_2
  if("FRR" %in% colnames(df)){
    df$FRR <- 100*df$FRR}
  if("FRR_1" %in% colnames(df)){
    df$FRR_1 <- 100*df$FRR_1
  }
  if("FRR_2" %in% colnames(df)){
    df$FRR_2 <- 100*df$FRR_2
  }
  colnames(df)[which(colnames(df) == "FRR")] <- "FRR"
  colnames(df)[which(colnames(df) == "FRR_1")] <- "FRR_1"
  colnames(df)[which(colnames(df) == "FRR_2")] <- "FRR_2"
  colnames(df)[which(colnames(df) == "MDRI")] <- "MDRI"
  colnames(df)[which(colnames(df) == "MDRI_1")] <- "MDRI_1"
  colnames(df)[which(colnames(df) == "MDRI_2")] <- "MDRI_2"
  colnames(df)[which(colnames(df) == "V1")] <- "N"
  colnames(df)[which(colnames(df) == "I_1 (%)")] <- "I_1"
  colnames(df)[which(colnames(df) == "I_2 (%)")] <- "I_2"
  colnames(df)[which(colnames(df) == "PrevH_1 (%)")] <- "PrevH_2"
  colnames(df)[which(colnames(df) == "PrevH_2 (%)")] <- "PrevH_2"
  df$N <- as.numeric(df$N)
  df$N[which(df$N<0)] <- NA
  return(df)
}

#' Creates the title plot based on different cases
plot_title <- function(x_var = "MDRI", Power = 0.8, alpha = 0.1, I_1 = 5, I_2 = 2.5, BigT = 720,
                       PrevH_1 = 15, PrevH_2 = 15, FRR = 1, FRR_1 = 1.5, FRR_2 = 2, MDRI = 360, scenario_case = 1) {
  if("MDRI" == x_var & 1 == scenario_case) {
    title <- paste0("Required sample size \n",
                    " Power: ", Power,",",  " alpha: ", alpha, "%,", "  FRR  : ", FRR, "%, ",
                    " BigT: ", BigT, " days, \n",
                    " I_1: ", I_1, "%,  I_2: ", I_2, "%,  PrevH_1: ", PrevH_1, "%,  PrevH_2: ", PrevH_2, "% \n")}

  if("MDRI" == x_var & 2 == scenario_case) {
    title <- paste0("Required sample size \n",
                    "Power: ", Power,",",  " alpha: ", alpha,"%,",  " FRR_1: ", FRR_1, "%,  FRR_2: ", FRR_2,  "%, BigT: ", BigT, " days, \n",
                    "I_1: ", I_1, "%,  I_2: ", I_2, "%,  PrevH_1: ", PrevH_1, "%,  PrevH_2: ", PrevH_2, "%\n")}

  if("FRR" == x_var & 1 == scenario_case) {
    title <-  paste0("Required sample size \n",
                     "Power: ", Power,",",  " alpha: ", alpha,",",  " MDRI: ", MDRI, " days, BigT: ", BigT, " days, \n",
                     "I_1: ", I_1, "%,  I_2: ", I_2, "%,  PrevH_1: ", PrevH_1, "%,  PrevH_2: ", PrevH_2, "%\n")}

  if(any("simulate" == x_var, "simulate_FRR" == x_var) & 1 == scenario_case) {
    title <- paste0("Required sample size \n",
                    "Power: ", Power,",",  " alpha: ", alpha, "%, BigT: ", BigT, " days, \n",
                    "I_1: ", I_1, "%,  I_2: ", I_2, "%,  PrevH_1: ", PrevH_1, "%,  PrevH_2: ", PrevH_2, "%\n")}

  return(title)
}

#' Creates additional text to display on the plot
legendbox_text <- function(scenario_case = 1, x_var = "MDRI",
                           RSE_FRR = 0.2 , RSE_FRR_1 = 0.3, RSE_FRR_2 = 0.4,
                           RSE_MDRI = 0.05, RSE_MDRI_1 = 0.06, RSE_MDRI_2 = 0.07,
                           CR_1 = 100, CR_2 = 100,
                           DE_H_1 = 1, DE_H_2 = 1, DE_R_1 = 1, DE_R_2 = 1) {
  if(1 == scenario_case) {
    legend <- paste0("Relative standard error of FRR: ", RSE_FRR, "%",
                     "\nRelative standard error of MDRI:", RSE_MDRI, "%",
                     "\nRecency test coverage in survey 1: ", CR_1, "%\n",
                     "Recency test coverage in survey 2: ", CR_2, "%\n",
                     "DE_1 prevalence: ", DE_H_1,
                     "\nDE_2 prevalence: ", DE_H_2,
                     "\nDE_1 recent prevalence among positives: ", DE_R_1,
                     "\nDE_2 recent prevalence among positives: ", DE_R_2)
  }
  if(2 == scenario_case) {
    legend <- paste0("FRR_1 covariance: ", RSE_FRR_1, "% ",
                     "FRR_2 covariance: ", RSE_FRR_2, "%",
                     "\nRelative standard error of MDRI:", RSE_MDRI, "%",
                     "\nRecency test coverage in survey 1: ", CR_1, "%\n",
                     "Recency test coverage in survey 2: ", CR_2, "%\n",
                     "DE_1 prevalence: ", DE_H_1,
                     "\nDE_2 prevalence: ", DE_H_2,
                     "\nDE_1 recent prevalence among positives: ", DE_R_1,
                     "\nDE_2 recent prevalence among positives: ", DE_R_2)
  }
  return(legend)
}


