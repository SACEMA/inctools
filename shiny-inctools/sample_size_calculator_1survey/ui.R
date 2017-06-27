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

#ui.R
#sample size table
library(shiny)
#library(shinyURL)
# Define UI
shinyUI(fluidPage(
  titlePanel("HIV incidence assays sample size calculator: one survey"),
  img(src='FIND_logo.jpg', align = "right", height = "75px"),
  img(src='SACEMA_logo.jpg', align = "right", height = "75px"),
  #br(),
  p("This tool allows you to upload a file with incidence/prevalence and MDRI/FRR [days/%] data, and obtain a table with the corresponding
    sample size required to detect a specified incidence in a cross-sectional survey, and the 95% confidence interval associated with the
    estimate at a given value alpha."),
  p("The data must be in .csv format, in one of these two formats:"),
  p("5 columns named: inc, prev, MDRI, FRR, scenario (the relative order is not important)"),
  p("5 columns containing incidence, prevalence, MDRI, FRR and scenario name (the column name is not important but the order is)"),
  p("The output will be displayed below the main input table, and the user can download the table by clicking on the button below"),
  p("Incidence and Prevalence data must be entered in % (e.g. 1.3 corresponds to 1.3%), the MDRI in days, the FRR in %"),
  p("Additional parameters include the RSE for the MDRI and FRR estimates, the RSE for the incidence estimate, design effect corrections,
    the coverage for recent testing and T."),
  tags$a('Download example file here', href = 'Table_example_one_survey.csv'),
  sidebarLayout(
    sidebarPanel(

      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      tags$hr(),
  downloadButton('downloadData', 'Download output table'),
  #numericInput('percent_reduction', label = h3("Incidence % reduction"), 50, step = 5),
  numericInput('DE_prev', label = h3("Design effect prevalence"), 1.3, step = 0.1),
  numericInput('DE_RgivenTested', label = h3("Design effect prevalence recent"), 1.3, step = 0.1),
  numericInput('inc_cov',label = h3("RSE of the incidence estimate [%]"), 20, step = 0.5),
  numericInput('rec_cov',label = h3("Coverage of recent testing [%]"), 100, step = 0.5),
  numericInput('TIME', label = h3("T [days]"), 730, step = 0.5),
  numericInput('alpha', label = h3("Alpha"), 0.05, step = 0.01),
  numericInput('mdrihatcov', label = h3("MDRI RSE [%]"), 5, step = .5),
  numericInput('frrhatcov', label = h3("FRR RSE [%]"), 30, step = .5)
  ),

    mainPanel(
      p("Input table"),
      tabPanel('contents', tableOutput('contents'), value = 1),
      p("Results table"),
      tabPanel('tab', tableOutput('tab'), value = 2)
      #tableOutput('tab')
    )
  )
))

