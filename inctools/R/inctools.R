# Incidence Estimation Tools Copyright (C) 2015-2016, DST/NRF Centre of
# Excellence in Epidemiological Modelling and Analysis (SACEMA) and individual
# contributors.  This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.  This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the GNU
# General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

#' inctools: a package for incidence estimation
#'
#' The inctools package uses methods described by Kassanjee, et al. 'A new general biomarker-based incidence estimator,' \emph{Epidemiology} (2012), to implement functions to calculate incidence and tests of incidence difference between two populations, as well as power and sample size constraints for differenct study scenarios.
#' abie also provides functions for calculation of mean duration of recent infection and false recency rates from assays for recent infection.
#'
#' @section inctools Functions:
#' inctools has functions \emph{frrcal} to estimate false recency rate; \emph{mdrical} to estimate mean duration or recent infection; \emph{incprops} and \emph{inccounts} to calculate estimates and confidence intervals for incidence and incidence difference, as well as other summary and inferential statistics related to the survey; \emph{sspower} to calculate sample size needed for a given power in a test of incidence difference, or vice versa; and \emph{ssprecision}, which gives sample size needed for a given precision in the incidence estimate.
#'
#' For a longer introduction, see the introductory vignette for this package. Use vignette('Introduction',package='inctools') to access the document.
#'
#' @docType package
#' @name inctools
NULL 
