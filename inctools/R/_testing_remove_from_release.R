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


### ss_precision


# Note to self: 'step' is number of steps between minimum I and maximum I in the
# calculation of a range of output. So supply a vector or max/min theoretical
# incidences, and the function gives a range of values (step number of values)
# for the output. This can be done for all input variables. DO WE WANT THIS
# OPTION??

# Only gives option of n or RSE_I to be 'out' or not. So only gives precision or
# n. Makes sense. Either give function n to get RSE_I or give it RSE_I to get n.

# This function covers what was made in spreadsheets
# ABIE_v3_Test_Performance_Calculator (which takes SS and other variables of
# test, hypothetical data, and gives precision), and what was done in sheet
# ABIE_v3_Sample_Size_Calculator, which gives SS for a given precision level.


# Examples of function use: default of spreadsheet
# ABIE_v3_Sample_Size_Calculator.xlsx

ssprecision(I = 0.015, RSE_I = 0.25, PrevH = 0.2, CR = 1, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = 1, DE_R = 1, n = "out", step = 5)

##################################################################################################################### 
ssprecision(I = 0.015, RSE_I = 0.25, PrevH = 0.2, CR = 1, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = c(530, 730), DE_H = 1, DE_R = 1, n = "out", 
    step = 5)


ssprecision(I = c(0.015, 0.02), RSE_I = 0.25, PrevH = 0.2, CR = 1, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = c(530, 730), DE_H = 1, DE_R = 1, n = "out", 
    step = 5)



# doesn't work when FRR goes above 3.7% for these values I comment out because
# the package will not build with error messages present ssprecision ( I =0.015,
# RSE_I =0.25, PrevH =0.20, CR =1, MDRI =200, RSE_MDRI =0.05, FRR =0.039, RSE_FRR
# =0.2, BigT = 730, DE_H = 1, DE_R = 1, n = 'out', step = 5)

##################################################################################################################### 

# default of spreadsheet ABIE_v3_Test_Performance_Calculator
ssprecision(I = 0.015, RSE_I = "out", PrevH = 0.2, CR = 0.7, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = 1, DE_R = 1, n = c(5000, 5500), 
    step = 5)


##################################################################################################################### 
ssprecision(I = 0.015, RSE_I = "out", PrevH = 0.2, CR = 0.7, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = 1, DE_R = 1, n = 3622, step = 5)



##################################################################################################################### 
ssprecision(I = 0.015, RSE_I = "out", PrevH = c(0.2, 0.22), CR = 1, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = c(1, 1.1), DE_R = 1, n = 3622, 
    step = 5)




ssprecision(I = 0.015, RSE_I = "out", PrevH = 0.2, CR = 0.7, MDRI = 200, RSE_MDRI = 0.05, 
    FRR = 0.01, RSE_FRR = 0.2, BigT = 730, DE_H = 1, DE_R = 1, n = 3622, step = 5)


### ss_power

##################### --- Test values against spreadsheets (Find Power)---- #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.2
PrevH2 <- 0.2
Power = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- 0.01
RSE_FRR <- 0.2
BigT <- 730
CR <- 1
DE_H <- 1
DE_R <- 1
n1 <- 5000
n2 <- 5000
alpha = 0.05
BMest = "same.test"

sspower(I1, I2, PrevH1, PrevH2, n1, n2, alpha = 0.05, Power = "out", SS = NULL, CR = 1, 
    DE_H = 1, DE_R = 1, BMest = "same.test", MDRI, RSE_MDRI, FRR, RSE_FRR, BigT = 730)


my.data <- sspower(I1, I2, PrevH1, PrevH2, n1, n2, alpha = 0.05, Power = "out", SS = NULL, 
    CR = 1, DE_H = 1, DE_R = 1, BMest = "same.test", MDRI, RSE_MDRI, FRR, RSE_FRR, 
    BigT = 730)
my.data$Inc.Difference.Statistics
my.data[[1]]
my.data[[1]][4]
my.data[[1]]$Power
names(my.data)
str(my.data)



##################### --- Test values against spreadsheets (Find Power)---- #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.1
PrevH2 <- 0.1
Power = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- c(0.01, 0.03)
RSE_FRR <- c(0.2, 0.21)
BigT <- 730
CR <- c(1, 0.9)
DE_H <- c(1, 1.2)
DE_R <- c(1.1, 1)
n1 <- 4000
n2 <- 4500
alpha = 0.05
BMest = "FRR.indep"

sspower(I1, I2, PrevH1, PrevH2, n1, n2, alpha = 0.05, Power = "out", SS = NULL, CR = 1, 
    DE_H = DE_H, DE_R = DE_H, BMest = "FRR.indep", MDRI, RSE_MDRI, FRR, RSE_FRR, 
    BigT = 730)





##################### --- Test values against spreadsheets (Find Power)---- #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.1
PrevH2 <- 0.1
Power = "out"
CR = c(1, 0.9)
MDRI <- c(200, 210)
RSE_MDRI <- c(0.05, 0.04)
FRR <- c(0.01, 0.03)
RSE_FRR <- c(0.2, 0.21)
BigT <- c(730, 720)
CR <- c(1, 0.9)
DE_H <- c(1, 1.2)
DE_R <- c(1.1, 1)
n1 <- 4000
n2 <- 4500
alpha = 0.05
BMest = "MDRI.FRR.indep"

sspower(I1 = I1, I2 = I2, PrevH1 = PrevH1, PrevH2 = PrevH2, n1 = n1, n2 = n2, alpha = alpha, 
    Power = "out", CR = CR, DE_H = DE_H, DE_R = DE_R, BMest = "MDRI.FRR.indep", MDRI = MDRI, 
    RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT)






##################### --- Test values against spreadsheets (Find SS)---- #######################

I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.2
PrevH2 <- 0.15
Power = 0.8
SS = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- 0.01
RSE_FRR <- 0.2
BigT <- 730
CR <- 1
DE_H <- 1
DE_R <- 1
n1 <- "both"
n2 <- "both"
alpha = 0.05
BMest = "same.test"



test.SS <- sspower(I1, I2, PrevH1, PrevH2, n1, n2, alpha = 0.05, Power = 0.8, SS = "out", 
    CR = 1, DE_H = DE_H, DE_R = DE_R, BMest = "same.test", MDRI, RSE_MDRI, FRR, RSE_FRR, 
    BigT = 730)

test.SS




##################### --- Test values against spreadsheets (Find SS)---- #######################


I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.2
PrevH2 <- 0.15
Power = 0.8
SS = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- c(0.01, 0.02)
RSE_FRR <- c(0.2, 0.19)
BigT <- 730
CR <- 1
DE_H <- c(1.1, 1.2)
DE_R <- 1
n1 <- "both"
n2 <- "both"
alpha = 0.05
BMest = "FRR.indep"

sspower(I1, I2, PrevH1, PrevH2, n1, n2, alpha = 0.05, Power = 0.8, SS = "out", CR = 1, 
    DE_H = DE_H, DE_R = DE_R, BMest = "FRR.indep", MDRI, RSE_MDRI, FRR, RSE_FRR, 
    BigT = 730)




##################### --- Test values against spreadsheets (Find SS)---- #######################


I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.2
PrevH2 <- 0.15
Power = 0.8
SS = "out"
MDRI <- 200
RSE_MDRI <- c(0.05, 0.04)
FRR <- c(0.01, 0.02)
RSE_FRR <- c(0.2, 0.19)
BigT <- c(730, 725)
CR <- 1
DE_H <- c(1.1, 1.2)
DE_R <- 1
n1 <- "both"
n2 <- "both"
alpha = 0.05
BMest = "MDRI.FRR.indep"

sspower(I1 = I1, I2 = I2, PrevH1 = PrevH1, PrevH2 = PrevH2, n1 = n1, n2 = n2, alpha = 0.05, 
    Power = 0.8, SS = "out", CR = CR, DE_H = DE_H, DE_R = DE_R, BMest = "MDRI.FRR.indep", 
    MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT)





##################### --- Test values against spreadsheets (Find SS--one value)----
##################### #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.2
PrevH2 <- 0.15
Power = 0.8
SS = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- c(0.01, 0.02)
RSE_FRR <- c(0.2, 0.19)
BigT <- 730
CR <- 1
DE_H <- 1
DE_R <- 1
n1 <- "out"
n2 <- 5000
alpha = 0.05
BMest = "FRR.indep"

sspower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.2, PrevH2 = 0.15, n1 = "out", n2 = 5000, 
    alpha = 0.05, Power = 0.8, SS = "out", CR = 1, DE_H = 1, DE_R = 1, BMest = "FRR.indep", 
    MDRI = 200, RSE_MDRI = 0.05, FRR = FRR, RSE_FRR = RSE_FRR, BigT = 730)




##################### --- Test values against spreadsheets (Find SS)---- #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.1
PrevH2 <- 0.1
Power = 0.8
SS = "out"
CR = c(1, 0.9)
MDRI <- c(200, 210)
RSE_MDRI <- c(0.05, 0.04)
FRR <- c(0.01, 0.03)
RSE_FRR <- c(0.2, 0.21)
BigT <- c(730, 720)
CR <- c(1, 0.9)
DE_H <- c(1, 1.2)
DE_R <- c(1.1, 1)
n1 <- "both"
n2 <- "both"
alpha = 0.05
BMest = "MDRI.FRR.indep"

sspower(I1 = I1, I2 = I2, PrevH1 = PrevH1, PrevH2 = PrevH1, n1 = n1, n2 = n2, alpha = 0.05, 
    Power = Power, SS = "out", CR = CR, DE_H = DE_H, DE_R = DE_R, BMest = "MDRI.FRR.indep", 
    MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT)




##################### --- Test values against spreadsheets (Find SS--one value)----
##################### #######################
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.1
PrevH2 <- 0.1
Power = 0.8
SS = "out"
CR = c(1, 0.9)
MDRI <- c(200, 210)
RSE_MDRI <- c(0.05, 0.04)
FRR <- c(0.01, 0.03)
RSE_FRR <- c(0.2, 0.21)
BigT <- c(730, 720)
CR <- c(1, 0.9)
DE_H <- c(1, 1.2)
DE_R <- c(1.1, 1)
n1 <- 4500
n2 <- "out"
alpha = 0.05
BMest = "MDRI.FRR.indep"

sspower(I1 = I1, I2 = I2, PrevH1 = PrevH1, PrevH2 = PrevH1, n1 = n1, n2 = n2, alpha = 0.05, 
    Power = Power, SS = "out", CR = CR, DE_H = DE_H, DE_R = DE_R, BMest = "MDRI.FRR.indep", 
    MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT)






##################### --- Test values against spreadsheets (Find Sample Size)----
##################### ####################### See how to break the thing
I1 <- 0.01
I2 <- 0.011
PrevH1 <- 0.2
PrevH2 <- 0.15
Power = 0.99
SS = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- 0.01
RSE_FRR <- 0.2
BigT <- 730
CR <- 1
DE_H <- 1
DE_R <- 1
n1 <- "both"
n2 <- "both"
alpha = 0.05
BMest = "same.test"



test.SS <- sspower(I1 = I1, I2 = I2, PrevH1 = PrevH1, PrevH2 = PrevH2, n1 = n1, n2 = n2, 
    alpha = 0.05, Power = Power, SS = "out", CR = CR, DE_H = DE_H, DE_R = DE_R, BMest = "same.test", 
    MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT)

test.SS




# another error message test: both SS and both n1 and n2 are labeled 'out'
I1 <- 0.05
I2 <- 0.03
PrevH1 <- 0.2
PrevH2 <- 0.15
Power = 0.8
SS = "out"
MDRI <- 200
RSE_MDRI <- 0.05
FRR <- c(0.01, 0.02)
RSE_FRR <- c(0.2, 0.19)
BigT <- 730
CR <- 1
DE_H <- 1
DE_R <- 1
n1 <- "out"
n2 <- "out"
alpha = 0.05
BMest = "FRR.indep"

sspower(I1 = 0.05, I2 = 0.03, PrevH1 = 0.2, PrevH2 = 0.15, n1 = "out", n2 = "out", 
    alpha = 0.05, Power = 0.8, SS = "out", CR = 1, DE_H = 1, DE_R = 1, BMest = "FRR.indep", 
    MDRI = 200, RSE_MDRI = 0.05, FRR = FRR, RSE_FRR = RSE_FRR, BigT = 730)




### incidence_estimation

# Begin Examples ----------------------------------------------------------
# ######## -- The rest of this code gives example to show the function agrees
# with spreadsheet answers, and to see when and how ######## the function breaks,
# what error messages it outputs. Final R package must omit this code
probs <- prevcounts(N = c(5000, 5000), N_H = c(1000, 1000), N_testR = c(1000, 970), 
    N_R = c(100, 70))
# ################### == Call - DM only
# ==######################################################
incprops(BS_Count = 10000, Boot = FALSE, BMest = "same.test", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730)
# # (single survey)
incprops(BS_Count = 10000, Boot = FALSE, BMest = "same.test", PrevH = probs[1, 1], 
    RSE_PrevH = probs[1, 3], PrevR = probs[1, 2], RSE_PrevR = probs[1, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730)



############################################################################################## 

################### == Call - BS only ==######################################################
incprops(BS_Count = 10000, Boot = TRUE, BMest = "same.test", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730)

# boot with a larger sample
incprops(BS_Count = 1e+05, Boot = TRUE, BMest = "same.test", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730, Covar_HR = 2e-05)

# boot with a larger sample and different testing scheme
incprops(BS_Count = 1e+05, Boot = TRUE, BMest = "FRR.indep", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = c(0.01, 0.008), RSE_FRR = c(0.2, 0.21), BigT = 730, Covar_HR = 2e-05)
############################################################################################## 



######## == Check with spread sheets ''''''same.test, three surveys''''''==###########
probs <- prevcounts(N = c(5000, 5000, 3000), N_H = c(1000, 1000, 1000), N_testR = c(1000, 
    1000, 960), N_R = c(100, 70, 120))

incprops(BS_Count = 10000, Boot = FALSE, BMest = "same.test", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730)

######## == Check with spread sheets ''''''FRR.indep''''''==###########
probs <- prevcounts(N = c(5000, 5000, 3000), N_H = c(1000, 1000, 1000), N_testR = c(1000, 
    1000, 955), N_R = c(100, 70, 120))

incprops(BS_Count = 10000, Boot = FALSE, BMest = "FRR.indep", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730)

######## == Check with spread sheets ''''''MDRI.FRR.indep''''''==###########
probs <- prevcounts(N = c(5000, 5000, 3000), N_H = c(1000, 1000, 1000), N_testR = c(1000, 
    1000, 980), N_R = c(100, 70, 120))

incprops(BS_Count = 10000, Boot = FALSE, BMest = "MDRI.FRR.indep", PrevH = probs[, 
    1], RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = c(0.01, 0.02, 0.001), RSE_FRR = 0.2, BigT = 730)


####### == bootstrap with warining from FRR=0 ==###########
probs <- prevcounts(N = c(5000, 5000, 3000), N_H = c(1000, 1000, 1000), N_testR = c(1000, 
    1000, 990), N_R = c(100, 70, 120))

incprops(BS_Count = 10000, Boot = TRUE, BMest = "MDRI.FRR.indep", PrevH = probs[, 
    1], RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, BigT = 730)

incprops(BS_Count = 10000, Boot = FALSE, BMest = "same.test", PrevH = probs[, 1], 
    RSE_PrevH = probs[, 3], PrevR = probs[, 2], RSE_PrevR = probs[, 4], MDRI = 200, 
    RSE_MDRI = 0.05, FRR = 0, RSE_FRR = 0, BigT = 730)


######## == incidence by counts, single survey ==###########
inccounts(N = 5000, N_H = 1000, N_testR = 1000, N_R = 100, DE_H = 1, DE_R = 1, BS_Count = 10000, 
    Boot = FALSE, BMest = "same.test", MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.2, 
    BigT = 730, Covar_HR = 0)

######## == incidence by counts, two surveys ==###########
inccounts(N = c(5000, 5000), N_H = c(1000, 1000), N_testR = c(1000, 1000), N_R = c(100, 
    70), DE_H = 1, DE_R = 1, BS_Count = 10000, Boot = FALSE, BMest = "same.test", 
    MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.06, BigT = 730, Covar_HR = 0)

inccounts(N = c(5000, 5000), N_H = c(1000, 1000), N_testR = c(1000, 1000), N_R = c(100, 
    70), DE_H = 1, DE_R = 1, BS_Count = 10000, Boot = TRUE, BMest = "same.test", 
    MDRI = 200, RSE_MDRI = 0.05, FRR = 0.01, RSE_FRR = 0.06, BigT = 730, Covar_HR = 0)



######## == incidence by counts, two surveys, FRR independent ==########### == incidence
######## by counts, two surveys, FRR independent ==###########
inccounts(N = c(5000, 5000), N_H = c(1000, 1000), N_testR = c(1000, 950), N_R = c(100, 
    70), DE_H = 1, DE_R = 1, BS_Count = 10000, Boot = FALSE, BMest = "MDRI.FRR.indep", 
    MDRI = c(200, 190), RSE_MDRI = c(0.05, 0.07), FRR = c(0.01, 0.02), RSE_FRR = 0.05, 
    BigT = c(730, 735), Covar_HR = 0)

inccounts(N = c(5000, 5000), N_H = c(1000, 1000), N_testR = c(1000, 950), N_R = c(100, 
    70), DE_H = 1, DE_R = 1, BS_Count = 10000, Boot = FALSE, BMest = "FRR.indep", 
    MDRI = 200, RSE_MDRI = 0.05, FRR = 0, RSE_FRR = 0.05, BigT = 730, Covar_HR = 0)



# End Examples ------------------------------------------------------------


## sspower

################################## Break ######################################### This section was an attempt to
################################## use the same machinery as recnecyI() wrt to first order terms. It was quite
################################## troublesome, and I wasn't able to get it to work, so I scrapped it and
################################## explicitly calculated the formulas

# DM_Var_PrevH<-PrevH*(1-PrevH)*DE_H/N DM_Var_PrevR<-PrevR*(1-PrevR)*DE_R/N


# # #next few lines make delta method matrix fot_Mat <- matrix (nrow=no_s,
# ncol=4) DM_Var_I <- vector(length=no_s) for (i in c(1:no_s)) { fot_Mat[i,] <-
# DM_FirstOrderTerms (prevH=PrevH[i], prevR=PrevR[i], mdri=MDRI[i], frr=FRR[i],
# bigt=BigT) #DM_FirstOrderTerms gives fot_prevH, fot_prevR, fot_mdri, fot_frr,
# so FOT for each prev. survey input #and each row of fot_Mat has the FOT for
# prev. survey i.  DM_Var_I[i] <- (fot_Mat[i,1]^2)*DM_Var_PrevH[i] +
# (fot_Mat[i,2]^2)*DM_Var_PrevR[i] + (fot_Mat[i,3]^2)*DM_Var_MDRI[i] +
# (fot_Mat[i,4]^2)*DM_Var_FRR[i] } #I tested this function explicitly. It works,
# and works in the other function scripts as well.

# # DM_Var_deltaI <- DM_VAR_deltaI (BMest=BMest, fot_prevH1=fot_Mat[1,1],
# fot_prevH2=fot_Mat[2,1], # fot_prevR1=fot_Mat[1,2], fot_prevR2=fot_Mat[2,2], #
# fot_mdri1=fot_Mat[1,3], fot_mdri2=fot_Mat[2,3], # fot_frr1=fot_Mat[1,4],
# fot_frr2=fot_Mat[2,4], # dm_var_prevH1=DM_Var_PrevH[1],
# dm_var_prevH2=DM_Var_PrevH[2], # dm_var_prevR1=DM_Var_PrevR[1],
# dm_var_prevR2=DM_Var_PrevR[2], # dm_var_mdri1=DM_Var_MDRI[1],
# dm_var_mdri2=DM_Var_MDRI[2], # dm_var_frr1=DM_Var_FRR[1],
# dm_var_frr2=DM_Var_FRR[2]) # # # DM_Var_deltaI.infSS <- DM_VAR_deltaI.infSS
# (BMest=BMest, # fot_mdri1=fot_Mat[1,3], fot_mdri2=fot_Mat[2,3], #
# fot_frr1=fot_Mat[1,4], fot_frr2=fot_Mat[2,4], # dm_var_mdri1=DM_Var_MDRI[1],
# dm_var_mdri2=DM_Var_MDRI[2], # dm_var_frr1=DM_Var_FRR[1],
# dm_var_frr2=DM_Var_FRR[2]) # # DM_SD_I <- sqrt(DM_Var_I) # DM_SD_deltaI <-
# sqrt(DM_Var_deltaI) # DM_SD_deltaI.infSS <- sqrt(DM_Var_deltaI.infSS) # Var_I
# <- DM_Var_I # RSE_I <- sqrt(Var_I)/I # SD_I <- sqrt(Var_I) # Var_deltaI <-
# DM_Var_deltaI #+ BS_Var_deltaI # RSE_deltaI <- sqrt(Var_deltaI)/abs(deltaI_Est)
# #sqrt(Var_deltaI)/abs(deltaI_Est_Vec) # RSE_deltaI.infSS <-
# sqrt(DM_Var_deltaI.infSS)/abs(deltaI_Est) # SD_deltaI <- sqrt(Var_deltaI) End
# Break ######################################### 
