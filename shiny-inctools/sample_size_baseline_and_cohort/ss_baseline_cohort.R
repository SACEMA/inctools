power <- function(X1, sigma1, X2, sigma2, alpha = 0.05) {
  DeltaTrue <- X2-X1
  SigmaDelta <- sqrt(sigma1^2 + sigma2^2)
  DeltaCrit <- qnorm(p = alpha/2, mean = 0, sd = SigmaDelta, lower.tail = DeltaTrue<0)
  power <- pnorm(DeltaCrit, mean = DeltaTrue, sd = SigmaDelta, lower.tail = DeltaTrue<0)
  return(power)
}

sigma1 <- function(n, Inc, Prev, CR, MDRI, RSE_MDRI, FRR, RSE_FRR, BigT, DE_H, DE_R) {
  RSE <- suppressWarnings(incprecision(RSE_I = "out", I = Inc, PrevH = Prev, CR = CR, MDRI = MDRI, 
                                       RSE_MDRI = RSE_MDRI, FRR = FRR, 
                                       RSE_FRR = RSE_FRR, BigT = BigT,
                                       DE_H = DE_H, DE_R = DE_R, n = n))$RSE_I[[1]]
  sigma1 <- RSE*Inc
  return(sigma1)
}

sigma2 <- function(n, Inc, Prev, FracIncRed, FUT, CohortCR, DE_C) {
  binomialP <- (1 - exp(-Inc*(1 - FracIncRed)*FUT))
  n_cohort <- (1 - Prev)*CohortCR*n
  sigma_SCs <- sqrt(binomialP*(1 - binomialP)*n_cohort)
  RSE <- sigma_SCs/(binomialP*n_cohort)
  sigma2 <- RSE*Inc*(1 - FracIncRed)*sqrt(DE_C)
  return(sigma2)
}

powerdif <- function(n, Inc, Prev, CR, MDRI, RSE_MDRI, FRR, RSE_FRR, 
                     DE_H, DE_R, BigT, FracIncRed, FUT, CohortCR, Power, alpha, DE_C) {
  sigma1 <- sigma1(n = n, Inc = Inc, Prev = Prev, CR = CR, 
                   MDRI = MDRI, RSE_MDRI = RSE_MDRI, 
                   FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT,
                   DE_H = DE_H, DE_R = DE_R)
  sigma2 <- sigma2(n = n, Inc = Inc, Prev = Prev, FracIncRed = FracIncRed, 
                   FUT = FUT, CohortCR = CohortCR, DE_C = DE_C)
  PowerN <- power(X1 = Inc, sigma1 = sigma1, X2 = Inc*(1 - FracIncRed), sigma2 = sigma2, alpha = alpha) 
  return(Power - PowerN)
}

#generates a table with all sample sizes as function of different power values
do_table_data_ss_bl_cohort <- function(Inc, Prev, FracIncRed, 
                                       Power = 0.8, alpha = 0.05, 
                                       MDRI, RSE_MDRI, FRR, RSE_FRR, 
                                       CR = 1, DE_H = 1, DE_R = 1, DE_C = 1, 
                                       BigT = 730, CohortCR = 1, FUT = 1) {
  #browser()
  Ns <- seq(0, 50000, length = 250)
  powers <- vector(length = length(Ns))
  for (i in 1:length(Ns)) {
    N_tmp <- Ns[i]
    sigma1 <- sigma1(n = N_tmp, Inc = Inc, Prev = Prev, 
                     CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, 
                     FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT,
                     DE_H = DE_H, DE_R = DE_R)
    sigma2 <- sigma2(n = N_tmp, Inc = Inc, Prev = Prev, 
                     FracIncRed = FracIncRed, FUT = FUT, 
                     CohortCR = CohortCR, DE_C = DE_C)
    powers[i] <- power(X1 = Inc, sigma1 = sigma1, X2 = Inc*(1-FracIncRed), 
                       sigma2 = sigma2, alpha = alpha) }
  df <- data.frame(Ns, powers)
  return(df)
}

#calculates the sample size as function of power and creates the table for the plot
#and plots the result (should be split)
ss_baseline_cohort <- function(Inc, Prev, FracIncRed, Power = 0.8, alpha = 0.05, 
                               MDRI, RSE_MDRI, FRR, RSE_FRR, 
                               CR = 1, DE_H = 1, DE_R = 1, DE_C = 1, 
                               BigT = 730, CohortCR = 1, FUT = 1, produce_plot = TRUE) {
  N <- tryCatch({
    stats::uniroot(powerdif, Inc = Inc, Prev = Prev, CR = CR, MDRI = MDRI, RSE_MDRI = RSE_MDRI, 
                   FRR = FRR, RSE_FRR = RSE_FRR, BigT = BigT, DE_H = DE_H, DE_R = DE_R, 
                   FracIncRed = FracIncRed, FUT = FUT, CohortCR = CohortCR, alpha = alpha, 
                   Power = Power, DE_C = DE_C, interval = c(1,1000000))$root
  },
  error = function(err) {
    print("Desired power level cannot be achieved!")
    return(NA)
  })
  if(produce_plot) {
    data <- do_table_data_ss_bl_cohort(Inc = Inc, Prev = Prev, FracIncRed = FracIncRed, 
                                       Power = Power, alpha = alpha, 
                                       MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR = FRR, RSE_FRR = RSE_FRR, 
                                       CR = CR, DE_H = DE_H, DE_R = DE_R, DE_C = DE_C, 
                                       BigT = BigT, CohortCR = CohortCR, FUT = FUT)
    plot <- do_power_plot(data, N = N, Power = Power)
  } else {plot <- NULL}
  return(list(RequiredN = round(N), plot = plot))
}

#wrapper to produce the incidence reduction plot
incred_plot <- function(data){
  plot <- do_inc_red_plot(data)
  return(plot)
}

#generates the data for the inc red plot
incred_data <- function(Inc, Prev, IncRedRange = c(0.25, 0.75), 
                        steps = 25, Power = 0.8, alpha = 0.05, 
                        MDRI, RSE_MDRI, FRR, RSE_FRR, CR = 1, DE_H = 1, DE_R = 1, DE_C = 1, 
                        BigT = 730, CohortCR = 1, FUT = 1) {
 
  IncReds <- seq(IncRedRange[1], IncRedRange[2], length.out = steps)
  Ns <- vector(length = length(IncReds))
  for (i in 1:length(IncReds)) {
    Ns[i] <- ss_baseline_cohort(Inc = Inc,
                                Prev = Prev,
                                FracIncRed = IncReds[i],
                                Power = Power,
                                alpha = alpha,
                                MDRI = MDRI,
                                RSE_MDRI = RSE_MDRI,
                                FRR = FRR,
                                RSE_FRR = RSE_FRR,
                                CR=CR,
                                DE_H = DE_H,
                                DE_R = DE_R,
                                DE_C = DE_C,
                                BigT = BigT,
                                CohortCR = CohortCR,
                                FUT = FUT,
                                produce_plot = FALSE)$RequiredN
    #if (is.na(Ns[i])) {stop("Desired power cannot be achieved for all values in specified range.")}
  }
  IncReds <- IncReds*100
  data <- data.frame(IncReds, Ns)
  #plot <- do_inc_red_plot(data)
  return(data)
}


#creates the power plot for one case
do_power_plot <- function(data, N = 10000, Power = 0.8) {
  ggplot(data, aes(x = Ns, y = powers)) + 
    geom_line(size = 1.5) + 
    xlab("Sample size") + 
    ylab("Power") +
    theme_bw() +
    theme(text = element_text(size = 20)) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(labels = comma) +
    geom_vline(xintercept = N, colour = "red") +
    geom_hline(yintercept = Power, colour = "red")
}

#creates the incidence reduction plot for one case
do_inc_red_plot <- function(data) {
  ggplot(data, aes(y = Ns, x = IncReds)) + 
    geom_line(size = 1.5) + 
    xlab("Incidence reduction (%)") + 
    ylab("Sample size") +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    scale_y_continuous(labels = comma)
}
# 
# r1 <- ss_baseline_cohort(Inc = 0.05, Prev = 0.2,
#                    FracIncRed = 0.5, Power = 0.8, alpha = 0.05,
#                    MDRI = 200, RSE_MDRI = 0.05,
#                    FRR = 0.01, RSE_FRR = 0.05,
#                    CR = 1, DE_H = 1, DE_R = 1, DE_C = 1,
#                   BigT = 730, CohortCR = 1, FUT = 1, produce_plot = TRUE)
# 

# do_table_data_ss_bl_cohort(Inc = 0.05, Prev = 0.2,
#                            FracIncRed = 0.5, Power = 0.8, alpha = 0.05,
#                            MDRI = 200, RSE_MDRI = 0.05,
#                            FRR = 0.01, RSE_FRR = 0.05,
#                            CR = 1, DE_H = 1, DE_R = 1, DE_C = 1,
#                            BigT = 730, CohortCR = 1, FUT = 1)

#creates a grid with parameters and then applies do_table_data_ss_bl_cohort() 
#to get all sample sizes and powers in a single table
do_big_grid_power_ss <- function(Inc = 0.05, Prev = 0.2,
                                 FracIncRed = seq(0.1, 0.8, 0.05), Power = 0.8, alpha = 0.05,
                                 MDRI = 200, RSE_MDRI = 0.05,
                                 FRR = 0.01, RSE_FRR = 0.05,
                                 CR = 1, DE_H = 1, DE_R = 1, DE_C = 1,
                                 BigT = 730, CohortCR = 1, FUT = 1) {
  
  grid <- expand.grid(Inc = Inc, Prev = Prev,
                      FracIncRed = FracIncRed, Power = Power, alpha = alpha,
                      MDRI = MDRI, RSE_MDRI = RSE_MDRI,
                      FRR = FRR, RSE_FRR = RSE_FRR,
                      CR = CR, DE_H = DE_H, DE_R = DE_R, DE_C = DE_C,
                      BigT = BigT, CohortCR = CohortCR, FUT = FUT)
  
  res <- grid %>% ddply(.(Inc, Prev, FracIncRed, Power, alpha, MDRI, RSE_MDRI, FRR, RSE_FRR, CR, 
                          DE_H, DE_R, DE_C, BigT, CohortCR, FUT), 
                        transform, z = do_table_data_ss_bl_cohort(Inc, Prev, FracIncRed, 
                                                                  Power, alpha, MDRI, 
                                                                  RSE_MDRI, FRR, RSE_FRR,
                                                                  CR, DE_H, DE_R, DE_C, BigT, 
                                                                  CohortCR, FUT))
  return(res)}
#do_power_plot_many_red(temp_res)

#plots the power vs sample size in the entire range of fractional increment
do_power_plot_many_red <- function(df) {
  df$inc_red <- as.character(100*df$FracIncRed)
  plot <- ggplot(df, aes(x = z.Ns, y = z.powers)) + 
    geom_line(size = 1.5, aes(colour = inc_red)) + 
    xlab("Sample size") + 
    ylab("Power") +
    theme_bw() +
    theme(text = element_text(size = 20)) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(labels = comma) +
    scale_colour_discrete(name = "Incidence\nreduction (%)")
  return(plot)
}

#wrapper to generate the grid and produce the plot for many inc red and power
do_inc_red_many <- function(Inc = 0.05, Prev = 0.2,
                            FracIncRed = seq(0.1, 0.8, 0.05), Power = 0.8, alpha = 0.05,
                            MDRI = 200, RSE_MDRI = 0.05,
                            FRR = 0.01, RSE_FRR = 0.05,
                            CR = 1, DE_H = 1, DE_R = 1, DE_C = 1,
                            BigT = 730, CohortCR = 1, FUT = 1) {
  temp <- do_big_grid_power_ss(Inc = Inc, Prev = Prev,
                               FracIncRed = FracIncRed, Power = Power, alpha = alpha,
                               MDRI = MDRI, RSE_MDRI = RSE_MDRI,
                               FRR = FRR, RSE_FRR = RSE_FRR,
                               CR = CR, DE_H = DE_H, DE_R = DE_R, DE_C = DE_C,
                               BigT = BigT, CohortCR = CohortCR, FUT = FUT)
  plot <- do_power_plot_many_red(temp)
  return(plot)
}

# temp <- do_big_grid_power_ss(Inc = 0.05, Prev = 0.2,
#                              FracIncRed = seq(0.1, 0.8, 0.05), Power = 0.8, alpha = 0.05,
#                              MDRI = 200, RSE_MDRI = 0.05,
#                              FRR = 0.01, RSE_FRR = 0.05,
#                              CR = 1, DE_H = 1, DE_R = 1, DE_C = 1,
#                              BigT = 730, CohortCR = 1, FUT = 1)
# do_power_plot_many_red(temp)
# head(temp)
# filter(temp, z.powers > 0.79 & z.powers < 0.81)


#creates a table with sample size and power closest to the specified 
#threshold based on the list of results provided
find_power_ss <- function(df, power_threshold = 0.8) {
  
  idx <- which(abs(df$z.powers - power_threshold) == min(abs(df$z.powers - power_threshold), na.rm = TRUE))
  res <- data.frame(N = ceiling(df$z.Ns[idx]), Power = df$z.powers[idx])
  return(res)
}

#applies find_power_ss to a whole data frame by incidence reduction
generate_table_power_ss_incred <- function(Inc = 0.05, Prev = 0.2,
                                           FracIncRed = seq(0.1, 0.8, 0.05), Power = 0.8, alpha = 0.05,
                                           MDRI = 200, RSE_MDRI = 0.05,
                                           FRR = 0.01, RSE_FRR = 0.05,
                                           CR = 1, DE_H = 1, DE_R = 1, DE_C = 1,
                                           BigT = 730, CohortCR = 1, FUT = 1, 
                                           power_threshold = 0.8) {
  
  df <- do_big_grid_power_ss(Inc = Inc, Prev = Prev,
                             FracIncRed = FracIncRed, Power = Power, alpha = alpha,
                             MDRI = MDRI, RSE_MDRI = RSE_MDRI,
                             FRR = FRR, RSE_FRR = RSE_FRR,
                             CR = CR, DE_H = DE_H, DE_R = DE_R, DE_C = DE_C,
                             BigT = BigT, CohortCR = CohortCR, FUT = FUT)
  
  res <- ddply(df, "FracIncRed", find_power_ss, power_threshold = power_threshold)
  return(res)
}