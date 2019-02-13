# -------------------------
# Install and load packages
# -------------------------

library(devtools)
# https://github.com/r-lib/devtools/issues/1772#issuecomment-388639815
assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
find_rtools()

# Install pmontman's M4metalearning package
devtools::install_github("pmontman/tsfeatures")
devtools::install_github("robjhyndman/M4metalearning")
devtools::install_github("pmontman/customxgboost")
devtools::install_github("carlanetto/M4comp2018")
install.packages("https://github.com/pmontman/M4metaresults/releases/download/v0.0.0.9000/M4metaresults_0.0.0.9000.tar.gz", 
                 repos = NULL, type="source")
devtools::install_github("M4Competition/M4-methods")
# Source code for all teams: https://github.com/M4Competition/M4-methods.
devtools::install_github("robjhyndman/anomalous")

library(tsfeatures)  # Create features for xgboost.
library(tidyverse)
library(anomalous)
library(purrr)
library(M4metalearning)
library(xgboost)
library(M4metaresults)  # where pmontman's forecasts are
library(M4comp2018)  # where the time series data themselves are
library(forecast)


# -------------------------------
# Evaluation metric and benchmark
# -------------------------------

RAAE_cal <- function(insample, outsample, forecasts){
  benchmarks <- as.numeric(naive_2(insample, length(outsample)))
  RAAE_vec <- rep(NA, length(outsample))
  for (i in 1:length(outsample)){
    benchmark <- benchmarks[i]
    pred <- as.numeric(forecasts[i])
    real <- as.numeric(outsample[i])
    forecast_AE <- abs(real - pred)
    benchmark_AE <- abs(real - benchmark)
    temp_RAAE <- 2 * forecast_AE/(forecast_AE + benchmark_AE)
    RAAE_vec[i] <- ifelse(is.nan(temp_RAAE), 1, temp_RAAE)
  }
  return(RAAE_vec)
}

SeasonalityTest <- function(input, ppy){
  # Used to determine whether a time series is seasonal
  tcrit <- 1.645
  if (length(input)<3*ppy){
    test_seasonal <- FALSE
  }else{
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- ( abs(xacf[ppy]) > clim[ppy] )
    
    if (is.na(test_seasonal)==TRUE){ test_seasonal <- FALSE }
  }
  return(test_seasonal)
}

naive_2 <- function(input, fh){
  # Estimate seasonaly adjusted time series
  ppy <- frequency(input); ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input,ppy) }
  if (ST==T){
    Dec <- decompose(input,type="multiplicative")
    des_input <- input/Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  } else {
    des_input <- input; SIout <- rep(1, fh)
  }
  return(naive(des_input, h=fh)$mean*SIout)
}


# ------------------------------------------
# Load M4 time series data and rearrange it.
# ------------------------------------------

data(M4)  # attach M4 data.

M4[[1]]
no_of_ts <- 100000
M4_letters_str <- c()
for (i in 1:no_of_ts) {
  M4_letters_str[i] <- substr(M4[[i]]$st, 0, 1)
}

submission_M4[[1]]
submission_M4_letters_str <- c()
for (i in 1:no_of_ts) {
  submission_M4_letters_str[i] <- substr(submission_M4[[i]]$st, 0, 1)
}

unique(M4_letters_str)
# Original order: D, H, M, Q, W, Y
unique(submission_M4_letters_str)
# New order: Y, Q, M, W, D, H

yearly_M4 <- Filter(function(l) l$period == "Yearly", M4)
y_len <- length(yearly_M4)  # 23,000
quarterly_M4 <- Filter(function(l) l$period == "Quarterly", M4)
q_len <-length(quarterly_M4)  # 24,000
monthly_M4 <- Filter(function(l) l$period == "Monthly", M4)
m_len <-length(monthly_M4)  # 48,000
weekly_M4 <- Filter(function(l) l$period == "Weekly", M4)
w_len <-length(weekly_M4)  # 359
daily_M4 <- Filter(function(l) l$period == "Daily", M4)
d_len <-length(daily_M4)  # 4,227
hourly_M4 <- Filter(function(l) l$period == "Hourly", M4)
h_len <-length(hourly_M4)  # 414

rearranged_M4 <- append(yearly_M4, quarterly_M4)
rearranged_M4 <- append(rearranged_M4, monthly_M4)
rearranged_M4 <- append(rearranged_M4, weekly_M4)
rearranged_M4 <- append(rearranged_M4, daily_M4)
rearranged_M4 <- append(rearranged_M4, hourly_M4)


# -------------------------------------
# Data from the screeing in L&W (2019).
# -------------------------------------

# New order: Y, Q, M, W, D, H
series_types <- c('Y', 'Q', 'M', 'W', 'D', 'H')

y_len <- 23000
q_len <- 24000
m_len <- 48000
w_len <- 359
d_len <- 4227
h_len <- 414

Y_final_models <- c(1, 2, 3, 4, 6)
Y_trim <- 0.2
Q_final_models <- c(1, 2, 4, 7)
Q_trim <- 0.25
M_final_models <- c(1, 2, 4, 7)
M_trim <- 0.25
W_final_models <- c(1, 2, 4, 6, 7, 8)
W_trim <- 1/6
D_final_models <- 8
D_trim <- 0
H_final_models <- c(1, 3, 4, 5, 9)
H_trim <- 0.2


# -------------------------------------------------------
# Evaluate accuracies of combinations in the testing set.
# -------------------------------------------------------

# https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
# https://github.com/M4Competition/M4-methods/blob/master/245%20-%20pmontman/M4-Method-Description.pdf

no_of_ts <- length(rearranged_M4)
no_of_models <- 4
RAAE_df <- as.data.frame(matrix(NA, no_of_ts, no_of_models))
colnames(RAAE_df) <- c('pmontman', 'sm', 'tm', 'screened')

# This loop takes a while.
ptm <- proc.time()
for (i in 1:no_of_ts) {
  # pmontman
  RAAE_df$pmontman[i] <- mean(RAAE_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(submission_M4[[i]]$y_hat)))
  # simple mean
  simple_mean <- apply(submission_M4[[i]]$ff, 2, mean)
  RAAE_df$sm[i] <- mean(RAAE_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(simple_mean)))
  # trimmed mean (lightly trimmed)
  trimmed_mean <- apply(submission_M4[[i]]$ff, 2, mean, trim = 1/9)
  RAAE_df$tm[i] <- mean(RAAE_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(trimmed_mean)))
}
proc.time() - ptm

# Calculate the screened combination's errors. Takes about 3 minutes.
ptm <- proc.time()
for (j in 1:6) {
  this_series_type <- series_types[j]
  if (this_series_type  == 'Y') starting_indices <- (1):(y_len)
  if (this_series_type  == 'Q') starting_indices <- (y_len + 1):(y_len + q_len)
  if (this_series_type  == 'M') starting_indices <- (y_len + q_len + 1):(y_len + q_len + m_len)
  if (this_series_type == 'W') starting_indices <- (y_len + q_len + m_len + 1):(y_len + q_len + m_len + w_len)
  if (this_series_type == 'D') starting_indices <- (y_len + q_len + m_len + w_len + 1):(y_len + q_len + m_len + w_len + d_len)
  if (this_series_type == 'H') starting_indices <- (y_len + q_len + m_len + w_len + d_len + 1):(y_len + q_len + m_len + w_len + d_len + h_len)
  
  if (this_series_type  == 'Y') final_models <- Y_final_models 
  if (this_series_type  == 'Q') final_models <- Q_final_models  
  if (this_series_type  == 'M') final_models <- M_final_models  
  if (this_series_type == 'W') final_models <- W_final_models  
  if (this_series_type == 'D') final_models <- D_final_models  
  if (this_series_type == 'H') final_models <- H_final_models 
  
  if (this_series_type  == 'Y') trim = Y_trim
  if (this_series_type  == 'Q') trim = Q_trim
  if (this_series_type  == 'M') trim = M_trim
  if (this_series_type == 'W') trim = W_trim 
  if (this_series_type == 'D') trim = D_trim  
  if (this_series_type == 'H') trim = H_trim 
  
  for (i in starting_indices) {
    forecasts <- submission_M4[[i]]$ff
    if (length(final_models) == 1) screened_combination <- forecasts[final_models, ]
    else screened_combination <- apply(forecasts[final_models, ], 2, mean, trim = trim)
    RAAE_df$screened[i] <- mean(RAAE_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, screened_combination))
  }
}
proc.time() - ptm

write.csv(RAAE_df, 'RAAE_df.csv')


# --------------
# Accuracy risk.
# --------------

RAAE_df <- read.csv('RAAE_df.csv')

# MRAAE in the testing set.
colMeans(RAAE_df)
#  pmontman        sm        tm  screened 
# 0.9125542 0.9727307 0.9375568 0.9143438 

# SDRAAE in the testing set.
apply(RAAE_df, 2, sd)
#pmontman        sm        tm  screened 
#0.2559562 0.1971782 0.1847634 0.2422236 

# Plot the ecdf's for RAAE.
par(mfrow=c(1,1))
plot(sort(RAAE_df$pmontman), ppoints(RAAE_df$pmontman), xlim=c(0,3), xlab = 'Series-Level OWA', ylab = 'Empirical CDF', main = '', bty='n', type ='l', lwd=2)
lines(sort(RAAE_df$screened), ppoints(RAAE_df$screened), col = 'green', lwd=2, lty=2)
legend(0.75, 0.25, legend=c('Montero-Manso et al. Combination', 'Trimmed Mean of Screened Methods'), col=c('black', 'green'), 
       box.lty=0, lty=c(1,2), lwd=c(2,2))

# Plot histograms of RAAE.
hist(RAAE_df$pmontman, breaks=50, xlim=c(0,2), ylim=c(0,20000), xlab = 'RAAE', ylab = 'Frequency', main = '')
hist(RAAE_df$screened, breaks=50, xlim=c(0,2), ylim=c(0,20000), xlab = 'RAAE', ylab = 'Frequency', main = '')

summary(RAAE_df$pmontman)
summary(RAAE_df$screened)

