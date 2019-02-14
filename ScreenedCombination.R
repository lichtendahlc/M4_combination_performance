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


# --------------------------------
# Evaluation metrics and benchmark
# --------------------------------

smape_cal <- function(outsample, forecasts){
  # Used to estimate sMAPE
  outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
  smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
  return(smape)
}

mase_cal <- function(insample, outsample, forecasts){
  # Used to estimate MASE
  frq <- frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  
  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
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


# ------------------------------------------
# Load validation set data and rearrange it.
# ------------------------------------------

# Load the base model's forecasts into the validation.
# This Rdata file is created by script in "BaseForecastsInValidationSet.R"

load('M4_valid_forec_100k.RData') 
baseForecasts_valid[[1]]

no_of_ts <- 100000
valid_M4_letters_str <- c()
for (i in 1:no_of_ts) {
  valid_M4_letters_str[i] <- substr(baseForecasts_valid[[i]]$st, 0, 1)
}
unique(valid_M4_letters_str)

yearly_valid_M4 <- Filter(function(l) l$period == "Yearly", baseForecasts_valid)
quarterly_valid_M4 <- Filter(function(l) l$period == "Quarterly", baseForecasts_valid)
monthly_valid_M4 <- Filter(function(l) l$period == "Monthly", baseForecasts_valid)
weekly_valid_M4 <- Filter(function(l) l$period == "Weekly", baseForecasts_valid)
daily_valid_M4 <- Filter(function(l) l$period == "Daily", baseForecasts_valid)
hourly_valid_M4 <- Filter(function(l) l$period == "Hourly", baseForecasts_valid)

rearranged_valid_M4 <- append(yearly_valid_M4, quarterly_valid_M4)
rearranged_valid_M4 <- append(rearranged_valid_M4, monthly_valid_M4)
rearranged_valid_M4 <- append(rearranged_valid_M4, weekly_valid_M4)
rearranged_valid_M4 <- append(rearranged_valid_M4, daily_valid_M4)
rearranged_valid_M4 <- append(rearranged_valid_M4, hourly_valid_M4)

rearranged_valid_M4_letters_str <- c()
for (i in 1:no_of_ts) {
  rearranged_valid_M4_letters_str[i] <- substr(rearranged_valid_M4[[i]]$st, 0, 1)
}
unique(rearranged_valid_M4_letters_str)

valid_M4_errors <- list()
valid_M4_errors[1] <- calc_errors(rearranged_valid_M4[1])
valid_M4_errors[[1]]
for (i in 2:no_of_ts) {
  valid_M4_errors <- append(valid_M4_errors, calc_errors(rearranged_valid_M4[i]))
}

valid_M4_errors_letters_str <- c()
for (i in 1:no_of_ts) {
  valid_M4_errors_letters_str[i] <- substr(valid_M4_errors[[i]]$st, 0, 1)
}
unique(valid_M4_errors_letters_str)


# -------------------------------------------
# Evaluate base models in the validation set.
# -------------------------------------------

valid_smape_df <- as.data.frame(matrix(NA, no_of_ts, 9))
colnames(valid_smape_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                              'theta', 'naive', 'snaive')
valid_mase_df <- as.data.frame(matrix(NA, no_of_ts, 9))
colnames(valid_mase_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                             'theta', 'naive', 'snaive')
ptm <- proc.time()
for (i in 1:no_of_ts) {
  for (j in 1:9) {
    valid_smape_df[i, j] <- mean(valid_M4_errors[[i]]$smape_err[j, ])
    valid_mase_df[i, j] <- mean(valid_M4_errors[[i]]$mase_err[j, ])
  }
}
proc.time() - ptm

ptm <- proc.time()
for (i in 1:no_of_ts) {
  naive_2_forecasts <- naive_2(rearranged_valid_M4[[i]]$x, rearranged_valid_M4[[i]]$h)
  valid_smape_df$naive_2[i] <- mean(smape_cal(rearranged_valid_M4[[i]]$xx, naive_2_forecasts))
  valid_mase_df$naive_2[i] <- mean(mase_cal(rearranged_valid_M4[[i]]$x, rearranged_valid_M4[[i]]$xx, naive_2_forecasts))
}
proc.time() - ptm

valid_long_range_percent_errors_df <- as.data.frame(matrix(NA, no_of_ts, 9))
colnames(valid_long_range_percent_errors_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                                                  'theta', 'naive', 'snaive')
ptm <- proc.time()
for (i in 1:no_of_ts) {
  h <- rearranged_valid_M4[[i]]$h
  for (j in 1:9) {
    valid_long_range_percent_errors_df[i, j] <- (as.numeric(rearranged_valid_M4[[i]]$xx[h] - rearranged_valid_M4[[i]]$ff[j, h]))/as.numeric(rearranged_valid_M4[[i]]$xx[h])
  }
}
proc.time() - ptm

write.csv(valid_smape_df, 'valid_smape.csv', row.names=FALSE)
write.csv(valid_mase_df, 'valid_mase.csv', row.names=FALSE)
write.csv(valid_long_range_percent_errors_df, 'valid_long_range_percent_errors_df.csv', row.names=FALSE)

class(valid_mase_df)


# ----------------------------------------------------
# Calculate expertise and diversity in validation set.
# ----------------------------------------------------

valid_smape_df <- read.csv('valid_smape.csv')
valid_mase_df <- read.csv('valid_mase.csv')
valid_long_range_percent_errors_df <- read.csv('valid_long_range_percent_errors_df.csv')

y_len <- 23000
q_len <- 24000
m_len <- 48000
w_len <- 359
d_len <- 4227
h_len <- 414

# New order: Y, Q, M, W, D, H
series_types <- c('Y', 'Q', 'M', 'W', 'D', 'H')

# Some series have problematic MASEs: either NaN or Inf. It's all the same problem, no change in the smaller training
# set (the larger training set minus the validation set).
which(is.na(valid_mase_df$ets))
# Three series have MASE  = NaN.
#  Y12146 Y21168 D2085
#   12146  21168  97444
which(is.infinite(valid_mase_df$naive_2))  
# One series has MASE = Inf.
#   Q5619
#   28619
bad_index <- 97444
rearranged_valid_M4[[bad_index]]$st
rearranged_valid_M4[[bad_index]]$x
rearranged_valid_M4[[bad_index]]$xx
(naive_2_forecasts <- naive_2(rearranged_valid_M4[[bad_index]]$x, rearranged_valid_M4[[bad_index]]$h))
mase_cal(rearranged_valid_M4[[bad_index]]$x, rearranged_valid_M4[[bad_index]]$xx, naive_2_forecasts)

# Remove series with problematic MASEs. 
bad_indices <- c(which(is.na(valid_mase_df)), which(is.infinite(valid_mase_df$naive_2)))

# Other series. These series produce zero sOWA, which is why we propose MsOWA, rather than a geometric mean of sOWA_i's.
valid_sOWA_df <- as.data.frame(matrix(NA, 7, 10))
colnames(valid_sOWA_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                            'theta', 'naive', 'snaive', 'naive2')
for (i in 1:100000) {
  valid_sOWA_df[i, ] = (valid_smape_df[i, ] / valid_smape_df$naive_2[i] + valid_mase_df[i, ] / valid_mase_df$naive_2[i]) / 2
}

which(valid_sOWA_df$ets == 0)
example_index <- 28662
rearranged_valid_M4[[example_index]]$st
rearranged_valid_M4[[example_index]]$x
rearranged_valid_M4[[example_index]]$xx
rearranged_M4[[example_index]]$x
rearranged_M4[[example_index]]$xx

# Run a loop to calculate: 
# (1) the OWA by base model and series type.
# (2) the correlations among the base model's long-range percentage errors.
valid_OWA_df <- as.data.frame(matrix(NA, 7, 10))
colnames(valid_OWA_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                            'theta', 'naive', 'snaive', 'naive2')
row.names(valid_OWA_df) <- c(series_types, 'Total')

cor_df <- as.data.frame(matrix(NA, 6 * 9, 9))
colnames(cor_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                      'theta', 'naive', 'snaive')
trim = 0.05
for (i in 1:6) {
  this_series_type <- series_types[i]
  if (this_series_type  == 'Y') starting_indices <- (1):(y_len)
  if (this_series_type  == 'Q') starting_indices <- (y_len + 1):(y_len + q_len)
  if (this_series_type  == 'M') starting_indices <- (y_len + q_len + 1):(y_len + q_len + m_len)
  if (this_series_type == 'W') starting_indices <- (y_len + q_len + m_len + 1):(y_len + q_len + m_len + w_len)
  if (this_series_type == 'D') starting_indices <- (y_len + q_len + m_len + w_len + 1):(y_len + q_len + m_len + w_len + d_len)
  if (this_series_type == 'H') starting_indices <- (y_len + q_len + m_len + w_len + d_len + 1):(y_len + q_len + m_len + w_len + d_len + h_len)
  
  indices <- setdiff(starting_indices, bad_indices) 
  # Calculate the OWA by base model and series type.
  valid_OWA_df[i, ] <- (apply(valid_smape_df[indices, ], 2, mean, trim=trim)/mean(valid_smape_df$naive_2[indices], trim=trim) 
                        + apply(valid_mase_df[indices, ], 2, mean, trim=trim)/mean(valid_mase_df$naive_2[indices], trim=trim))/2
  # Calculate the correlations among the base model's long-range percentage errors.
  cor_df[(9 * (i - 1) + 1):(9 * i), ] <- cor(na.omit(valid_long_range_percent_errors_df[starting_indices, ]))
}

# Calculate the overall OWA by base model across the series types.
starting_indices <- 1:100000
indices <- setdiff(starting_indices, bad_indices) 
valid_OWA_df[7, ] <- (apply(valid_smape_df[indices, ], 2, mean, trim=trim)/mean(valid_smape_df$naive_2[indices], trim=trim) 
                      + apply(valid_mase_df[indices, ], 2, mean, trim=trim)/mean(valid_mase_df$naive_2[indices], trim=trim))/2

t(valid_OWA_df)  # L&W paper's Table 1.
# Apply expertise screen. 
Y_expert_models <- c(1, 2, 3, 4, 6, 7, 8, 9)
Q_expert_models <- c(1, 2, 4, 7)
M_expert_models <- c(1, 2, 4, 7)
W_expert_models <- c(1, 2, 4, 6, 7, 8, 9)
D_expert_models <- c(8, 9)
H_expert_models <- c(1, 2, 3, 4, 5, 9)

# Apply diversity screen.
# New order: Y, Q, M, W, D, H
model_names <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 'theta', 'naive', 'snaive')

i = 1 # Y
Y_cor_df <- cor_df[(9 * (i - 1) + 1):(9 * i), ]
row.names(Y_cor_df) <- model_names
Y_cor_df[Y_expert_models, Y_expert_models]  # L&W paper's Table 2.
Y_final_models <- c(1, 2, 3, 4, 6)
Y_trim <- 0.2

i = 2 # Q
Q_cor_df <- cor_df[(9 * (i - 1) + 1):(9 * i), ]
row.names(Y_cor_df) <- model_names
Q_cor_df[Q_expert_models, Q_expert_models]
Q_final_models <- c(1, 2, 4, 7)
Q_trim <- 0.25

i = 3 # M
M_cor_df <- cor_df[(9 * (i - 1) + 1):(9 * i), ]
row.names(Y_cor_df) <- model_names
M_cor_df[M_expert_models, M_expert_models]
M_final_models <- c(1, 2, 4, 7)
M_trim <- 0.25

i = 4 # W
W_cor_df <- cor_df[(9 * (i - 1) + 1):(9 * i), ]
row.names(Y_cor_df) <- model_names
W_cor_df[W_expert_models, W_expert_models]
W_final_models <- c(1, 2, 4, 6, 7, 8)
W_trim <- 1/6

i = 5 # D
D_cor_df <- cor_df[(9 * (i - 1) + 1):(9 * i), ]
row.names(Y_cor_df) <- model_names
D_cor_df[D_expert_models, D_expert_models]
D_final_models <- 8
D_trim <- 0

i = 6 # H
H_cor_df <- cor_df[(9 * (i - 1) + 1):(9 * i), ]
row.names(Y_cor_df) <- model_names
H_cor_df[H_expert_models, H_expert_models]
H_final_models <- c(1, 3, 4, 5, 9)
H_trim <- 0.2


# -------------------------------------------------------
# Evaluate accuracies of combinations in the testing set.
# -------------------------------------------------------

# https://github.com/M4Competition/M4-methods/blob/master/Benchmarks%20and%20Evaluation.R
# https://github.com/M4Competition/M4-methods/blob/master/245%20-%20pmontman/M4-Method-Description.pdf

no_of_ts <- length(rearranged_M4)
no_of_models <- 14  # pmontman's 9 base models, naive_2, pmontman's combination, 
                    # simple mean, trimmed mean, and proposed screened combination.

smape_df <- as.data.frame(matrix(NA, no_of_ts, no_of_models))
colnames(smape_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                        'theta', 'naive', 'snaive', 'naive_2', 'pmontman', 'sm', 'tm', 'screened')
mase_df <- as.data.frame(matrix(NA, no_of_ts, no_of_models))
colnames(mase_df) <- c('arima', 'ets', 'nnetar', 'tbats', 'stlm', 'rw', 
                       'theta', 'naive', 'snaive', 'naive_2', 'pmontman', 'sm', 'tm', 'screened')

# This loop takes about an hour.
ptm <- proc.time()
for (i in 1:no_of_ts) {
  # Nine base models
  for (j in 1:9) {
    smape_df[i, j] <- mean(smape_cal(rearranged_M4[[i]]$xx, as.numeric(submission_M4[[i]]$ff[j, ])))
    mase_df[i,j] <- mean(mase_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(submission_M4[[i]]$ff[j, ])))
  }
  # naive2
  naive_2_forecasts <- naive_2(rearranged_M4[[i]]$x, rearranged_M4[[i]]$h)
  smape_df$naive_2[i] <- mean(smape_cal(rearranged_M4[[i]]$xx, naive_2_forecasts))
  mase_df$naive_2[i] <- mean(mase_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, naive_2_forecasts))
  # pmontman
  smape_df$pmontman[i] <- mean(smape_cal(rearranged_M4[[i]]$xx, as.numeric(submission_M4[[i]]$y_hat)))
  mase_df$pmontman[i] <- mean(mase_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(submission_M4[[i]]$y_hat)))
  # simple mean
  simple_mean <- apply(submission_M4[[i]]$ff, 2, mean)
  smape_df$sm[i] <- mean(smape_cal(rearranged_M4[[i]]$xx, as.numeric(simple_mean)))
  mase_df$sm[i] <- mean(mase_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(simple_mean)))
  # trimmed mean (lightly trimmed)
  trimmed_mean <- apply(submission_M4[[i]]$ff, 2, mean, trim = 1/9)
  smape_df$tm[i] <- mean(smape_cal(rearranged_M4[[i]]$xx, as.numeric(trimmed_mean)))
  mase_df$tm[i] <- mean(mase_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, as.numeric(trimmed_mean)))
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
    smape_df$screened[i] <- mean(smape_cal(rearranged_M4[[i]]$xx, screened_combination))
    mase_df$screened[i] <- mean(mase_cal(rearranged_M4[[i]]$x, rearranged_M4[[i]]$xx, screened_combination))
  }
}
proc.time() - ptm

write.csv(smape_df, 'smape.csv', row.names=FALSE)
write.csv(mase_df, 'mase.csv', row.names=FALSE)

smape_df <- read.csv('smape.csv')
mase_df <- read.csv('mase.csv')

# OWA in the testing set (L&W paper's Table 3, column 5).
(colMeans(smape_df)/mean(smape_df$naive_2) + colMeans(mase_df)/mean(mase_df$naive_2))/2

# OWA in the testing set, by series type (L&W paper's Table 3, columns 1-4).
this_series_type <- 'O'
if (this_series_type  == 'Y') starting_indices <- (1):(y_len)
if (this_series_type  == 'Q') starting_indices <- (y_len + 1):(y_len + q_len)
if (this_series_type  == 'M') starting_indices <- (y_len + q_len + 1):(y_len + q_len + m_len)
if (this_series_type == 'O') starting_indices <- (y_len + q_len + m_len + 1):(y_len + q_len + m_len + w_len + d_len + h_len)
(colMeans(smape_df[starting_indices, c(11, 14)])/mean(smape_df$naive_2[starting_indices]) + colMeans(mase_df[starting_indices, c(11, 14)])/mean(mase_df$naive_2[starting_indices]))/2


# --------------
# Accuracy risk.
# --------------

# Series-by-series OWA in the testing set (L&W paper's Table 4)
sOWA_df <- (smape_df/smape_df$naive_2 + mase_df/mase_df$naive_2)/2
colMeans(sOWA_df)
apply(sOWA_df, 2, sd)

# Create L&W paper's Figure 1.
par(mfrow=c(1,1))
plot(sort(sOWA_df$pmontman), ppoints(sOWA_df$pmontman), xlim=c(0,3), xlab = 'Series-Level OWA', ylab = 'Empirical CDF', main = '', bty='n', type ='l', lwd=2)
lines(sort(sOWA_df$screened), ppoints(sOWA_df$screened), col = 'green', lwd=2, lty=2)
lines(sort(sOWA_df$tm), ppoints(sOWA_df$tm), col = 'blue', lwd=2, lty=3)
lines(sort(sOWA_df$sm), ppoints(sOWA_df$sm), col = 'red', lwd=2, lty=4)
legend(0.75, 0.25, legend=c('Montero-Manso et al. Combination', 'Trimmed Mean of Screened Methods', 
                            'Trimmed Mean of All Nine Methods', 'Simple Mean of All Nine Methods'), col=c('black', 'green', 'blue', 'red'), 
       box.lty=0, lty=c(1,2,3, 4), lwd=c(2,2,2,2))

