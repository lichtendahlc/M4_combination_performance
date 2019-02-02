library(devtools)
devtools::install_github("carlanetto/M4comp2018")
# install Rob's M4metalearning package
devtools::install_github("robjhyndman/M4metalearning")

library(M4comp2018)
library(M4metalearning)


# ----------------------------------------
# Base model forecasts into validation set
# ----------------------------------------

data(M4)
ts_validation <- temp_holdout(M4)
ts_validation[1]

# Save chunks of the time series into multiple Rdata files. Later we aggregate these files into one Rdata file.
# Iterate through the chunks. When chunk = 5000 and iterations = 20 (recommended), the loop below will produce forecasts from the base nine models into the validation set.

chunk <- 5000 
n_iterations <- 20

# Test one chunk
i = 1
ptm <- proc.time() 
M4_valid_forec <- calc_forecasts(ts_validation[(chunk*(i-1)+1):(chunk*i)], forec_methods(), n.cores = 40)
proc.time() - ptm

# Run the loop
for (i in 1:n_iterations) {
  ptm <- proc.time()
  M4_valid_forec <- calc_forecasts(ts_validation[(chunk*(i-1)+1):(chunk*i)], forec_methods(), n.cores=40)
  proc.time() - ptm
  save(M4_valid_forec, file = paste("M4_valid_baseModelForecasts_", i, ".RData", sep = ""))
  rm(M4_valid_forec)
}

# Read the files back in and bind them together.
baseForecasts_valid <- c()
for (i in 1:20) {
  load(paste('M4_valid_baseModelForecasts_',i,'.Rdata',sep=''))
  baseForecasts_valid <- c(baseForecasts_valid, M4_valid_forec)
}

# Write the validation set data out.
save(baseForecasts_valid, file = "M4_valid_forec_100k.RData")