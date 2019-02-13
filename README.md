This repository contains code to reproduce the results in Lichtendahl and Winkler's (2019) paper "Why Do Some Combinations Perform Better Than Others?" This paper proposes a simple, robust combination of a subset of the nine methods used in the M4 competition’s best combination by Montero-Manso, Talagala, Hyndman, and Athanasopoulos (2018). The proposed combination performs almost as well as that combination and is easier to implement. We screened out methods with low accuracy or highly correlated errors and combined the remaining methods with a trimmed mean. We also investigated accuracy risk (the risk of a bad forecast), proposing a series-level accuracy measure. Our trimmed mean and the trimmed mean of all nine methods had less accuracy risk than the best combination in the M4 competition and the simple mean of the nine methods. 

Run the code in the BaseForecastsInValidationSet.R file first. This code will produce the nine base models's forecasts into the validation set. These forecasts are saved in an Rdata file for use with the code in ScreenedCombination.R. The nine base time series models are auto_arima, ets, nnetar, tbats, stlm_ar, rw_drift, thetaf, naive, and snaive from the forecast package in R (Hyndman 2018).

Run the code in ScreenedCombination.R. This code reproduce the results in Tables 1-4 and Figure 1 in Lichtendahl and Winkler (2019). Run the code in NewRelativeError.R. This code reproduces the results in Table 5 and Figures 2 and 3. 
