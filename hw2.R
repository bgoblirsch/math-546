#################
##  Problem 1  ##
#################

# Set seed for reproducibility
set.seed(123)

# Parameters for the ARMA process
ar_params <- c(0.2, 0.08) # AR coefficients
ma_params <- c(0.7)       # MA coefficient

# Number of simulations
num_simulations <- 10000

# Time series lengths
time_series_lengths <- c(100, 1000, 10000)

# True autocorrelation at lag 3 for the given process
true_acf_lag3 <- ARMAacf(ar = ar_params, ma = ma_params, lag.max = 3)[4] 

# Initialize list to store ACF estimates
acf_lag3_list <- list()

# Perform simulations for each time series length
for (n in time_series_lengths) {
  acf_lag3 <- numeric(num_simulations)
  for (i in 1:num_simulations) {
    # Simulate the ARMA process
    sim_data <- arima.sim(model = list(ar = ar_params, ma = ma_params), n = n)
    # Compute the sample ACF at lag 3
    acf_vals <- acf(sim_data, lag.max = 3, plot = FALSE)$acf
    acf_lag3[i] <- acf_vals[4] # ACF at lag 3
  }
  acf_lag3_list[[as.character(n)]] <- acf_lag3
}

# Find the global range of the data across all time series lengths
global_min <- min(unlist(acf_lag3_list))
global_max <- max(unlist(acf_lag3_list))

# Define consistent breakpoints, x-axis, and y-axis ranges
breakpoints <- seq(global_min - 0.01, global_max + 0.01, by = 0.01)
x_limits <- c(global_min, global_max)
y_limits <- c(0, 3500) 

# Plot histograms horizontally
par(mfrow = c(1, 3))

for (n in time_series_lengths) {
  hist(acf_lag3_list[[as.character(n)]], breaks = breakpoints,
       main = paste("10,000 simulations of rho.hat(2)\nTime Series of Length", n),
       xlab = "Sample ACF at Lag 3", ylab = "Frequency",
       col = "lightblue", border = "black", xlim = x_limits, ylim = y_limits)
  abline(v = true_acf_lag3, col = "blue", lwd = 2, lty = 2) # True value
}

###################
##  Problem 2.A  ##
###################

set.seed(4)

# Define parameters
phi1 <- 0.7  # coefficient for Xt-1
phi2 <- -0.6  # coefficient for Xt-2
sigma_w <- 1  # variance of white noise
n <- 100  # number of timesteps

# Initialize process
Xt <- numeric(n)
Wt <- rnorm(n, mean = 0, sd = sqrt(sigma_w))  # white noise

# Simulate AR(2) process
for (t in 3:n) {
  Xt[t] <- phi1 * Xt[t - 1] + phi2 * Xt[t - 2] + Wt[t]
}

# Compute sample ACF
sample_acf <- acf(Xt, lag.max = 10, plot = FALSE)$acf

# Compute sample ACVF at lag 0 (variance of Xt)
sample_acvf_lag0 <- var(Xt)

# Sample ACF for first few lags
print(sample_acf)
# Sample ACVF at lag 0
print(sample_acvf_lag0)

###################
##  Problem 2.B  ##
###################

# get ACF values for Yule-Walker equations
acf_1 <- sample_acf[2]  # ACF at lag 1
acf_2 <- sample_acf[3]  # ACF at lag 2

#Yule-Walker equations
R_matrix <- matrix(c(1, acf_1, acf_1, 1), nrow = 2, ncol = 2)  # Autocorrelation matrix
r_vector <- c(acf_1, acf_2)

# AR coefficients
phi_hat <- solve(R_matrix, r_vector)

# compute variance of white noise
sigma_w_hat <- sample_acvf_lag0 * (1 - sum(phi_hat * r_vector))

# Yule-Walker Estimated Coefficients
print(phi_hat)

#########################
##  Problem 2.C & 2.D  ##
#########################
library(astsa)

# Maximum Likelihood Estimation
model <- sarima(Xt, p = 2, d = 0, q = 0, no.constant = TRUE)

# Estimated coefficients
mle_phi <- model$ttable[, 1]  
# Standard errors
mle_se <- model$ttable[, 2]  

relative_differences <- abs((phi_hat - mle_phi) / mle_se)

# Maximum Likelihood Estimators
print(mle_phi)

# Standard Errors
print(mle_se)

# Yule-Walker Estimators
print(phi_hat)

# Relative Differences (Yule-Walker vs MLE)
print(relative_differences)

###################
##  Problem 2.E  ##
###################

WN.var = NULL; AIC = NULL; AICC = NULL

for(p in 1:3){
  X.model <- sarima(Xt, p,0,0, no.constant = TRUE, details = FALSE)
  WN.var[p] <- X.model$fit$sigma2
  AIC[p] <- X.model$fit$aic
  AICC[p] = -2*X.model$fit$loglik + 2*(p+0+1)*n/(n-p-0-2)
}

#print(WN.var)
print(AIC)
print(AICC)

###################
##  Problem 3.A  ##
###################

set.seed(4)

# Simulate an AR(1) process with non-zero mean
phi1 <- 0.5  # AR(1) coefficient
mu <- 2      # Non-zero mean
sigma_w <- 1 # Variance of white noise
n <- 100     # Number of timesteps

Xt <- numeric(n)
Wt <- rnorm(n, mean = 0, sd = sqrt(sigma_w))  # White noise

for (t in 2:n) {
  Xt[t] <- mu + phi1 * (Xt[t - 1] - mu) + Wt[t]
}

# Fit AR(1) model with non-zero mean using sarima
sarima_fit <- sarima(Xt, p = 1, d = 0, q = 0, no.constant = FALSE, details = FALSE)

# Extract log-likelihood from sarima output
loglik <- sarima_fit$fit$loglik


# Compute AIC manually, treating the mean as an additional parameter
p <- 1  # number of AR terms
q <- 0  # number of MR terms
aic_manual <- -2 * loglik + 2 * (p + q + 2) 


# AIC (computed manually)
print(aic_manual)
# AIC (from sarima)
print(sarima_fit$fit$aic)

aicc_manual <- -2 * loglik + 2 * (p + q + 2) * ((n-3) / (n - p - q - 2))

# AICc (computed manually)
print(aicc_manual)

# AICc (from sarima)
print(sarima_fit$ICs[1]*100)


