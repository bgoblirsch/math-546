library(astsa)
library(tseries)
library(forecast)
library(ggplot2)
library(tidyverse)

chicken_df <- tibble(
  date = time(chicken),
  value = as.numeric(chicken)
)

# plot data
ggplot(chicken_df, aes(x = date, y = value)) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(floor(min(chicken_df$date)), ceiling(max(chicken_df$date)), by = 1)
  ) +
  labs(
    title = "Monthly Chicken Prices\n(Aug. 2001 - Jul. 2016)",
    x = NULL,
    y = "Price (cents per pound)"
  )
)


# check variance levels after spotting suspected increase in variance over time
window_size <- 3

# vector to store rolling st. dev.
rolling_stdev <- rep(NA, length(chicken))

# Compute rolling st. dev.
for (i in seq(window_size, length(chicken))) {
  rolling_stdev[i] <- sd(chicken[(i - window_size + 1):i])
}
# Convert rolling st. dev. to a time series object with the same frequency and start as chicken
rolling_stdev_ts <- ts(rolling_stdev, start = start(chicken), frequency = frequency(chicken))
plot(rolling_stdev_ts, main = "Rolling Standard Deviation (3-month window)", ylab = "Standard Deviation", xlab = NULL)

chicken_df <- chicken_df %>%
  mutate(log_price = log(value))

ggplot(chicken_df, aes(x = date)) +  
  geom_line(aes(y = value, color = "Raw Data"), size = 1) +
  geom_line(aes(y = log_price, color = "Log-Transformed Data"), size = 1) +
  scale_color_manual(values = c("Raw Data" = "blue", "Log-Transformed Data" = "red")) +
  labs(title = "Comparison of Raw and Log-Transformed Prices",
       x = "Time",
       y = "Price",
       color = "Legend") +  # Legend title
  theme_minimal() +
  theme(legend.position = c(0.8, 0.3),  # Position the legend (x, y) within the plot
        legend.background = element_rect(fill = "white", color = "black"),  # Add background to the legend
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

# first-order difference
chicken_df <- chicken_df %>%
  mutate(diff_log_price = log_price - lag(log_price))

# seasonal difference - lag 12
chicken_df <- chicken_df %>%
  mutate(seasonal_diff_log_price = diff_log_price - lag(diff_log_price, 12))  
ggplot(chicken_df, aes(x = date)) +
  geom_line(aes(y = seasonal_diff_log_price), color = "blue", size = 1) +
  labs(title = "Differenced Log-Transformed Chicken Prices \n(First-Order + Seasonal [lag = 12])",
       x = "Time",
       y = "Seasonally Differenced Log(Price)") +
  theme_minimal()

# augmented Dickey–Fuller test
adf_result <- adf.test(na.omit(chicken_df$seasonal_diff_log_price), alternative = "stationary")
adf_result$p.value


par(mfrow = c(1, 2))
acf(na.omit(chicken_df$seasonal_diff_log_price), 
    main = "ACF of Seasonally Differenced\nLog-Transformed Prices")
pacf(na.omit(chicken_df$seasonal_diff_log_price), 
     main = "PACF of Seasonally Differenced\nLog-Transformed Prices")

# Identify optimal ARIMA model
auto_model <- auto.arima(chicken_df$seasonal_diff_log_price, seasonal = TRUE, trace = TRUE)
auto_model

# Fit the selected ARIMA(1,0,2) with zero mean
best_model <- Arima(chicken_df$seasonal_diff_log_price, 
                    order = c(2, 0, 2), include.mean = FALSE)

log_price_ts <- ts(chicken_df$log_price, start = c(2001, 8), frequency = 12)  # Adjust start date as needed

model_sma2 <- Arima(log_price_ts, 
                    order = c(1, 1, 1), 
                    seasonal = c(1, 1, 2))


# Summary of the model
summary(best_model)

# Check residual diagnostics
checkresiduals(best_model)

# Ljung-Box test for residual autocorrelation at lag 12
ljung_box_test <- Box.test(residuals(best_model), lag = 12, type = "Ljung-Box")

# Print the test result
print(ljung_box_test)

arima_model <- sarima(chicken, p = 1, d = 0, q = 2)

arima_model <- sarima(chicken, p = 1, d = 0, q = 3)

arima_model <- sarima(chicken, p = 2, d = 0, q = 2)


