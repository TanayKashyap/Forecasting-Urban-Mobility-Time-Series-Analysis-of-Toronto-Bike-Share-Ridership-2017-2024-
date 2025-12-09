library(tidyverse)
library(janitor)
library(lubridate)
library(MASS)
library(car)
library(randtests)
library(forecast)
library(tseries)
library(glmnet)
library(patchwork)

## -------------------------------------------------------------------------
## DATA WRANGLING: Professor, you will not be able to run this portion as
## the raw data files were on our local computers - this is just to demo how 
## wrangling was done. 
## -------------------------------------------------------------------------

library(tidyverse)
library(janitor)

data_folder <- "all_data_files"
file_list <- list.files(path = data_folder, pattern = "*.csv", full.names = TRUE)

process_bikeshare_file <- function(file_path) {
  df <- read_csv(file_path, col_types = cols(.default = "c")) %>% 
    # clean the messed up col names
    clean_names()
  if ("trip_start_time" %in% names(df)) {
    df <- df %>% rename(start_time = trip_start_time)
  }
  if ("trip_stop_time" %in% names(df)) {
    df <- df %>% rename(end_time = trip_stop_time)
  }
  # Select only the columns we want
  df_clean <- df %>%
    select(any_of(c("trip_id", "start_time", "end_time")))
  
  return(df_clean)
}

all_rides_data <- map_dfr(file_list, process_bikeshare_file)
print(paste("Total rides compiled:", nrow(all_rides_data)))

library(lubridate)

all_rides_data <- all_rides_data %>%
  mutate(
    # Parse timestamps in the common format
    start_time_parsed = mdy_hm(start_time, quiet = TRUE),
    
    # Extract date portion
    trip_date = as.Date(start_time_parsed)
  )

#aggregate the trips for daily 
trips_per_day <- all_rides_data %>%
  group_by(trip_date) %>%
  summarise(
    n_trips = n(),
    .groups = "drop"
  ) %>%
  arrange(trip_date)

# remove all NA dates
all_rides_data <- all_rides_data %>%
  filter(!is.na(trip_date))

# pull in 2017-2018 dates from a diff file
old_trips <- read_csv("2017_2018_trips.csv",
                      show_col_types = FALSE)

trips_per_day <- old_trips %>%
  bind_rows(trips_per_day) %>%
  arrange(trip_date)

write_csv(trips_per_day, "trips_per_day.csv")


## -------------------------------------------------------------------------
## PROJECT BEGINS HERE: trips_per_day.csv is not the original data
## file we sent you (it is the wrangled version). If you'd like to run this 
## code we can attach that updated datafile to the submission 
## -------------------------------------------------------------------------


# Load Data 
bike <- read.csv("trips_per_day.csv")
bike <- bike %>% filter(!is.na(trip_date))
bike$trip_date <- as.Date(bike$trip_date)

# Filter for 2017-2024 
bike <- bike[bike$trip_date >= as.Date("2017-01-01"), ]

# Create Time Series Object w/ Daily Frequency
y <- bike$n_trips
bike_ts <- ts(y, start = c(2017, 1), frequency = 365)

# TS Plot
plot(bike_ts, main = "Daily BikeShare Trips (2017-2024)", ylab = "Number of Trips")

bike$trip_date <- as.Date(bike$trip_date)
bike$n_month <- as.numeric(format(bike$trip_date, "%m"))
bike$month <- factor(bike$n_month, levels = 1:12, labels = month.abb)
boxplot(n_trips ~ month, data = bike,
        main = "Boxplot of Annual BikeShare Trips", ylab = "Number of Trips",
        xlab = "Month", col = "steelblue")

# Check ACF with high lag max to get period
acf(y, main = "ACF of Raw Data")
spec <- spectrum(y, spans = 5, main = "Spectral Density", plot = FALSE)

# Check dominant period
print(1 / spec$freq[which.max(spec$spec)])
acf(y, lag.max = 400, main = "ACF of Raw Data")

# Box-Cox 
bc <- boxcox(lm(y ~ 1), lambda = seq(-2, 2, 0.1))
lam <- bc$x[which.max(bc$y)]
print(paste("Optimal Lambda:", lam))

if (lam == 0) {
  y_trans <- log(bike_ts)
} else if (lam > 0) {
  y_trans <- (bike_ts^lam - 1) / lam
} else { 
  y_trans <- -(bike_ts^lam)
}

plot(y_trans, main = "Transformed Daily Trips", ylab = "Transformed trips")
bike$y_trans <- y_trans # Add to dataframe for regression

# Classical Decomposition
decomp <- stl(ts(y_trans, frequency=365), s.window="periodic")
plot(decomp, main="Decomposition of Transformed Data")

# Add indicators 
bike <- bike %>%
  mutate(
    # 1. Time Index
    time_index = 1:n(),
    
    # 2. Factor-based Seasonality
    month_fac = as.factor(month(trip_date)),
    weekday_fac = factor(wday(trip_date, label = TRUE), ordered = FALSE),
    is_weekend = as.factor(ifelse(wday(trip_date) %in% c(1, 7), "Yes", "No")),
    
    # 3. Fourier Term
    # Using period = 360
    sin_year = sin((2 * pi * time_index / 360) + 0.75),
  )

# Train on 2017-2022, Test on 2023-2024
train_seas <- subset(bike, trip_date < "2023-01-01")
valid_seas <- subset(bike, trip_date >= "2023-01-01")

# --- Model Fitting ---

# 1. Trend + Fourier Only
mod1 <- lm(y_trans ~ time_index + sin_year , data = train_seas)

# 2. Trend + Fourier + Month + Day of Week
mod2 <- lm(y_trans ~ time_index + sin_year + weekday_fac, data = train_seas)

# 3. Trend + Fourier + Month
mod3 <- lm(y_trans ~ time_index + sin_year  + month_fac, data = train_seas)

# 4. Trend + Fourier + Month + Weekend/Weekday
mod4 <- lm(y_trans ~ time_index + sin_year  + month_fac + is_weekend, data = train_seas)

# 5. Trend + Fourier + Month + Day of Week
mod5 <- lm(y_trans ~ time_index + sin_year  + month_fac + weekday_fac, data = train_seas)


# invert the boxcox

inv_boxcox <- function(z, lambda) {
  if (lambda == 0) {
    # inverse of log
    return(exp(z))
  } else if (lambda > 0) {
    # inverse of (y^λ - 1)/λ
    return(((lambda * z) + 1)^(1 / lambda))
  } else {
    # you used: y_trans = -(y^λ)
    # so y^λ = -y_trans  →  y = (-y_trans)^(1/λ)
    return((-z)^(1 / lambda))
  }
}


# --- Model Evaluation ---

# Function to calculate evaluation metrics
evaluate_model <- function(model, train_df, valid_df, lambda) {
  # Predictions on validation set (still on transformed scale)
  preds_trans <- predict(model, newdata = valid_df)
  
  # Back-transform predictions AND actuals to raw scale
  preds_raw   <- inv_boxcox(preds_trans, lambda)
  actual_raw  <- valid_df$n_trips   # original outcome
  
  # APSE on raw scale
  apse <- mean((actual_raw - preds_raw)^2)
  
  # Number of parameters (k): coefficients (beta) + 1 (sigma)
  num_params <- length(coef(model)) + 1 
  
  n <- nrow(train_df)
  aicc <- AIC(model) + (2 * num_params * (num_params + 1)) / (n - num_params - 1)
  adj_r2 <- summary(model)$adj.r.squared
  n_coeffs <- length(coef(model))
  
  return(c(APSE = apse, AICc = aicc, Adj_R2 = adj_r2, Num_Params = n_coeffs))
}


results <- rbind(
  "Trend + Fourier" = evaluate_model(mod1, train_seas, valid_seas, lam),
  "Trend + Fourier + DayOfWeek" = evaluate_model(mod2, train_seas, valid_seas, lam),
  "Trend + Fourier + Month" = evaluate_model(mod3, train_seas, valid_seas, lam),
  "Trend + Fourier + Month + Weekend" = evaluate_model(mod4, train_seas, valid_seas, lam),
  "Trend + Fourier + Month + DayOfWeek" = evaluate_model(mod5, train_seas, valid_seas, lam)
)

print("Model Evaluation Results:")
print(results)


# Visualizing the 5 fits
plot_data <- bike %>%
  dplyr::select(trip_date, y_trans, time_index, month_fac, weekday_fac, is_weekend, sin_year)
plot_data$Pred_Mod1 <- predict(mod1, newdata = plot_data)
plot_data$Pred_Mod2 <- predict(mod2, newdata = plot_data)
plot_data$Pred_Mod3 <- predict(mod3, newdata = plot_data)
plot_data$Pred_Mod4 <- predict(mod4, newdata = plot_data)
plot_data$Pred_Mod5 <- predict(mod5, newdata = plot_data)

plot_long <- plot_data %>%
  dplyr::select(trip_date, y_trans, starts_with("Pred")) %>%  # <--- ADD dplyr:: HERE
  pivot_longer(cols = c(y_trans, starts_with("Pred")), 
               names_to = "Series", 
               values_to = "Value")

# Plot
ggplot(plot_long, aes(x = trip_date, y = Value, color = Series)) +
  geom_line(alpha = 0.6) +
  geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed", color = "black") +
  labs(title = "Model Fits vs. Actual Data (2017-2024)",
       subtitle = "Models trained on data before 2023 (left of dashed line)",
       y = "Transformed Trips", x = "Date") +
  theme_minimal() +
  scale_color_manual(values = c("y_trans" = "gray", 
                                "Pred_Mod1" = "orange", 
                                "Pred_Mod2" = "blue", 
                                "Pred_Mod3" = "green", 
                                "Pred_Mod4" = "red",
                                "Pred_Mod5" = "purple")) +
  theme(legend.position = "bottom")

plot_one_model <- function(pred_col, model_name) {
  ggplot(plot_data, aes(x = trip_date)) +
    geom_line(aes(y = y_trans), color = "gray", alpha = 0.5) +
    geom_line(aes(y = .data[[pred_col]]), color = "blue") +
    geom_vline(xintercept = as.Date("2023-01-01"), linetype = "dashed") +
    labs(title = paste("Fit:", model_name), y = "Transformed Trips", x = "Date") +
    theme_minimal()
}

# plot each one by one
plot_one_model("Pred_Mod1", "Trend + Fourier")
plot_one_model("Pred_Mod2", "Trend + Fourier + DayOfWeek")
plot_one_model("Pred_Mod3", "Trend + Fourier + Month")
plot_one_model("Pred_Mod4", "Trend + Fourier + Month + Weekend")
plot_one_model("Pred_Mod5", "Trend + Fourier + Month + DayOfWeek")

# CV to select polynomial degree

# Setup
max_degree <- 8
results_poly <- data.frame(Degree = 1:max_degree, APSE = NA)
min_train_size <- 730 # first 2 years minimum

print("Running Rolling Origin CV for Trend with Best Seasonal Structure")

# Limit CV to pre-2024 data to avoid looking at the final test set
cv_data <- subset(bike, trip_date < "2024-01-01")
n_cv <- nrow(cv_data)

for (p in 1:max_degree) {
  errors <- numeric()
  poly_basis <- poly(cv_data$time_index, p, raw = TRUE)
  
  for (i in seq(min_train_size, n_cv - 30, by = 30)) {
    
    train_idx <- 1:i
    test_idx <- (i + 1):(i + 30) 
    x_train <- data.frame(
      y = cv_data$y_trans[train_idx],
      poly = poly_basis[train_idx, , drop=FALSE],
      month = cv_data$month_fac[train_idx],
      weekday = cv_data$weekday_fac[train_idx],
      sin_year = cv_data$sin_year[train_idx]    )
    
    x_test <- data.frame(
      poly = poly_basis[test_idx, , drop=FALSE],
      month = cv_data$month_fac[test_idx],
      weekday = cv_data$weekday_fac[test_idx],
      sin_year = cv_data$sin_year[test_idx]    )
    # Fit Model
    fit <- lm(y ~ ., data = x_train)
    
    # Predict on transformed scale
    preds <- predict(fit, newdata = x_test)
    
    preds_raw  <- inv_boxcox(preds, lam)
    actual_raw <- cv_data$n_trips[test_idx]
    errors <- c(errors, (actual_raw - preds_raw)^2)
  }
  
  results_poly$APSE[p] <- mean(errors)
  print(paste("Degree", p, "APSE:", round(results_poly$APSE[p], 4)))
}

# Plot Results
plot(results_poly$Degree, results_poly$APSE, type='b', main="Polynomial Selection (with Best Seasonality)", ylab = "APSE",
     xlab = "Degree of Polynomial$")
best_p <- which.min(results_poly$APSE)
points(best_p, results_poly$APSE[best_p], col = "red", cex = 2, pch = 1)
text(best_p, results_poly$APSE[best_p], paste("Best:", best_p), pos = 3, col = "red")

print(paste("Best Polynomial Degree:", best_p))


# Fit final Regression Model
# Define Training Data (2017-2023)
train_full <- subset(bike, trip_date < "2024-01-01")

# Create basis for the trend
poly_basis_final <- poly(train_full$time_index, 1, raw = TRUE)

# Build the dataframe for lm
final_train_data <- data.frame(
  y = train_full$y_trans,
  poly = poly_basis_final,
  month = train_full$month_fac,
  weekday = train_full$weekday_fac,
  sin_year = train_full$sin_year
)

# Fit the model
final_reg <- lm(y ~ ., data = final_train_data)

summary(final_reg)


# check VIF for multicollinearity
library(car)
vif_values <- vif(final_reg)
print(vif_values)
high_vif <- vif_values[,"GVIF^(1/(2*Df))"] > 2
print("Variables with potentially problematic collinearity:")
print(vif_values[high_vif, ])

# Trying LASSO

library(glmnet)
X_train <- model.matrix(
  y_trans ~ poly(time_index, 1, raw = TRUE) + month_fac + weekday_fac + sin_year,
  data = train_full
)[, -1]   # remove intercept

Y_train <- train_full$y_trans

set.seed(443)
cv_lasso <- cv.glmnet(
  X_train,
  Y_train,
  alpha = 1,
  standardize = TRUE,
  intercept = TRUE,
  family = "gaussian"
)

best_lambda <- cv_lasso$lambda.min
best_lambda

lasso_model <- glmnet(
  X_train, Y_train,
  alpha = 1,
  lambda = best_lambda,
  standardize = TRUE,
  intercept = TRUE,
  family = "gaussian"
)

coef(lasso_model)


# AIC/BIC for OLS vs Lasso using Gaussian likelihood 
calc_AIC <- function(rss, n, k) {
  sigma2 <- rss / n
  loglik <- -n/2 * (log(2*pi*sigma2) + 1)
  AIC <- -2*loglik + 2*k
  return(AIC)
}

calc_BIC <- function(rss, n, k) {
  sigma2 <- rss / n
  loglik <- -n/2 * (log(2*pi*sigma2) + 1)
  BIC <- -2*loglik + log(n)*k
  return(BIC)
}

calc_AICc <- function(aic, n, k) {
  aic + (2*k*(k+1)) / (n - k - 1)
}

# 1. OLS model (final_reg)
rss_ols <- sum(residuals(final_reg)^2)
n_ols <- nobs(final_reg)
k_ols <- length(coef(final_reg))

AIC_ols  <- calc_AIC(rss_ols, n_ols, k_ols)
BIC_ols  <- calc_BIC(rss_ols, n_ols, k_ols)
AICc_ols <- calc_AICc(AIC_ols, n_ols, k_ols)
adjR2_ols <- summary(final_reg)$adj.r.squared

# 2. Lasso model
lasso_pred_train <- predict(lasso_model, newx = X_train)

rss_lasso <- sum((Y_train - lasso_pred_train)^2)
n_lasso <- length(Y_train)
k_lasso <- sum(coef(lasso_model) != 0)  # df = non-zero coefficients

AIC_lasso  <- calc_AIC(rss_lasso, n_lasso, k_lasso)
BIC_lasso  <- calc_BIC(rss_lasso, n_lasso, k_lasso)
AICc_lasso <- calc_AICc(AIC_lasso, n_lasso, k_lasso)

# Adjusted R2
tss <- sum((Y_train - mean(Y_train))^2)
adjR2_lasso <- 1 - ((rss_lasso/(n_lasso - k_lasso - 1)) /
                      (tss/(n_lasso - 1)))

# 3. Final comparison table
comparison <- data.frame(
  Model = c("OLS", "Lasso"),
  AIC = c(AIC_ols, AIC_lasso),
  AICc = c(AICc_ols, AICc_lasso),
  BIC = c(BIC_ols, BIC_lasso),
  Adjusted_R2 = c(adjR2_ols, adjR2_lasso)
)

comparison

# Residual Diagnostics on Selected Model

RegressionDiagnosicsPlots <- function(model, my_residuals=NULL, my_fitted=NULL, title_suffix="") {
  if(is.null(my_residuals)) {
    my_residuals <- residuals(model)
    my_fitted <- fitted(model)
    model_name <- deparse(substitute(model))
  } else {
    model_name <- title_suffix
  }
  
  # Histogram
  hist(my_residuals, xlab = "Residuals", main = paste("Hist:", model_name)) 
  
  # QQ Plot
  car::qqPlot(my_residuals, pch = 16, col = adjustcolor("black", 0.7),
              xlab = "Theoretical Quantiles (Normal)", ylab = "Sample Quantiles",
              main = "Normal Q-Q Plot")
  
  # Fitted vs Residuals
  if(!is.null(my_fitted)){
    plot(my_fitted, my_residuals, pch = 16, col = adjustcolor("black", 0.5),
         xlab = "Fitted Values", ylab = "Residuals", main = "Fitted vs Residuals")
    abline(h = 0, lty = 2, col = 'red')
  } else {
    plot.new() # Empty plot if no fitted values
  }
  
  # Time vs Residuals
  plot(my_residuals, pch = 16, col = adjustcolor("black", 0.5),
       xlab = "Time", ylab = "Residuals", main = "Residuals vs Time")
  abline(h = 0, lty = 2, col = 'red')
  
  # ACF
  acf(my_residuals, main = "ACF of Residuals")
  
  par(mfrow = c(1, 1)) # Reset layout
}

RegressionDiagnosicsTests <- function(model, my_residuals=NULL, segments=6) {
  
  if(is.null(my_residuals)) {
    my_residuals <- residuals(model)
  }
  
  print("--- Shapiro-Wilk Test (Normality) ---")
  if(length(my_residuals) > 5000) {
    print(shapiro.test(sample(my_residuals, 5000)))
  } else {
    print(shapiro.test(my_residuals))
  }
  
  print("--- Kolmogorov-Smirnov Test (Normality) ---")
  print(ks.test(scale(my_residuals), "pnorm"))
  
  print("--- Fligner-Killeen Test (Homoscedasticity) ---")
  n <- length(my_residuals)
  seg <- factor(cut(1:n, breaks = segments, labels = FALSE))
  print(fligner.test(my_residuals, seg)) 
  
  print("--- Runs Test (Randomness) ---")
  par(mfrow = c(1, 1))
  print(randtests::runs.test(my_residuals))
}

# Generate Plots
RegressionDiagnosicsPlots(final_reg)

# Generate Formal Tests
RegressionDiagnosicsTests(final_reg)


## ------------------------------------------------------------------------- ##
## Regression + SARIMA Model                                                 ##
## ------------------------------------------------------------------------- ##
print("Hybrid SARIMA analysis on Regression Residuals")

# Extract residuals from the OLS model
resid_reg <- residuals(final_reg)

# Convert to TS with weekly frequency (Seasonality = 7 days)
resid_ts <- ts(resid_reg, frequency = 7)

# Visualizing Residuals ACF/PACF
p_ts <- autoplot(resid_ts) + ggtitle("Residuals")
p_acf <- ggAcf(resid_ts) + ggtitle("ACF")
p_pacf <- ggPacf(resid_ts) + ggtitle("PACF")

print((p_ts / (p_acf | p_pacf)) + plot_annotation(title = "ACF & PACF of Regression Residuals"))

# Determine Order of Seasonal and Non-Seasonal Lags
trend_lag <- ndiffs(resid_ts)
seasonal_lag <- ndiffs(resid_ts, m=7)
print(paste("Recommended trend diffs:", trend_lag, "| Seasonal diffs:", seasonal_lag))

# Differencing analysis

## Trend Differencing
resid_diff1 <- diff(resid_ts, differences = trend_lag)

## Seasonal Differencing
resid_diff2 <- diff(resid_diff1, lag = 7, differences = max(1, seasonal_lag))

ts_diff_plot   <- autoplot(resid_diff2) + ggtitle("Seasonally Differenced Residuals")
acf_diff_plot  <- ggAcf(resid_diff2) + ggtitle("ACF")
pacf_diff_plot <- ggPacf(resid_diff2) + ggtitle("PACF")

# Plot PACF + ACF of Differenced Residuals
print((ts_diff_plot / (acf_diff_plot | pacf_diff_plot)) +
        plot_annotation(title = "ACF & PACF of Seasonal & Non-Seasonal Differenced Residuals"))

# Fit SARIMA Models to Residuals
print("Fitting SARIMA models to residuals...")

# Model 1: SARIMA(0,1,1)(0,1,1)[7]
model1 <- Arima(resid_ts,
                order=c(0,1,1),
                seasonal=list(order=c(0,1,1), period=7),
                method="ML")

# Model 2: SARIMA(1,1,1)(0,1,1)[7]
model2 <- Arima(resid_ts,
                order=c(1,1,1),
                seasonal=list(order=c(0,1,1), period=7),
                method="ML")

# Model 3: Auto ARIMA
model_auto <- auto.arima(resid_ts, seasonal=TRUE, stepwise=FALSE, approximation=FALSE)

# Diagnostic Function from 443.Rmd
diagnose_model <- function(model, model_name = "Model") {
  cat("\n==============================\n")
  cat(" Diagnostics for:", model_name, "\n")
  cat("==============================\n\n")
  
  res <- residuals(model)
  
  cat("AIC :", AIC(model), "\n")
  cat("BIC :", BIC(model), "\n")
  rmse_in <- sqrt(mean(res^2, na.rm = TRUE))
  cat("In-sample RMSE:", rmse_in, "\n")
  
  print(checkresiduals(model))
  
  return(data.frame(
    Model = model_name,
    AIC = AIC(model),
    BIC = BIC(model),
    RMSE_in = rmse_in
  ))
}

diag1 <- diagnose_model(model1, "SARIMA(0,1,1)(0,1,1)[7]")
diag2 <- diagnose_model(model2, "SARIMA(1,1,1)(0,1,1)[7]")
diag_auto <- diagnose_model(model_auto, "auto.arima()")

comparison_table <- bind_rows(diag1, diag2, diag_auto)
print(comparison_table)

print("Selected Model for Hybrid Forecasting: SARIMA(1,1,1)(0,1,1)[7] (Model 2)")

## -------------------------------------------------------------------------
## HYBRID FORECASTING (2024 Test Set)
## -------------------------------------------------------------------------

# Create Test Data
test_data_2024 <- subset(bike, trip_date >= "2024-01-01")

# Create polynomial basis like for Train Data
poly_basis_test <- poly(test_data_2024$time_index, 1, raw = TRUE)

final_test_data <- data.frame(
  y = test_data_2024$y_trans, # Actual transformed values
  poly = poly_basis_test,
  month = test_data_2024$month_fac,
  weekday = test_data_2024$weekday_fac,
  sin_year = test_data_2024$sin_year
)

# Forecast Regression Component
reg_forecast_trans <- predict(final_reg, newdata = final_test_data)

# Forecast Residual Component (SARIMA)
h <- nrow(final_test_data)
sarima_fc <- forecast(model2, h = h)
sarima_resid_forecast <- sarima_fc$mean

# Combine Components
hybrid_forecast_trans <- reg_forecast_trans + sarima_resid_forecast

# Inverse Transformation (Box-Cox)
hybrid_forecast <- InvBoxCox(hybrid_forecast_trans, lambda = lam)
actual_2024 <- InvBoxCox(final_test_data$y, lambda = lam)

# Evaluate Accuracy
accuracy_hybrid <- accuracy(hybrid_forecast, actual_2024)
apse_test <- mean((actual_2024 - hybrid_forecast)^2)

accuracy_hybrid_df <- as.data.frame(accuracy_hybrid)
accuracy_hybrid_df$APSE <- apse_test

print("Hybrid Model Accuracy (2024):")
print(accuracy_hybrid_df)

# Final Plot: 2024 Forecast vs Actuals
plot(test_data_2024$trip_date, actual_2024, type = "l", col = "gray", lwd = 2,
     main = "2024 Forecast vs Actuals (Hybrid OLS + SARIMA)", 
     ylab = "Number of Trips", xlab = "Date")
lines(test_data_2024$trip_date, hybrid_forecast, col = "red", lwd = 2)
legend("topright", legend = c("Actual 2024 Data", "Hybrid Forecast"), 
       col = c("gray", "red"), lty = 1, lwd = 2)



# Holt Winters 

data <- read.csv("trips_per_day.csv")
data$trip_date <- as.Date(data$trip_date)
data <- data[data$trip_date >= as.Date("2017-01-01"), ]

y <- data$n_trips
y <- y[!is.na(y)]
data_ts <- ts(y, start = c(2017, 1), frequency = 365)


#multiplicative model
hw.multiplicative <- HoltWinters(y_trans, seasonal = "multiplicative")
pred.multiplicative <- predict(hw.multiplicative, n.ahead=365)

#invert predictions
y_mult_inv <- InvBoxCox(pred.multiplicative, lam)
pred_ts <- ts(y_mult_inv,
              start = end(data_ts) + c(0, 1),
              frequency = frequency(data_ts))


#plot of transformed predictions
plot(hw.multiplicative, predict(hw.multiplicative, n.ahead=365), main = "Holt Winters (Multiplicative)", ylab = "Transformed Values")
legend("bottomright",
       legend = c("Observed", "Last Observed", "Forecast"),
       col = c("black", "black", "red"),
       lty = c(1, 2, 1),
       bty = "o")

#plot of inverted predictions
plot(bike_ts, main = "Holt Winters Multiplicative", ylab = "Trips", xlim = c(start(bike_ts)[1] + (start(bike_ts)[2]-1)/frequency(bike_ts),
                                                                             2026 + 0))
lines(pred_ts, col="red",lwd = 1)
legend("topleft",
       legend = c("Observed", "Inverted Forecast"),
       col = c("black", "red"),
       lty = c(1, 1),
       bty = "o")


#additive model
hw.additive <- HoltWinters(y_trans, seasonal = "additive")
pred.additive = predict(hw.additive, n.ahead=365)

plot(hw.additive, predict(hw.additive, n.ahead=365), main = "Holt Winters (Additive)", ylab = "Transformed Values")
legend("topleft",
       legend = c("Observed", "Last Observed", "Forecast"),
       col = c("black", "black", "red"),
       lty = c(1, 2, 1),
       bty = "o")


pred_ts <- ts(pred.additive,
              start = end(data_ts) + c(0, 1),
              frequency = frequency(data_ts))


y_add_inv <- InvBoxCox(pred.additive, lam)
pred_ts <- ts(y_add_inv,
              start = end(data_ts) + c(0, 1),
              frequency = frequency(data_ts))

plot(bike_ts, main = "Holt Winters Additive", ylab = "Trips", xlim = c(start(bike_ts)[1] + (start(bike_ts)[2]-1)/frequency(bike_ts),
                                                                       2026 + 0))

lines(pred_ts, col="red",lwd = 1)
legend("topleft",
       legend = c("Observed", "Inverted Forecast"),
       col = c("black", "red"),
       lty = c(1, 1),
       bty = "o")


# inverted fitted values
y_fitted_add_inv <- InvBoxCox(hw.additive$fitted[, 1], lam)
y_fitted_mult_inv <- InvBoxCox(hw.multiplicative$fitted[, 1], lam)
y_trimmed <- y[-(1:365)]

#residuals
resid_add <- y_trimmed - y_fitted_add_inv
resid_mult <- y_trimmed - y_fitted_mult_inv

#transformed residuals
y_add <- hw.additive$fitted[, 1]
y_mult <- hw.multiplicative$fitted[, 1]
y_tr <- y_trans[-(1:365)]

res_add <- y_tr - y_add
res_mult <- y_tr - y_mult

#acf
acf(resid_add, main="ACF of HW Residuals for Additive Model")
acf(resid_mult, main="ACF of HW Residuals for Multiplicative Model")

#residual vs time
plot(resid_add)
plot(resid_mult)

period <- 365   # 12 months of daily data
n <- length(y)

y_test <- y[(n - period + 1):n]

#transformed test and trains sets
y_train_trans <- y_trans[1:(n - period)]
y_test_trans  <- y_trans[(n - period + 1):n]
y_train_ts <- ts(y_train_trans,
                 start = start(data_ts),
                 frequency = frequency(data_ts))


#additive model on transformed data
hw_train_add <- HoltWinters(y_train_ts, seasonal = "additive")
hw_test_add = predict(hw_train_add, n.ahead=365)

#multiplicative model on transformed data
hw_train_mult <- HoltWinters(y_train_ts, seasonal = "multiplicative")
hw_test_mult = predict(hw_train_mult, n.ahead=365)

y_mult_inv <- InvBoxCox(hw_test_mult, lam)
y_add_inv <- InvBoxCox(hw_test_add, lam)


#APSE
MSE.add <- mean((y_test - y_add_inv)^2)
MSE.mult <- mean((y_test - y_mult_inv)^2)
MSE.add
MSE.mult

