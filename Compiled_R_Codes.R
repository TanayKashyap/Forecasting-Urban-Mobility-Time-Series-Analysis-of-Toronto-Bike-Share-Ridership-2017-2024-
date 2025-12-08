library(tidyverse)
library(janitor)
library(lubridate)
library(MASS)
library(car)
library(randtests)
library(forecast)

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

# Model Fitting

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


# Model Eval

# 
evaluate_model <- function(model, train_df, valid_df) {
  # predictions on validation set
  preds <- predict(model, newdata = valid_df)
  
  # APSE 
  apse <- mean((valid_df$y_trans - preds)^2)
  
  # Number of parameters
  num_params <- length(coef(model)) + 1 
  
  # AICc
  n <- nrow(train_df)
  aicc <- AIC(model) + (2 * num_params * (num_params + 1)) / (n - num_params - 1)
  
  # Adjusted R-squared
  adj_r2 <- summary(model)$adj.r.squared
  
  n_coeffs <- length(coef(model))
  
  return(c(APSE = apse, AICc = aicc, Adj_R2 = adj_r2, Num_Params = n_coeffs))
}

# Collect results
results <- rbind(
  "Trend + Fourier" = evaluate_model(mod1, train_seas, valid_seas),
  "Trend + Fourier + DayOfWeek" = evaluate_model(mod2, train_seas, valid_seas),
  "Trend + Fourier + Month" = evaluate_model(mod3, train_seas, valid_seas),
  "Trend + Fourier + Month + Weekend" = evaluate_model(mod4, train_seas, valid_seas),
  "Trend + Fourier + Month + DayOfWeek" = evaluate_model(mod5, train_seas, valid_seas)
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
    fit <- lm(y ~ ., data = x_train)
    
    preds <- predict(fit, newdata = x_test)
    
    errors <- c(errors, (cv_data$y_trans[test_idx] - preds)^2)
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
saveRDS(final_reg, "final_reg.rds")

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