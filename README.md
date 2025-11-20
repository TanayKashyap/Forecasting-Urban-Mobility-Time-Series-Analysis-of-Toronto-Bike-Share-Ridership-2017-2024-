# Forecasting-Urban-Mobility-Time-Series-Analysis-of-Toronto-Bike-Share-Ridership-2017-2024-
Final group project for 443 Fall 

hello 



# Paste at the top of files:

# Packages 
library(MASS)      # boxcox
library(car)       # qqPlot
library(randtests) # runs.test

bike <- read.csv("trips_per_day.csv")
bike$trip_date <- as.Date(bike$trip_date)

str(bike)
head(bike)
range(bike$trip_date)

bike <- bike[order(bike$trip_date), ]
y <- bike$n_trips

bike_ts <- ts(
  y,
  start = c(as.numeric(format(min(bike$trip_date), "%Y")), 
            as.numeric(format(min(bike$trip_date), "%j"))),
  frequency = 365
)

tim <- time(bike_ts)  # continuous time index

month <- factor(format(bike$trip_date, "%m"))

reg_for_bc <- lm(bike_ts ~ tim + month)
boxcox_mod <- MASS::boxcox(reg_for_bc, lambda = seq(-2, 2, 0.1))
(lambda_opt_mod <- boxcox_mod$x[which.max(boxcox_mod$y)])

lam <- lambda_opt_mod   # keep this for later

if (lam == 0) {
  y_trans <- log(bike_ts)
} else if (lam > 0) {
  y_trans <- (bike_ts^lam - 1) / lam
} else { 
  # negative lambda → use minus sign trick like in lectures
  y_trans <- -(bike_ts^lam)
}

# The transformed data is y_trans
