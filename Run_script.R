###
# PRELIMINARY SETTINGS
###

# Set parameters
r <- 0.05

tranche1 <- list(Ku=0.1, Kd=0, n=100, upfront=TRUE)
tranche2 <- list(Ku=.15, Kd=.1, n=100, upfront=TRUE)
tranche3 <- list(Ku=.25, Kd=.15, n=100, upfront=FALSE)
tranche4 <- list(Ku=.35, Kd=.25, n=100, upfront=FALSE)
tranche5 <- list(Ku=1, Kd=.35, n=100, upfront=FALSE)

tranche <- list(tranche1, tranche2, tranche3, tranche4, tranche5)
  
method <- "fourier"
theta0 <- c(1, 0.1, 10, 1)

#Advanced settings
truncate_inverse_transform <- TRUE
print_distribution <- FALSE
show_distribution_plot <- FALSE
number_points_transform <- 200   ####### TODO: THINK ABOUT THIS
n_points_integral <- 12

#Run model (for debugging/test purposes)
model(tranche = tranche[[1]], theta=theta0, method=method)

###
# IMPORT AND PREPROCESS DATA
###
setwd("/Users/bern/Documents/Courses/MS&E 347/Project/Git Project")
library(gdata)
data <- read.table("data_to_import.csv", header=TRUE, sep=",")

dates <- unique(data$Date)
data_list <- list()
keep_columns <- c("T1", "T2", "T3", "T4", "T5", "Maturity")
for (d in dates) {
  data_list[[d]] <- data[data$Date==d, keep_columns]
}

###
# RUN OPTIMIZATION
###

M <- 100 #The number of initial points
# We now create the data frame containing results
na_rep <- rep(NA, length(dates))
opt_results <- data.frame(date=dates, c=na_rep, kappa=na_rep, delta=na_rep, lambda0=na_rep, obj_value=na_rep, AAPE=na_rep)
# We run the optimization for each date available
for (d in dates) {
  mid <- data_list[[d]]
  # not generating min or max
  min_so_far <- +Inf
  params <- c(0,0,0,0)
  
  for (i in 1:M){
    # Generating theta_0
    c <- runif(1,0,8)
    kappa <- runif(1,0,10)
    delta <- runif(1,0,20)
    lambda0 <- runif(1,0,30)
    
    theta0 <- c(c, kappa, delta, lambda0)
    
    # Optimization
    res <- optim(theta0, objective)
    val <- res$value
    
    if (val<min_so_far){
      min_so_far <- val
      params <- res$par
    }
  }  
  aape <- AAPE(params)
  opt_results[d, c("c", "kappa", "delta", "opt", "obj_value", "AAPE")] <- c(params, min_so_far, aape)
}

###  
# COMPARISON
###

# Time comparison between inversion methods


# Compare results with the plots on the paper

