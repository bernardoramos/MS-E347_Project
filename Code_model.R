

model <- function(tranche = list(Ku=0.1, Kd=0, n=100, upfront=TRUE), theta, T, method="fourier") {
  # This function computes the rate (either upfront or running) implied by the parameters in theta.
  # PARAMETERS:
  # Tranche is a list containing the following elements: 
  #   Ku: upper tranche limit
  #   Kd: lower tranche limit
  #   T: the maturity of the CDS
  #   n: the number of contracts
  #   upfront: "upfront" (TRUE) or "running" (FALSE), to indicate whether the desired quantity is the upfront rate or the quarterly CDS spread
  # Theta is a 4-dimensional vector containing parameters (c, kappa, delta, lambda0)
  
  
  Ku <- tranche$Ku
  Kd <- tranche$Kd
  n <- tranche$n
  is_upfront <- tranche$upfront
  #We first compute the distribution of N_t and L_t based on the transforms
  adv_settings <- mget(c("print_distribution", "show_distribution_plot", "print_spreads"), ifnotfound = c(FALSE, FALSE, FALSE), envir=.GlobalEnv)
  N_distribution <- get_distribution(theta, T, method="fourier", distr="N")
  if (adv_settings$show_distribution_plot) {
    pllot <- ggplot(data=data.frame(x=seq(length(N_distribution)),y=N_distribution), aes(x=x,y=y)) + geom_line()
    plot(pllot)
  }
  if (adv_settings$print_distribution) {
    print(N_distribution)
  }
  #Compute the value of payments
  D <- PVPayments(tranche, theta, T, method = method)
  # We now compute the 
  if (is_upfront) {
    K <- n * (Ku - Kd)
    F <- D / K
    if (print_spreads) {print(c("Upfront rate:", D/K))}
    return(F)
  } else {
    sum <- 0
    premium_dates <- seq(0, T - 1/4, 1/4)
    cm <- 1/4
    for (tm in premium_dates) {
      sum <- sum + exp(-r*tm) * cm * premium_notional(tm, N_distribution)
    }
    S <- D / sum
    if (print_spreads) {print(c("Running spread:", S))}
    return(S)
  }
}

PVPayments <- function(tranche, theta, T, method) {
  # This function computes the present value of premium payments, D_t given in formula 34
  Ku <- tranche$Ku
  Kd <- tranche$Kd
  F <- tranche$F 
  n <- tranche$n
  distribution <- get_distribution(theta, T, method=method, distr="L") #######
  func <- function(s) {
    distribution <- get_distribution(theta=theta, s, method=method, distr="L")
    exp(-r * s) * expectation_U(distribution, s, Ku, Kd, n, method = method)
  }
  adv_settings <- mget(c("n_points_integral"), ifnotfound = c(30), envir=.GlobalEnv)
  D <- exp(-r * T) * expectation_U(distribution, T, Ku, Kd, n, method = method) + r * integr(func, 0, T, delta = 1/adv_settings$n_points_integral)
  return(D)
}

expectation_U <- function(distribution, s, Ku, Kd, n=10, method) {
  # This is the expectation E(U_t), computed using the fact that U_t = (L_t - Kd * n)_+  - (L_t - Ku * n)_+
  # distribution must be a matrix 
  ### TODO: For now, we assume L = 0.6 * N. Change this code along with the computation of a and b when accounting for general Ls
  values <- 0.6 * seq(0, length(distribution)-1)
  expectation <- sum((pmax(values - Kd * n, 0) - pmax(values - Ku * n, 0))* distribution)
  return(expectation)
}

integr <- function(f, lower, upper, delta=1/12) {
  val <- 0
  if (lower<upper) {
    n <- length(seq(lower+delta, upper, by=delta))
    for (s in seq(lower+delta, upper, by=delta)) { val <- val + f(s) / n }
  }
  return(val)
}

premium_notional <- function(n, distribution) {
  # This function computes the expectation of the premium notional, E(I_t)
  # Note: verification that sum(distribution) = 1?
  values <- seq(0, length(distribution)-1) 
  expectation <- n - sum(pmin(n, values) * distribution)
  return(expectation)
}

get_distribution <- function(theta, T, method="fourier", distr="N") {
  if (distr == "N") {
    u <- c(0, 1)
  } else {
    u <- c(1, 0)
  }
  if (method == "fourier") {
    adv_settings <- mget(c("truncate_inverse_transform", "number_points_transform"), ifnotfound = c(FALSE, n*2), envir=.GlobalEnv)
    N <- adv_settings$number_points_transform
    if (distr == "N") {
      u <- matrix(c(0,1), nrow=N, ncol=2, byrow=TRUE)
    } else {
      u <- matrix(c(1,0), nrow=N, ncol=2, byrow=TRUE)
    }
    u <- u * 1i*2*pi*(0:(N-1))/N
    alp <- a(u, theta, T, N=100) 
    bet <- b(u, theta, T, N=100) 
    lambda0 <- theta[4]
    transform <- exp(alp + bet * lambda0)
    distr <- fft(transform, inverse=TRUE)
    distr <- distr/length(distr) #We divide to normalize the transform
    distr <- Mod(distr)
    if (adv_settings$truncate_inverse_transform) {
      distr <- distr[1:floor(N*0.6)]/sum(distr[1:floor(N*0.6)]) #We normalize the distribution ## TODO: what to do about this?
    } else {
      distr <- distr[1:N]/sum(distr[1:N]) #We normalize the distribution
    }
  } else if (method == "derivative") {
    if (method != "N") {warning("Using the derivative method for the distribution of L")}
  }
  return(distr)
}


# params = c, kappa, delta, lambda0

# ImplÈmentÈ ‡ la main

# N number of steps
# T maturity
# Initial condition : 0

partial_b <- function(y, u, params){
  c <- params[1]
  kappa <- params[2]
  delta <- params[3]
  lambda0 <- params[4]
  
  theta <- exp(.6*(delta*y+u[,1]))
  res <- kappa*y+1-theta*exp(u[,2])
  return(-res)
}

# returns a vector of size N+1 from t=0 to t=T que l'on renverse

b <- function(u, params, T, N, full=FALSE){
  h <- T/N
  
  y <- matrix(0, ncol=N+1, nrow=nrow(u))
  res <- matrix(0, ncol=N+1, nrow=nrow(u))
  
  for (i in 1:N){
    k1 <- partial_b(y[, i], u, params)
    k2 <- partial_b(y[, i]+h/2*k1, u, params)
    k3 <- partial_b(y[, i]+h/2*k2, u, params)
    k4 <- partial_b(y[, i]+h*k3, u, params)
    
    y[, i+1] <- y[, i] +h/6*(k1+2*k2+2*k3+k4)
    res[, N+1-i] <- y[, i+1]
  }
  if (full) {
    return(res)
  } else {
    return(res[, 1])
  }
}

a <- function(u, params, T, N){
  c <- params[1]
  kappa <- params[2]
  b_t <- b(u, params, T, 2*N, full=TRUE)
  
  h <- T/N
  
  y <- 0
  res <- 0
  
  for (i in 1:N){
    k1 <- -kappa*c*b_t[, 2*i-1]
    k2 <- -kappa*c*b_t[, 2*i]
    k3 <- -kappa*c*b_t[, 2*i]
    k4 <- -kappa*c*b_t[, 2*i+1]
    
    y <- y +h/6*(k1+2*k2+2*k3+k4)
    res <- y
  }
  return(res)
}

