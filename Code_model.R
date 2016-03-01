
r <- 0.05
K <- 0

model <- function(tranche = list(Ku=1, Kd=0, F=0, T=3, n=100), theta, method="fourier") {
  # Tranche is a list containing the following elements: 
  # Ku: upper tranche limit
  # Kd: lower tranche limit
  # F: the upfront payment rate for the tranche
  # T: the maturity of the CDS
  # n: the number of contracts
  # Theta is a 4-dimensional vector containing parameters (c, kappa, delta, lambda0)

  Ku <- tranche$Ku
  Kd <- tranche$Kd
  T <- tranche$T
  F <- tranche$F 
  n <- tranche$n
  #We first compute the distribution of N_t and L_t based on the transforms
  N_distribution <- get_distribution(theta, T, method="fourier", distr="N")
  print(N_distribution)
  #Compute the value of payments
  D <- PVPayments(tranche, theta, method = method)
  sum <- 0
  premium_dates <- seq(0, T - 1/4, 1/4)
  cm <- 1/4
  for (tm in premium_dates) {
    sum <- sum + exp(-r*tm) * cm * premium_notional(tm, N_distribution)
  }
  print(sum)
  print(D)
  print(F)
  print(K)
  S <- (D - F * K) / sum
  return(S)
}

PVPayments <- function(tranche, theta, method) {
  # This function computes the present value of premium payments, D_t given in formula 34
  Ku <- tranche$Ku
  Kd <- tranche$Kd
  T <- tranche$T
  F <- tranche$F 
  n <- tranche$n
  distribution <- get_distribution(theta, T, method=method, distr="L") #######
  func <- function(s) {
    distribution <- get_distribution(theta=theta, s, method=method, distr="L")
    exp(-r * s) * expectation_U(distribution, s, Ku, Kd, n, method = method)
  }
  D <- exp(-r * T) * expectation_U(distribution, T, Ku, Kd, n, method = method) + r * integr(func, 0, T, delta = 1/100)
  print("D")
  print(D)
  print("expectation_U")
  print(expectation_U(distribution, T, Ku, Kd, n, method = method))
  print("distribution")
  print(distribution)
  return(D)
}

expectation_U <- function(distribution, s, Ku, Kd, n=10, method) {
  # This is the expectation E(U_t), computed using the fact that U_t = (L_t - Kd * n)_+  - (L_t - Ku * n)_+
  # distribution must be a matrix 
  ### For now, we assume L = 0.6 * N
  values <- 0.6 * seq(0, length(distribution)) 
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
  values <- seq(0, length(distribution))
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
    N <- 120 ### OR N = 1100? ####### THINK ABOUT THIS
    if (distr == "N") {
      u <- matrix(c(0,1), nrow=N, ncol=2, byrow=TRUE)
    } else {
      u <- matrix(c(0,1), nrow=N, ncol=2, byrow=TRUE)
    }
    u <- u * 1i*2*pi*(1:(N))/N
    alp <- a(u, theta, T, N=100) 
    bet <- b(u, theta, T, N=100) 
    print("alp")
    print(alp)
    print("bet")
    print(bet)
    print("transform")
    print(transform)
    lambda0 <- theta[4]
    transform <- exp(alp + bet * lambda0)
    distr <- fft(transform, inverse=TRUE)
    distr <- distr/length(distr) #We divide to normalize the transform
    distr <- Mod(distr)
  } else if (method == "derivative") {
    if (method != "N") {warning("Using the derivative method for the distribution of L")}
  }
  print("distributoon in get_distr")
  print(distr)
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
  
  y <- matrix(0, ncol=N+1, nrow=length(u))
  res <- matrix(0, ncol=N+1, nrow=length(u))
  
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
