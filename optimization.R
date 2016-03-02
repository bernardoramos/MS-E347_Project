library(stats)



m<-5
mid <- rep(0.7,5)

objective <- function(theta){
  # mid must be a data frame where rows correspond to the different maturities, and whose first five columns correspond to the spreads of the tranches in increasing order.
  tot <- 0
  for (j in 1:nrow(mid)) {
    for (i in 1:5){
      tran <- tranche[[i]]
      T <- mid$Maturity[j]
      model_value <- model(tranche=tranche[[i]], theta, T, method=method)
      tot <- tot + (mid[j, i]-model_value)^2/mid[j, i]
    }
  }
  print(c("Objective:", tot))
  return(tot)
}


AAPE <- function(theta){
  # mid must be a data frame where rows correspond to the different maturities, and whose first five columns correspond to the spreads of the tranches in increasing order.
  tot <- 0
  for (j in 1:nrow(mid)) {
    for (i in 1:5){
      tran <- tranche[[i]]
      T <- mid$Maturity[j]
      model_value <- model(tranche=tranche[[i]], theta, T, method=method)
      tot <- tot + abs(mid[j, i]-model_value)/mid[j, i] # Usual square objective
    }
  }
  return(tot)
}
