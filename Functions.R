# Simple 2 arm trial
# Common variance
# 2 subpopulations

Trial_Std  <- function(n, lambda, mu_A_S1, mu_A_S2, mu_B_S1, mu_B_S2, sd) {
  # n: total sample size
  n1 <- n * lambda
  n2 <- n * (1 - lambda)
  
  # Realisations for subpopulation 1
  A_S1 <- rnorm(n1/2, mean = mu_A_S1, sd = sd)
  B_S1 <- rnorm(n1/2, mean = mu_B_S1, sd = sd)
  
  # Realisations for subpopulation 2
  A_S2 <- rnorm(n2/2, mean = mu_A_S2, sd = sd)
  B_S2 <- rnorm(n2/2, mean = mu_B_S2, sd = sd)
  
  return(list(A_S1 = A_S1, B_S1 = B_S1, A_S2 = A_S2, B_S2 = B_S2))
}



# Calculate trt effects
calculate_thetas <- function(A_S1, B_S1, A_S2, B_S2, lambda) {
  # Means for subpopulation 1
  mu_hat_A_S1 <- mean(A_S1)
  mu_hat_B_S1 <- mean(B_S1)
  
  # Means for subpopulation 2
  mu_hat_A_S2 <- mean(A_S2)
  mu_hat_B_S2 <- mean(B_S2)
  
  
  # Treatment effects
  theta1 <- mu_hat_A_S1 - mu_hat_B_S1
  theta2 <- mu_hat_A_S2 - mu_hat_B_S2
  
  theta3 <- lambda * theta1 + (1 - lambda) * theta2
  
  return(list(theta1 = theta1, theta2 = theta2, theta3 = theta3))
}


# Calculate Z and P values
calculate_zvals <- function(theta1, theta3, n, lambda, sd) {
  z1 <- (sqrt(lambda * n) * theta1) / (2 * sd)
  z3 <- (sqrt(n) * theta3) / (2 * sd)
  
  return(list(z1 = z1, z3 = z3))
}

calculate_p_vals <- function(z1, z3) {
  p1 <- 1 - pnorm(z1)
  p3 <- 1 - pnorm(z3)
  
  return(list(p1 = p1, p3 = p3))
}


# Simes method

simes_method <- function(p1, p3) {
  p13 <- min(2*min(p1,p3),max(p1, p3))
  
  return(p13)
}

# Weighted inverse normal method

win <- function(z_1, z_2){
  # CHANGE WEIGHTS LATER
  w1 <- 1/sqrt(2)
  w2 <- 1/sqrt(2)
  
  zc <- w1*z_1 + w2*z_2
  
  pc <- 1 - pnorm(zc)
  
  return(list(zc = zc, pc = pc))
}
