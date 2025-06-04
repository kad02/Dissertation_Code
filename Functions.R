# Pass treatment effects for subs
# Variance, total sample size n, and fraction in S1

# Will need to update with tau for multiple stages

treatment_effect <- function(theta1, theta2, sigma, lambda, n) {
  hat_theta1 <- rnorm(1, mean = theta1, sd = (4 * sigma^2) / (lambda * n))
  hat_theta2 <- rnorm(1, mean = theta2, sd = (4 * sigma^2) / ((1 - lambda) * n))
  theta3 <- lambda * theta1 + (1 - lambda) * theta2
  hat_theta3 <- rnorm(1, mean = theta3, sd = (4 * sigma^2) / n)
  
  return(list(hat_theta1 = hat_theta1, hat_theta2 = hat_theta2, hat_theta3 = hat_theta3))
}











#############################

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
