

StageOne <- function(n, lambda, tau, sigma, theta1, theta2){
  # Calculate theta3
  theta3 <- lambda * theta1 + (1 - lambda) * theta2
  
  # Simulate the trial
  hat_theta1_1 <- rnorm(1, mean = theta1, sd = sqrt((4*sigma^2)/(tau*lambda*n)))
  hat_theta2_1 <- rnorm(1, mean = theta2, sd = sqrt((4*sigma^2)/(tau*(1-lambda)*n)))
  hat_theta3_1 <- rnorm(1, mean = theta3, sd = sqrt((4*sigma^2)/(tau*n)))
  
  return(list(
    hat_theta1_1 = hat_theta1_1,
    hat_theta2_1 = hat_theta2_1,
    hat_theta3_1 = hat_theta3_1
  ))
}

StageOne_Z_P <- function(n, lambda, tau, sigma, hat_theta1_1, hat_theta2_1, hat_theta3_1) {
  # Calculate z-values
 Z1_1 <- (sqrt(lambda*tau*n)*hat_theta1_1)/(2*sigma)
 Z3_1 <- (sqrt(tau*n)*hat_theta3_1)/(2*sigma)
 
  # Calculate p-values
  P1_1 <- 1 - pnorm(Z1_1)
  P3_1 <- 1 - pnorm(Z3_1)
  
  return(list(
    Z1_1 = Z1_1,
    Z3_1 = Z3_1,
    P1_1 = P1_1,
    P3_1 = P3_1
  ))
}

StageTwo_NoEnrich <- function(n, lambda, tau, sigma, theta1, theta2) {
  # Calculate theta3
  theta3 <- lambda * theta1 + (1 - lambda) * theta2
  
  # Simulate the trial
  hat_theta1_2 <- rnorm(1, mean = theta1, sd = sqrt((4*sigma^2)/((1-tau)*lambda*n)))
  hat_theta2_2 <- rnorm(1, mean = theta2, sd = sqrt((4*sigma^2)/((1-tau)*(1-lambda)*n)))
  hat_theta3_2 <- rnorm(1, mean = theta3, sd = sqrt((4*sigma^2)/((1-tau)*n)))
  
  return(list(
    hat_theta1_2 = hat_theta1_2,
    hat_theta2_2 = hat_theta2_2,
    hat_theta3_2 = hat_theta3_2
  ))
}

StageTwo_NoEnrich_Z_P <- function(n, lambda, tau, sigma, hat_theta1_2, hat_theta2_2, hat_theta3_2) {
  # Calculate z-values
  Z1_2 <- (sqrt((1-tau)*lambda*n)*hat_theta1_2)/(2*sigma)
  Z3_2 <- (sqrt((1-tau)*n)*hat_theta3_2)/(2*sigma)
  
  # Calculate p-values
  P1_2 <- 1 - pnorm(Z1_2)
  P3_2 <- 1 - pnorm(Z3_2)
  
  return(list(
    Z1_2 = Z1_2,
    Z3_2 = Z3_2,
    P1_2 = P1_2,
    P3_2 = P3_2
  ))
}

StageTwo_Enrich <- function(n, tau, sigma, theta1) {
  # Simulate the trial
  hat_theta1_2_enrich <- rnorm(1, mean = theta1, sd = sqrt((4*sigma^2)/((1-tau)*n)))
  
  return(list(
    hat_theta1_2_enrich = hat_theta1_2_enrich
  ))
}

StageTwo_Enrich_Z_P <- function(n, tau, sigma, hat_theta1_2_enrich) {
  # Calculate z-values
  Z1_2_enrich <- (sqrt((1-tau)*n)*hat_theta1_2_enrich)/(2*sigma)
  
  # Calculate p-values
  P1_2_enrich <- 1 - pnorm(Z1_2_enrich)
  
  return(list(
    Z1_2_enrich = Z1_2_enrich,
    P1_2_enrich = P1_2_enrich
  ))
}


# Intersection test
simes_method <- function(p1, p3) {
  p13 <- min(2*min(p1,p3),max(p1, p3))
  
  return(p13)
}


# Combination
win <- function(z_1, z_2, w1, w2){
  # Z value
  zc <- w1*z_1 + w2*z_2
  
  # P value
  pc <- 1 - pnorm(zc)
  
  return(list(zc = zc, pc = pc))
}


weights <- function(Enrich, lambda, tau){
  if (Enrich) {
    w1 <- sqrt((tau*lambda)/(tau*lambda-tau+1))
    w2 <- sqrt((1-tau)/(tau*lambda-tau+1))
  } else {
    w1 <- sqrt(tau)
    w2 <- sqrt(1 - tau)
  }
  
  return(list(w1 = w1, w2 = w2))
}



# Decision
Decision_identity <- function(hat_theta1, hat_theta3){
  enrich <- FALSE
  if (hat_theta1 > hat_theta3) {
    enrich <- TRUE
  } else {
    enrich <- FALSE
  }
  return(enrich)
}




