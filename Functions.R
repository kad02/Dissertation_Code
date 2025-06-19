
# STAGE 1
StageOne <- function(n_sim, n, lambda, tau, sigma, theta1, theta2){
  # Calculate theta3
  theta3 <- lambda * theta1 + (1 - lambda) * theta2
  
  # Simulate the trial
  hat_theta1_1 <- rnorm(n_sim, mean = theta1, sd = sqrt((4*sigma^2)/(tau*lambda*n)))
  hat_theta2_1 <- rnorm(n_sim, mean = theta2, sd = sqrt((4*sigma^2)/(tau*(1-lambda)*n)))
  #hat_theta3_1 <- rnorm(n_sim, mean = theta3, sd = sqrt((4*sigma^2)/(tau*n)))
  
  hat_theta3_1 <- lambda * hat_theta1_1 + (1 - lambda) * hat_theta2_1
  
  dat <- data.frame(
    hat_theta1_1 = hat_theta1_1,
    hat_theta2_1 = hat_theta2_1,
    hat_theta3_1 = hat_theta3_1
  ) %>%
    mutate(
      z1_1 = (sqrt(lambda*tau*n)*hat_theta1_1)/(2*sigma),
      z3_1 = (sqrt(tau*n)*hat_theta3_1)/(2*sigma),
      p1_1 = 1 - pnorm(z1_1),
      p3_1 = 1 - pnorm(z3_1),
      p_13_1 = simes_method(p1_1, p3_1),
      enrich = Decision_identity(hat_theta1_1, hat_theta3_1)
    )%>%
    arrange(
      enrich
    )
  
  return(dat)
}


# STAGE 2
StageTwo_NoEnrich <- function(n_not, n, lambda, tau, sigma, theta1, theta2) {
  # Calculate theta3
  theta3 <- lambda * theta1 + (1 - lambda) * theta2
  
  #hat_theta3_2 = rnorm(n_not, mean = theta3, sd = sqrt((4*sigma^2)/((1-tau)*n))),
  
  # Simulate the trial
  no_en <- data.frame(
    hat_theta1_2 = rnorm(n_not, mean = theta1, sd = sqrt((4*sigma^2)/((1-tau)*lambda*n))),
    hat_theta2_2 = rnorm(n_not, mean = theta2, sd = sqrt((4*sigma^2)/((1-tau)*(1-lambda)*n))),
    hat_theta3_2 = lambda * hat_theta1_2 + (1 - lambda) * hat_theta2_2
  ) %>%
    mutate(
      z1_2 = (sqrt((1-tau)*lambda*n)*hat_theta1_2)/(2*sigma),
      z3_2 = (sqrt((1-tau)*n)*hat_theta3_2)/(2*sigma),
      p1_2 = 1 - pnorm(z1_2),
      p3_2 = 1 - pnorm(z3_2),
      p_13_2 = simes_method(p1_2, p3_2)
    )
  return(no_en)
}


StageTwo_Enrich <- function(n_en, n, tau, sigma, theta1) {
  # Simulate the trial
  en <- data.frame(
    hat_theta1_2 = rnorm(n_en, mean = theta1, sd = sqrt((4*sigma^2)/((1-tau)*n)))
  ) %>%
    mutate(
      z1_2 = (sqrt((1-tau)*n)*hat_theta1_2)/(2*sigma),
      z3_2 = -Inf,
      p1_2 = 1 - pnorm(z1_2),
      p3_2 = 1,  # p3 is not used in this case
      p_13_2 = p1_2
    )
  return(en)
}




# Intersection test
simes_method <- function(p1, p3) {
  p13 <- pmin(2*pmin(p1,p3),pmax(p1, p3))
  
  return(p13)
}


# Combination
win <- function(z_1, z_2, w1, w2){
  # Z value
  zc <- w1*z_1 + w2*z_2
  
  return(zc)
}


win_p <- function(p_1, p_2, w1, w2){
  # Z value
  zc <- w1*qnorm(1 - p_1) + w2*qnorm(1 - p_2)
  return(zc)
}


weights <- function(lambda, tau){
  #if (Enrich) {
  #  w1 <- sqrt((tau*lambda)/(tau*lambda-tau+1))
  #  w2 <- sqrt((1-tau)/(tau*lambda-tau+1))
  #} else {
  #  w1 <- sqrt(tau)
  #  w2 <- sqrt(1 - tau)
  #}
  
  w1 <- sqrt(tau)
  w2 <- sqrt(1 - tau)
  
  return(list(w1 = w1, w2 = w2))
}



# Decision
Decision_identity <- function(hat_theta1, hat_theta3){
  enrich <- hat_theta1 > hat_theta3
  return(enrich)
}


# Combination stuff
Comb <- function(data, lambda, tau){
  ws <- weights(lambda, tau)
  
  data <- data %>%
    mutate(
      z1_c = win(z1_1, z1_2, ws$w1, ws$w2),
      z3_c = win(z3_1, z3_2, ws$w1, ws$w2),
      z13_c = win_p(p_13_1, p_13_2, ws$w1, ws$w2),
      p1_c = 1 - pnorm(z1_c),
      p3_c = 1 - pnorm(z3_c),
      p13_c = 1 - pnorm(z13_c)
    )
  return(data)
}





Trial <- function(n_sim, n, lambda, tau, sigma, theta1, theta2){
  stage1 <- StageOne(n_sim, n, lambda, tau, sigma, theta1, theta2)
  
  n_en <- sum(stage1$enrich)
  n_not <- n_sim - n_en
  
  stage2_no_enrich <- StageTwo_NoEnrich(n_not, n, lambda, tau, sigma, theta1, theta2)
  stage2_enrich <- StageTwo_Enrich(n_en, n, tau, sigma, theta1)
  
  stage2 <- bind_rows(stage2_no_enrich, stage2_enrich)
  
  dat_full <- bind_cols(stage1, stage2)
  
  dat_full_combined <- Comb(dat_full, lambda, tau)
  
  return(dat_full_combined)
}
