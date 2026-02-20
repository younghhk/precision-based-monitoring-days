############################################################
# precision_functions.R
#
# Core simulation functions for evaluating the precision-based design rule
# and generating simulation results reported in Figure 1 and Table 1.
#   1. compute_kstar_indep()
#   2. compute_kstar_ar1()
#   3. evaluate_design()
############################################################


############################################################
# 1. Compute k* under independence
############################################################
compute_kstar_indep <- function(N, sigma_b, sigma_w,
                                h, alpha = 0.05,
                                k_max = 200) {

  z <- qnorm(1 - alpha/2)

  # Lower bound when k -> infinity
  min_var <- sigma_b^2 / N

  # If infinite monitoring cannot meet precision
  if (z * sqrt(min_var) > h) {
    return(Inf)
  }

  for (k in 1:k_max) {

    var_mu <- sigma_b^2 / N +
              sigma_w^2 / (N * k)

    if (z * sqrt(var_mu) <= h) {
      return(k)
    }
  }

  return(NA)
}


############################################################
# 2. Compute k* under AR(1)
############################################################
compute_kstar_ar1 <- function(N, sigma_b, sigma_w,
                              phi,
                              h,
                              alpha = 0.05,
                              k_max = 200) {

  z <- qnorm(1 - alpha/2)

  min_var <- sigma_b^2 / N

  if (z * sqrt(min_var) > h) {
    return(Inf)
  }

  for (k in 1:k_max) {

    # AR(1) correlation matrix
    R <- outer(1:k, 1:k, function(i, j) phi^abs(i - j))

    # Design effect
    DE_k <- sum(R) / k

    # Effective number of days
    k_eff <- k / DE_k

    var_mu <- sigma_b^2 / N +
              sigma_w^2 / (N * k_eff)

    if (z * sqrt(var_mu) <= h) {
      return(k)
    }
  }

  return(NA)
}


############################################################
# 3. Monte Carlo evaluation function
############################################################
evaluate_design <- function(nsim,
                            N,
                            k,
                            mu,
                            sigma_b,
                            sigma_w,
                            phi = 0,
                            alpha = 0.05) {
  
  # Safety check
  if (is.infinite(k) || is.na(k)) {
    stop("k is Inf or NA — precision target cannot be achieved.")
  }
  
  z <- qnorm(1 - alpha/2)
  est <- numeric(nsim)
  
  ##########################################################
  # Monte Carlo simulation
  ##########################################################
  
  for (s in 1:nsim) {
    
    Ybar <- numeric(N)
    
    for (i in 1:N) {
      
      b_i <- rnorm(1, 0, sigma_b)
      
      if (phi == 0) {
        eps <- rnorm(k, 0, sigma_w)
      } else {
        eps <- as.numeric(
          arima.sim(model = list(ar = phi),
                    n = k,
                    sd = sigma_w)
        )
      }
      
      Y_i <- mu + b_i + eps
      Ybar[i] <- mean(Y_i)
    }
    
    est[s] <- mean(Ybar)
  }
  
  ##########################################################
  # Empirical quantities
  ##########################################################
  
  emp_var <- var(est)
  emp_sd  <- sd(est)
  emp_halfwidth <- z * emp_sd
  
  ##########################################################
  # Theoretical variance (analytic CI formula)
  ##########################################################
  
  if (phi == 0) {
    
    theoretical_var <- sigma_b^2 / N +
      sigma_w^2 / (N * k)
    
  } else {
    
    R <- outer(1:k, 1:k, function(i, j) phi^abs(i - j))
    DE_k <- sum(R) / k
    k_eff <- k / DE_k
    
    theoretical_var <- sigma_b^2 / N +
      sigma_w^2 / (N * k_eff)
  }
  
  theoretical_sd <- sqrt(theoretical_var)
  theoretical_halfwidth <- z * theoretical_sd
  
  ##########################################################
  # Coverage using analytic variance
  ##########################################################
  
  lower <- est - z * theoretical_sd
  upper <- est + z * theoretical_sd
  coverage <- mean(lower <= mu & mu <= upper)
  
  return(list(
    emp_sd = emp_sd,
    emp_var = emp_var,
    emp_halfwidth = emp_halfwidth,
    theoretical_sd = theoretical_sd,
    theoretical_var = theoretical_var,
    theoretical_halfwidth = theoretical_halfwidth,
    coverage = coverage
  ))
}
