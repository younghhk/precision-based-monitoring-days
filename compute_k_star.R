# ==============================================================
# Precision-Based Monitoring Days Calculator
# Direct computation with known parameters
# ==============================================================

# --------------------------------------------------------------
# Helper: AR(1) design effect
# --------------------------------------------------------------

design_effect_ar1 <- function(k, phi) {
  if (k <= 1) return(1)
  h <- seq_len(k - 1)
  1 + 2 * sum((1 - h / k) * phi^h)
}

k_eff_ar1 <- function(k, phi) {
  k / design_effect_ar1(k, phi)
}

# --------------------------------------------------------------
# Main Function for Continuous Outcomes
# --------------------------------------------------------------

compute_k_star_continuous <- function(
    N,                    # Sample size
    h,                    # Tolerance (desired margin of error)
    sigma_b2,            # Between-person variance
    sigma_w2,            # Within-person variance
    phi = 0,             # Lag-1 autocorrelation (default = 0)
    k_max = 14,          # Maximum days to check
    z = 1.96             # Z-value for CI (default 95% CI)
) {
  
  k_vals <- seq_len(k_max)
  
  # Adjust for autocorrelation if phi != 0
  k_eff <- if (phi != 0) {
    sapply(k_vals, k_eff_ar1, phi = phi)
  } else {
    k_vals
  }
  
  # Compute margin of error for each k
  margin_of_error <- z * sqrt((sigma_b2 + sigma_w2 / k_eff) / N)
  
  # Find minimum k that meets tolerance
  valid_k <- k_vals[margin_of_error <= h]
  k_star  <- if (length(valid_k)) min(valid_k) else NA_integer_
  
  # Create results table
  results_table <- data.frame(
    k = k_vals,
    margin_of_error = margin_of_error,
    meets_tolerance = margin_of_error <= h
  )
  
  # Print summary
  cat("\n=== Precision-Based Monitoring Days (Continuous Outcome) ===\n")
  cat("Sample size (N):", N, "\n")
  cat("Desired precision (h):", h, "\n")
  cat("Between-person variance (σ²_b):", sigma_b2, "\n")
  cat("Within-person variance (σ²_w):", sigma_w2, "\n")
  if (phi != 0) cat("Lag-1 autocorrelation (φ):", phi, "\n")
  cat("\nRequired monitoring days (k*):", 
      ifelse(is.na(k_star), "Not attainable within k_max", k_star), "\n")
  cat("\nNote: margin_of_error is the ±CI half-width achieved at each k.\n")
  cat("      k* is the minimum k where margin_of_error ≤ h.\n\n")
  
  print(results_table)
  
  invisible(list(
    k_star = k_star,
    outcome_type = "continuous",
    N = N,
    h = h,
    sigma_b2 = sigma_b2,
    sigma_w2 = sigma_w2,
    phi = phi,
    results = results_table
  ))
}

# --------------------------------------------------------------
# Main Function for Binary Outcomes (Updated)
# --------------------------------------------------------------

compute_k_star_binary <- function(
    N,                    # Sample size
    h,                    # Tolerance (desired margin of error)
    mu,                   # Marginal prevalence/mean
    rho,                  # Intraclass correlation
    k_max = 14,           # Maximum days to check
    z = 1.96              # Z-value for CI (default 95% CI)
) {
  
  k_vals <- seq_len(k_max)
  
  # Compute margin of error for each k 
  margin_of_error <- z * sqrt(
    (mu * (1 - mu) / N) * (rho + (1 - rho) / k_vals)
  )
  
  # Find minimum k that meets tolerance
  valid_k <- k_vals[margin_of_error <= h]
  k_star  <- if (length(valid_k)) min(valid_k) else NA_integer_
  
  # Create results table
  results_table <- data.frame(
    k = k_vals,
    margin_of_error = margin_of_error,
    meets_tolerance = margin_of_error <= h
  )
  
  # Print summary
  cat("\n=== Precision-Based Monitoring Days (Binary Outcome) ===\n")
  cat("Sample size (N):", N, "\n")
  cat("Desired precision (h):", h, "\n")
  cat("Marginal prevalence (μ):", mu, "\n")
  cat("Intraclass correlation (ρ):", rho, "\n")
  cat("\nRequired monitoring days (k*):", 
      ifelse(is.na(k_star), "Not attainable within k_max", k_star), "\n")
  cat("\nNote: margin_of_error is the ±CI half-width achieved at each k.\n")
  cat("      k* is the minimum k where margin_of_error ≤ h.\n\n")
  
  print(results_table)
  
  invisible(list(
    k_star = k_star,
    outcome_type = "binary",
    N = N,
    h = h,
    mu = mu,
    rho = rho,
    results = results_table
  ))
}


