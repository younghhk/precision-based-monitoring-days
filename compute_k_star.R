# ==============================================================
# Precision-Based Monitoring Day Selection
# Hong & Matthews (2026)
# A Precision-Based Approach to Determining Accelerometer
# Monitoring Duration for Population Surveillance
# ==============================================================


#' Compute Required Number of Monitoring Days (k*)
#'
#' Implements the precision-based decision rule described in:
#' Hong & Matthews (2026), "A Precision-Based Approach to Determining
#' Accelerometer Monitoring Duration for Population Surveillance".
#'
#' The function selects the smallest integer number of monitoring days (k*)
#' such that the half-width of a nominal 95% confidence interval for the
#' population mean does not exceed a prespecified tolerance h.
#'
#' Variance formula:
#'   Var(mu_hat) = (1 / N) * (sigma_b2 + sigma_w2 / k)
#'
#' When monitoring days vary across participants, the harmonic mean
#' of observed days is used to evaluate the observed design.
#'
#' @param sigma_b2 Estimated between-person variance (σ_b^2).
#' @param sigma_w2 Estimated within-person variance (σ_w^2).
#' @param N Number of participants contributing at least one observed day.
#' @param k_obs Numeric vector of observed monitoring days per participant.
#' @param h Precision tolerance (maximum acceptable half-width).
#' @param k_max Maximum candidate number of monitoring days to evaluate.
#'        Default is 14.
#'
#' @return An object of class "precision_k_star" containing:
#'   \describe{
#'     \item{k_harmonic}{Harmonic mean of observed monitoring days.}
#'     \item{half_width_observed}{Achieved half-width under observed design.}
#'     \item{meets_precision}{Logical; does observed design meet tolerance?}
#'     \item{required_k}{Smallest integer k satisfying precision rule.}
#'     \item{tolerance}{Specified precision tolerance (h).}
#'   }
#'
#' @export
compute_k_star <- function(
    sigma_b2,
    sigma_w2,
    N,
    k_obs,
    h,
    k_max = 14
) {
  
  # ---------------------------
  # Input validation
  # ---------------------------
  if (any(c(sigma_b2, sigma_w2, N, h) <= 0)) {
    stop("sigma_b2, sigma_w2, N, and h must be positive.")
  }
  
  if (!is.numeric(k_obs) || any(k_obs <= 0)) {
    stop("k_obs must be a numeric vector of positive values.")
  }
  
  if (!is.numeric(k_max) || k_max < 1) {
    stop("k_max must be a positive integer.")
  }
  
  # ---------------------------
  # 95% confidence multiplier
  # ---------------------------
  z <- 1.96
  
  # ---------------------------
  # Harmonic mean of observed days
  # ---------------------------
  k_harm <- 1 / mean(1 / k_obs)
  
  # ---------------------------
  # Precision under observed design
  # ---------------------------
  var_obs <- (sigma_b2 + sigma_w2 / k_harm) / N
  half_width_obs <- z * sqrt(var_obs)
  
  # ---------------------------
  # Evaluate candidate k values
  # ---------------------------
  k_candidates <- seq_len(k_max)
  
  half_widths <- z * sqrt(
    (sigma_b2 + sigma_w2 / k_candidates) / N
  )
  
  valid_k <- k_candidates[half_widths <= h]
  
  required_k <- if (length(valid_k) == 0) NA_integer_ else min(valid_k)
  
  # ---------------------------
  # Output object
  # ---------------------------
  result <- list(
    k_harmonic = k_harm,
    half_width_observed = half_width_obs,
    meets_precision = half_width_obs <= h,
    required_k = required_k,
    tolerance = h
  )
  
  class(result) <- "precision_k_star"
  
  return(result)
}


# ==============================================================
# Print
# ==============================================================

#' @export
print.precision_k_star <- function(x, ...) {
  
  cat("\nPrecision-Based Monitoring Day Evaluation\n")
  cat("------------------------------------------------\n")
  cat(sprintf("Harmonic mean of observed days: %.2f\n", x$k_harmonic))
  cat(sprintf("Observed half-width: %.6f\n", x$half_width_observed))
  cat(sprintf("Tolerance (h): %.6f\n", x$tolerance))
  cat(sprintf("Meets precision target: %s\n",
              ifelse(x$meets_precision, "Yes", "No")))
  
  if (is.na(x$required_k)) {
    cat("Required number of days (k*): Not achievable within k_max\n")
  } else {
    cat(sprintf("Required number of days (k*): %d\n", x$required_k))
  }
  
  invisible(x)
}
