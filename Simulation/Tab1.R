############################################################
# Script to run Monte Carlo simulations for Table 1
############################################################

# Load required functions
source("precision_functions.R")

# Set random seed for reproducibility
set.seed(123)

############################################################
# 1. Define scenario parameters
############################################################

mu <- 7
sigma2_total <- 2
ICC <- 0.4

sigma_b <- sqrt(ICC * sigma2_total)
sigma_w <- sqrt((1 - ICC) * sigma2_total)

N <- 2000
h <- 0.05
nsim <- 2000


############################################################
# 2. Independence scenario
############################################################

k_indep <- compute_kstar_indep(N, sigma_b, sigma_w, h)

result_indep <- evaluate_design(
  nsim = nsim,
  N = N,
  k = k_indep,
  mu = mu,
  sigma_b = sigma_b,
  sigma_w = sigma_w,
  phi = 0
)


############################################################
# 3. AR(1) misspecified scenario
############################################################

phi_true <- 0.5

result_misspec <- evaluate_design(
  nsim = nsim,
  N = N,
  k = k_indep,
  mu = mu,
  sigma_b = sigma_b,
  sigma_w = sigma_w,
  phi = phi_true
)


############################################################
# 4. AR(1) corrected scenario
############################################################

k_corrected <- compute_kstar_ar1(
  N = N,
  sigma_b = sigma_b,
  sigma_w = sigma_w,
  phi = phi_true,
  h = h
)

result_corrected <- evaluate_design(
  nsim = nsim,
  N = N,
  k = k_corrected,
  mu = mu,
  sigma_b = sigma_b,
  sigma_w = sigma_w,
  phi = phi_true
)


############################################################
# 5. Small-sample scenario
############################################################

N_small <- 500
k_test <- 30

result_smallN <- evaluate_design(
  nsim = nsim,
  N = N_small,
  k = k_test,
  mu = mu,
  sigma_b = sigma_b,
  sigma_w = sigma_w,
  phi = 0
)


############################################################
# 6. Summary table
############################################################

sim_table <- data.frame(
  Scenario = c("Correct model (phi=0)",
               "AR1 true, independence assumed",
               "AR1 corrected",
               "Small N (k=30)"),
  Halfwidth = c(result_indep$emp_halfwidth,
                result_misspec$emp_halfwidth,
                result_corrected$emp_halfwidth,
                result_smallN$emp_halfwidth),
  Coverage = c(result_indep$coverage,
               result_misspec$coverage,
               result_corrected$coverage,
               result_smallN$coverage)
)

print(sim_table)




######################################
#Optional
cat("\n============================\n")
cat("Variance Check: Independence\n")
cat("============================\n")
cat("Empirical variance   :", result_indep$emp_var, "\n")
cat("Theoretical variance :", result_indep$theoretical_var, "\n")
cat("Difference           :", 
    result_indep$emp_var - result_indep$theoretical_var, "\n")
cat("Relative difference  :", 
    (result_indep$emp_var - result_indep$theoretical_var) /
      result_indep$theoretical_var, "\n")
############################################################
# 7. Variance Comparison Table
############################################################

variance_check <- data.frame(
  Scenario = c("Independence",
               "AR1 misspecified",
               "AR1 corrected",
               "Small N"),
  Empirical_Var = c(result_indep$emp_var,
                    result_misspec$emp_var,
                    result_corrected$emp_var,
                    result_smallN$emp_var),
  Theoretical_Var = c(result_indep$theoretical_var,
                      result_misspec$theoretical_var,
                      result_corrected$theoretical_var,
                      result_smallN$theoretical_var)
)

variance_check$Difference <-
  variance_check$Empirical_Var -
  variance_check$Theoretical_Var

variance_check$Relative_Diff <-
  variance_check$Difference /
  variance_check$Theoretical_Var

print(variance_check)
