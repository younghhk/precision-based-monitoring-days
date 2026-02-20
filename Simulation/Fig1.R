############################################################
# Script to generate design-space results for Figure 1
############################################################

############################################################
# 0. Setup
############################################################

rm(list = ls())

source("precision_functions.R")
library(ggplot2)

set.seed(123)

# Precision settings
h <- 0.05
alpha <- 0.05

# Total variance used for design exploration
sigma2_total <- 2


############################################################
# 1. Illustrative Scenario Demonstration
############################################################

cat("\n============================================\n")
cat("Illustrative Scenario Results (ICC = 0.5)\n")
cat("============================================\n")

illustrative_scenarios <- list(
  small_N  = list(N = 500,  ICC = 0.5),
  medium_N = list(N = 2000, ICC = 0.5),
  large_N  = list(N = 5000, ICC = 0.5)
)

for (scenario_name in names(illustrative_scenarios)) {
  
  s <- illustrative_scenarios[[scenario_name]]
  
  sigma_b <- sqrt(s$ICC * sigma2_total)
  sigma_w <- sqrt((1 - s$ICC) * sigma2_total)
  
  kstar <- compute_kstar_indep(
    N = s$N,
    sigma_b = sigma_b,
    sigma_w = sigma_w,
    h = h,
    alpha = alpha
  )
  
  if (is.infinite(kstar)) {
    cat(scenario_name, ": Precision target NOT achievable (increase N)\n")
  } else {
    cat(scenario_name, ": k* =", kstar, "\n")
  }
}


############################################################
# 2. Full Design Grid (ICC × N × phi)
############################################################

ICC_vals <- seq(0.1, 0.9, by = 0.1)
N_vals   <- c(500, 1000, 2000, 5000, 10000)
phi_vals <- c(0, 0.5)

grid_list <- list()
counter <- 1

for (phi in phi_vals) {
  for (N in N_vals) {
    for (ICC in ICC_vals) {
      
      sigma_b <- sqrt(ICC * sigma2_total)
      sigma_w <- sqrt((1 - ICC) * sigma2_total)
      
      if (phi == 0) {
        kstar <- compute_kstar_indep(
          N = N,
          sigma_b = sigma_b,
          sigma_w = sigma_w,
          h = h,
          alpha = alpha
        )
      } else {
        kstar <- compute_kstar_ar1(
          N = N,
          sigma_b = sigma_b,
          sigma_w = sigma_w,
          phi = phi,
          h = h,
          alpha = alpha
        )
      }
      
      grid_list[[counter]] <- data.frame(
        ICC      = ICC,
        N        = N,
        phi      = phi,
        kstar    = ifelse(is.infinite(kstar), NA, kstar),
        feasible = !is.infinite(kstar)
      )
      
      counter <- counter + 1
    }
  }
}

grid_results <- do.call(rbind, grid_list)


############################################################
# 3. Print Results Separately for phi = 0 and phi = 0.5
############################################################

cat("\n============================================\n")
cat("Design Grid: Conditional Independence (phi = 0)\n")
cat("============================================\n")
print(subset(grid_results, phi == 0))

cat("\n============================================\n")
cat("Design Grid: AR(1) Serial Correlation (phi = 0.5)\n")
cat("============================================\n")
print(subset(grid_results, phi == 0.5))



############################################################
# 4. Generate Publication Figure
############################################################

pdf("kstar_vs_ICC_vs_N.pdf", width = 8, height = 6)

ggplot(grid_results,
       aes(x = ICC,
           y = kstar,
           color = factor(N),
           group = factor(N))) +
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ phi,
             labeller = labeller(
               phi = c(
                 "0"   = "Conditional independence (phi = 0)",
                 "0.5" = "AR(1) serial correlation (phi = 0.5)"
               )
             )) +
  labs(
    x = "Intraclass Correlation (ICC, rho)",
    y = expression(k^"* (Required Monitoring Days)"),
    color = "Sample Size (N)"
  ) +
  theme_minimal(base_size = 14)

dev.off()

cat("\nFigure saved as: kstar_vs_ICC_vs_N.pdf\n")

