This repository provides an R implementation of a **precision-based design criterion** for selecting the number of monitoring days required to estimate population-level outcomes with a prespecified level of accuracy. The method is described in detail in the accompanying manuscript.

The goal is to support **study design and planning** for accelerometer-based (and similar) studies by translating empirically estimated variability into guidance on how many days of monitoring are needed.

---

## Overview

The approach focuses on the precision of the estimated **population mean**. Specifically, it evaluates how the variance of the estimated mean depends on:

- between-person variability,
- within-person (day-to-day) variability,
- the number of participants, and
- the number of observed days per participant.

Under a simple variance-components framework, the variance of the estimated population mean decreases as more days are averaged per participant. The required number of monitoring days is defined as the smallest value that achieves a user-specified precision target.

When the number of observed days varies across participants, the method uses the **harmonic mean** of observed days to account for unequal contributions.

---

## Function

### `compute_k_star()`

This function evaluates whether the observed monitoring-day structure achieves a desired level of precision and returns the effective number of days implied by the data.

#### Arguments

- `sigma_b2`  
  Estimated between-person variance.

- `sigma_w2`  
  Estimated within-person (day-to-day) variance.

- `N`  
  Number of participants contributing at least one observed day.

- `k_obs`  
  A numeric vector giving the number of observed days for each participant.  
  When days vary across participants, the harmonic mean of this vector is used.

- `h`  
  Precision tolerance, defined as the maximum acceptable half-width of a nominal 95% confidence interval for the population mean.

- `alpha` (optional, default = 0.05)  
  Significance level used to define the confidence interval.

---

## Example Usage

```r
# Example inputs
sigma_b2 <- 0.60     # estimated between-person variance
sigma_w2 <- 0.40     # estimated within-person variance
N <- 2500            # number of participants
k_obs <- c(5, 6, 7, 7, 6, 4, 5)  # observed days vary by participant
h <- 0.25            # precision tolerance

compute_k_star(
  sigma_b2 = sigma_b2,
  sigma_w2 = sigma_w2,
  N = N,
  k_obs = k_obs,
  h = h
)

