
This repository provides an R implementation of a **precision-based design criterion** for selecting the number of monitoring days required to estimate population-level outcomes with a prespecified level of accuracy. The method is described in detail in the accompanying manuscript.

The primary goal is to support **study design and planning** for accelerometer-based (and similar) studies by translating empirically estimated variability into guidance on how many days of monitoring are needed.

---

## Overview

The approach targets the precision of the estimated **population mean** by evaluating the variance
\[
\mathrm{Var}(\hat{\mu}) = \frac{1}{N}\left(\sigma_b^2 + \frac{\sigma_w^2}{k}\right),
\]
where:
- \(\sigma_b^2\) is the between-person variance,
- \(\sigma_w^2\) is the within-person (day-to-day) variance,
- \(N\) is the number of participants, and
- \(k\) is the number of observed days per participant.

When the number of observed days varies across participants, the method uses the **harmonic mean** of observed days to account for unequal contributions.

The required number of monitoring days, denoted \(k^\star\), is defined as the smallest value that satisfies a prespecified precision criterion:
\[
1.96 \sqrt{\mathrm{Var}(\hat{\mu})} \le h,
\]
where \(h\) is a user-defined tolerance for the half-width of a nominal 95% confidence interval.

---

## Function

### `compute_k_star()`

This function evaluates whether a given set of observed monitoring days achieves the desired precision and returns the effective number of days implied by the data.

#### Arguments

- `sigma_b2`  
  Estimated between-person variance (\(\sigma_b^2\)).

- `sigma_w2`  
  Estimated within-person variance (\(\sigma_w^2\)).

- `N`  
  Number of participants contributing at least one observed day.

- `k_obs`  
  A vector giving the number of observed days for each participant.  
  When days vary across participants, the harmonic mean of `k_obs` is used.

- `h`  
  Precision tolerance (half-width of a nominal 95% confidence interval).

- `alpha` (optional, default = 0.05)  
  Significance level for the confidence interval.

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
