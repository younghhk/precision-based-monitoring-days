# Precision-Based Monitoring Days

This repository provides an R implementation of a **precision-based design criterion** for selecting the number of monitoring days required to estimate population-level outcomes with a prespecified level of accuracy. The method is described in detail in the accompanying manuscript.

The goal is to support **study design and planning** for accelerometer-based (and similar repeated-measures) studies by translating empirically estimated variability into guidance on how many days of monitoring are needed.

---

## Overview

The approach focuses on the precision of the estimated **population mean**. Specifically, it evaluates how the variance of the estimated mean depends on:

* between-person variability,
* within-person (day-to-day) variability,
* the number of participants, and
* the number of observed days per participant.

Under a variance-components framework, the variance of the estimated population mean decreases as more days are averaged per participant. The required number of monitoring days is defined as the smallest value that achieves a user-specified precision target (i.e., a maximum acceptable half-width of a 95% confidence interval).

### Key Features

* Works for **continuous outcomes** (linear mixed-effects models via `lmer`)
* Works for **binary outcomes** (logistic mixed-effects models via `glmer`)
* Handles unequal monitoring days across participants (via harmonic mean logic in the manuscript)
* Optional sensitivity adjustment for **within-person lag-1 autocorrelation** in consecutive monitoring days

---

## Primary Function

### `compute_k_star()`

This function computes the required number of monitoring days (`k*`) based on a fitted mixed-effects model and a user-specified precision tolerance.

### Arguments

* `model`
  A fitted mixed-effects model (`lmerMod` for continuous outcomes or `glmerMod` for binary outcomes).

* `outcome_type`
  `"continuous"` or `"binary"`.

* `N`
  Number of participants contributing at least one observed day.

* `h`
  Precision tolerance. Defined as the maximum acceptable half-width of a nominal 95% confidence interval for the population mean.

* `id_var`
  Name of the participant ID variable (character string).

* `day_var`
  Name of the day variable used to determine ordering (character string).

* `k_max` (default = 14)
  Maximum number of monitoring days to evaluate.

* `adjust_autocorrelation` (default = `FALSE`)
  If `TRUE`, applies an AR(1)-based effective-days adjustment using estimated lag-1 residual autocorrelation.
  This is intended as a **sensitivity analysis**.

* `n_sim` (binary models only; default = 5000)
  Number of Monte Carlo samples used to estimate marginal prevalence.

---

## Modeling Assumptions

The primary method assumes:

* A random-intercept mixed-effects model
* Conditional independence of daily residuals after accounting for:

  * participant-level random effects
  * fixed effects (e.g., season, day-of-week)

Daily observations are **not assumed independent** overall — correlation across days arises naturally through the participant random intercept.

If consecutive monitoring days exhibit additional short-term persistence (e.g., lag-1 autocorrelation), the optional adjustment estimates an effective number of monitoring days:

```
k_eff = k / design_effect
```

In short monitoring windows (e.g., 7 consecutive days), this adjustment typically has modest impact but provides a principled sensitivity analysis.

---

## Example Usage

### Example 1: Continuous Outcome (Primary Method)

```r
fit_cont <- lmer(
  sleep_duration ~ season + (1 | id),
  data = mydata1
)

result_cont <- compute_k_star(
  model = fit_cont,
  outcome_type = "continuous",
  N = 1000,
  h = 0.1,
  id_var = "id",
  day_var = "day"
)


result_cont

Precision-Based Monitoring Evaluation
---------------------------------------
Outcome type: continuous 
Sample size (N): 1000 
Tolerance (h): 0.1 
Required monitoring days (k*): 3 
```

---

### Example 2: Continuous Outcome with Autocorrelation Sensitivity

```r
fit_cont <- lmer(
  sleep_duration ~ season + (1 | id),
  data = mydata2
)

result_cont <- compute_k_star(
  model = fit_cont2,
  outcome_type = "continuous",
  N = 1000,
  h = 0.1,
  id_var = "id",
  day_var = "day",
  data = mydata2,
  adjust_autocorrelation = TRUE
)

result_cont

Precision-Based Monitoring Evaluation
---------------------------------------
Outcome type: continuous 
Sample size (N): 1000 
Tolerance (h): 0.1 
Lag-1 residual correlation: 0.127 
Required monitoring days (k*): 4 
```

---

### Example 3: Binary Outcome (Primary Method)

```r
fit_bin <- glmer(
  slept_7h ~ season + (1 | id),
  data = mydata3,
  family = binomial
)

result_bin <- compute_k_star(
  model = fit_bin,
  outcome_type = "binary",
  N = length(unique(mydata3$id)),
  h = 0.02,
  id_var = "id",
  day_var = "day"
)

result_bin

Precision-Based Monitoring Evaluation
---------------------------------------
Outcome type: binary 
Sample size (N): 1000 
Tolerance (h): 0.02 
Required monitoring days (k*): 3 
```

---

## Interpretation of Output

The function returns:

* Estimated variance components (continuous) or prevalence/ICC (binary)
* Optional lag-1 autocorrelation estimate (if requested)
* Effective number of monitoring days (if adjusted)
* Required number of monitoring days (`k*`)

If no value ≤ `k_max` satisfies the precision criterion, `k*` is returned as `NA`.

---

## Supplementary Materials

The accompanying manuscript provides full derivations for:

* Continuous outcomes under a variance-components model
* Binary outcomes under a logistic mixed-effects model
* Harmonic mean adjustment for unequal monitoring days
* Effective-days adjustment for serial correlation

---

## Citation

If you use this code, please cite:

> Hong, H.G. & Matthews, C. (2026).
> *A Precision-Based Approach to Determining Accelerometer Monitoring Duration for Population Surveillance.*

