# Precision-Based Monitoring Days


This repository provides an R implementation of a precision-based rule for determining the number of monitoring days required per participant to achieve a desired level of accuracy when estimating a population mean or prevalence.

Using variance components estimated from existing data (e.g., pilot studies or large observational cohorts), the method addresses a practical design question:

Given the observed between-person and day-to-day variability, how many monitoring days are needed per participant so that the 95% confidence interval for the population estimate is sufficiently narrow?

This tool is intended to support study design and protocol decisions in accelerometer-based and other repeated-measures research settings.

---

## Overview

The method focuses on the precision of the estimated **population mean** (or prevalence for binary outcomes). Specifically, it determines the smallest number of monitoring days per participant required to achieve a user-specified confidence interval half-width.

Under a variance-components framework, the variance of the estimated population mean depends on:

* Between-person variability
* Within-person (day-to-day) variability
* Total sample size (number of participants)
* Number of monitoring days per participant

The required monitoring duration (k*) is the smallest number of monitoring days per participant needed so that the half-width of a 95% confidence interval does not exceed a user-specified tolerance (h).

---



## Modeling Framework

The primary method assumes a random-intercept mixed-effects model:

* Daily observations are conditionally independent given participant-level random effects.
* Correlation across days arises naturally through the participant random intercept.

The design rule assumes a planned study with equal monitoring days per participant. Estimated variance components from existing data are used to determine the required monitoring duration.

### Optional Autocorrelation Adjustment

If short-term persistence remains after accounting for random intercepts (e.g., AR(1) structure), an optional sensitivity analysis estimates an effective number of monitoring days:

```
k_eff = k / design_effect
```

This adjustment is most relevant for consecutive-day monitoring protocols.

---

## Primary Function

### `compute_k_star()`

Computes the required number of monitoring days (`k*`) to achieve a prespecified precision target.

### Arguments

* `model`
  A fitted mixed-effects model (`lmerMod` or `glmerMod`).

* `outcome_type`
  `"continuous"` or `"binary"`.

* `N`
  Number of participants contributing at least one observed day.

* `h`
  Precision tolerance. Defined as the maximum acceptable half-width of a 95% confidence interval.

* `id_var`
  Name of the participant ID variable (character string).

* `day_var`
  Name of the day variable (used for ordering residuals when estimating autocorrelation).

* `k_max` (default = 14)
  Maximum number of monitoring days evaluated.

* `adjust_autocorrelation` (default = `FALSE`)
  If `TRUE`, applies an AR(1)-based effective-days adjustment as a sensitivity analysis.

* `n_sim` (binary models only; default = 5000)
  Number of Monte Carlo draws used to estimate marginal prevalence.

---

## Example Usage

### Example 1: Continuous Outcome

```r
library(lme4)

fit_cont <- lmer(
  sleep_duration ~ season + (1 | id),
  data = mydata
)

result_cont <- compute_k_star(
  model = fit_cont,
  outcome_type = "continuous",
  N = length(unique(mydata$id)),
  h = 0.10,
  id_var = "id",
  day_var = "day"
)

result_cont
```

Example output:

```
Precision-Based Monitoring Evaluation
---------------------------------------
Outcome type: continuous
Sample size (N): 1000
Tolerance (h): 0.10
Required monitoring days (k*): 3
```

---

### Example 2: Continuous Outcome with Autocorrelation Sensitivity

```r
result_cont_adj <- compute_k_star(
  model = fit_cont,
  outcome_type = "continuous",
  N = length(unique(mydata$id)),
  h = 0.10,
  id_var = "id",
  day_var = "day",
  adjust_autocorrelation = TRUE
)

result_cont_adj
```

Example output:

```
Precision-Based Monitoring Evaluation
---------------------------------------
Outcome type: continuous
Sample size (N): 1000
Tolerance (h): 0.10
Lag-1 residual correlation: 0.127
Required monitoring days (k*): 4
```

---

### Example 3: Binary Outcome

```r
fit_bin <- glmer(
  slept_7h ~ season + (1 | id),
  data = mydata,
  family = binomial
)

result_bin <- compute_k_star(
  model = fit_bin,
  outcome_type = "binary",
  N = length(unique(mydata$id)),
  h = 0.02,
  id_var = "id",
  day_var = "day"
)

result_bin
```

---

## Interpretation of Output

The function returns:

* Estimated variance components (continuous outcome)
* Estimated marginal prevalence and ICC (binary outcome)
* Optional lag-1 residual correlation (if requested)
* Required number of monitoring days (`k*`)

If no value less than or equal to `k_max` satisfies the precision constraint, `k*` is returned as `NA`.


---

## Citation

If you use this code, please cite:

Hong, H.G. & Matthews, C. (2026).
*A Statistical Framework for Determining Monitoring Duration in Population Surveillance*

