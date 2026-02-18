# Precision-Based Monitoring Days

This repository provides an R implementation of a precision-based rule for determining the number of monitoring days required per participant to achieve a desired level of accuracy when estimating a population mean or prevalence.

Using variance components estimated from existing data (e.g., pilot studies or large observational cohorts), the method addresses a practical design question:

> Given the observed between-person and day-to-day variability, how many monitoring days are needed per participant so that the 95% confidence interval for the population estimate is sufficiently narrow?

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

## Primary Functions

### `compute_k_star_continuous()`
Computes the required number of monitoring days (k*) for **continuous outcomes**.

### `compute_k_star_binary()`
Computes the required number of monitoring days (k*) for **binary outcomes**.

These simplified functions allow direct computation of k* when variance components are known, without requiring a fitted model or original dataset.

---

## Arguments

### For Continuous Outcomes
* **`N`**  
  Number of participants.
* **`h`**  
  Precision tolerance. Maximum acceptable half-width (margin of error) of a 95% confidence interval.
* **`sigma_b2`**  
  Between-person variance (σ²_b).
* **`sigma_w2`**  
  Within-person variance (σ²_w).
* **`phi`** (default = 0)  
  Lag-1 autocorrelation. Set to 0 for independence, or provide estimated value from prior data/literature (typically 0.2-0.4 for EMA studies).
* **`k_max`** (default = 14)  
  Maximum number of monitoring days to evaluate.
* **`z`** (default = 1.96)  
  Z-value for confidence interval (1.96 for 95% CI).

### For Binary Outcomes
* **`N`**  
  Number of participants.
* **`h`**  
  Precision tolerance. Maximum acceptable half-width (margin of error) of a 95% confidence interval.
* **`mu`**  
  Marginal prevalence/mean probability.
* **`rho`**  
  Intraclass correlation (ICC).
* **`phi`** (default = 0)  
  Lag-1 autocorrelation. Set to 0 for independence, or provide estimated value from prior data/literature.
* **`k_max`** (default = 14)  
  Maximum number of monitoring days to evaluate.
* **`z`** (default = 1.96)  
  Z-value for confidence interval (1.96 for 95% CI).

---

## Example Usage

### Example 1: Continuous Outcome (No Autocorrelation)

```r
result_cont <- compute_k_star_continuous(
  N = 100,
  h = 0.10,
  sigma_b2 = 0.5,
  sigma_w2 = 1.0,
  phi = 0,
  k_max = 14
)
```

**Example output:**
```
=== Precision-Based Monitoring Days (Continuous Outcome) ===
Sample size (N): 100
Desired precision (h): 0.1
Between-person variance (σ²_b): 0.5
Within-person variance (σ²_w): 1

Required monitoring days (k*): 7

Note: margin_of_error is the ±CI half-width achieved at each k.
      k* is the minimum k where margin_of_error ≤ h.

    k margin_of_error meets_tolerance
1   1       0.2400500          FALSE
2   2       0.1697056          FALSE
3   3       0.1385641          FALSE
4   4       0.1200250          FALSE
5   5       0.1073313          FALSE
6   6       0.0980297          TRUE
7   7       0.0907376          TRUE
...
```

---

### Example 2: Continuous Outcome with Autocorrelation

**Note:** `phi` can be estimated from pilot data using the `estimate_lag1_acf()` helper function (requires fitted model and data), obtained from published literature, or assumed based on typical values in EMA/accelerometry studies (0.2-0.4).

```r
result_cont_adj <- compute_k_star_continuous(
  N = 100,
  h = 0.10,
  sigma_b2 = 0.5,
  sigma_w2 = 1.0,
  phi = 0.3,  # Estimated from pilot data or literature
  k_max = 14
)
```

**Example output:**
```
=== Precision-Based Monitoring Days (Continuous Outcome) ===
Sample size (N): 100
Desired precision (h): 0.1
Between-person variance (σ²_b): 0.5
Within-person variance (σ²_w): 1
Lag-1 autocorrelation (φ): 0.3

Required monitoring days (k*): 9

Note: margin_of_error is the ±CI half-width achieved at each k.
      k* is the minimum k where margin_of_error ≤ h.

    k margin_of_error meets_tolerance
1   1       0.2400500          FALSE
2   2       0.1886953          FALSE
3   3       0.1685423          FALSE
...
```

---

### Example 3: Binary Outcome

```r
result_bin <- compute_k_star_binary(
  N = 100,
  h = 0.05,
  mu = 0.3,
  rho = 0.25,
  phi = 0,
  k_max = 14
)
```

**Example output:**
```
=== Precision-Based Monitoring Days (Binary Outcome) ===
Sample size (N): 100
Desired precision (h): 0.05
Marginal prevalence (μ): 0.3
Intraclass correlation (ρ): 0.25

Required monitoring days (k*): 5

Note: margin_of_error is the ±CI half-width achieved at each k.
      k* is the minimum k where margin_of_error ≤ h.

    k margin_of_error meets_tolerance
1   1      0.08977724          FALSE
2   2      0.06347467          FALSE
3   3      0.05182954          FALSE
4   4      0.04488862          TRUE
5   5      0.04016623          TRUE
...
```

---

## Interpretation of Output

Each function returns a list containing:
* **`k_star`** - Required number of monitoring days (or `NA` if not attainable within `k_max`)
* Input parameters (N, h, variance components, etc.)
* **`results`** - Data frame showing margin of error and whether tolerance is met for each k value

The function automatically prints:
* Summary of input parameters
* Required monitoring days (k*)
* Table showing margin of error at each k value
* Indicator of which k values meet the precision tolerance

---

## Notes on Autocorrelation (phi)

### What is phi?
* Lag-1 autocorrelation between consecutive day residuals
* Represents short-term persistence not captured by random intercepts
* Typical range in EMA/accelerometry studies: 0.2 to 0.4
* phi = 0 assumes conditional independence given random effects

### How to specify phi

**Option 1: Assume independence**  
Set `phi = 0` (simplest approach, but may underestimate required days if autocorrelation exists)

**Option 2: Use literature values**  
Set `phi = 0.3` (reasonable middle-ground for most EMA/accelerometry data)

**Option 3: Estimate from pilot data**  
If you have existing data with a fitted model, use the `estimate_lag1_acf()` helper function:
```r
# Fit model to pilot data
fit <- lmer(outcome ~ (1 | id), data = pilot_data)

# Estimate phi
phi_est <- estimate_lag1_acf(
  model = fit,
  data = pilot_data,
  id_var = "id",
  day_var = "day"
)

# Use in k* calculation
result <- compute_k_star_continuous(
  N = 100,
  h = 0.10,
  sigma_b2 = 0.5,
  sigma_w2 = 1.0,
  phi = phi_est
)
```

**Option 4: Use domain knowledge**  
Consult published studies in your research area for typical autocorrelation values

### When autocorrelation matters
* Higher phi values require more monitoring days to achieve the same precision
* The effect is accounted for through AR(1) design effect adjustment
* Most relevant for consecutive-day monitoring protocols
* Less relevant for sampling designs with gaps between measurements



## Installation

```r
# Source the functions directly
source("compute_k_star_functions.R")
```

---

## Citation

If you use this code, please cite:

Hong, H.G. & Matthews, C. (2026).  
*A Statistical Framework for Determining Monitoring Duration in Population Surveillance*

---

## License





## Contact

grace.hong@nih.gov
