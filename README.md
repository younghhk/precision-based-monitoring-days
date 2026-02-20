
[![Back to Hub](https://img.shields.io/badge/⬅️%20Back%20to%20Hub-2962FF?style=for-the-badge)](https://github.com/younghhk/NCI)

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

The primary method assumes a random-intercept mixed-effects model for estimating the population mean:

**y = μ + b + e**

where:
* **y** is the daily observation
* **μ** is the population mean (fixed effect)
* **b ~ N(0, σ²_b)** is the subject-specific random effect (between-person variability)
* **e ~ N(0, σ²_w)** is the day-level random effect (within-person variability)

Key assumptions:
* Both b and e are random effects; b is subject-specific and constant across days within a participant
* Daily observations are conditionally independent given the subject-specific random effect b
* Correlation across days arises naturally through the shared subject-specific random effect
* The model focuses on estimating the population mean without additional covariates

The design rule assumes a planned study with equal monitoring days per participant. Estimated variance components from existing data are used to determine the required monitoring duration.

### Optional Autocorrelation Adjustment

If short-term persistence remains after accounting for random intercepts (e.g., AR(1) structure), an optional sensitivity analysis estimates an effective number of monitoring days:

```
k_eff = k / design_effect
```

This adjustment is most relevant for consecutive-day monitoring protocols.

---

### For Binary Outcomes

Estimate the marginal prevalence and intraclass correlation:

```r
library(lme4)

# Fit logistic random-intercept model
fit_binary <- glmer(binary_outcome ~ (1 | participant_id), 
                    data = pilot_data, 
                    family = binomial)

# Extract variance components
sigma_b2_logit <- as.data.frame(VarCorr(fit_binary))$vcov[1]

# Calculate ICC on probability scale
# For logistic model: ρ ≈ σ²_b / (σ²_b + π²/3)
rho <- sigma_b2_logit / (sigma_b2_logit + pi^2/3)

# Calculate marginal prevalence
mu <- mean(pilot_data$binary_outcome, na.rm = TRUE)

# Display estimates
print(paste("Marginal prevalence (μ):", round(mu, 3)))
print(paste("Intraclass correlation (ρ):", round(rho, 3)))
```


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
  Between-person variance (σ²_b). Extract from `VarCorr()` output after fitting `lmer()` model.
* **`sigma_w2`**  
  Within-person variance (σ²_w). Extract from `VarCorr()` output after fitting `lmer()` model.
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
  Marginal prevalence/mean probability. Calculate as `mean(binary_outcome)` from pilot data.
* **`rho`**  
  Intraclass correlation (ICC). Calculate from random-intercept logistic model (see example above).
* **`phi`** (default = 0)  
  Lag-1 autocorrelation. Set to 0 for independence, or provide estimated value from prior data/literature.
* **`k_max`** (default = 14)  
  Maximum number of monitoring days to evaluate.
* **`z`** (default = 1.96)  
  Z-value for confidence interval (1.96 for 95% CI).

---

## Example Usage

### Example 1: Complete Workflow for Continuous Outcome (No Autocorrelation)

```r
library(lme4)

# Step 1 (if needed): Estimate variance components from pilot data
# Only required when sigma_b2 and sigma_w2 are not available
# from prior studies or existing surveillance data.

# library(lme4)
# fit <- lmer(outcome ~ 1 + (1 | id), data = pilot_data)
# vc <- VarCorr(fit)
# sigma_b2 <- as.numeric(vc$id)        # Between-person variance
# sigma_w2 <- attr(vc, "sc")^2         # Within-person variance

# Step 2: Determine required monitoring days for new study
result_cont <- compute_k_star_continuous(
  N = 100,
  h = 0.2,
  sigma_b2 = 0.5,
  sigma_w2 = 1,
  phi = 0,
  k_max = 14
)
```

**Example output:**
```
=== Precision-Based Monitoring Days (Continuous Outcome) ===
Sample size (N): 100 
Desired precision (h): 0.2 
Between-person variance (σ²_b): 0.5 
Within-person variance (σ²_w): 1 

Required monitoring days (k*): 2 

Note: margin_of_error is the ±CI half-width achieved at each k.
      k* is the minimum k where margin_of_error ≤ h.

    k margin_of_error meets_tolerance
1   1       0.2400500           FALSE
2   2       0.1960000            TRUE
3   3       0.1789227            TRUE
4   4       0.1697410            TRUE
5   5       0.1639854            TRUE
6   6       0.1600333            TRUE
7   7       0.1571496            TRUE
...
```

---

### Example 2: Continuous Outcome with Autocorrelation

**Note:** `phi` can be estimated from pilot data using the `estimate_lag1_acf()` helper function (requires fitted model and data), obtained from published literature, or assumed based on typical values in EMA/accelerometry studies (0.2-0.4).

```r
# Option A: Estimate autocorrelation from pilot data
phi_est <- estimate_lag1_acf(
  model = fit,
  data = pilot_data,
  id_var = "id",
  day_var = "day"
)

# Option B: Or use literature value
result_cont_adj <- compute_k_star_continuous(
  N = 200,
  h = 0.10,
  sigma_b2 = 0.3,
  sigma_w2 = 0.8,
  phi = 0.3,  # Estimated from pilot data or literature
  k_max = 14
)
```

**Example output:**
```
=== Precision-Based Monitoring Days (Continuous Outcome) ===
Sample size (N): 200 
Desired precision (h): 0.1 
Between-person variance (σ²_b): 0.3 
Within-person variance (σ²_w): 0.8 
Lag-1 autocorrelation (φ): 0.3 

Required monitoring days (k*): 6 

Note: margin_of_error is the ±CI half-width achieved at each k.
      k* is the minimum k where margin_of_error ≤ h.

    k margin_of_error meets_tolerance
1   1      0.14535749           FALSE
2   2      0.12550124           FALSE
3   3      0.11506830           FALSE
4   4      0.10830663           FALSE
5   5      0.10353313           FALSE
6   6      0.09998190            TRUE
...
```

---

### Example 3: Complete Workflow for Binary Outcome

```r
library(lme4)

# Step 1: Estimate parameters from pilot data
fit_binary <- glmer(binary_outcome ~ (1 | id), 
                    data = pilot_data, 
                    family = binomial)

sigma_b2_logit <- as.data.frame(VarCorr(fit_binary))$vcov[1]
rho <- sigma_b2_logit / (sigma_b2_logit + pi^2/3)
mu <- mean(pilot_data$binary_outcome, na.rm = TRUE)

# Step 2: Determine required monitoring days
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

Required monitoring days (k*): 13 

Note: margin_of_error is the ±CI half-width achieved at each k.
      k* is the minimum k where margin_of_error ≤ h.

    k margin_of_error meets_tolerance
1   1      0.08981848           FALSE
2   2      0.07100775           FALSE
3   3      0.06351126           FALSE
...
13 13      0.04982233            TRUE
14 14      0.04948757            TRUE
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


## Installation

```r
# Source the functions directly
source("compute_k_star_functions.R")
```

---

## Citation

If you use this code, please cite:

Hong, H.G. & Matthews, C. (in preparation).  
*A Statistical Framework for Determining Monitoring Duration in Population Surveillance*

Or cite this repository directly:
```
Hong, H.G. & Matthews, C. (2026). Precision-Based Monitoring Days [R code]. 
[https://github.com/younghhk/precision-based-monitoring-days](https://github.com/younghhk/precision-based-monitoring-days)
```

---

## License

[To be added]

---

## Contact

grace.hong@nih.gov
