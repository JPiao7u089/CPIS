# CPIS

CPIS implements the proposed conformal prediction procedure for right-censored
survival outcomes with coverage calibration across prespecified, potentially
overlapping subgroups. Three working survival models are available:

- `prop.weibull`: Weibull accelerated failure-time model
- `prop.cox`: Cox proportional-hazards model
- `prop.rf`: random survival forest

This repository contains only the proposed CPIS methods. It does not implement
comparator methods from the manuscript.

## Requirements

R packages `survival` and `quantreg` are required. Package `randomForestSRC` is
required only for the random-survival-forest working model.

```r
install.packages(c("survival", "quantreg"))
install.packages("randomForestSRC") # optional
```

From the repository directory, load all functions with:

```r
source("load_cpis.R")
```

## Data format

By default, the training data must contain:

- `y`: observed time, equal to the minimum of failure and censoring times;
- `delta`: failure indicator (`1` = failure, `0` = right-censored);
- `x1`, ..., `x4`: working-model predictors; and
- `eta1`, ..., `eta3`: prespecified 0/1 subgroup indicators.

The predictor, subgroup, time, and status column names can be changed through
the corresponding function arguments. Subgroups may overlap.

## Fit and predict

The recommended interface separates model fitting from prediction. Prediction
on a new subject requires covariates and subgroup indicators, but not the
subject's failure time.

```r
dat_train <- read.csv("dat_train.csv")
dat_test <- read.csv("dat_test.csv")

fit <- cpis_fit(
  dat_train = dat_train,
  alpha = 0.10,
  side = "two",
  split = 0.50,
  subgroup = c("eta1", "eta2", "eta3"),
  methods = c("prop.weibull", "prop.cox"),
  seed = 2026
)

intervals <- predict(fit, newdata = dat_test)
head(intervals$prop.weibull)
```

Each prediction data frame contains subject-level lower and upper prediction
limits and the calibrated lower- and upper-tail probability levels.

To use one working model directly:

```r
set.seed(2026)
proper_index <- sample.int(nrow(dat_train), floor(0.5 * nrow(dat_train)))

weibull_fit <- fit_cpis_weibull(
  proper_train = dat_train[proper_index, ],
  calibration = dat_train[-proper_index, ],
  alpha = 0.10,
  side = "two",
  subgroup_cols = c("eta1", "eta2", "eta3")
)

weibull_intervals <- predict(weibull_fit, newdata = dat_test)
```

## Simulation evaluation

`main_fun()` is retained for simulation studies in which the true test failure
time is available in column `t`. The test outcome is used only to evaluate the
completed intervals and is never used during model fitting or prediction.

```r
results <- main_fun(
  dat_train = dat_train,
  dat_test = dat_test,
  alpha = 0.10,
  side = "two",
  split = 0.50,
  subgroup = 3,
  methods = c("prop.weibull", "prop.cox"),
  seed = 2026
)

results$prop.weibull
```


