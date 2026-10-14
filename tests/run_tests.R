source("weight.R")
source("prop_weibull.R")
source("prop_cox.R")
source("prop_rf.R")
source("main.R")

dat_train <- read.csv("dat_train.csv")
dat_test <- read.csv("dat_test.csv")
newdata <- dat_test[, setdiff(names(dat_test), c("t", "c", "y", "delta"))]

# The censoring distribution is estimated from proper training, not calibration.
proper_check <- data.frame(y = c(1, 2, 3), delta = c(1, 1, 1))
calibration_check <- data.frame(y = c(1.5, 2.5), delta = c(1, 0))
stopifnot(identical(
  as.numeric(weight.fun(proper_check, calibration_check)),
  c(1, 0)
))

# Fit/predict does not require the unobserved test outcome.
fit <- cpis_fit(
  dat_train,
  alpha = 0.10,
  side = "two",
  split = 0.50,
  subgroup = 3,
  methods = c("prop.weibull", "prop.cox"),
  seed = 2026
)
intervals <- predict(fit, newdata = newdata)
stopifnot(
  identical(names(intervals), c("prop.weibull", "prop.cox")),
  all(vapply(intervals, nrow, integer(1)) == nrow(newdata)),
  all(vapply(intervals, function(x) all(c("lower", "upper") %in% names(x)), logical(1)))
)

# Method order is preserved and results are labeled correctly.
ordered_results <- main_fun(
  dat_train,
  dat_test,
  alpha = 0.10,
  side = "two",
  split = 0.50,
  subgroup = 3,
  methods = c("prop.cox", "prop.weibull"),
  seed = 2026
)
stopifnot(identical(names(ordered_results), c("prop.cox", "prop.weibull")))
standard_results <- main_fun(
  dat_train,
  dat_test,
  alpha = 0.10,
  side = "two",
  split = 0.50,
  subgroup = 3,
  methods = c("prop.weibull", "prop.cox"),
  seed = 2026
)
stopifnot(
  isTRUE(all.equal(ordered_results$prop.cox, standard_results$prop.cox)),
  isTRUE(all.equal(ordered_results$prop.weibull, standard_results$prop.weibull))
)

# One-sided prediction returns a finite lower bound and infinite upper bound.
one_sided_fit <- cpis_fit(
  dat_train,
  alpha = 0.10,
  side = "one",
  split = 0.50,
  subgroup = 3,
  methods = "prop.weibull",
  seed = 7
)
one_sided_intervals <- predict(one_sided_fit, newdata = newdata)$prop.weibull
stopifnot(
  all(is.finite(one_sided_intervals$lower)),
  all(is.infinite(one_sided_intervals$upper))
)

# Run the RSF path when its optional dependency is installed.
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  rsf_fit <- cpis_fit(
    dat_train,
    alpha = 0.10,
    side = "two",
    split = 0.50,
    subgroup = 3,
    methods = "prop.rf",
    seed = 2026
  )
  rsf_intervals <- predict(rsf_fit, newdata = newdata)$prop.rf
  rsf_single_interval <- predict(rsf_fit, newdata = newdata[1, , drop = FALSE])$prop.rf
  stopifnot(
    nrow(rsf_intervals) == nrow(newdata),
    nrow(rsf_single_interval) == 1L
  )
}

message("All available CPIS tests passed.")
