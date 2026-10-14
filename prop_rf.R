# CPIS with a random survival forest working model.

cpis_rsf_cdf_at_time <- function(prediction, times) {
  survival <- prediction$survival
  if (is.null(dim(survival))) {
    survival <- matrix(survival, nrow = 1L)
  } else {
    survival <- as.matrix(survival)
  }
  time_grid <- prediction$time.interest
  if (nrow(survival) != length(times)) {
    stop("The RSF prediction does not match the requested observations.", call. = FALSE)
  }
  index <- findInterval(times, time_grid)
  cdf <- numeric(length(times))
  after_first_time <- index > 0L
  rows <- which(after_first_time)
  cdf[rows] <- 1 - survival[cbind(rows, index[rows])]
  cdf
}

cpis_rsf_quantile <- function(prediction, probabilities) {
  survival <- prediction$survival
  if (is.null(dim(survival))) {
    survival <- matrix(survival, nrow = 1L)
  } else {
    survival <- as.matrix(survival)
  }
  time_grid <- prediction$time.interest
  if (nrow(survival) != length(probabilities)) {
    stop("The RSF prediction does not match the requested probabilities.",
         call. = FALSE)
  }

  vapply(seq_along(probabilities), function(i) {
    probability <- probabilities[i]
    if (is.na(probability)) {
      return(NA_real_)
    }
    if (probability <= 0) {
      return(0)
    }
    cdf <- 1 - survival[i, ]
    index <- which(cdf >= probability)[1]
    if (is.na(index)) Inf else time_grid[index]
  }, numeric(1))
}

fit_cpis_rsf <- function(proper_train,
                         calibration,
                         alpha = 0.10,
                         side = c("two", "one"),
                         predictor_cols = paste0("x", 1:4),
                         subgroup_cols = paste0("eta", 1:3),
                         time_col = "y",
                         status_col = "delta",
                         nsplit = 10,
                         ...) {
  side <- cpis_validate_alpha_side(alpha, side)
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop("Package 'randomForestSRC' is required for method 'prop.rf'.",
         call. = FALSE)
  }
  cpis_validate_survival_data(proper_train, time_col, status_col, "proper_train")
  cpis_validate_survival_data(calibration, time_col, status_col, "calibration")
  cpis_check_columns(proper_train, predictor_cols, "proper_train")
  cpis_check_columns(calibration, c(predictor_cols, subgroup_cols), "calibration")
  cpis_validate_subgroups(calibration, subgroup_cols, "calibration")

  formula <- stats::as.formula(
    paste0(
      "Surv(", time_col, ", ", status_col, ") ~ ",
      paste(predictor_cols, collapse = " + ")
    )
  )
  working_model <- randomForestSRC::rfsrc(
    formula,
    data = proper_train,
    nsplit = nsplit,
    ntime = 0,
    na.action = "na.impute",
    ...
  )
  calibration_prediction <- stats::predict(
    working_model,
    newdata = calibration[predictor_cols]
  )
  probability <- cpis_rsf_cdf_at_time(
    calibration_prediction,
    calibration[[time_col]]
  )
  score <- stats::qnorm(cpis_clip_probability(probability))
  weights <- weight.fun(proper_train, calibration, time_col, status_col)
  probabilities <- if (side == "two") c(alpha / 2, 1 - alpha / 2) else alpha
  coefficients <- cpis_fit_score_quantiles(
    calibration,
    score,
    weights,
    subgroup_cols,
    probabilities
  )

  object <- list(
    method = "prop.rf",
    working_model = working_model,
    alpha = alpha,
    side = side,
    predictor_cols = predictor_cols,
    subgroup_cols = subgroup_cols,
    time_col = time_col,
    status_col = status_col,
    beta_lower = coefficients[[1]],
    beta_upper = if (side == "two") coefficients[[2]] else NULL,
    censoring_fit = attr(weights, "censoring_fit"),
    calibration_failures = sum(weights > 0)
  )
  class(object) <- c("cpis_rsf", "cpis_fit")
  object
}

predict_cpis_rsf <- function(object, newdata) {
  cpis_check_columns(
    newdata,
    c(object$predictor_cols, object$subgroup_cols),
    "newdata"
  )
  design <- cpis_design_matrix(newdata, object$subgroup_cols)
  prediction <- stats::predict(
    object$working_model,
    newdata = newdata[object$predictor_cols]
  )
  lower_probability <- stats::pnorm(
    as.numeric(design %*% object$beta_lower)
  )
  lower <- cpis_rsf_quantile(prediction, lower_probability)

  if (object$side == "two") {
    upper_probability <- stats::pnorm(
      as.numeric(design %*% object$beta_upper)
    )
    upper <- cpis_rsf_quantile(prediction, upper_probability)
  } else {
    upper_probability <- rep(1, nrow(newdata))
    upper <- rep(Inf, nrow(newdata))
  }

  data.frame(
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    lower_probability = as.numeric(lower_probability),
    upper_probability = as.numeric(upper_probability),
    row.names = row.names(newdata)
  )
}

# Backward-compatible evaluation wrapper. New applications should use
# fit_cpis_rsf() followed by predict().
prop_method_rf <- function(dat_train_I1,
                           dat_train_I2,
                           dat_test,
                           alpha = 0.10,
                           side = c("two", "one"),
                           subgroup = 3,
                           predictor_cols = paste0("x", 1:4),
                           time_col = "y",
                           status_col = "delta",
                           truth_col = "t",
                           ...) {
  side <- match.arg(side)
  subgroup_cols <- cpis_subgroup_columns(subgroup)
  fit <- fit_cpis_rsf(
    dat_train_I1,
    dat_train_I2,
    alpha,
    side,
    predictor_cols,
    subgroup_cols,
    time_col,
    status_col,
    ...
  )
  intervals <- predict_cpis_rsf(fit, dat_test)
  evaluation <- evaluate_cpis(intervals, dat_test, subgroup_cols, truth_col, side)
  metrics <- if (side == "two") c("Coverage", "Lower", "Upper", "Length") else "Coverage"
  values <- as.numeric(t(as.matrix(evaluation[metrics])))
  coefficients <- if (side == "two") c(fit$beta_upper, fit$beta_lower) else fit$beta_lower
  c(values, coefficients)
}
