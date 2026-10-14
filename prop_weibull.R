# CPIS with a Weibull accelerated failure-time working model.

fit_cpis_weibull <- function(proper_train,
                             calibration,
                             alpha = 0.10,
                             side = c("two", "one"),
                             predictor_cols = paste0("x", 1:4),
                             subgroup_cols = paste0("eta", 1:3),
                             time_col = "y",
                             status_col = "delta") {
  side <- cpis_validate_alpha_side(alpha, side)
  cpis_validate_survival_data(proper_train, time_col, status_col, "proper_train")
  cpis_validate_survival_data(calibration, time_col, status_col, "calibration")
  cpis_check_columns(proper_train, predictor_cols, "proper_train")
  cpis_check_columns(calibration, c(predictor_cols, subgroup_cols), "calibration")
  cpis_validate_subgroups(calibration, subgroup_cols, "calibration")

  formula <- stats::as.formula(
    paste0(
      "survival::Surv(", time_col, ", ", status_col, ") ~ ",
      paste(predictor_cols, collapse = " + ")
    )
  )
  working_model <- survival::survreg(
    formula,
    dist = "weibull",
    data = proper_train
  )

  linear_predictor <- stats::predict(
    working_model,
    newdata = calibration,
    type = "lp"
  )
  probability <- stats::pweibull(
    calibration[[time_col]],
    shape = 1 / working_model$scale,
    scale = exp(linear_predictor)
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
    method = "prop.weibull",
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
  class(object) <- c("cpis_weibull", "cpis_fit")
  object
}

predict_cpis_weibull <- function(object, newdata) {
  cpis_check_columns(
    newdata,
    c(object$predictor_cols, object$subgroup_cols),
    "newdata"
  )
  design <- cpis_design_matrix(newdata, object$subgroup_cols)
  linear_predictor <- stats::predict(
    object$working_model,
    newdata = newdata,
    type = "lp"
  )
  lower_probability <- stats::pnorm(
    as.numeric(design %*% object$beta_lower)
  )
  lower <- stats::qweibull(
    lower_probability,
    shape = 1 / object$working_model$scale,
    scale = exp(linear_predictor)
  )

  if (object$side == "two") {
    upper_probability <- stats::pnorm(
      as.numeric(design %*% object$beta_upper)
    )
    upper <- stats::qweibull(
      upper_probability,
      shape = 1 / object$working_model$scale,
      scale = exp(linear_predictor)
    )
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
# fit_cpis_weibull() followed by predict().
prop_method_weibull <- function(dat_train_I1,
                                dat_train_I2,
                                dat_test,
                                alpha = 0.10,
                                side = c("two", "one"),
                                subgroup = 3,
                                predictor_cols = paste0("x", 1:4),
                                time_col = "y",
                                status_col = "delta",
                                truth_col = "t") {
  side <- match.arg(side)
  subgroup_cols <- cpis_subgroup_columns(subgroup)
  fit <- fit_cpis_weibull(
    dat_train_I1,
    dat_train_I2,
    alpha,
    side,
    predictor_cols,
    subgroup_cols,
    time_col,
    status_col
  )
  intervals <- predict_cpis_weibull(fit, dat_test)
  evaluation <- evaluate_cpis(intervals, dat_test, subgroup_cols, truth_col, side)
  metrics <- if (side == "two") c("Coverage", "Lower", "Upper", "Length") else "Coverage"
  values <- as.numeric(t(as.matrix(evaluation[metrics])))
  coefficients <- if (side == "two") c(fit$beta_upper, fit$beta_lower) else fit$beta_lower
  c(values, coefficients)
}
