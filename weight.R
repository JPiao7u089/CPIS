# Utilities shared by the CPIS implementations.

cpis_check_columns <- function(data, columns, data_name = "data") {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        data_name,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

cpis_validate_alpha_side <- function(alpha, side) {
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  match.arg(side, c("two", "one"))
}

cpis_validate_survival_data <- function(data, time_col, status_col, data_name) {
  cpis_check_columns(data, c(time_col, status_col), data_name)

  time <- data[[time_col]]
  status <- data[[status_col]]
  if (!is.numeric(time) || any(!is.finite(time)) || any(time < 0)) {
    stop(sprintf("%s$%s must contain finite, nonnegative times.", data_name, time_col),
         call. = FALSE)
  }
  if (any(is.na(status)) || !all(status %in% c(0, 1))) {
    stop(sprintf("%s$%s must contain only 0 (censored) and 1 (failure).",
                 data_name, status_col), call. = FALSE)
  }
  invisible(TRUE)
}

cpis_validate_subgroups <- function(data, subgroup_cols, data_name = "data") {
  if (length(subgroup_cols) == 0L) {
    return(invisible(TRUE))
  }
  cpis_check_columns(data, subgroup_cols, data_name)
  valid <- vapply(data[subgroup_cols], function(x) {
    (is.numeric(x) || is.logical(x)) && !anyNA(x) && all(x %in% c(0, 1))
  }, logical(1))
  if (!all(valid)) {
    stop(
      sprintf(
        "Subgroup columns must be complete 0/1 indicators. Invalid column(s): %s",
        paste(subgroup_cols[!valid], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

cpis_subgroup_columns <- function(subgroup) {
  if (is.character(subgroup)) {
    if (length(subgroup) == 0L || anyNA(subgroup) || any(!nzchar(subgroup))) {
      stop("subgroup must contain one or more valid column names.", call. = FALSE)
    }
    return(subgroup)
  }
  if (length(subgroup) != 1L || is.na(subgroup) || subgroup < 1 ||
      subgroup != as.integer(subgroup)) {
    stop("subgroup must be a positive integer or a character vector of subgroup columns.",
         call. = FALSE)
  }
  paste0("eta", seq_len(as.integer(subgroup)))
}

cpis_design_matrix <- function(data, subgroup_cols) {
  cpis_validate_subgroups(data, subgroup_cols)
  design <- cbind("(Intercept)" = 1, as.matrix(data[subgroup_cols]))
  storage.mode(design) <- "double"
  design
}

cpis_clip_probability <- function(probability) {
  epsilon <- .Machine$double.eps
  pmin(pmax(probability, epsilon), 1 - epsilon)
}

# Fit the censoring survival distribution G on the proper-training sample.
cpis_fit_censoring_km <- function(proper_train,
                                  time_col = "y",
                                  status_col = "delta") {
  cpis_validate_survival_data(proper_train, time_col, status_col, "proper_train")

  if (all(proper_train[[status_col]] == 1)) {
    return(structure(list(all_observed = TRUE), class = "cpis_censoring_km"))
  }

  censoring_data <- data.frame(
    time = proper_train[[time_col]],
    censoring_event = 1 - proper_train[[status_col]]
  )
  fit <- survival::survfit(
    survival::Surv(time, censoring_event) ~ 1,
    data = censoring_data
  )

  structure(
    list(
      all_observed = FALSE,
      time = fit$time,
      survival = fit$surv
    ),
    class = "cpis_censoring_km"
  )
}

# Evaluate the right-continuous Kaplan--Meier estimate G(t).
cpis_censoring_survival <- function(censoring_fit, times) {
  if (isTRUE(censoring_fit$all_observed)) {
    return(rep(1, length(times)))
  }

  index <- findInterval(times, censoring_fit$time)
  result <- rep(1, length(times))
  after_first_time <- index > 0L
  result[after_first_time] <- censoring_fit$survival[index[after_first_time]]
  result
}

# IPCW weights for calibration observations. G is estimated exclusively from
# the proper-training sample; weights are not truncated or bounded.
weight.fun <- function(dat_train_I1,
                       dat_train_I2,
                       time_col = "y",
                       status_col = "delta") {
  cpis_validate_survival_data(dat_train_I1, time_col, status_col, "dat_train_I1")
  cpis_validate_survival_data(dat_train_I2, time_col, status_col, "dat_train_I2")

  censoring_fit <- cpis_fit_censoring_km(dat_train_I1, time_col, status_col)
  censoring_survival <- cpis_censoring_survival(
    censoring_fit,
    dat_train_I2[[time_col]]
  )

  failures <- dat_train_I2[[status_col]] == 1
  invalid <- failures & (!is.finite(censoring_survival) | censoring_survival <= 0)
  if (any(invalid)) {
    stop(
      paste0(
        "The estimated censoring survival probability is zero or undefined for ",
        sum(invalid), " uncensored calibration observation(s). The censoring ",
        "positivity condition is not satisfied; CPIS cannot be fitted without ",
        "changing the target population or imposing an explicit horizon."
      ),
      call. = FALSE
    )
  }

  weights <- numeric(nrow(dat_train_I2))
  weights[failures] <- 1 / censoring_survival[failures]
  attr(weights, "censoring_fit") <- censoring_fit
  weights
}

cpis_fit_score_quantiles <- function(calibration,
                                     score,
                                     weights,
                                     subgroup_cols,
                                     probabilities) {
  if (length(score) != nrow(calibration) || length(weights) != nrow(calibration)) {
    stop("score and weights must have one value per calibration observation.",
         call. = FALSE)
  }
  keep <- is.finite(score) & is.finite(weights) & weights > 0
  if (sum(keep) <= length(subgroup_cols) + 1L) {
    stop("Too few uncensored calibration observations to fit subgroup quantiles.",
         call. = FALSE)
  }

  regression_data <- data.frame(score = score[keep])
  if (length(subgroup_cols) > 0L) {
    regression_data[subgroup_cols] <- calibration[keep, subgroup_cols, drop = FALSE]
  }
  formula <- if (length(subgroup_cols) == 0L) {
    score ~ 1
  } else {
    stats::reformulate(subgroup_cols, response = "score")
  }

  fits <- lapply(probabilities, function(tau) {
    fit <- quantreg::rq(
      formula,
      tau = tau,
      data = regression_data,
      weights = weights[keep]
    )
    coefficients <- stats::coef(fit)
    expected <- c("(Intercept)", subgroup_cols)
    coefficients <- coefficients[expected]
    if (any(!is.finite(coefficients))) {
      stop(
        paste0(
          "The weighted subgroup quantile regression is rank deficient. ",
          "Check for empty, duplicate, or linearly dependent subgroup indicators."
        ),
        call. = FALSE
      )
    }
    coefficients
  })
  names(fits) <- as.character(probabilities)
  fits
}
