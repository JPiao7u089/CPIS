# Main fit, predict, and evaluation functions for CPIS.

cpis_validate_options <- function(alpha, side, split) {
  side <- cpis_validate_alpha_side(alpha, side)
  if (length(split) != 1L || !is.finite(split) || split <= 0 || split >= 1) {
    stop("split must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  side
}

cpis_normalize_methods <- function(methods) {
  proposed_methods <- c("prop.weibull", "prop.cox", "prop.rf")
  if (length(methods) == 1L && identical(methods, "all")) {
    return(proposed_methods)
  }
  methods <- unique(methods)
  invalid <- setdiff(methods, proposed_methods)
  if (length(methods) == 0L || length(invalid) > 0L) {
    stop(
      sprintf(
        "methods must contain only: %s",
        paste(proposed_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  methods
}

# Fit one or more proposed CPIS methods using a shared random split.
cpis_fit <- function(dat_train,
                     alpha = 0.10,
                     side = c("two", "one"),
                     split = 0.50,
                     subgroup = 3,
                     methods = c("prop.weibull", "prop.cox", "prop.rf"),
                     predictor_cols = paste0("x", 1:4),
                     time_col = "y",
                     status_col = "delta",
                     seed = NULL) {
  side <- cpis_validate_options(alpha, side, split)
  methods <- cpis_normalize_methods(methods)
  subgroup_cols <- cpis_subgroup_columns(subgroup)
  cpis_validate_survival_data(dat_train, time_col, status_col, "dat_train")
  cpis_check_columns(dat_train, c(predictor_cols, subgroup_cols), "dat_train")
  cpis_validate_subgroups(dat_train, subgroup_cols, "dat_train")

  if (nrow(dat_train) < 2L * (length(subgroup_cols) + 2L)) {
    stop("dat_train is too small for the requested split and subgroup model.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or one finite number.", call. = FALSE)
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  proper_size <- floor(split * nrow(dat_train))
  proper_index <- sample.int(nrow(dat_train), size = proper_size, replace = FALSE)
  calibration_index <- setdiff(seq_len(nrow(dat_train)), proper_index)
  proper_train <- dat_train[proper_index, , drop = FALSE]
  calibration <- dat_train[calibration_index, , drop = FALSE]

  fits <- list()
  for (method in methods) {
    fits[[method]] <- switch(
      method,
      "prop.weibull" = fit_cpis_weibull(
        proper_train, calibration, alpha, side, predictor_cols,
        subgroup_cols, time_col, status_col
      ),
      "prop.cox" = fit_cpis_cox(
        proper_train, calibration, alpha, side, predictor_cols,
        subgroup_cols, time_col, status_col
      ),
      "prop.rf" = fit_cpis_rsf(
        proper_train, calibration, alpha, side, predictor_cols,
        subgroup_cols, time_col, status_col
      )
    )
  }

  result <- list(
    fits = fits,
    methods = methods,
    alpha = alpha,
    side = side,
    split = split,
    subgroup_cols = subgroup_cols,
    predictor_cols = predictor_cols,
    time_col = time_col,
    status_col = status_col,
    proper_index = proper_index,
    calibration_index = calibration_index
  )
  class(result) <- "cpis_collection"
  result
}

predict.cpis_fit <- function(object, newdata, ...) {
  switch(
    object$method,
    "prop.weibull" = predict_cpis_weibull(object, newdata),
    "prop.cox" = predict_cpis_cox(object, newdata),
    "prop.rf" = predict_cpis_rsf(object, newdata),
    stop("Unknown fitted CPIS method.", call. = FALSE)
  )
}

predict.cpis_collection <- function(object, newdata, ...) {
  predictions <- lapply(object$fits, stats::predict, newdata = newdata, ...)
  names(predictions) <- object$methods
  predictions
}

print.cpis_collection <- function(x, ...) {
  cat("CPIS fit\n")
  cat("  Methods:", paste(x$methods, collapse = ", "), "\n")
  cat("  Coverage target:", 1 - x$alpha, "\n")
  cat("  Interval:", if (x$side == "two") "two-sided" else "one-sided lower", "\n")
  cat("  Proper-training/calibration sizes:",
      length(x$proper_index), "/", length(x$calibration_index), "\n")
  invisible(x)
}

# Evaluate already-constructed intervals. The outcome is used only here, never
# in fitting or interval construction.
evaluate_cpis <- function(intervals,
                          test_data,
                          subgroup_cols,
                          truth_col = "t",
                          side = c("two", "one")) {
  side <- match.arg(side)
  cpis_check_columns(test_data, c(truth_col, subgroup_cols), "test_data")
  cpis_validate_subgroups(test_data, subgroup_cols, "test_data")
  if (nrow(intervals) != nrow(test_data)) {
    stop("intervals and test_data must have the same number of rows.",
         call. = FALSE)
  }
  truth <- test_data[[truth_col]]
  if (!is.numeric(truth) || any(!is.finite(truth)) || any(truth < 0)) {
    stop("The test failure times must be finite and nonnegative.", call. = FALSE)
  }

  summarize_group <- function(index, label) {
    if (!any(index)) {
      warning(sprintf("%s contains no test observations.", label), call. = FALSE)
      return(data.frame(
        Group = label,
        N = 0L,
        Coverage = NA_real_,
        Lower = NA_real_,
        Upper = NA_real_,
        Length = NA_real_
      ))
    }
    covered <- truth[index] >= intervals$lower[index]
    if (side == "two") {
      covered <- covered & truth[index] <= intervals$upper[index]
    }
    data.frame(
      Group = label,
      N = sum(index),
      Coverage = mean(covered),
      Lower = mean(intervals$lower[index]),
      Upper = mean(intervals$upper[index]),
      Length = mean(intervals$upper[index] - intervals$lower[index])
    )
  }

  output <- list(summarize_group(rep(TRUE, nrow(test_data)), "Marginal"))
  for (j in seq_along(subgroup_cols)) {
    output[[j + 1L]] <- summarize_group(
      test_data[[subgroup_cols[j]]] == 1,
      paste("Subgroup", j)
    )
  }
  do.call(rbind, output)
}

# Backward-compatible simulation evaluator. Only the proposed CPIS methods are
# available. For prediction on new patients, use cpis_fit() and predict().
main_fun <- function(dat_train,
                     dat_test,
                     alpha = 0.10,
                     side = c("two", "one"),
                     split = 0.50,
                     subgroup = 3,
                     methods = c("prop.weibull", "prop.cox", "prop.rf"),
                     predictor_cols = paste0("x", 1:4),
                     time_col = "y",
                     status_col = "delta",
                     truth_col = "t",
                     seed = NULL) {
  side <- match.arg(side)
  fitted <- cpis_fit(
    dat_train = dat_train,
    alpha = alpha,
    side = side,
    split = split,
    subgroup = subgroup,
    methods = methods,
    predictor_cols = predictor_cols,
    time_col = time_col,
    status_col = status_col,
    seed = seed
  )
  intervals <- stats::predict(fitted, newdata = dat_test)
  evaluations <- lapply(intervals, function(x) {
    evaluate_cpis(x, dat_test, fitted$subgroup_cols, truth_col, side)
  })

  if (side == "one") {
    output <- do.call(rbind, lapply(names(evaluations), function(method) {
      data.frame(
        Method = method,
        Group = evaluations[[method]]$Group,
        Coverage = evaluations[[method]]$Coverage,
        row.names = NULL
      )
    }))
    return(output)
  }

  output <- lapply(names(evaluations), function(method) {
    data.frame(
      Method = method,
      Subgroup = evaluations[[method]]$Group,
      evaluations[[method]][c("Coverage", "Lower", "Upper", "Length")],
      row.names = NULL
    )
  })
  names(output) <- names(evaluations)
  output
}
