#' Main Function for Conformal Inference
#'
#' Performs conformal inference on training and test data using various user-specified methods.
#' Supports both marginal and subgroup coverage guarantees, with options for one-sided or two-sided prediction intervals.
#'
#' @param dat_train A data frame or matrix containing the training data. Should include predictors and outcomes used to build the model.
#' @param dat_test A data frame or matrix containing the test data to which the conformal method will be applied.
#' @param max.t.adjust Logical. If TRUE, coverage is calculated using only observations from the test data with failure times less than or equal to the maximum observed failure time in the training data.
#' @param alpha The miscoverage level (e.g., 0.1 corresponds to 90% coverage).
#' @param side Character string indicating the type of interval. Options are "two" for two-sided or "one" for one-sided.
#' @param split Numeric. The proportion of training data used to fit the working model. The remaining (1 - split) portion is reserved for quantile regression.
#' @param subgroup A vector or factor indicating subgroup membership for the test observations.
#' @param methods A character vector specifying the conformal inference methods to apply.
#' Options include proposed CPI methods that ensure both marginal and subgroup coverage ("prop.weibull", "prop.cox", "prop.rf"),
#' and naive methods that ensure only marginal coverage but not subgroup coverage ("marginal.weibull", "marginal.cox", "marginal.rf").
#' The suffix (after the dot) indicates the working model used.
#'
#' @return A list or data frame containing coverage rates for the marginal group and subgroups (for one-sided intervals),
#' or coverage, lower bound, upper bound, and interval length (for two-sided intervals).
#'
#' @export
#' @name main_fun


main_fun = function(dat_train, dat_test, max.t.adjust=TRUE,
                   alpha,side,
                   split,
                   subgroup,
                   methods=methods){


  trainRows<-runif(nrow(dat_train))<split #randomly put aside 1-split of the data
  dat_train_I1 <- dat_train[trainRows,]
  dat_train_I2 <- dat_train[!trainRows,]

  res <- vector(mode='list', length=6)

  ########################################
  #### CPI weibull
  ########################################
  if ("prop.weibull" %in% methods | "all" %in% methods){
    res.prop.weibull = prop_method_weibull(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=max.t.adjust,
                                           alpha=alpha,side=side,subgroup=subgroup)[1:(ifelse(side=="two",4,1)*(subgroup+1))]

    res[[1]]=res.prop.weibull
  }


  ########################################
  #### CPI cox
  ########################################
  if ("prop.cox" %in% methods | "all" %in% methods){
    res.prop.cox = prop_method_cox(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=max.t.adjust,
                                   alpha=alpha,side=side,subgroup=subgroup)[1:(ifelse(side=="two",4,1)*(subgroup+1))]
    res[[2]]=res.prop.cox
  }

  ########################################
  #### CPI random forest
  ########################################
  if ("prop.rf" %in% methods | "all" %in% methods){
    res.prop.rf = prop_method_rf(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=max.t.adjust,
                                 alpha=alpha,side=side,subgroup=subgroup)[1:(ifelse(side=="two",4,1)*(subgroup+1))]
    res[[3]]=res.prop.rf
  }


  ########################################
  #### marginal weibull
  ########################################
  if ("marginal.weibull" %in% methods | "all" %in% methods){
    res.marginal.weibull = naive_method_weibull(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=max.t.adjust,
                                                alpha=alpha,side=side,subgroup=subgroup)[1:(ifelse(side=="two",4,1)*(subgroup+1))]
    res[[4]]=res.marginal.weibull
  }

  ########################################
  #### marginal cox
  ########################################
  if ("marginal.cox" %in% methods | "all" %in% methods){
    res.marginal.cox = naive_method_cox(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=max.t.adjust,
                                        alpha=alpha,side=side,subgroup=subgroup)[1:(ifelse(side=="two",4,1)*(subgroup+1))]
    res[[5]]=res.marginal.cox
  }

  ########################################
  ####  marginal random forest
  ########################################
  if ("marginal.rf" %in% methods | "all" %in% methods){
    res.marginal.rf = naive_method_rf(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=max.t.adjust,
                                      alpha=alpha,side=side,subgroup=subgroup)[1:(ifelse(side=="two",4,1)*(subgroup+1))]
    res[[6]]=res.marginal.rf
  }



  res = res[lengths(res) != 0]
  names(res) = methods

  if("all" %in% methods) {names(res) = c("prop.weibull", "prop.cox","prop.rf",
                                         "marginal.weibull","marginal.cox","marginal.rf")}else{
                                           names(res) = methods
                                         }

  generate_group_names <- function(subgroup) {
    c("Marginal", paste0("Subgroup ", seq_len(subgroup)))
  }
  group_names <- generate_group_names(subgroup)

  if(side=="two") {
    metrics<- c("Coverage", "Lower", "Upper", "Length")
    # Convert
    reshape_model <- function(model_name, values) {
      matrix_data <- matrix(values, ncol = 4, byrow = TRUE)
      df <- as.data.frame(matrix_data)
      colnames(df) <- metrics
      df <- cbind(Method = model_name, Subgroup = group_names, df)
      return(df)
    }
    output <- mapply(reshape_model, names(res), res, SIMPLIFY = FALSE)
    }
  else{
    # Convert to data frame
    output <- as.data.frame(do.call(rbind, res))
    colnames(output) <- group_names
  }

  return(output)
}

