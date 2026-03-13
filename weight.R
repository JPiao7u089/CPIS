weight.fun <- function(dat_train_I2){
  # estimate weights for the failed subjects in the reserved dataset
  if (sum(dat_train_I2$delta)==nrow(dat_train_I2)) {wt=rep(1, nrow(dat_train_I2))} else{
    cox.censor=survfit(Surv(y,(1-delta))~1, data=dat_train_I2)
    censor_surv = rep(NA,nrow(dat_train_I2))
    for (i in 1:nrow(dat_train_I2)){
      censor_surv[i] = summary(cox.censor, times=dat_train_I2$y[i])$surv
    }
    
    wt = dat_train_I2$delta/censor_surv
    summary(wt)
    wt[wt==Inf]=0; wt[is.na(wt)]=0
  }
  return(wt)
}
