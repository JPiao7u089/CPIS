naive_method_rf = function(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=FALSE,
                           alpha,side,subgroup){
  if (max.t.adjust==TRUE){
    eta = max(dat_train_I1$y[dat_train_I1$delta==1])
    dat_test =dat_test[dat_test$t<=eta,]    
  }
  
  ## working model: random forest
  xnames <- paste0("x",1:4)
  fmla2 = as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse= "+")))
  rf <- rfsrc(fmla2, data = dat_train_I1,nsplit = 10,ntime=0, na.action = "na.impute") # ntime=0 uses all observed event times.
  
  pred = predict(rf, dat_train_I2)
  
  # get cdf of a given time 
  ind1 = as.numeric(sapply(dat_train_I2$y,function(time) sum(pred$time.interest<=time)))
  Sy <- as.numeric(sapply(1:nrow(pred$survival), function(i) pred$survival[i, ind1[i]]))
  Sy[ind1==0] <- max(Sy, na.rm=TRUE)
  
  Fy = 1-Sy
  Fy[Fy==1]=0.9999
  Fy[Fy==0]=0.0001
  
  Fyy = qnorm(Fy); 
  wt = weight.fun(dat_train_I2)
  dat_train_I2 = cbind(dat_train_I2, Fyy, wt)
  
  
  
  if (side=="one"){
    
    fit_quantiles_marginal_naive <- rq(Fyy ~ 1, tau = alpha, data = dat_train_I2,weights=wt)
    bhat_marginal_naive = summary(fit_quantiles_marginal_naive, se = "boot")$coef[,1]
    
    cutoff_naive =   bhat_marginal_naive
    
    pred2 = predict(rf, dat_test[,xnames])
    ind2 = as.numeric(sapply(dat_test$t,function(time) sum(pred2$time.interest<=time)))
    temp <- 1-as.numeric(sapply(1:nrow(pred2$survival), function(i) pred2$survival[i, ind2[i]]))
    
    res_naive = mean(temp>=pnorm(cutoff_naive),na.rm=TRUE)
    for (j in 1:subgroup){
      subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
      res_naive = c(res_naive,mean(temp[subgoup.ind]>=pnorm(cutoff_naive),na.rm=TRUE))
    }
  }
  
  if (side=="two"){
    
    fit_quantiles_marginal_1_naive <- rq(Fyy ~ 1, tau = 1-alpha/2, data = dat_train_I2,weights=wt)
    bhat_marginal_1_naive = summary(fit_quantiles_marginal_1_naive, se = "boot")$coef[,1]
    
    fit_quantiles_marginal_2_naive <- rq(Fyy ~ 1, tau = alpha/2, data = dat_train_I2,weights=wt)
    bhat_marginal_2_naive = summary(fit_quantiles_marginal_2_naive, se = "boot")$coef[,1]
    
    pred2 = predict(rf, dat_test[,xnames])
    ind2 = as.numeric(sapply(dat_test$t,function(time) sum(pred2$time.interest<=time)))
    temp <- 1-as.numeric(sapply(1:nrow(pred2$survival), function(i) pred2$survival[i, ind2[i]]))
    
    cutoff_1 = bhat_marginal_1_naive
    cutoff_2 = bhat_marginal_2_naive
    
    u=pnorm(cutoff_1)
    l=pnorm(cutoff_2)
    
    # bound of failure times
    uuu = as.numeric(sapply(1:nrow(pred2$survival), function(i) {
      sfun1 = stepfun(1-pred2$survival[i,], c(pred2$time.interest,max(pred2$time.interest)))
      sfun1(u)   } ))
    lll = as.numeric(sapply(1:nrow(pred2$survival), function(i) {
      sfun1 = stepfun(1-pred2$survival[i,], c(pred2$time.interest,max(pred2$time.interest)))
      sfun1(l) } ))  # note: removed i from the function 
    
    
    res_naive = c(mean(temp<=pnorm(cutoff_1) & temp>= pnorm(cutoff_2),na.rm=TRUE),
                  mean(lll,na.rm=TRUE),mean(uuu,na.rm=TRUE),mean(uuu-lll,na.rm=TRUE))
    for (j in 1:subgroup){
      subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
      res_naive = c(res_naive,mean(temp[subgoup.ind]>=pnorm(cutoff_2) &
                                     temp[subgoup.ind]<=pnorm(cutoff_1),na.rm=TRUE),
                    mean(lll[subgoup.ind],na.rm=TRUE),
                    mean(uuu[subgoup.ind],na.rm=TRUE),
                    mean(uuu[subgoup.ind]-lll[subgoup.ind],na.rm=TRUE))
    }
    
    cutoff_naive=c(cutoff_1,cutoff_2)
    
  }
  
  
  return(c(res_naive,cutoff_naive))
}