naive_method_weibull = function(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=FALSE,
                                alpha,side,subgroup){
  if (max.t.adjust==TRUE){
    eta = max(dat_train_I1$y[dat_train_I1$delta==1])
    dat_test =dat_test[dat_test$t<=eta,]    
  }
  
  ## working model: Weibull
  fit.weibull = survreg(Surv(y, delta) ~ x1+x2+x3+x4, dist="weibull", data=dat_train_I1) 
  train.xb = predict(fit.weibull, newdata=dat_train_I2, type='lp')
  Fy = pweibull(dat_train_I2$y,shape=1/fit.weibull$scale,scale=exp(train.xb))
  Fy[Fy==1]=0.999
  Fy[Fy==0]=0.001

  Fyy = qnorm(Fy)
  
  wt = weight.fun(dat_train_I2)
  dat_train_I2 = cbind(dat_train_I2, Fyy, wt)
  
  if (side=="one"){
    
    fit_quantiles_marginal_naive <- rq(Fyy ~ 1, tau = alpha, data = dat_train_I2,weights=wt)
    bhat_marginal_naive = summary(fit_quantiles_marginal_naive, se = "boot")$coef[,1]
    
    cutoff_naive =   bhat_marginal_naive
    
    test.xb = predict(fit.weibull, newdata=dat_test, type='lp')
    temp=pweibull(dat_test$t,shape=1/fit.weibull$scale,scale=exp(test.xb))
    
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
    
    test.xb = predict(fit.weibull, newdata=dat_test, type='lp')
    temp=pweibull(dat_test$t,shape=1/fit.weibull$scale,scale=exp(test.xb))
    
    
    ## add lower and upper bound to get interval on 09_29
    cutoff_1 = bhat_marginal_1_naive
    cutoff_2 = bhat_marginal_2_naive
    
    uu=pnorm(cutoff_1)
    ll=pnorm(cutoff_2)
    
    uuu = as.numeric((-log(1-uu))^fit.weibull$scale*exp(test.xb))
    lll = as.numeric((-log(1-ll))^fit.weibull$scale*exp(test.xb))
    
    res_naive = c(mean(temp<=pnorm(cutoff_1) & temp>= pnorm(cutoff_2),na.rm=TRUE),
                  mean(lll,na.rm=TRUE),mean(uuu,na.rm=TRUE),mean(uuu-lll,na.rm=TRUE))
    for (j in 1:subgroup){
      subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
      res_naive = c(res_naive,mean(temp[subgoup.ind]>=pnorm(cutoff_2) &
                                     temp[subgoup.ind]<=pnorm(cutoff_1),na.rm=TRUE),
                    mean(lll[subgoup.ind],na.rm=TRUE),mean(uuu[subgoup.ind],na.rm=TRUE),
                    mean(uuu[subgoup.ind]-lll[subgoup.ind],na.rm=TRUE))
    }
    
    cutoff_naive=c(cutoff_1,cutoff_2)
    
  }
  
  
  return(c(res_naive,cutoff_naive))
}

