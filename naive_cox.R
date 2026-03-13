naive_method_cox = function(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=FALSE,
                        alpha,side,subgroup){
  if (max.t.adjust==TRUE){
    eta = max(dat_train_I1$y[dat_train_I1$delta==1])
    dat_test =dat_test[dat_test$t<=eta,]    
  }
  
  ## working model: Cox model
  xnames <- paste0("x",1:4)
  cox1 = coxph(as.formula(paste("Surv(y, delta) ~ ", paste(xnames, collapse= "+"))), data=dat_train_I1)
  beta.1=cox1$coefficients
  bh <- basehaz(cox1, centered=FALSE)
  
  sfun0 = stepfun(bh$time, c(0,bh$hazard))
  sfun1_upper = stepfun(bh$hazard,c(bh$time,max(bh$time))) # take the upper bound of a range, max of observed censroed time if infnity
  
  train.xb = as.numeric(t(beta.1) %*% t(dat_train_I2[,xnames]))  
  Fy = 1-exp(-sfun0(dat_train_I2$y)*exp(train.xb)) 
  
  Fy[Fy==1]=0.999
  Fy[Fy==0]=0.001
  
  Fyy = qnorm(Fy); 
  wt = weight.fun(dat_train_I2)
  dat_train_I2 = cbind(dat_train_I2, Fyy, wt)
  
  if (side=="one"){
    
    fit_quantiles_marginal_naive <- rq(Fyy ~ 1, tau = alpha, data = dat_train_I2,weights=wt)
    bhat_marginal_naive = summary(fit_quantiles_marginal_naive, se = "boot")$coef[,1]
    
    cutoff_naive =   bhat_marginal_naive
    
    test.xb = as.numeric(t(beta.1) %*% t(dat_test[,xnames]) )
    temp = 1-exp(-sfun0(dat_test$t)*exp(test.xb)) 
    
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
    
    test.xb = as.numeric(t(beta.1) %*% t(dat_test[,xnames]) )
    temp = 1-exp(-sfun0(dat_test$t)*exp(test.xb)) 
    
    cutoff_1 = bhat_marginal_1_naive
    cutoff_2 = bhat_marginal_2_naive
    
    u=pnorm(cutoff_1)
    l=pnorm(cutoff_2)
    
    uu = -log(1-u)*exp(- test.xb)
    ll = -log(1-l)*exp(- test.xb) 
    
    uuu = sfun1_upper(uu); 
    lll = sfun1_upper(ll)
    
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
