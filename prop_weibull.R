prop_method_weibull = function(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=FALSE,
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

  
  Fyy = qnorm(Fy); 
  wt = weight.fun(dat_train_I2)
  dat_train_I2 = cbind(dat_train_I2, Fyy, wt)
  
  etanames <- paste0("eta",1:subgroup)
  fmla <- as.formula(paste("Fyy ~ ", paste(etanames, collapse= "+")))
  
  if (side=="one"){
    fit_quantiles_marginal <- rq(fmla, tau = alpha, data = dat_train_I2,weights=wt)
    bhat_marginal = summary(fit_quantiles_marginal, se = "boot")$coef[,1]
    
    test.xb = predict(fit.weibull, newdata=dat_test, type='lp')
    temp=pweibull(dat_test$t,shape=1/fit.weibull$scale,scale=exp(test.xb))
    
    design.eta =  as.matrix(cbind(1,dat_test[,etanames]))
    cutoff = design.eta %*% matrix(bhat_marginal,nc=1)
    
    res = mean(temp>=pnorm(cutoff),na.rm=TRUE)
    for (j in 1:subgroup){
      subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
      res = c(res,mean(temp[subgoup.ind]>=pnorm(cutoff)[subgoup.ind],na.rm=TRUE))
    }
    
    bhat=bhat_marginal
    
  }
  
  if (side=="two"){
    
    fit_quantiles_marginal_1 <- rq(fmla, tau = 1-alpha/2, data = dat_train_I2,weights=wt)
    bhat_marginal_1 = summary(fit_quantiles_marginal_1, se = "boot")$coef[,1]
    
    fit_quantiles_marginal_2 <- rq(fmla, tau = alpha/2, data = dat_train_I2,weights=wt)
    bhat_marginal_2 = summary(fit_quantiles_marginal_2, se = "boot")$coef[,1]
    
    test.xb = predict(fit.weibull, newdata=dat_test, type='lp')
    temp=pweibull(dat_test$t,shape=1/fit.weibull$scale,scale=exp(test.xb))   ### note this is cdf
    
    design.eta =  as.matrix(cbind(1,dat_test[,etanames]))
    
    cutoff_1 =  design.eta %*% matrix(bhat_marginal_1,nc=1)
    cutoff_2 =  design.eta %*% matrix(bhat_marginal_2,nc=1)
    
    ## add lower and upper bound to get interval on 09_29
    uu=pnorm(cutoff_1)
    ll=pnorm(cutoff_2)
    
    uuu = as.numeric((-log(1-uu))^fit.weibull$scale*exp(test.xb))
    lll = as.numeric((-log(1-ll))^fit.weibull$scale*exp(test.xb))
    
    ## test lll and uuu
    # piece1 =  ifelse(temp<=uu & temp>= ll,1,0)
    # piece2 =  ifelse(dat_test$t<=uuu &dat_test$t>=lll,1,0)
    # sum(piece1-piece2)
    
    res = c(mean(temp<=pnorm(cutoff_1) & temp>= pnorm(cutoff_2),na.rm=TRUE),
            mean(lll,na.rm=TRUE),mean(uuu,na.rm=TRUE),mean(uuu-lll,na.rm=TRUE))
    for (j in 1:subgroup){
      subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
      res = c(res,mean(temp[subgoup.ind]>=pnorm(cutoff_2)[subgoup.ind] &
                         temp[subgoup.ind]<=pnorm(cutoff_1)[subgoup.ind],na.rm=TRUE),
              mean(lll[subgoup.ind],na.rm=TRUE),mean(uuu[subgoup.ind],na.rm=TRUE),
              mean(uuu[subgoup.ind]-lll[subgoup.ind],na.rm=TRUE))
    }
    
    # res = mean(dat_test$t>=lll &dat_test$t<=uuu)
    # for (j in 1:subgroup){
    #   subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
    #   res = c(res,mean(dat_test$t[subgoup.ind]>=lll[subgoup.ind] &
    #                      dat_test$t[subgoup.ind]<=uuu[subgoup.ind]))
    # }
    
    bhat=c(bhat_marginal_1,bhat_marginal_2)
    
    
  }
  
  return(c(res,bhat))
}