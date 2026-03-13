prop_method_rf = function(dat_train_I1,dat_train_I2,dat_test,max.t.adjust=FALSE,
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
  
  etanames <- paste0("eta",1:subgroup)
  fmla <- as.formula(paste("Fyy ~ ", paste(etanames, collapse= "+")))
  
  if (side=="one"){
    fit_quantiles_marginal <- rq(fmla, tau = alpha, data = dat_train_I2,weights=wt)
    bhat_marginal = summary(fit_quantiles_marginal, se = "boot")$coef[,1]
    
    pred2 = predict(rf, dat_test[,xnames])
    ind2 = as.numeric(sapply(dat_test$t,function(time) sum(pred2$time.interest<=time)))
    temp <- 1-as.numeric(sapply(1:nrow(pred2$survival), function(i) pred2$survival[i, ind2[i]]))
    
    
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
    
    pred2 = predict(rf, dat_test[,xnames])
    ind2 = as.numeric(sapply(dat_test$t,function(time) sum(pred2$time.interest<=time)))
    temp <- 1-as.numeric(sapply(1:nrow(pred2$survival), function(i) pred2$survival[i, ind2[i]]))
    
    design.eta =  as.matrix(cbind(1,dat_test[,etanames]))
    cutoff_1 =  design.eta %*% matrix(bhat_marginal_1,nc=1)
    cutoff_2 =  design.eta %*% matrix(bhat_marginal_2,nc=1)
    
    u=pnorm(cutoff_1)
    l=pnorm(cutoff_2)
    
    # coverage 
    # res = mean(temp<=pnorm(cutoff_1) & temp>= pnorm(cutoff_2))
    # for (j in 1:subgroup){
    #   subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
    #   res = c(res,mean(temp[subgoup.ind]>=pnorm(cutoff_2)[subgoup.ind] &
    #                      temp[subgoup.ind]<=pnorm(cutoff_1)[subgoup.ind]))
    # }
    # res
    
    # bound of failure times
    uuu = as.numeric(sapply(1:nrow(pred2$survival), function(i) {
      sfun1 = stepfun(1-pred2$survival[i,], c(pred2$time.interest,max(pred2$time.interest)))
      sfun1(u[i])   } ))
    lll = as.numeric(sapply(1:nrow(pred2$survival), function(i) {
      sfun1 = stepfun(1-pred2$survival[i,], c(pred2$time.interest,max(pred2$time.interest)))
      sfun1(l[i]) } ))
    
    # verify it 
    # res = mean(dat_test$t>=lll &dat_test$t<=uuu)
    # for (j in 1:subgroup){
    #   subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
    #   res = c(res,mean(dat_test$t[subgoup.ind]>=lll[subgoup.ind] &
    #                      dat_test$t[subgoup.ind]<=uuu[subgoup.ind]))
    # }
    # res
    # 
    
    
    res = c(mean(temp<=pnorm(cutoff_1) & temp>= pnorm(cutoff_2),na.rm=TRUE),
            mean(lll,na.rm=TRUE),mean(uuu,na.rm=TRUE),mean(uuu-lll,na.rm=TRUE))
    for (j in 1:subgroup){
      subgoup.ind = eval(parse(text = paste('dat_test$eta', j, sep='')))==1
      res = c(res,mean(temp[subgoup.ind]>=pnorm(cutoff_2)[subgoup.ind] &
                         temp[subgoup.ind]<=pnorm(cutoff_1)[subgoup.ind],na.rm=TRUE),
              mean(lll[subgoup.ind],na.rm=TRUE),
              mean(uuu[subgoup.ind],na.rm=TRUE),
              mean(uuu[subgoup.ind]-lll[subgoup.ind],na.rm=TRUE))
    }
    bhat=c(bhat_marginal_1,bhat_marginal_2)
    
    
  }
  
  return(c(res,bhat))
}