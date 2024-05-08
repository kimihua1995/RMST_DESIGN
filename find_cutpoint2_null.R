rmst.reg <- function(dat, tau, trt){
  dat$y <- pmin(dat$time,tau)
  dat$d <- dat$status
  dat$d[dat$y == tau] <- 1
  id=order(dat$y)
  dat=dat[id,]
  fitc=my.func_surv(dat$y, 1-dat$d)
  dat$weight=dat$d/rep(pmax(fitc$surv,0.001), table(dat$y))
  if (trt == 1){
    fit <- lm(y ~ 1 + X, weights=weight, data = dat)
  }else if (trt == 0){
    fit <- lm(y ~ 1, weights=weight, data = dat)
  }
  
  return(fit)
}



fit.rmst.reg <- function(dat0,dat1,tau){
  tau.temp <- min(max(dat0$time), max(dat1$time))
  tau.temp <- min(tau.temp, tau)
  
  dat0$A=0
  dat1$A=1
  dat <- rbind(dat0,dat1)
  cov <- data.frame(A=dat$A, AX=dat$A*dat$X)
  rmst_fit <- my.rmst2reg(y = dat$time,
                          delta = dat$status,
                          x = cov,
                          arm = dat$A,
                          tau = tau.temp)
  return(rmst_fit)
}




find.cutpoint <- function(rmst_fit){
  
  gamma <- coef(rmst_fit)
  cutpoint <- -gamma[2]/gamma[3]
  cutpoint <- min(max(cutpoint,0),1)
  sign <- gamma[3]
  
  
  if ((cutpoint == 0 & sign < 0) | (cutpoint == 1 & sign > 0)){
    direction <- "all negative"
  }else if (sign > 0){
    direction <- "right"
  }else if (sign < 0){
    direction <- "left"
  }
  

  return(c(cutpoint, direction))
}








true.cut.piece.null <- function(tau,lambda0,lambda1,
                                beta0,beta1,t.change0,t.change1){
  mu <- function(x,t,t.change,lambda,beta){
    if (t <= t.change){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change -
                                   lambda[2]*exp(beta*(1-x))*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  
  
  if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))*
      (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) < 0){
    cut_t <- uniroot(function(x){mu(x,tau,t.change1,lambda1,beta1) - mu(x,tau,t.change0,lambda0,beta0)},
                     interval = c(0,1), tol = 1e-10)$root
    if (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0) > 0){
      direction <- "left"}else {
        direction <- "right"}
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0)) >= 0 & 
      (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) >= 
      (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))){
    cut_t <- 0
    direction <- "right"
  }else if ((mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) >= 0 & 
      (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) <= 
      (mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0))){
    cut_t <- 1
    direction <- "left"
  }else if ((mu(0,tau,t.change1,lambda1,beta1) - mu(0,tau,t.change0,lambda0,beta0)) < 0 & 
      (mu(1,tau,t.change1,lambda1,beta1) - mu(1,tau,t.change0,lambda0,beta0)) < 0){
    cut_t <- NA
    direction <- "all negative"
  }
  
  
  
  return(c(cut_t, direction))
}




fit.piece.null.1 <- function(dat,k,t.change,t.analysis){
  dat$X1 <- 1 - dat$X
  if (k == 0){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ X1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }else if (k >= 1){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.change,t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ -1 + interval + X1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }
  return(fit)
}


fit.piece.null.0 <- function(dat,k,t.change,t.analysis){
  if (k == 0){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ 1 + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }else if (k >= 1){
    dat_split <- survSplit(formula = Surv(time, status) ~ ., data = dat, 
                           cut = c(t.change,t.analysis), start="start",end="stop",episode="episode") %>%
      mutate(interval = factor(start),interval_length = stop - start)
    
    fit <- glm(status ~ -1 + interval + offset(log(interval_length)),
               data = dat_split, family = poisson())
  }
  return(fit)
}







# prediction
cut.pred1.null <- function(dat0, dat1, tau, t.analysis){
  # when the number of cutpoint and location are known
  fit0 <- fit.piece.null.0(dat0,0,tau,t.analysis)
  fit1 <- fit.piece.null.1(dat1,0,tau,t.analysis)
  #fit0 <- fit.piece.null(dat0,0,tau,t.analysis)
  #fit1 <- fit.piece.null(dat1,0,tau,t.analysis)
  
  beta1 <- coef(fit1)[2]
  if (abs(beta1/sqrt(vcov(fit1)[2,2])) <= qnorm(0.975)){
    beta1 <- 0
  }
  cutpoint0 <- true.cut.piece.null(tau,
                                  lambda0 = rep(exp(coef(fit0)[1]),2),
                                  lambda1 = rep(exp(coef(fit1)[1]),2),
                                  beta0 = 0,
                                  beta1 = beta1,
                                  t.change0 = tau,
                                  t.change1 = tau)
  return(list(cutpoint0 = cutpoint0, beta1 = beta1))
}







cut.pred3.null <- function(dat0, dat1, tau, t.analysis){
  # when the number of cutpoint is known, but location is unknown
  rej0 <- 1
  k0 <- 0
  while (rej0 == 1 & k0 < 2){
    k0 <- k0 + 1
    fit_piece0 <- piecewiseExp_MLE(dat0$time, dat0$status, k0)
    test0 <- piecewiseExp_test_changepoint(fit_piece0)
    rej0 <- prod(test0$reject)
  }
  k0 <- k0 - 1
  fit_piece0 <- piecewiseExp_MLE(dat0$time, dat0$status, k0)
  t.change0 <- fit_piece0$tau
  fit0 <- fit.piece.null.0(dat0,k0,t.change0,t.analysis)
  
  if (k0 == 0){
    lambda0 <- rep(exp(coef(fit0)[1]),2)
    beta0 <- 0
    t.change0 <- tau
  }else if (k0 == 1){
    lambda0 <- exp(coef(fit0)[1:2])
    beta0 <- 0
    t.change0 <- pmin(t.change0, tau)
  }
  
  
  rej1 <- 1
  k1 <- 0
  while (rej1 == 1 & k1 < 2){
    k1 <- k1 + 1
    fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, k1)
    test1 <- piecewiseExp_test_changepoint(fit_piece1)
    rej1 <- prod(test1$reject)
  }
  k1 <- k1 - 1
  fit_piece1 <- piecewiseExp_MLE(dat1$time, dat1$status, k1)
  t.change1 <- fit_piece1$tau
  fit1 <- fit.piece.null.1(dat1,k1,t.change1,t.analysis)
  
  if (k1 == 0){
    lambda1 <- rep(exp(coef(fit1)[1]),2)
    beta1 <- ifelse(abs(coef(fit1)[2]/sqrt(vcov(fit1)[2,2])) <= qnorm(0.975),
                    0, coef(fit1)[2])
    t.change1 <- tau
  }else if (k1 == 1){
    lambda1 <- exp(coef(fit1)[1:2])
    beta1 <- ifelse(abs(coef(fit1)[3]/sqrt(vcov(fit1)[3,3])) <= qnorm(0.975),
                    0, coef(fit1)[3])
    t.change1 <- pmin(t.change1, tau)
  }
  
  
  
  cutpoint0 <- true.cut.piece.null(tau,
                              lambda0, lambda1,
                              beta0, beta1,
                              t.change0, t.change1)
  
  
  return(list(cutpoint0 = cutpoint0, 
              t.change0 = t.change0, 
              t.change1 = t.change1, 
              k0 = k0, k1 = k1,
              beta1 = beta1))
}









