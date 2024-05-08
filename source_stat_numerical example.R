my.func_surv <- function(y, d){
  #--input--
  #y=time
  #d=status
  
  #--
  id=order(y)
  y=y[id]
  d=d[id]
  
  #--
  t_idx = unique(c(0,y))
  ny = length(y)
  
  #--
  Y = N = C = S = H = D = E = rep(0,length(t_idx))
  
  #i=1
  Y[1] = ny
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0
  
  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]
    N[i] = ifelse(sum(y==t_idx[i] & d==1)>0, sum(y==t_idx[i] & d==1), 0)
    C[i] = ifelse(sum(y==t_idx[i] & d==0)>0, sum(y==t_idx[i] & d==0), 0)
    
    if(Y[i]<0){Y[i] = 0}
    
    S[i] = ifelse(Y[i]==0, S[i-1], S[i-1]*(1-(N[i]/Y[i])))
    H[i] = ifelse(Y[i]*(Y[i]-N[i])==0, 0, N[i]/(Y[i]*(Y[i]-N[i])))
    
    if(S[i]<0){S[i] = 0}
    
    D[i] = sum(H[2:i])
    E[i] = sqrt((S[i]**2)*D[i])
    
    if(is.na(S[i])){S[i] = 0}
    if(is.na(E[i])){E[i] = 0}
  }
  
  #--output--
  out           = as.data.frame(cbind(t_idx, Y, N, C, S, E))
  colnames(out) = c("t_idx", "n_risk", "n_event", "n_censor", "surv", "se")
  
  #--to match the output of survfit--
  out2 = out[t_idx!=0,]
  
  #--
  Z2 = list()
  Z2$out      = out2
  Z2$t_idx    = out2[,"t_idx"]
  Z2$n_risk   = out2[,"n_risk"]
  Z2$n_event  = out2[,"n_event"]
  Z2$n_censor = out2[,"n_censor"]
  Z2$surv     = out2[,"surv"]
  Z2$se       = out2[,"se"]
  
  return(Z2)
}



my.rmst2reg=function(y, delta, arm, x, tau, w=rep(1,length(y))){
  
  n=length(y)
  x=as.matrix(cbind(1, x))
  p=length(x[1,])
  
  y0=pmin(y, tau)
  d0=delta
  d0[y0==tau]=1
  
  d10=d0[arm==1]
  d00=d0[arm==0]
  y10=y0[arm==1]
  y00=y0[arm==0]
  x1=x[arm==1,]
  x0=x[arm==0,]
  n1=length(d10)
  n0=length(d00)
  
  id1=order(y10)
  y10=y10[id1]
  d10=d10[id1]
  x1=x1[id1,]
  
  id0=order(y00)
  y00=y00[id0]
  d00=d00[id0]
  x0=x0[id0,]
  
  fitc1=my.func_surv(y10, 1-d10)
  fitc0=my.func_surv(y00, 1-d00)
  
  weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
  weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))
  
  w1=w[arm==1]
  w0=w[arm==0]
  w1=w1[id1]
  w0=w0[id0]
  weights=c(weights1, weights0)*c(w1,w0)
  
  
  fitt=lm(c(y10,y00)~ rbind(x1, x0)-1, weights=weights)
  
  return(fitt)
}




fit.rmst.reg <- function(dat0,dat1,tau){
  tau.temp <- min(max(dat0$time), max(dat1$time))
  tau.temp <- min(tau.temp, tau)
  
  dat0$A=0
  dat1$A=1
  dat <- rbind(dat0,dat1)
  cov <- data.frame(A=dat$A, X=dat$X, AX=dat$A*dat$X)
  rmst_fit <- my.rmst2reg(y = dat$time,
                          delta = dat$status,
                          x = cov,
                          arm = dat$A,
                          tau = tau.temp)
  return(rmst_fit)
}


find.cutpoint <- function(rmst_fit){
  
  gamma <- coef(rmst_fit)
  cutpoint <- -gamma[2]/gamma[4]
  cutpoint <- min(max(cutpoint,0),1)
  sign <- gamma[4]
  
  
  if ((cutpoint == 0 & sign < 0) | (cutpoint == 1 & sign > 0)){
    direction <- "all negative"
  }else if (sign > 0){
    direction <- "right"
  }else if (sign < 0){
    direction <- "left"
  }
  
  
  
  return(c(cutpoint, direction))
}



stat1 <- function(dat0,dat1,tau){
  time_pos=c(dat0$time, dat1$time)
  status_pos=c(dat0$status, dat1$status)
  arm_pos=c(dat0$A, dat1$A)
  
  fit_pos=rmst2(time_pos, status_pos, arm_pos, tau=tau)
  rmst1est_pos=fit_pos$RMST.arm1$rmst[1:2]
  rmst0est_pos=fit_pos$RMST.arm0$rmst[1:2]
  d_pos=(rmst1est_pos[1]-rmst0est_pos[1])
  se_pos=sqrt(rmst1est_pos[2]^2+rmst0est_pos[2]^2)
  z_pos=d_pos/se_pos
  
  return(as.numeric(c(z_pos, d_pos, se_pos)))
}



stat3 <- function(rmst_fit,dat){
  gamma <- coef(rmst_fit)
  vgamma <- vcov(rmst_fit)
  
  d_pos <- mean(gamma[2] + gamma[4]*dat$X)
  J <- matrix(c(1, mean(dat$X)), nrow = 2, ncol = 1)
  se_pos <- as.numeric(sqrt(t(J) %*% vgamma[c(2,4),c(2,4)] %*% J))
  
  z_pos <- d_pos/se_pos
  
  
  return(as.numeric(c(z_pos,d_pos,se_pos)))
}


dat_modify <- function(dat, rmst_fit, tau){
  A <- dat$A
  y0=pmin(dat$time, tau)
  d0=dat$status
  d0[y0==tau]=1
  
  d10=d0[A==1]; d00=d0[A==0]
  y10=y0[A==1]; y00=y0[A==0]
  
  id1=order(y10); y10=y10[id1]; d10=d10[id1]
  id0=order(y00); y00=y00[id0]; d00=d00[id0]
  
  fitc1=my.func_surv(y10, 1-d10)
  fitc0=my.func_surv(y00, 1-d00)
  
  weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
  weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))
  
  gamma <- coef(rmst_fit)
  mu1 <- gamma[1] + gamma[2] + (gamma[3]+gamma[4])*dat$X
  mu0 <- gamma[1] + gamma[3]*dat$X
  
  mu11 <- mu1[A==1]; mu11 <- mu11[id1]
  mu10 <- mu1[A==0]; mu10 <- mu10[id0]
  mu01 <- mu0[A==1]; mu01 <- mu01[id1]
  mu00 <- mu0[A==0]; mu00 <- mu00[id0]
  
  dat_mod <- data.frame(weights = c(weights1, weights0),
                        Y = c(y10, y00),
                        A = c(rep(1,sum(A==1)), rep(0,sum(A==0))),
                        mu1 = c(mu11, mu10),
                        mu0 = c(mu01, mu00))
  return(dat_mod)
}



stat4 <- function(dat_mod){
  
  m_fun <- function(data){
    A <- data$A
    Y <- data$Y
    w <- data$weights
    function(theta){
      c(A*w*(Y-theta[1]),
        (1-A)*w*(Y-theta[2]))
    }
  }
  
  res_RMSTD <- geex::m_estimate(
    estFUN = m_fun,
    data = dat_mod,
    root_control = setup_root_control(start = c(0.5,0.5))
  )
  
  c1 <- matrix(c(1,-1), nrow = 1, ncol = 2)
  d_pos <- as.numeric(c1 %*% coef(res_RMSTD))
  se_pos <- as.numeric(sqrt(c1%*%vcov(res_RMSTD)%*%t(c1)))
  z_pos <- d_pos/se_pos
  
  return(as.numeric(c(z_pos, d_pos, se_pos)))
  
}






stat5 <- function(dat_mod){
  
  m_fun <- function(data){
    A <- data$A
    Y <- data$Y
    mu1 <- data$mu1
    mu0 <- data$mu0
    w <- data$weights
    function(theta){
      c(A*w*(Y-mu1-theta[1]),
        (1-A)*w*(Y-mu0-theta[2]),
        mu1-mu0-theta[3])
    }
  }
  
  res_RMSTD <- geex::m_estimate(
    estFUN = m_fun,
    data = dat_mod,
    root_control = setup_root_control(start = c(0.5,0.5,0.5))
  )
  
  c2 <- matrix(c(1,-1,1), nrow = 1, ncol = 3)
  d_pos <- as.numeric(c2 %*% coef(res_RMSTD))
  se_pos <- as.numeric(sqrt(c2%*%vcov(res_RMSTD)%*%t(c2)))
  z_pos <- d_pos/se_pos
  
  return(as.numeric(c(z_pos, d_pos, se_pos)))
  
}