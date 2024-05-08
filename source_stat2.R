# RMST estimate in positive group
stat1 <- function(dat0,dat1,cutpoint,pos,tau){
  if (pos == "right"){
    dat0_pos <- dat0[dat0$X >= cutpoint,]
    dat1_pos <- dat1[dat1$X >= cutpoint,]
  }else if (pos == "left"){
    dat0_pos <- dat0[dat0$X <= cutpoint,]
    dat1_pos <- dat1[dat1$X <= cutpoint,]
  }
  
  time_pos=c(dat0_pos$time, dat1_pos$time)
  status_pos=c(dat0_pos$status, dat1_pos$status)
  arm_pos=c(rep(0,nrow(dat0_pos)), rep(1,nrow(dat1_pos)))
  
  if (min(max(dat0_pos$time), max(dat1_pos$time)) < tau) {
    tau.temp <- min(max(dat0_pos$time), max(dat1_pos$time))
  }else {tau.temp <- tau}
  
  fit_pos=rmst2(time_pos, status_pos, arm_pos, tau=tau.temp)
  rmst1est_pos=fit_pos$RMST.arm1$rmst[1:2]
  rmst0est_pos=fit_pos$RMST.arm0$rmst[1:2]
  d_pos=(rmst1est_pos[1]-rmst0est_pos[1])
  se_pos=sqrt(rmst1est_pos[2]^2+rmst0est_pos[2]^2)
  z_pos=d_pos/se_pos
  

  
  return(as.numeric(c(z_pos,d_pos,se_pos)))
}





cw <- function(dat,cutpoint,pos){
  cutpoint <- as.numeric(cutpoint)
  
  if (pos == "right"){
    g <- c((1+cutpoint)/2, (1-cutpoint)^2/12 + (1+cutpoint)^2/4)
    dat_pos <- dat[dat$X >= cutpoint,]
  }else if (pos == "left"){
    g <- c(cutpoint/2, cutpoint^2/12 + cutpoint^2/4)
    dat_pos <- dat[dat$X <= cutpoint,]
  }
  
  dat_pos <- dat_pos %>% arrange(time)
  
  fn <- function(x){
    y1 <- sum(exp(x[1]*dat_pos$X + x[2]*(dat_pos$X)^2) * 
                (dat_pos$X - g[1]))
    y2 <- sum(exp(x[1]*dat_pos$X + x[2]*(dat_pos$X)^2) * 
                ((dat_pos$X)^2 - g[2]))
    return(c(y1, y2)) 
  }
  
  sol_lambda <- nleqslv(c(1,1), fn,
                        control = list(btol = 1e-5))
  lambda <- sol_lambda$x
  p <- exp(lambda[1]*dat_pos$X + lambda[2]*(dat_pos$X)^2)
  p <- p/sum(p)
  
  dat_pos$p <- p
  
  return(dat_pos)
}



cw2 <- function(dat,cutpoint,pos){
  cutpoint <- as.numeric(cutpoint)
  dat_1 <- dat[dat$S==1,]
  if (pos == "right"){
    dat_pos <- dat[dat$X >= cutpoint,]
    dat_1_pos <- dat_1[dat_1$X >= cutpoint,]
  }else if (pos == "left"){
    dat_pos <- dat[dat$X <= cutpoint,]
    dat_1_pos <- dat_1[dat_1$X <= cutpoint,]
  }
  
  dat_pos <- dat_pos %>% arrange(time)
  g <- c(mean(dat_1_pos$X), mean(dat_1_pos$X^2))
  
  
  fn <- function(x){
    y1 <- sum(exp(x[1]*dat_pos$X + x[2]*(dat_pos$X)^2) * 
                (dat_pos$X - g[1]))
    y2 <- sum(exp(x[1]*dat_pos$X + x[2]*(dat_pos$X)^2) * 
                ((dat_pos$X)^2 - g[2]))
    return(c(y1, y2)) 
  }
  
  sol_lambda <- nleqslv(c(1,1), fn,
                        control = list(btol = 1e-5))
  lambda <- sol_lambda$x
  p <- exp(lambda[1]*dat_pos$X + lambda[2]*(dat_pos$X)^2)
  p <- p/sum(p)
  
  dat_pos$p <- p
  
  return(dat_pos)
}




my_akm_rmst <- function(time, status, weight=NULL, tau=NULL){
  
  data <- data.frame(time, status, weight) %>% arrange(time)
  #--- AKM ---
  # Based on 'adjusted.KM' function from {IPWsurvival} package
  # Author: F. Le Borgne and Y. Foucher
  tj <- c(0,sort(unique(data$time[data$status==1])))
  dj <- sapply(tj, function(x){sum(data$weight[data$time==x & data$status==1])})
  yj <- sapply(tj, function(x){sum(data$weight[data$time>=x])})
  st <- cumprod(1-(dj/yj))
  m <- sapply(tj, function(x){sum((data$weight[data$time>=x])^2)})
  mj <- ((yj^2)/m)
  #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
  ft <- data.frame(tj, yj, dj, st, mj)
  
  #--- RMST ---
  # Based on 'rmst1 function' from {survRM2} package
  # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
  rtime <- ft$tj<=tau
  tj_r <- sort(c(ft$tj[rtime],tau))
  st_r <- ft$st[rtime]
  yj_r <- ft$yj[rtime]
  dj_r <- ft$dj[rtime]
  time_diff <- diff(c(0, tj_r))
  areas <- time_diff * c(1, st_r)
  rmst <- sum(areas)
  
  #--- Variance ---
  mj_r <- ft$mj[rtime]
  var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
  #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
  var_r <- c(var_r,0)
  rmst_var <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
  
  return(data.frame(mu = rmst, V = rmst_var))     
}



stat2 <- function(dat_pos, tau){
  p <- dat_pos$p
  A <- dat_pos$A
  
  mu1 <- my_akm_rmst(dat_pos$time[A==1], 
                     dat_pos$status[A==1], 
                     p[A==1], tau)
  mu0 <- my_akm_rmst(dat_pos$time[A==0], 
                     dat_pos$status[A==0], 
                     p[A==0], tau)
  
  d_pos <- mu1$mu - mu0$mu
  se_pos <- sqrt(mu1$V + mu0$V)
  z_pos <- d_pos/se_pos
  
  return(as.numeric(c(z_pos,d_pos,se_pos)))
  
}



stat3 <- function(rmst_fit,dat_pos){
  gamma <- coef(rmst_fit)
  vgamma <- vcov(rmst_fit)
  p <- dat_pos$p
  
  d_pos <- sum(p*(gamma[2] + gamma[4]*dat_pos$X))
  J <- matrix(c(1, sum(p*dat_pos$X)), nrow = 2, ncol = 1)
  se_pos <- as.numeric(sqrt(t(J) %*% vgamma[c(2,4),c(2,4)] %*% J))
  
  z_pos <- d_pos/se_pos
  
  
  return(as.numeric(c(z_pos,d_pos,se_pos)))
}





dat_modify <- function(dat_pos, rmst_fit, tau){
  A <- dat_pos$A
  y0=pmin(dat_pos$time, tau)
  d0=dat_pos$status
  d0[y0==tau]=1
  
  d10=d0[A==1]; d00=d0[A==0]
  y10=y0[A==1]; y00=y0[A==0]
  
  id1=order(y10); y10=y10[id1]; d10=d10[id1]
  id0=order(y00); y00=y00[id0]; d00=d00[id0]
  
  fitc1=my.func_surv(y10, 1-d10)
  fitc0=my.func_surv(y00, 1-d00)
  
  weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
  weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))
  
  p1 <- dat_pos$p[A==1]; p1 <- p1[id1]
  p0 <- dat_pos$p[A==0]; p0 <- p0[id0]
  
  
  gamma <- coef(rmst_fit)
  mu1 <- gamma[1] + gamma[2] + (gamma[3]+gamma[4])*dat_pos$X
  mu0 <- gamma[1] + gamma[3]*dat_pos$X
  
  mu11 <- mu1[A==1]; mu11 <- mu11[id1]
  mu10 <- mu1[A==0]; mu10 <- mu10[id0]
  mu01 <- mu0[A==1]; mu01 <- mu01[id1]
  mu00 <- mu0[A==0]; mu00 <- mu00[id0]
  
  dat_mod <- data.frame(p = c(p1, p0),
                        weights = c(weights1, weights0),
                        Y = c(y10, y00),
                        A = c(rep(1,length(p1)), rep(0,length(p0))),
                        mu1 = c(mu11, mu10),
                        mu0 = c(mu01, mu00))
  return(dat_mod)
}






stat4 <- function(dat_mod){
  
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$weights
    function(theta){
      c(p*A*w*(Y-theta[1]),
        p*(1-A)*w*(Y-theta[2]))
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
    p <- data$p
    A <- data$A
    Y <- data$Y
    mu1 <- data$mu1
    mu0 <- data$mu0
    w <- data$weights
    function(theta){
      c(p*A*w*(Y-mu1-theta[1]),
        p*(1-A)*w*(Y-mu0-theta[2]),
        p*(mu1-mu0-theta[3]))
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




