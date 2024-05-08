sim.data.piece <- function(X,lambda,beta,t.change,t.enroll.start,t.enroll.end,t.IA){
  # administrative censoring only
  n <- length(X)
  u <- runif(n,0.00001,1)
  time.event <- -log(u)*exp(-(1-X)*beta)/lambda[1]
  pos <- time.event > t.change
  #time.event[pos] <- (-log(u[pos])-lambda[1]*exp(beta*X[pos])*t.change+lambda[2]*exp(beta*X[pos])*t.change)*
  #  exp(-X[pos]*beta)/lambda[2]
  time.event[pos] <- -log(u[pos])*exp(-(1-X[pos])*beta)/lambda[2] + (1 - lambda[1]/lambda[2])*t.change
  time.enroll <- runif(n, t.enroll.start, t.enroll.end)
  time.loss <- rexp(n)/0.12
  #time.loss <- rep(10000,n)
  time.censor <- pmin(time.loss, t.IA - time.enroll)
  time <- pmin(time.event, time.censor)
  status <- as.numeric(time.event < time.censor)
  
  dat <- data.frame(X = X, time.event = time.event, time.enroll = time.enroll,
                    time.loss = time.loss, time.censor = time.censor, time = time, status = status)
  return(dat)
}




true.cut.piece0 <- function(tau,lambda0,lambda1,beta0,beta1,t.change){
  mu0 <- function(x,t){
    if (t <= t.change){
      integrate(function(y){exp(-lambda0[1]*exp(beta0*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda0[1]*exp(beta0*(1-x))*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda0[1]*exp(beta0*(1-x))*t.change-lambda0[2]*exp(beta0*(1-x))*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  mu1 <- function(x,t){
    if (t <= t.change){
      integrate(function(y){exp(-lambda1[1]*exp(beta1*(1-x))*y)},
                lower = 0, upper = t, rel.tol = 1e-10)$value
    }else if (t > t.change){
      integrate(function(y){exp(-lambda1[1]*exp(beta1*(1-x))*y)},
                lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
        integrate(function(y){exp(-lambda1[1]*exp(beta1*(1-x))*t.change-lambda1[2]*exp(beta1*(1-x))*(y-t.change))},
                  lower = t.change, upper = t, rel.tol = 1e-10)$value
    }
  }
  
  if ((mu1(0,tau) - mu0(0,tau))*(mu1(1,tau) - mu0(1,tau)) < 0){
    cut_t <- uniroot(function(x){mu1(x,tau)-mu0(x,tau)},
                     interval = c(0,1), tol = 1e-10)$root
    if ((mu1(0,tau) - mu0(0,tau)) > 0){
      direction <- "left"}else {
        direction <- "right"}
  }
  
  if ((mu1(0,tau) - mu0(0,tau)) > 0 & 
      (mu1(1,tau) - mu0(1,tau)) > (mu1(0,tau) - mu0(0,tau))){
    cut_t <- 0
    direction <- "right"
  }
  if ((mu1(1,tau) - mu0(1,tau)) > 0 & 
            (mu1(1,tau) - mu0(1,tau)) < (mu1(0,tau) - mu0(0,tau))){
    cut_t <- 1
    direction <- "left"
  }
  if ((mu1(0,tau) - mu0(0,tau)) < 0 & 
     (mu1(1,tau) - mu0(1,tau)) < 0){
    cut_t <- NA
    direction <- "all negative"
  }
  
  
  
  #skip_to_next <- FALSE
  #tryCatch(cut_t <- uniroot(function(x){mu1(x,t)-mu0(x,t)},
  #                          interval = c(0,1), tol = 1e-10)$root,
  #         error = function(e){skip_to_next <<- TRUE})
  #if (skip_to_next) {cut_t <- NA}
  return(c(cut_t, direction))
}






rmst.true <- function(lambda, beta, t.change, tau, lower, upper){
  mu <- function(x, tau){
    ifelse(tau <= t.change,
          integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                     lower = 0, upper = tau, rel.tol = 1e-10)$value,
          integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*y)},
                     lower = 0, upper = t.change, rel.tol = 1e-10)$value + 
          integrate(function(y){exp(-lambda[1]*exp(beta*(1-x))*t.change-lambda[2]*exp(beta*(1-x))*(y-t.change))},
                       lower = t.change, upper = tau, rel.tol = 1e-10)$value)
  }
  
  integrate(Vectorize(function(x){mu(x,tau)/(upper-lower)}), lower = lower, upper = upper, rel.tol = 1e-10)$value
}








