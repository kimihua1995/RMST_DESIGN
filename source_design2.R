# function for modifying censoring time
dat.modify <- function(dat,t.IA){
  dat$time.censor <- pmin(dat$time.loss, t.IA - dat$time.enroll)
  dat$time <- pmin(dat$time.event, dat$time.censor)
  dat$status <- as.numeric(dat$time.event < dat$time.censor)
  return(dat)
}





####################################################################################
my.design <- function(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                    t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                    tau, pred.method=1, design=1,
                    alpha1, alpha2, true.cutoff){
  
  stat.FA1 <- stat.FA2 <- stat.FA3 <- stat.FA4 <- stat.FA5 <- 
    data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  stat.FA1.true <- stat.FA2.true <- stat.FA3.true <- stat.FA4.true <- stat.FA5.true <- 
    data.frame(z_pos = rep(NA,B), d_pos = rep(NA,B), se_pos = rep(NA,B))
  
  
  cutpoint0.list <- cutpoint.FA.list <- data.frame(cutpoint = rep(NA,B), direction = rep(NA,B))

  change.point <- rep(NA,B)
  change.point0 <- change.point1 <- matrix(NA, B, 2)
  k0.list <- k1.list <- rep(NA,B)
  
  rej1 <- rej2 <- rej3 <- rej4 <- rej5 <- rep(FALSE,B)
  rej1.true <- rej2.true <- rej3.true <- rej4.true <- rej5.true <- rep(FALSE,B)
  gamma3.list <- beta1.list <- rep(NA,B)
  N.neg <- rep(NA,B)
  
  set.seed(seed)
  for (b in 1:B){
    if (!b%%1000) {print(b)}
    # Stage I
    X0_1 <- runif(n1, 0, 1)
    X1_1 <- runif(n1, 0, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda1,beta1,t.change,0,t.enroll1,t.enrich)
    dat0_1$S = 1
    dat1_1$S = 1
    
    if (design == 1){
      X0_2 <- runif(n2, 0, 1)
      X1_2 <- runif(n2, 0, 1)
    }else if (design == 2){
      if (pred.method == 1){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred1(dat0_1, dat1_1, tau, t.change, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {next}
      }else if (pred.method == 2){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred2(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {
          next
        }else {
          change.point[b] <- cut_pred$t.change
        }
      }else if (pred.method == 3){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred3(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {
          next
        }else {
          change.point0[b,] <- cut_pred$t.change0
          change.point1[b,] <- cut_pred$t.change1
          k0.list[b] <- cut_pred$k0
          k1.list[b] <- cut_pred$k1
        }
      }

      cutpoint0 <- cut_pred$cutpoint0
      cutpoint0.list[b,] <- cutpoint0
      beta1.list[b] <- cut_pred$beta1
    
      # Stage II
      if (cutpoint0[2] == "all negative") {
        cutpoint0.list[b,] <- c(NA, "all negative")
        X0_2 <- runif(n2, 0, 1)
        X1_2 <- runif(n2, 0, 1)
      }else if (cutpoint0[2] == "left"){
        X0_2 <- runif(n2, 0, as.numeric(cutpoint0[1]))
        X1_2 <- runif(n2, 0, as.numeric(cutpoint0[1]))
      }else if (cutpoint0[2] == "right"){
        X0_2 <- runif(n2, as.numeric(cutpoint0[1]), 1)
        X1_2 <- runif(n2, as.numeric(cutpoint0[1]), 1)
      }
    }
    
    
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda1,beta1,t.change,t.enrich,t.enroll2,t.FA)
    dat0_2$S = 2
    dat1_2$S = 2
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
    
  
    
    #####################################
    ## FA
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    dat <- rbind(dat0, dat1)
    
    ## estimate the cutoff point at FA
    rmst_fit_FA <- fit.rmst.reg(dat0, dat1, tau)
    gamma3 <- coef(rmst_fit_FA)[4]/sqrt(vcov(rmst_fit_FA)[4,4])
    gamma3.list[b] <- gamma3
    cutpoint.FA <- find.cutpoint(rmst_fit_FA)
    cutpoint.FA.list[b,] <- cutpoint.FA
    
    
    if (cutpoint.FA[2] == "all negative"){
      cutpoint.FA.list[b,] <- c(NA, NA, "all negative")
    }else if (gamma3 <= qnorm(1 - alpha1)){
      stat.FA1[b,] <- stat1(dat0,dat1,0,"right",tau)
      rej1[b] <- stat.FA1[b,1] > qnorm(1 - alpha2)
      
      skip_to_next <- FALSE
      tryCatch(dat_pos <- cw(dat,0,"right"),
               error = function(e){skip_to_next <<- TRUE})
      if (!skip_to_next) {
             stat.FA2[b,] <- stat2(dat_pos,tau)
             rej2[b] <- stat.FA2[b,1] > qnorm(1 - alpha2)
             stat.FA3[b,] <- stat3(rmst_fit_FA, dat_pos)
             rej3[b] <- stat.FA3[b,1] > qnorm(1 - alpha2)
             
             dat_mod <- dat_modify(dat_pos, rmst_fit_FA, tau)
             skip_to_next <- FALSE
             tryCatch(stat.FA4[b,] <- stat4(dat_mod),
                      error = function(e){skip_to_next <<- TRUE})
             if (!skip_to_next) {rej4[b] <- stat.FA4[b,1] > qnorm(1 - alpha2)}
             skip_to_next <- FALSE
             tryCatch(stat.FA5[b,] <- stat5(dat_mod),
                      error = function(e){skip_to_next <<- TRUE})
             if (!skip_to_next) {rej5[b] <- stat.FA5[b,1] > qnorm(1 - alpha2)}
             }
    }else{
      stat.FA1[b,] <- stat1(dat0,dat1,cutpoint.FA[1],cutpoint.FA[2],tau)
      rej1[b] <- stat.FA1[b,1] > qnorm(1 - alpha2)
      
      skip_to_next <- FALSE
      tryCatch(dat_pos <- cw(dat,cutpoint.FA[1],cutpoint.FA[2]),
               error = function(e){skip_to_next <<- TRUE})
      if (!skip_to_next) {
             stat.FA2[b,] <- stat2(dat_pos,tau)
             rej2[b] <- stat.FA2[b,1] > qnorm(1 - alpha2)
             stat.FA3[b,] <- stat3(rmst_fit_FA, dat_pos)
             rej3[b] <- stat.FA3[b,1] > qnorm(1 - alpha2)
             
             dat_mod <- dat_modify(dat_pos, rmst_fit_FA, tau)
             skip_to_next <- FALSE
             tryCatch(stat.FA4[b,] <- stat4(dat_mod),
                      error = function(e){skip_to_next <<- TRUE})
             if (!skip_to_next) {rej4[b] <- stat.FA4[b,1] > qnorm(1 - alpha2)}
             skip_to_next <- FALSE
             tryCatch(stat.FA5[b,] <- stat5(dat_mod),
                      error = function(e){skip_to_next <<- TRUE})
             if (!skip_to_next) {rej5[b] <- stat.FA5[b,1] > qnorm(1 - alpha2)}
             }
    }
    

    N.neg[b] <- sum(dat$X < true.cutoff)

    
    # true.cutoff
    stat.FA1.true[b,] <- stat1(dat0,dat1,true.cutoff,"right",tau)
    rej1.true[b] <- stat.FA1.true[b,1] > qnorm(1 - alpha2)
    
    skip_to_next <- FALSE
    tryCatch(dat_pos.true <- cw(dat,true.cutoff,"right"),
             error = function(e){skip_to_next <<- TRUE})
    if (!skip_to_next) {
      stat.FA2.true[b,] <- stat2(dat_pos.true,tau)
      rej2.true[b] <- stat.FA2.true[b,1] > qnorm(1 - alpha2)
      stat.FA3.true[b,] <- stat3(rmst_fit_FA, dat_pos.true)
      rej3.true[b] <- stat.FA3.true[b,1] > qnorm(1 - alpha2)
      
      dat_mod.true <- dat_modify(dat_pos.true, rmst_fit_FA, tau)
      skip_to_next <- FALSE
      tryCatch(stat.FA4.true[b,] <- stat4(dat_mod.true),
               error = function(e){skip_to_next <<- TRUE})
      if (!skip_to_next) {rej4.true[b] <- stat.FA4.true[b,1] > qnorm(1 - alpha2)}
      skip_to_next <- FALSE
      tryCatch(stat.FA5.true[b,] <- stat5(dat_mod.true),
               error = function(e){skip_to_next <<- TRUE})
      if (!skip_to_next) {rej5.true[b] <- stat.FA5.true[b,1] > qnorm(1 - alpha2)}
    }
    
    
    
  }
  
  return(list(stat.FA1 = stat.FA1, stat.FA2 = stat.FA2, stat.FA3 = stat.FA3, 
              stat.FA4 = stat.FA4, stat.FA5 = stat.FA5, 
              stat.FA1.true = stat.FA1.true, stat.FA2.true = stat.FA2.true, 
              stat.FA3.true = stat.FA3.true, stat.FA4.true = stat.FA4.true, 
              stat.FA5.true = stat.FA5.true, 
              cutpoint0 = cutpoint0.list, cutpoint.FA = cutpoint.FA.list,
              rej1 = rej1, rej2 = rej2, rej3 = rej3,
              rej4 = rej4, rej5 = rej5,
              rej1.true = rej1.true, rej2.true = rej2.true, rej3.true = rej3.true,
              rej4.true = rej4.true, rej5.true = rej5.true,
              gamma3 = gamma3.list, beta1 = beta1.list,
              N.neg = N.neg,
              change.point = change.point,
              change.point0 = change.point0, change.point1 = change.point1,
              k0.list = k0.list, k1.list = k1.list))
  
}







####################################################################################
my.design.null <- function(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                           t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, pred.method=1, design=1,
                           alpha1, alpha2, true.cutoff){
  
  cutpoint0.list <- cutpoint.FA.list <- data.frame(cutpoint = rep(NA,B), direction = rep(NA,B))
  
  change.point <- rep(NA,B)
  change.point0 <- change.point1 <- matrix(NA, B, 2)
  k0.list <- k1.list <- rep(NA,B)
  
  a1 <- length(alpha1); a2 <- length(alpha2)
  rej1 <- rej2 <- rej3 <- rej4 <- rej5 <- matrix(FALSE,B,a1*a2)
  rej1.true <- rej2.true <- rej3.true <- rej4.true <- rej5.true <- matrix(FALSE,B,a2)
  gamma3.list <- beta1.list <- rep(NA,B)
  
  set.seed(seed)
  for (b in 1:B){
    if (!b%%1000) {print(b)}
    #set.seed(b)
    # Stage I
    X0_1 <- runif(n1, 0, 1)
    X1_1 <- runif(n1, 0, 1)
    dat0_1 <- sim.data.piece(X0_1,lambda0,beta0,t.change,0,t.enroll1,t.enrich)
    dat1_1 <- sim.data.piece(X1_1,lambda1,beta1,t.change,0,t.enroll1,t.enrich)
    dat0_1$S = 1
    dat1_1$S = 1
    
    if (design == 1){
      X0_2 <- runif(n2, 0, 1)
      X1_2 <- runif(n2, 0, 1)
    }else if (design == 2){
      if (pred.method == 1){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred1.null(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {next}
      }else if (pred.method == 3){
        skip_to_next <- FALSE
        tryCatch(cut_pred <- cut.pred3.null(dat0_1, dat1_1, tau, t.enrich),
                 error = function(e){skip_to_next <<- TRUE})
        if (skip_to_next) {
          next
        }else {
          change.point0[b,] <- cut_pred$t.change0
          change.point1[b,] <- cut_pred$t.change1
          k0.list[b] <- cut_pred$k0
          k1.list[b] <- cut_pred$k1
        }
      }
      
      cutpoint0 <- cut_pred$cutpoint0
      cutpoint0.list[b,] <- cutpoint0
      beta1.list[b] <- cut_pred$beta1
      
      # Stage II
      if (cutpoint0[2] == "all negative") {
        cutpoint0.list[b,] <- c(NA, "all negative")
        #cutpoint.FA.list[b,] <- c(NA, NA, "all negative")
        #next
        X0_2 <- runif(n2, 0, 1)
        X1_2 <- runif(n2, 0, 1)
      }else if (cutpoint0[2] == "left"){
        X0_2 <- runif(n2, 0, as.numeric(cutpoint0[1]))
        X1_2 <- runif(n2, 0, as.numeric(cutpoint0[1]))
      }else if (cutpoint0[2] == "right"){
        X0_2 <- runif(n2, as.numeric(cutpoint0[1]), 1)
        X1_2 <- runif(n2, as.numeric(cutpoint0[1]), 1)
      }
    }
    
    
    
    dat0_2 <- sim.data.piece(X0_2,lambda0,beta0,t.change,t.enrich,t.enroll2,t.FA)
    dat1_2 <- sim.data.piece(X1_2,lambda1,beta1,t.change,t.enrich,t.enroll2,t.FA)
    dat0_2$S = 2
    dat1_2$S = 2
    dat0 <- rbind(dat0_1, dat0_2); dat0$A = 0
    dat1 <- rbind(dat1_1, dat1_2); dat1$A = 1
    
    
    
    #####################################
    ## FA
    dat0 <- dat.modify(dat0, t.FA)
    dat1 <- dat.modify(dat1, t.FA)
    dat <- rbind(dat0, dat1)
    
    ## estimate the cutoff point at FA
    rmst_fit_FA <- fit.rmst.reg(dat0, dat1, tau)
    gamma3 <- coef(rmst_fit_FA)[4]/sqrt(vcov(rmst_fit_FA)[4,4])
    gamma3.list[b] <- gamma3
    cutpoint.FA <- find.cutpoint(rmst_fit_FA)
    cutpoint.FA.list[b,] <- cutpoint.FA
    
    
    ## true cutoff
    stat.FA1.true <- stat1(dat0,dat1,true.cutoff,"right",tau)
    for (j in 1:a2){
      rej1.true[b,j] <- stat.FA1.true[1] > qnorm(1-alpha2[j])
    }
    skip_to_next <- FALSE
    tryCatch(dat_pos.true <- cw(dat,true.cutoff,"right"),
             error = function(e){skip_to_next <<- TRUE})
    if (!skip_to_next) {
      stat.FA2.true <- stat2(dat_pos.true,tau)
      stat.FA3.true <- stat3(rmst_fit_FA, dat_pos.true)
      dat_mod.true <- dat_modify(dat_pos.true, rmst_fit_FA, tau)
      skip_to_next <- FALSE
      tryCatch(stat.FA4.true <- stat4(dat_mod.true),
               error = function(e){skip_to_next <<- TRUE})
      if (skip_to_next) {stat.FA4.true <- c(NA,NA,NA)}
      skip_to_next <- FALSE
      tryCatch(stat.FA5.true <- stat5(dat_mod.true),
               error = function(e){skip_to_next <<- TRUE})
      if (skip_to_next) {stat.FA5.true <- c(NA,NA,NA)}
      
      for (j in 1:a2){
        rej2.true[b,j] <- stat.FA2.true[1] > qnorm(1-alpha2[j])
        rej3.true[b,j] <- stat.FA3.true[1] > qnorm(1-alpha2[j])
        rej4.true[b,j] <- stat.FA4.true[1] > qnorm(1-alpha2[j])
        rej5.true[b,j] <- stat.FA5.true[1] > qnorm(1-alpha2[j])
      }
    }
    
    
    
    
    #######################################
    stat.FA1.all <- stat1(dat0,dat1,0,"right",tau)
    dat_pos.all <- cw(dat,0,"right")
    stat.FA2.all <- stat2(dat_pos.all,tau)
    stat.FA3.all <- stat3(rmst_fit_FA, dat_pos.all)
    dat_mod.all <- dat_modify(dat_pos.all, rmst_fit_FA, tau)
    skip_to_next <- FALSE
    tryCatch(stat.FA4.all <- stat4(dat_mod.all),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {stat.FA4.all <- c(NA,NA,NA)}
    skip_to_next <- FALSE
    tryCatch(stat.FA5.all <- stat5(dat_mod.all),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {stat.FA5.all <- c(NA,NA,NA)}
    
    
    
    for (i in 1:a1){
      if (cutpoint.FA[2] == "all negative"){
        cutpoint.FA.list[b,] <- c(NA, "all negative")
      }else if (gamma3 <= qnorm(1 - alpha1[i])){
        for (j in 1:a2){
          rej1[b,a2*(i-1)+j] <- stat.FA1.all[1] > qnorm(1 - alpha2[j])
          rej2[b,a2*(i-1)+j] <- stat.FA2.all[1] > qnorm(1 - alpha2[j])
          rej3[b,a2*(i-1)+j] <- stat.FA3.all[1] > qnorm(1 - alpha2[j])
          rej4[b,a2*(i-1)+j] <- stat.FA4.all[1] > qnorm(1 - alpha2[j])
          rej5[b,a2*(i-1)+j] <- stat.FA5.all[1] > qnorm(1 - alpha2[j])
        }
      }else{
        stat.FA1 <- stat1(dat0,dat1,cutpoint.FA[1],cutpoint.FA[2],tau)
        for (j in 1:a2){
          rej1[b,a2*(i-1)+j] <- stat.FA1[1] > qnorm(1 - alpha2[j])
        }
        
        skip_to_next <- FALSE
        tryCatch(dat_pos <- cw(dat,cutpoint.FA[1],cutpoint.FA[2]),
                 error = function(e){skip_to_next <<- TRUE})
        if (!skip_to_next) {
          stat.FA2 <- stat2(dat_pos,tau)
          stat.FA3 <- stat3(rmst_fit_FA, dat_pos)
          dat_mod <- dat_modify(dat_pos, rmst_fit_FA, tau)
          skip_to_next <- FALSE
          tryCatch(stat.FA4 <- stat4(dat_mod),
                   error = function(e){skip_to_next <<- TRUE})
          if (skip_to_next) {stat.FA4 <- c(NA,NA,NA)}
          skip_to_next <- FALSE
          tryCatch(stat.FA5 <- stat5(dat_mod),
                   error = function(e){skip_to_next <<- TRUE})
          if (skip_to_next) {stat.FA5 <- c(NA,NA,NA)}
          
          for (j in 1:a2){
            rej2[b,a2*(i-1)+j] <- stat.FA2[1] > qnorm(1-alpha2[j])
            rej3[b,a2*(i-1)+j] <- stat.FA3[1] > qnorm(1-alpha2[j])
            rej4[b,a2*(i-1)+j] <- stat.FA4[1] > qnorm(1-alpha2[j])
            rej5[b,a2*(i-1)+j] <- stat.FA5[1] > qnorm(1-alpha2[j])
          }
          
        }
      }
    }
    
    
    
    
    
  }
  
  return(list(cutpoint0 = cutpoint0.list, cutpoint.FA = cutpoint.FA.list,
              rej1 = rej1, rej2 = rej2, rej3 = rej3,
              rej4 = rej4, rej5 = rej5,
              rej1.true = rej1.true, rej2.true = rej2.true, 
              rej3.true = rej3.true, rej4.true = rej4.true, 
              rej5.true = rej5.true,
              gamma3 = gamma3.list, beta1 = beta1.list,
              change.point = change.point,
              change.point0 = change.point0, change.point1 = change.point1,
              k0.list = k0.list, k1.list = k1.list))
  
}




