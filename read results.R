library(writexl)


# alternative setting
load("resa_pred1_design1_n1_505_n2_505_seed333.RData"); resa_ds1 <- res
load("resa_pred1_design2_n1_505_n2_505_seed333.RData"); resa_ds2 <- res
load("resa_pred2_design2_n1_505_n2_505_seed333.RData"); resa_ds3 <- res
load("resa_pred3_design2_n1_505_n2_505_seed333.RData"); resa_ds4 <- res

RMSTD_pos <- 0.1339246

res_alt <- function(res){
  alpha1 <- 0.025
  B <- nrow(res$cutpoint0)
  # all negative
  #all_neg <- data.frame(stage = c(1,2),
  #                all_neg = c(sum(res$cutpoint0[,2] == "all negative")/B,
  #                sum(res$cutpoint.FA[,3] == "all negative")/B))
  
  # biomarker cutpoint
  bio_cut <- data.frame(stage = c(1,2),
            mean = c(mean(as.numeric(res$cutpoint0[,1]), na.rm = T),
                     mean(as.numeric(res$cutpoint.FA[,1]), na.rm = T)),
            sd = c(sd(as.numeric(res$cutpoint0[,1]), na.rm = T),
                   sd(as.numeric(res$cutpoint.FA[,1]), na.rm = T)))

  
  # significant interaction
  #pos_stageII <- which(res$cutpoint0[,2] != "all negative")
  pos <- which(res$gamma3 > qnorm(1-alpha1))
  #sig_int <- data.frame(scenario = c("overall","not terminate"),
  #           int = c(length(pos)/B,
  #           length(intersect(pos,pos_stageII))/length(pos_stageII)))
  sig_int <- length(pos)/B 
  
  # power
  power <- data.frame(method = 1:5,
                      est.cut = c(mean(res$rej1, na.rm = T),
                                  mean(res$rej2, na.rm = T),
                                  mean(res$rej3, na.rm = T),
                                  mean(res$rej4, na.rm = T),
                                  mean(res$rej5, na.rm = T)),
                      true.cut = c(mean(res$rej1.true, na.rm = T),
                                   mean(res$rej2.true, na.rm = T),
                                   mean(res$rej3.true, na.rm = T),
                                   mean(res$rej4.true, na.rm = T),
                                   mean(res$rej5.true, na.rm = T)))
  
  
  # RMST difference in positive patients
  CP <- function(d_pos, se_pos){
    d_l <- d_pos - qnorm(0.975)*se_pos
    d_h <- d_pos + qnorm(0.975)*se_pos
    cp <- (d_l <= RMSTD_pos) * (d_h >= RMSTD_pos)
    return(mean(cp, na.rm=T))
  }
  RMSTD <- data.frame(method = 1:5,
           est.delta = c(mean(res$stat.FA1$d_pos, na.rm = T),
                         mean(res$stat.FA2$d_pos, na.rm = T),
                         mean(res$stat.FA3$d_pos, na.rm = T),
                         mean(res$stat.FA4$d_pos, na.rm = T),
                         mean(res$stat.FA5$d_pos, na.rm = T)),
           est.cp = c(CP(res$stat.FA1$d_pos, res$stat.FA1$se_pos),
                      CP(res$stat.FA2$d_pos, res$stat.FA2$se_pos),
                      CP(res$stat.FA3$d_pos, res$stat.FA3$se_pos),
                      CP(res$stat.FA4$d_pos, res$stat.FA4$se_pos),
                      CP(res$stat.FA5$d_pos, res$stat.FA5$se_pos)),
           true.delta = c(mean(res$stat.FA1.true$d_pos, na.rm = T),
                          mean(res$stat.FA2.true$d_pos, na.rm = T),
                          mean(res$stat.FA3.true$d_pos, na.rm = T),
                          mean(res$stat.FA4.true$d_pos, na.rm = T),
                          mean(res$stat.FA5.true$d_pos, na.rm = T)),
           true.sd.delta = c(mean(res$stat.FA1.true$se_pos, na.rm = T),
                             mean(res$stat.FA2.true$se_pos, na.rm = T),
                             mean(res$stat.FA3.true$se_pos, na.rm = T),
                             mean(res$stat.FA4.true$se_pos, na.rm = T),
                             mean(res$stat.FA5.true$se_pos, na.rm = T)),
           true.delta.se = c(sd(res$stat.FA1.true$d_pos, na.rm = T),
                             sd(res$stat.FA2.true$d_pos, na.rm = T),
                             sd(res$stat.FA3.true$d_pos, na.rm = T),
                             sd(res$stat.FA4.true$d_pos, na.rm = T),
                             sd(res$stat.FA5.true$d_pos, na.rm = T)),
           true.cp = c(CP(res$stat.FA1.true$d_pos, res$stat.FA1.true$se_pos),
                       CP(res$stat.FA2.true$d_pos, res$stat.FA2.true$se_pos),
                       CP(res$stat.FA3.true$d_pos, res$stat.FA3.true$se_pos),
                       CP(res$stat.FA4.true$d_pos, res$stat.FA4.true$se_pos),
                       CP(res$stat.FA5.true$d_pos, res$stat.FA5.true$se_pos))
           )
  
  # enrolled negative patients
  N.neg <- mean(res$N.neg, na.rm = T)
  
  change.point <- mean(res$change.point, na.rm = T)
  
  k.list <- rbind(1:3,table(res$k0.list),table(res$k1.list))
  
  return(list(bio_cut = bio_cut, 
              sig_int = sig_int,
              power = power, RMSTD = RMSTD, N.neg = N.neg,
              change.point = change.point,
              k.list = k.list))
}


res_alt(resa_ds1)
res_alt(resa_ds2)
res_alt(resa_ds3)
res_alt(resa_ds4)








# global null
load("resnull_pred1_design1_n1_505_n2_505_seed333.RData"); res0_ds1 <- res0
load("resnull_pred1_design2_n1_505_n2_505_seed333.RData"); res0_ds2 <- res0
load("resnull_pred3_design2_n1_505_n2_505_seed333.RData"); res0_ds4 <- res0

res_null <- function(res_0){
  alpha1 <- c(0.015,0.020,0.025,0.035)
  alpha2 <- c(0.025,0.023,0.020)
  alpha1.list <- rep(alpha1, each=length(alpha2))
  alpha2.list <- rep(alpha2,length(alpha1))
  B <- nrow(res_0$cutpoint0)
  
  tab20 <- data.frame(alpha1 = alpha1.list,
                      alpha2 = alpha2.list,
                      method1 = colMeans(res_0$rej1, na.rm = T),
                      method2 = colMeans(res_0$rej2, na.rm = T),
                      method3 = colMeans(res_0$rej3, na.rm = T),
                      method4 = colMeans(res_0$rej4, na.rm = T),
                      method5 = colMeans(res_0$rej5, na.rm = T))
  tab20_true <- data.frame(alpha2 = alpha2,
                      method1 = colMeans(res_0$rej1.true, na.rm = T),
                      method2 = colMeans(res_0$rej2.true, na.rm = T),
                      method3 = colMeans(res_0$rej3.true, na.rm = T),
                      method4 = colMeans(res_0$rej4.true, na.rm = T),
                      method5 = colMeans(res_0$rej5.true, na.rm = T))
  
  tab20_mix <- data.frame(alpha1 = alpha1.list,
                          alpha2 = alpha2.list,
                          method1 = NA,
                          method2 = NA,
                          method3 = NA,
                          method4 = NA,
                          method5 = NA)
  for (i in 1:length(alpha1.list)){
    pos <- which(res_0$gamma3 > qnorm(1-alpha1.list[i]))
    j <- which(alpha2 == alpha2.list[i])
    tab20_mix$method1[i] <- (sum(res_0$rej1[-pos,i], na.rm = T)+
                               sum(res_0$rej1.true[pos,j], na.rm = T))/B
    tab20_mix$method2[i] <- (sum(res_0$rej2[-pos,i], na.rm = T)+
                               sum(res_0$rej2.true[pos,j], na.rm = T))/B
    tab20_mix$method3[i] <- (sum(res_0$rej3[-pos,i], na.rm = T)+
                               sum(res_0$rej3.true[pos,j], na.rm = T))/B
    tab20_mix$method4[i] <- (sum(res_0$rej4[-pos,i], na.rm = T)+
                               sum(res_0$rej4.true[pos,j], na.rm = T))/B
    tab20_mix$method5[i] <- (sum(res_0$rej5[-pos,i], na.rm = T)+
                               sum(res_0$rej5.true[pos,j], na.rm = T))/B
  }
  
  
  sig_int <- rep(NA, length(alpha1))
  for (i in 1:length(alpha1)){
    sig_int[i] <- sum(res_0$gamma3 > qnorm(1-alpha1[i]), na.rm = T)/B
  }
  tab_sig_int <- data.frame(alpha1 = alpha1,
                      int = sig_int)
  return(list(tab20 = tab20, 
              tab20_true = tab20_true, 
              tab20_mix = tab20_mix, 
              tab_sig_int = tab_sig_int))
}

res_null(res0_ds1)
res_null(res0_ds2)
res_null(res0_ds4)


write.csv(res_null(res0_ds1)$tab20, "tab20_ds1.csv")
write.csv(res_null(res0_ds2)$tab20, "tab20_ds2.csv")
write.csv(res_null(res0_ds4)$tab20, "tab20_ds4.csv")
write.csv(res_null(res20_p3_d2)$tab20, "tab20.csv")





