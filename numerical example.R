library(survRM2)
library(survival)
library(survminer)
library(tidyverse)
library(OpenMx)
#library(twang)
#ibrary(PSW)
library(Rsolnp)
library(eventTrack)
library(nleqslv)
#library(pseudo)
#library(drgee)
library(geex)
library(ggpubr)
library(ggplot2)


source("source_numerical example.R")
source("source_stat_numerical example.R")



##########################################################
lambda0 <- c(log(2)*2.5,log(2)*2.5)
lambda1 <- c(log(2)*6,log(2)*2.0)
beta0 <- 0
beta1 <- -0.8
t.change <- 1/6
n <- 250
tau <- 1.5
seed <- 212




KM_plot_all <- KM_plot(seed, n, 
                       lambda0, lambda1, beta0, beta1, 
                       t.change, tau, 0.01, "PD-L1 >= 1%")

KM_plot_pos <- KM_plot(seed, n, 
                       lambda0, lambda1, beta0, beta1, 
                       t.change, tau, 0.5, "PD-L1 >= 50%")

pdf("KM_curves.pdf", width = 8, height = 6, onefile = T)
KM_plot_all$curve
KM_plot_pos$curve
dev.off()

true.cutoff <- true.cut.piece(tau,lambda0,lambda1,beta0,beta1,t.change)
(true.cutoff <- as.numeric(true.cutoff[1]))
(RMSTD_pos <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,true.cutoff,1))
(RMSTD_all <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,0.01,1))




true_RMSTreg <- true.RMSTreg(seed, n=5000, B=10000,
                             lambda0, lambda1, beta0, beta1, 
                             t.change, tau)
colMeans(true_RMSTreg$betas)
apply(true_RMSTreg$betas, 2, sd)
colMeans(true_RMSTreg$vcovs)
mean(true_RMSTreg$betas[,4]/true_RMSTreg$vcovs[,4] > qnorm(0.975))

#########################################################



t.enroll1 <- 0.5 # stage I
t.enrich <- 0.5
t.enroll2 <- 1 # stage II
t.FA <- 2.5 # follow-up
#alpha <- 0.025
#alpha0 <- 0.025
#power <- 0.9
#n_pos <- sample_size(seed, B=1000, lambda0, lambda1, beta0, beta1,
#            t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
#            tau, alpha, power)
#ceiling(ceiling(n_pos)/(1.99-true.cutoff)*1.99)


# calculate asymptotic variance

  true.cutoff <- true.cut.piece(tau,lambda0,lambda1,beta0,beta1,t.change)
  true.cutoff <- as.numeric(true.cutoff[1])
  RMSTD_pos <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,true.cutoff,1)
  RMSTD_all <- rmst.diff.true(tau,lambda0,lambda1,beta0,beta1,t.change,0.01,1)
  
  stat_pos <- est.sigma(seed, B=10000, lambda0, lambda1, beta0, beta1,
                   t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                   tau, true.cutoff)
  #I_pos_N <- sqrt(10000)*(stat_pos$Naive$d_pos - RMSTD_pos)
  #I_pos_GF <- sqrt(10000)*(stat_pos$GF$d_pos - RMSTD_pos)
  #I_pos_AG <- sqrt(10000)*(stat_pos$AG$d_pos - RMSTD_pos)
  #sigma_pos_N <- sqrt(sum(I_pos_N^2)/2000)
  #sigma_pos_GF <- sqrt(sum(I_pos_GF^2)/2000)
  #sigma_pos_AG <- sqrt(sum(I_pos_AG^2)/2000)
  sigma_pos_N <- mean(stat_pos$Naive$se_pos)*sqrt(10000)
  sigma_pos_GF <- mean(stat_pos$GF$se_pos)*sqrt(10000)
  sigma_pos_HJ <- mean(stat_pos$HJ$se_pos)*sqrt(10000)
  sigma_pos_AG <- mean(stat_pos$AG$se_pos)*sqrt(10000)
  
  
  stat_all <- est.sigma(seed, B=10000, lambda0, lambda1, beta0, beta1,
                       t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                       tau, 0.01)
  #I_all_N <- sqrt(10000)*(stat_all$Naive$d_pos - RMSTD_all)
  #I_all_GF <- sqrt(10000)*(stat_all$GF$d_pos - RMSTD_all)
  #I_all_AG <- sqrt(10000)*(stat_all$AG$d_pos - RMSTD_all)
  #sigma_all_N <- sqrt(sum(I_all_N^2)/2000)
  #sigma_all_GF <- sqrt(sum(I_all_GF^2)/2000)
  #sigma_all_AG <- sqrt(sum(I_all_AG^2)/2000)
  sigma_all_N <- mean(stat_all$Naive$se_pos)*sqrt(10000)
  sigma_all_GF <- mean(stat_all$GF$se_pos)*sqrt(10000)
  sigma_all_HJ <- mean(stat_all$HJ$se_pos)*sqrt(10000)
  sigma_all_AG <- mean(stat_all$AG$se_pos)*sqrt(10000)
  
  true_RMSTreg <- true.RMSTreg(seed, n=5000, B=10000,
                               lambda0, lambda1, beta0, beta1, 
                               t.change, tau)
  beta3 <- mean(true_RMSTreg$betas[,4])
  sigma_beta3 <- mean(true_RMSTreg$vcovs[,4])*sqrt(10000)
  
  
  
  
  
  
# critical values
B <- 10000
res.null <- my.design.null(seed, B, lambda0, beta0, t.change,
                           t.enroll1, t.enroll2, t.enrich, t.FA, 
                           tau, alpha0, true.cutoff)
  
alpha_list <- seq(0.015, 0.025, by = 0.001)
q_list <- qnorm(1 - alpha_list)
stat_null <- res.null$stat.FA$z_pos
rej <- rep(NA, length(alpha_list))
for (i in 1:length(q_list)){
  rej[i] <- sum(stat_null > q_list[i], na.rm = T)/B
}
  
  
  save(stat_pos, stat_all, true_RMSTreg, res.null, file = "res10000_212.RData")

  
  
  
  n_list <- seq(700,1000,by=5)
  q <- qnorm(1-0.023)
  q0 <- qnorm(1-0.025)
  n_pos_list1 <- ceiling(n_list*(1-true.cutoff)/0.99)
  n_pos_list2 <- ceiling(n_list/2*(1-true.cutoff)/0.99 + n_list/2)
  eta_list <- 1 - pnorm(q0 - beta3/sigma_beta3*sqrt(n_list))
  power1_N <- (1 - pnorm(q - RMSTD_pos/sigma_pos_N*sqrt(n_pos_list1)))*eta_list +
    (1 - pnorm(q - RMSTD_all/sigma_all_N*sqrt(n_list)))*(1-eta_list)
  power2_N <- (1 - pnorm(q - RMSTD_pos/sigma_pos_N*sqrt(n_pos_list2)))*eta_list +
    (1 - pnorm(q - RMSTD_all/sigma_all_N*sqrt(n_list)))*(1-eta_list)
  power2_GF <- (1 - pnorm(q - RMSTD_pos/sigma_pos_GF*sqrt(n_pos_list2)))*eta_list +
    (1 - pnorm(q - RMSTD_all/sigma_all_GF*sqrt(n_list)))*(1-eta_list)
  power2_HJ <- (1 - pnorm(q - RMSTD_pos/sigma_pos_HJ*sqrt(n_pos_list2)))*eta_list +
    (1 - pnorm(q - RMSTD_all/sigma_all_HJ*sqrt(n_list)))*(1-eta_list)
  power2_AG <- (1 - pnorm(q - RMSTD_pos/sigma_pos_AG*sqrt(n_pos_list2)))*eta_list +
    (1 - pnorm(q - RMSTD_all/sigma_all_AG*sqrt(n_list)))*(1-eta_list)
  
  
  plot_dat <- data.frame(n = rep(n_list,5),
                         power = c(power1_N, power2_N, power2_GF, power2_HJ, power2_AG),
                         Design = rep(c("All-Comer", "Enrich-KM", "Enrich-GF",
                                        "Enrich-HJ","Enrich-AG"), each = length(n_list)))
  plot_dat$design <- factor(plot_dat$Design,
                            levels = c("All-Comer","Enrich-KM","Enrich-GF","Enrich-HJ","Enrich-AG"),
                            labels = c("All-Comer","Enrich-KM","Enrich-GF","Enrich-HJ","Enrich-AG"))
  
  
  library(RColorBrewer)
  f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
  (cols <- f("Dark2"))
  (cols <- f("Set1"))
  
  p <- ggplot(data = plot_dat, aes(x = n, y = power, color = Design)) +
    geom_line(size=1) +
    #scale_color_brewer(palette = "PuOr") +
    scale_color_manual(values=c("#666666","#1B9E77","#984EA3","#FF7F00","#377EB8"),
                       breaks = c("All-Comer","Enrich-KM","Enrich-GF","Enrich-HJ","Enrich-AG"))+
    geom_hline(yintercept=0.9, linetype="dashed", color = "red", size = 1) +
    scale_x_continuous(breaks = seq(700,1000,50)) +
    labs(x = "Total Sample Size n", y = "Global Power") +
    theme_bw() +
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
  
  pdf("sample size.pdf", width = 8, height = 6, onefile = T)
  p
  dev.off()





