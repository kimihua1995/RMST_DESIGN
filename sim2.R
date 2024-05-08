library(survRM2)
library(survival)
library(tidyverse)
library(eventTrack)
library(nleqslv)
library(geex)


source("my_survRM2.R")
source("paper_simulation_source2.R")
source("source_design2.R")
source("find_cutpoint2.R")
source("source_stat2.R")

args=(commandArgs(TRUE))

if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}



t.change <- 0.25
t.enroll1 <- 1
t.enrich <- 1
t.enroll2 <- 2
t.FA <- 4
tau <- 2



true.cutoff <- true.cut.piece0(tau,lambda0,lambda1,beta0,beta1,t.change)
true.cutoff <- as.numeric(true.cutoff[1])


seed=333
B=10000
lambda0=c(0.9,0.9)
lambda1=c(0.9,0.45)
beta0=0
beta1=0.9
n1=505
n2=505
alpha1=0.025
alpha2=0.025
resa_ds1 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                 t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                 tau, pred.method=1, design=1, alpha1, alpha2, true.cutoff)
resa_ds2 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, pred.method=1, design=2, alpha1, alpha2, true.cutoff)
resa_ds3 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, pred.method=2, design=2, alpha1, alpha2, true.cutoff)
resa_ds4 <- my.design(seed, B, n1, n2, lambda0, lambda1, beta0, beta1,
                      t.change, t.enroll1, t.enroll2, t.enrich, t.FA, 
                      tau, pred.method=3, design=2, alpha1, alpha2, true.cutoff)





