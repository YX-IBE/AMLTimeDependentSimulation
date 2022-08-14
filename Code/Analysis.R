#
# Simulation study to explore the guarantee-time bias due to the wait for 
# stem cell transplant (SCT) regarding the treatment effect of SCT on prognostic outcomes:
# relapse (Rel) or overall survival (OS)
# comment out one of the two lines (21&22) to switch the setting: line 21:table 1 (A), line 21: table 1(B)
#
library(survival)
library(simsurv)
library(dynpred)
#
set.seed(11321)
setwd("~/AMLTimeDependentSimulation/Code")
source("Functions.R")
# 
my.simulation.rfc <- function(k = 1, n.pat = 600) {
  #
  # We analyse the data with three approaches: 
  # naive Cox, Cox with time-dependent SCT, and the super landmark model
  # The simulation study contains n.pat patients
  #
  #dat.simul <- simfun(n = n.pat, trt.lower = 0, trt.upper = 5)
  dat.simul <- simfun(n = n.pat, trt.lower = 0, trt.upper = 5, lambda.os = 0.15, lambda.rel = 0.05)
  long.format.OS <- long.format.OS.rfc(dat.simul)
  long.format.RFS <- long.format.RFS.rfc(dat.simul)
  LMdata.OS <- super.OS.rfc(dat.simul)
  LMdata.RFS <- super.RFS.rfc(dat.simul) 
  #
  # naive Cox: RFS and OS
  #
  cox.1.RFS <- coxph(Surv(t.RFS, s.RFS) ~ SCT + var1 + var2, data = dat.simul)
  
  cox.1.OS <- coxph(Surv(t.death, s.death) ~ SCT + s.relapse + var1 + var2, data = dat.simul)
  #
  # time-dependent Cox: RFS and OS
  #
  cox.2.RFS <- coxph(Surv(start, stop, s.RFS) ~ SCT + var1 + var2, data = long.format.RFS)
  
  cox.2.OS <- coxph(Surv(start, stop, s.death) ~ SCT + rel + var1 + var2, data = long.format.OS)
  #
  # super landmark model: RFS and OS
  #
  cox.3.RFS <- coxph(Surv(LM, t.RFS, s.RFS) ~ t.SCT + var1 + var2 + strata(LM.strat)
                     + cluster(id), data = LMdata.RFS, method = "breslow")
  
  cox.3.OS <- coxph(Surv(LM, t.death, s.death) ~ t.SCT + rel + var1 + var2 + strata(LM.strat)
                    + cluster(id), data = LMdata.OS, method = "breslow")
  #
  # Preparing the results: RFS
  #
  cc.1.RFS <- summary(cox.1.RFS)$coefficients[1, c(1,5)]
  cc.2.RFS <- summary(cox.2.RFS)$coefficients[1, c(1,5)]
  cc.3.RFS <- summary(cox.3.RFS)$coefficients[1, c(1,6)]
  #
  # Preparing the results: OS
  #
  cc.1.OS <- summary(cox.1.OS)$coefficients[1, c(1,5)]
  cc.2.OS <- summary(cox.2.OS)$coefficients[1, c(1,5)]
  cc.3.OS <- summary(cox.3.OS)$coefficients[1, c(1,6)]
  
  
  res <- c(cc.1.RFS, cc.2.RFS, cc.3.RFS,
           cc.1.OS, cc.2.OS, cc.3.OS, 
           n.pat)
  
  names(res)<-c("coef.RSF.naive", "p.RSF.naive",
                "coef.RSF.td", "p.RSF.td",
                "coef.RSF.LM", "p.RSF.LM",
                "coef.OS.naive", "p.OS.naive",
                "coef.OS.td", "p.OS.td",
                "coef.OS.LM", "p.OS.LM",
                "N")
  
  return(res)
}

res <- lapply(1:1000, my.simulation.rfc)

mean(unlist(lapply(res, function(x){return(x[2])})) < 0.05)
mean(unlist(lapply(res, function(x){return(x[4])})) < 0.05)
mean(unlist(lapply(res, function(x){return(x[6])})) < 0.05)
mean(unlist(lapply(res, function(x){return(x[8])})) < 0.05)
mean(unlist(lapply(res, function(x){return(x[10])})) < 0.05)
mean(unlist(lapply(res, function(x){return(x[12])})) < 0.05)

exp(mean(unlist(lapply(res, function(x){return(x[1])}))))
exp(mean(unlist(lapply(res, function(x){return(x[3])}))))
exp(mean(unlist(lapply(res, function(x){return(x[5])}))))
exp(mean(unlist(lapply(res, function(x){return(x[7])}))))
exp(mean(unlist(lapply(res, function(x){return(x[9])}))))
exp(mean(unlist(lapply(res, function(x){return(x[11])}))))

