#
# simfun creates survival data of the following setting: 
# patients could experience relapse, death without relapse, or death after relapse
# there are two fixed prognostic factors that are dichotomous: var1 and var2 
# patients can receive SCT any time within a certain period of time, if relapse/death have not occurred
# we assume no effect of SCT on Rel or OS
#
simfun <- function(n = 600, t.censor = 5,
                   trt.lower = 0, trt.upper = 5,
                   lambda.os = 0.05, gammas.os = 1.5,
                   lambda.rel = 0.15, gammas.rel = 1.5,
                   lambda.post = 0.3, gammas.post = 1.5) {
  
  df <- data.frame(id = 1:n,
                   trt = runif(n = n, trt.lower, trt.upper),
                   var1 = rbinom(n = n, size = 1, prob = 0.5),
                   var2 = rbinom(n = n, size = 1, prob = 0.4))
  
  os <- simsurv(lambdas = lambda.os , gammas = gammas.os,
                betas = c(var1 = 0.5, var2 = 0.75), x = df)
  df <- merge(df, os)
  names(df)[5:6] <- c("t.d", "s.d")
  
  rel <- simsurv(lambdas = lambda.rel, gammas = gammas.rel,
                 betas = c(var1 = 0.25, var2 = 0.5), x = df)
  df <- merge(df, rel)
  names(df)[7:8] <- c("t.r", "s.r")
  
  df$t.relapse <- ifelse(df$t.r < df$t.d, df$t.r, df$t.d)
  df$s.relapse <- ifelse(df$t.r < df$t.d, 1, 0)
  
  t.d.post <- simsurv(lambdas = lambda.post, gammas = gammas.post,
                      betas = c(var1 = 0.5, var2 = 0.75), x = df)
  df$death.new <- df$t.r + t.d.post$eventtime
  
  df$t.death <- ifelse(df$t.r < df$t.d, df$death.new, df$t.d)
  df$s.death <- rep(1, n)
  
  df$t.RFS <- ifelse(df$t.relapse < df$t.death, df$t.relapse, df$t.death)
  df$s.RFS <- ifelse(df$t.relapse < df$t.death, 1, df$s.death)
  #
  # administrative censoring at t.censor
  #
  df$s.relapse[df$t.relapse > t.censor] <- 0
  df$t.relapse[df$t.relapse > t.censor] <- t.censor
  
  df$s.death[df$t.death > t.censor] <- 0
  df$t.death[df$t.death > t.censor] <- t.censor
  
  df$s.RFS[df$t.RFS > t.censor] <- 0
  df$t.RFS[df$t.RFS > t.censor] <- t.censor
  
  df$SCT <- ifelse(df$trt < df$t.relapse, 1, 0)
  df$t.SCT <- ifelse(df$trt < df$t.relapse, df$trt, t.censor)
  
  return(df)
}
#
#
# create a long-format data set for OS
#
#
long.format.OS.rfc <- function(dat.simul) {
  
  dat.simul.OS.long <- NULL
  
  n.pat <- nrow(dat.simul)
  
  for (i in 1:n.pat) {
    
    ind.1 <- (dat.simul$SCT[i]==0 & dat.simul$s.relapse[i]==0)
    if(ind.1) 
      yy <- c(id = i, start = 0, stop = dat.simul$t.death[i],
              s.death = dat.simul$s.death[i],
              SCT = 0, rel = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i])
    
    ind.2 <- (dat.simul$SCT[i]==1 & dat.simul$s.relapse[i]==0)
    if (ind.2) {
      yy <- c(id = i, start = 0, stop = dat.simul$t.SCT[i], 
              s.death = 0, 
              SCT = 0, rel = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i])
      
      yy <- rbind(yy, 
                  c(id = i, start = dat.simul$t.SCT[i], stop = dat.simul$t.death[i], 
                    s.death = dat.simul$s.death[i],
                    SCT = 1, rel = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i]))
    }
    
    ind.3 <- (dat.simul$SCT[i]==0 & dat.simul$s.relapse[i]==1)
    if(ind.3) {
      yy <- c(id = i, start = 0, stop = dat.simul$t.relapse[i], 
              s.death = 0, 
              SCT = 0, rel = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i])
      
      yy <- rbind(yy, 
                  c(id = i, start = dat.simul$t.relapse[i], stop = dat.simul$t.death[i],
                    s.death = dat.simul$s.death[i],
                    SCT = 0, rel = 1, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i]))
    }
    
    ind.4 <- (dat.simul$SCT[i]==1 & dat.simul$s.relapse[i]==1)
    if(ind.4) {
      yy <- c(id = i, start = 0, stop = dat.simul$t.SCT[i], 
              s.death = 0, 
              SCT = 0, rel = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i])
      
      yy <- rbind(yy, 
                  c(id = i, start = dat.simul$t.SCT[i], stop = dat.simul$t.relapse[i], 
                    s.death = 0,
                    SCT = 1, rel = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i]))
      
      yy <- rbind(yy, 
                  c(id = i, start = dat.simul$t.relapse[i], stop = dat.simul$t.death[i], 
                    s.death = dat.simul$s.death[i],
                    SCT = 1, rel = 1, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i]))
    }
    
    dat.simul.OS.long <- rbind(dat.simul.OS.long, yy)
  }
  
  dat.simul.OS.long <- as.data.frame(dat.simul.OS.long)
  
  return(dat.simul.OS.long)
}
#
#
# create a long-format data set for RFS
#
#
long.format.RFS.rfc <- function(dat.simul) {
  
  dat.simul.RFS.long <- NULL
  
  n.pat <- nrow(dat.simul)
  
  for (i in 1:n.pat) {
    
    ind.1 <- (dat.simul$SCT[i]==0)
    
    if (ind.1) yy <- c(id = i, start = 0, stop = dat.simul$t.RFS[i], 
                       s.RFS = dat.simul$s.RFS[i],
                       SCT = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i])
    
    if (!ind.1) {
      yy <- c(id = i, start = 0, stop = dat.simul$t.SCT[i], s.RFS = 0, 
              SCT = 0, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i])
      yy <- rbind(yy, 
                  c(id = i, start = dat.simul$t.SCT[i], stop = dat.simul$t.RFS[i], s.RFS = dat.simul$s.RFS[i],
                    SCT = 1, var1 = dat.simul$var1[i], var2 = dat.simul$var2[i]))
    }
    
    dat.simul.RFS.long <- rbind(dat.simul.RFS.long, yy)
  }
  
  dat.simul.RFS.long <- as.data.frame(dat.simul.RFS.long)
  
  return(dat.simul.RFS.long)
}

#
#
# create a super landmark data set for OS
#
#
super.OS.rfc <- function(dat.simul = dat.simul,
                         LM.horizon = 2, LM.lower = 0,
                         LM.upper = 3, by.delta = 0.1) {
  
  LMs <- seq(LM.lower, LM.upper, by = by.delta)
  
  # landmark = 0
  LMdata <- cutLM(data = dat.simul, outcome = list(time = "t.death", status = "s.death"),
                  LM = 0, horizon = LM.horizon, 
                  covs = list(fixed=c("var1", "var2", "id"), varying = "t.SCT"))
  
  # sliding landmark
  for (i in 2:length(LMs))
    LMdata <- rbind(LMdata, cutLM(data = dat.simul, outcome = list(time = "t.death", status = "s.death"),
                                  LM = LMs[i], horizon = LMs[i] + LM.horizon,
                                  covs = list(fixed=c("var1", "var2", "id"), varying = "t.SCT")))
  LMdata$LM.strat<-LMdata$LM
  
  rel <- rep(0, nrow(LMdata))
  
  for (ii in dat.simul$id) {
    ind <- (LMdata$id == ii)
    xx <- LMdata[ind, ]
    t.relapse <- dat.simul$t.relapse[dat.simul$id == ii]
    s.relapse <- dat.simul$s.relapse[dat.simul$id == ii]
    res <- ifelse(xx$LM + LM.horizon >= t.relapse & s.relapse == 1, 1, 0)
    rel[ind] <- res
  }
  
  LMdata$rel <- rel
  LMdata.OS <- LMdata
  return(LMdata.OS)
}

#
#
# create a super landmark data set for RFS
#
#
super.RFS.rfc <- function(dat.simul = dat.simul,
                          LM.horizon = 2, LM.lower = 0,
                          LM.upper = 3, by.delta = 0.1) {
  
  LMs <- seq(LM.lower, LM.upper, by = by.delta)
  
  LMdata <- cutLM(data = dat.simul, outcome = list(time="t.RFS", status="s.RFS"),
                  LM = 0, horizon = LM.horizon, 
                  covs = list(fixed = c("var1", "var2", "id"), varying = "t.SCT"))
  
  for (i in 2:length(LMs))
    LMdata <- rbind(LMdata, cutLM(data = dat.simul, outcome = list(time = "t.RFS", status="s.RFS"),
                                  LM = LMs[i], horizon = LMs[i] + LM.horizon,
                                  covs = list(fixed = c("var1", "var2", "id"), varying = "t.SCT")))
  
  LMdata$LM.strat <- LMdata$LM
  
  LMdata.RFS<-LMdata
  return(LMdata.RFS)
}
