library(tidyverse) # data.frame manipulation
library(Hmisc) # cstat computation
library(magrittr) # cstat computation
library(doParallel) # parallel computing
library(gtools) # calibration computation
library(survival) # TTE outcome
library(coxme) # cox regression with mixed effects
source("coxvc.R") # rank-1 with TTE outcome
library(splines) # used in rank-1
library(MASS) # ginv in rank-1

setwd("C:/Users/Florie BRION-BOUVIER/Documents/These/ITE with IPD-MA/Fichiers R/Revised code")

mse = function(y_true,y_pred){mean((y_true - y_pred)^2)}

expit = function(x) exp(x)/(1+exp(x))

cstat4ben = function(outcome, treatment, score, conf_lvl = 0.95, nbr_rep = 100){
  # return the C statistic for benefit and it's confidence interval of level 'conf_level'
  # check if the required packages are loaded
  if (!all(c("Hmisc", "magrittr") %in% .packages()))
    stop("packages \"Hmisc\" and \"magrittr\" must be loaded")
  
  # check if the treatment is binary in c(0,1)
  if (!all(unique(treatment) %in% 0:1))
    stop("'treatment' must be in c(0,1)")
  
  # score has to be continuous (e.g. ITE)
  if (length(unique(score)) < length(score)/5)
    warning("'score' has to be continuous")
  
  y = as.numeric(outcome)
  n0 = sum(treatment == 0); n1 = sum(treatment == 1)
  n.pair = min(n0, n1)
  replicate(
    n = nbr_rep,
    expr = {
      # sample n.pair observations in the largest group
      ind.0 = which(treatment==0)[sample(n0, n.pair)]
      ind.1 = which(treatment==1)[sample(n1, n.pair)]
      
      # get the index of the observations ordered by their score
      ind.0 = ind.0[order(score[ind.0])]
      ind.1 = ind.1[order(score[ind.1])]
      
      # compute the predicted and observed benefit in matched pairs
      pred.ben.avg = (score[ind.1] + score[ind.0]) / 2
      obs.ben = y[ind.1] - y[ind.0]
      
      # Benefit c-statistic
      cindex = rcorr.cens(pred.ben.avg, obs.ben)
      
      return(c(cindex["C Index"][[1]], cindex["S.D."][[1]]/2, pred.ben.avg, obs.ben))
    }) %>%
    { .[ ,.[1, ] %in% apply(., 1, quantile, probs = 0.5, type = 1)[1]] } %>%
    as.matrix(.) %>% .[ ,1]
}

#### scenario 13: 3 covariates & tte outcome & total sample size = 2800 ####

#### s-learner ####

#### naive model ####
m <-1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.na[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.na),
                           coef_mat.1 = mean(coef_mat[1, ]),
                           coef_mat.2 = mean(coef_mat[2, ]),
                           coef_mat.3 = mean(coef_mat[3, ]),
                           coef_mat.4 = mean(coef_mat[4, ]),
                           coef_mat.5 = mean(coef_mat[5, ]),
                           coef_mat.6 = mean(coef_mat[6, ]),
                           coef_mat.7 = mean(coef_mat[7, ]),
                           se_mat.1 = mean(se_mat[1, ]),
                           se_mat.2 = mean(se_mat[2, ]),
                           se_mat.3 = mean(se_mat[3, ]),
                           se_mat.4 = mean(se_mat[4, ]),
                           se_mat.5 = mean(se_mat[5, ]),
                           se_mat.6 = mean(se_mat[6, ]),
                           se_mat.7 = mean(se_mat[7, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.sl3.tte.n2800 <- t(res.na.sl3[1:12, ])
se.na.sl3.tte.n2800 <- t(res.na.sl3[13:19, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 7, ncol = k)
  se_mat <- matrix(NA, nrow = 7, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+(1+treat|trial),
                 data = train)
    coef <- mod$coefficients
    coef_mat[,i] <- coef
    se_mat[,i] <- sqrt(diag(vcov(mod)))
    
    #recalibration of event rates adapted to train
    lp <- mod$linear.predictor
    temp <- coxph(Surv(train$time,train$status)~offset(lp))
    bh <- basehaz(temp, centered=F)
    
    #ite estimation
    pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
    pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    coef_mat.1 = mean(coef_mat[1,]),
    coef_mat.2 = mean(coef_mat[2,]),
    coef_mat.3 = mean(coef_mat[3,]),
    coef_mat.4 = mean(coef_mat[4,]),
    coef_mat.5 = mean(coef_mat[5,]),
    coef_mat.6 = mean(coef_mat[6,]),
    coef_mat.7 = mean(coef_mat[7,]),
    se_mat.1 = mean(se_mat[1,]),
    se_mat.2 = mean(se_mat[2,]),
    se_mat.3 = mean(se_mat[3,]),
    se_mat.4 = mean(se_mat[4,]),
    se_mat.5 = mean(se_mat[5,]),
    se_mat.6 = mean(se_mat[6,]),
    se_mat.7 = mean(se_mat[7,])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.sl3.tte.n2800 <- t(res.re.sl3[1:12, ])
se.re.sl3.tte.n2800 <- t(res.re.sl3[13:19, ])

#### stratified intercept ####
m <-1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 7, ncol = k)
  se_mat <- matrix(NA, nrow = 7, ncol = k)
  
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+strata(trial)+(0+treat|trial),
                 data = train)
    coef <- mod$coefficients
    coef_mat[,i] <- coef
    se_mat[,i] <- sqrt(diag(vcov(mod)))
    
    #recalibration of event rates adapted to train
    lp <- mod$linear.predictor
    temp <- coxph(Surv(train$time,train$status)~offset(lp))
    bh <- basehaz(temp, centered=F)
    
    #ite estimation
    pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
    pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[6, ]),
    coef_mat.7 = mean(coef_mat[7, ]),
    se_mat.1 = mean(se_mat[1, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[6, ]),
    se_mat.7 = mean(se_mat[7, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.sl3.tte.n2800 <- t(res.si.sl3[1:12, ])
se.si.sl3.tte.n2800 <- t(res.si.sl3[13:19, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  source("coxvc.R")
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    Ft <- cbind(rep(1,nrow(train)),bs(train$time))
    mod <- coxvc(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+trial, Ft, rank = 1, data = train)
    coef <- mod$b
    coef_mat[,i] <- coef
    
    #recalibration of event rates adapted to train
    lp <- unlist(as.matrix(train[2:5]) %*% coef[1:4] + train[5] * (as.matrix(train[2:4]) %*% coef[6:8]))
    temp <- coxph(Surv(train$time,train$status)~offset(lp))
    bh <- basehaz(temp, centered=F)
    
    #ite estimation
    pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[6:8])))
    pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[6:8])))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben <- cstat[1]
    c.ben.se <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[6, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl3.tte.n2800 <- t(res.r1.sl3[1:12, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat <- matrix(NA, nrow = 27, ncol = k)
  se_mat <- matrix(NA, nrow = 27, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial),
                 data = train)
    coef <- mod$coefficients
    coef_mat[,i] <- coef
    se_mat[,i] <- sqrt(diag(vcov(mod)))
    
    #recalibration of event rates adapted to train
    lp <- mod$linear.predictors
    temp <- coxph(Surv(train$time,train$status)~offset(lp))
    bh <- basehaz(temp, centered=F)
    
    #ite estimation
    pred0 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test0$trial+
             mean(coef[c(6,13:17)])*test0$age+mean(coef[c(7,18:22)])*test0$sex+
             mean(coef[c(8,23:27)])*test0$sbp+coef[9]*test0$treat+
             (coef[10]*test0$age+coef[11]*test0$sex+coef[12]*test0$sbp)*test0$treat))
    pred1 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test1$trial+
             mean(coef[c(6,13:17)])*test1$age+mean(coef[c(7,18:22)])*test1$sex+
             mean(coef[c(8,23:27)])*test1$sbp+coef[9]*test1$treat+
             (coef[10]*test1$age+coef[11]*test1$sex+coef[12]*test1$sbp)*test1$treat))
    ite <- unlist(pred1 - pred0)

    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
    coef_mat.1 = mean(c(coef_mat[6,]+coef_mat[13:17,])),
    coef_mat.2 = mean(c(coef_mat[7,]+coef_mat[18:22,])),
    coef_mat.3 = mean(c(coef_mat[8,]+coef_mat[23:27,])),
    coef_mat.4 = mean(coef_mat[9,]),
    coef_mat.5 = mean(coef_mat[10,]),
    coef_mat.6 = mean(coef_mat[11,]),
    coef_mat.7 = mean(coef_mat[12,]),
    se_mat.1 = mean(c(se_mat[6,]+se_mat[13:17,])),
    se_mat.2 = mean(c(se_mat[7,]+se_mat[18:22,])),
    se_mat.3 = mean(c(se_mat[8,]+se_mat[23:27,])),
    se_mat.4 = mean(se_mat[9,]),
    se_mat.5 = mean(se_mat[10,]),
    se_mat.6 = mean(se_mat[11,]),
    se_mat.7 = mean(se_mat[12,])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.fs.sl3.tte.n2800 <- t(res.fs.sl3[1:12, ])
se.fs.sl3.tte.n2800 <- t(res.fs.sl3[13:19, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train0)
                           mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.na[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.na),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.tl3.tte.n2800 <- t(res.na.tl3[1:12, ])
se.na.tl3.tte.n2800 <- t(res.na.tl3[13:18, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
  se_mat0 <- matrix(NA, nrow = 3, ncol = k)
  se_mat1 <- matrix(NA, nrow = 3, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- coxme(Surv(time,status)~age+factor(sex)+sbp+(1|trial),
                  data = train0)
    mod1 <- coxme(Surv(time,status)~age+factor(sex)+sbp+(1|trial),
                  data = train1)
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    se_mat0[,i] <- sqrt(diag(vcov(mod0)))
    se_mat1[,i] <- sqrt(diag(vcov(mod1)))
    
    #recalibration of event rates adapted to train
    lp0 <- mod0$linear.predictor
    lp1 <- mod1$linear.predictor
    temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
    temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
    bh0 <- basehaz(temp0, centered=F)
    bh1 <- basehaz(temp1, centered=F)
    int0[i] <- bh0[nrow(bh0),]$hazard
    int1[i] <- bh1[nrow(bh1),]$hazard
    
    #ite estimation
    pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
    pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(int1)-mean(int0),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.tl3.tte.n2800 <- t(res.re.tl3[1:12, ])
se.re.tl3.tte.n2800 <- t(res.re.tl3[13:18, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
  se_mat0 <- matrix(NA, nrow = 3, ncol = k)
  se_mat1 <- matrix(NA, nrow = 3, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+strata(trial),data = train0)
    mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+strata(trial),data = train1)
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    se_mat0[,i] <- sqrt(diag(vcov(mod0)))
    se_mat1[,i] <- sqrt(diag(vcov(mod1)))
    
    #recalibration of event rates adapted to train
    lp0 <- mod0$linear.predictors
    lp1 <- mod1$linear.predictors
    temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
    temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
    bh0 <- basehaz(temp0, centered=F)
    bh1 <- basehaz(temp1, centered=F)
    int0[i] <- bh0[nrow(bh0),]$hazard
    int1[i] <- bh1[nrow(bh1),]$hazard
    
    #ite estimation
    pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
    pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(int1)-mean(int0),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.tl3.tte.n2800 <- t(res.si.tl3[1:12, ])
se.si.tl3.tte.n2800 <- t(res.si.tl3[13:18, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  source("coxvc.R")
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    Ft0<-cbind(rep(1,nrow(train0)),bs(train0$time))
    mod0 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+trial, Ft0, rank = 1, data = train0)
    Ft1<-cbind(rep(1,nrow(train0)),bs(train1$time))
    mod1 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+trial, Ft1, rank = 1, data = train1)
    
    coef0 <- mod0$b
    coef1 <- mod1$b
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    
    #recalibration of event rates adapted to train
    lp0 <- unlist(as.matrix(train0[2:4]) %*% coef[1:3])
    lp1 <- unlist(as.matrix(train1[2:4]) %*% coef[1:3])
    temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
    temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
    bh0 <- basehaz(temp0, centered=F)
    bh1 <- basehaz(temp1, centered=F)
    int0[i] <- bh0[nrow(bh0),]$hazard
    int1[i] <- bh1[nrow(bh1),]$hazard
    
    #ite estimation
    pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
    pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(int1)-mean(int0),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ])),
    seed = j,
    success = TRUE)
  
}
stopCluster(cl)
res.r1.tl3.tte.n2800 <- t(res.r1.tl3[1:12, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df3 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82)
  sd_age <- c(4,2,1,3,4,6,2)
  df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df3$age <- df3$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
  df3$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197)
  sd_sbp <- c(13,16.5,9.4,12,9,10,21)
  df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df3$sbp <- df3$sbp-185
  
  b0 <- 50
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0.5, 10)
  c[c>5] <- 5 # censoring
  df3$time <- pmin(t,c)  # observed time is min of censored and true
  df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 23, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 23, ncol = k)
  se_mat0 <- matrix(NA, nrow = 23, ncol = k)
  se_mat1 <- matrix(NA, nrow = 23, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(trial), data = train0)
    mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(trial), data = train1)
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    se_mat0[,i] <- sqrt(diag(vcov(mod0)))
    se_mat1[,i] <- sqrt(diag(vcov(mod1)))
    
    #recalibration of event rates adapted to train
    lp0 <- mod0$linear.predictors
    lp1 <- mod1$linear.predictors
    temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
    temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
    bh0 <- basehaz(temp0, centered=F)
    bh1 <- basehaz(temp1, centered=F)
    int0[i] <- bh0[nrow(bh0),]$hazard
    int1[i] <- bh1[nrow(bh1),]$hazard
    
    #ite estimation
    pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(mean(coef0[4:8])*test$trial+
             mean(coef0[c(1,9:13)])*test$age+mean(coef0[c(2,14:18)])*test$sex+
             mean(coef0[c(3,19:23)])*test$sbp))
    pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(mean(coef1[4:8])*test$trial+
             mean(coef1[c(1,9:13)])*test$age+mean(coef1[c(2,14:18)])*test$sex+
             mean(coef1[c(3,19:23)])*test$sbp))
    ite <- unlist(pred1 - pred0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(1-ite, q=5)
    pben <-  tapply(1-ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(oben,pben)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[9:13, ]),
    coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[14:18, ]),
    coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[19:23, ]),
    coef_mat.4 = mean(int1)-mean(int0),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[9:13, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[9:13, ])),
    coef_mat.6 = (mean(coef_mat1[2, ]) + mean(coef_mat1[14:18, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[14:18, ])),
    coef_mat.7 = (mean(coef_mat1[3, ]) + mean(coef_mat1[19:23, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[19:23, ])),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[9:13, ]),
    se_mat.2 = mean(se_mat0[2, ]) + mean(se_mat0[14:18, ]),
    se_mat.3 = mean(se_mat0[3, ]) + mean(se_mat0[19:23, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[9:13, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[9:13, ])),
    se_mat.6 = (mean(se_mat1[2, ]) + mean(se_mat1[14:18, ])) - (mean(se_mat0[2, ]) + mean(se_mat0[14:18, ])),
    se_mat.7 = (mean(se_mat1[3, ]) + mean(se_mat1[19:23, ])) - (mean(se_mat0[3, ]) + mean(se_mat0[19:23, ]))),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.fs.tl3.tte.n2800 <- t(res.fs.tl3[1:12, ])
se.fs.tl3.tte.n2800 <- t(res.fs.tl3[13:18, ])

save(res.na.sl3.tte.n2800,se.na.sl3.tte.n2800,
     res.re.sl3.tte.n2800,se.re.sl3.tte.n2800,
     res.si.sl3.tte.n2800,se.si.sl3.tte.n2800,
     res.fs.sl3.tte.n2800,se.fs.sl3.tte.n2800,
     res.na.tl3.tte.n2800,se.na.tl3.tte.n2800,
     res.re.tl3.tte.n2800,se.re.tl3.tte.n2800,
     res.si.tl3.tte.n2800,se.si.tl3.tte.n2800,
     res.fs.tl3.tte.n2800,se.fs.tl3.tte.n2800,
     file = "res_scenario13.Rdata") #res.r1.sl3.tte.n2800, res.r1.tl3.tte.n2800,

#### scenario 14: 3 covariates & tte outcome & total sample size = 1400 ####

#### s-learner ####

#### naive model ####
m <-1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.na[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.na),
                           coef_mat.1 = mean(coef_mat[1, ]),
                           coef_mat.2 = mean(coef_mat[2, ]),
                           coef_mat.3 = mean(coef_mat[3, ]),
                           coef_mat.4 = mean(coef_mat[4, ]),
                           coef_mat.5 = mean(coef_mat[5, ]),
                           coef_mat.6 = mean(coef_mat[6, ]),
                           coef_mat.7 = mean(coef_mat[7, ]),
                           se_mat.1 = mean(se_mat[1, ]),
                           se_mat.2 = mean(se_mat[2, ]),
                           se_mat.3 = mean(se_mat[3, ]),
                           se_mat.4 = mean(se_mat[4, ]),
                           se_mat.5 = mean(se_mat[5, ]),
                           se_mat.6 = mean(se_mat[6, ]),
                           se_mat.7 = mean(se_mat[7, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.sl3.tte.n1400 <- t(res.na.sl3[1:12, ])
se.na.sl3.tte.n1400 <- t(res.na.sl3[13:19, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying the model to train
                           mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+(1+treat|trial),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0) # ite for tte outcome is E[Y(1)-Y(0)]
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.re[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.re),
                           coef_mat.1 = mean(coef_mat[1,]),
                           coef_mat.2 = mean(coef_mat[2,]),
                           coef_mat.3 = mean(coef_mat[3,]),
                           coef_mat.4 = mean(coef_mat[4,]),
                           coef_mat.5 = mean(coef_mat[5,]),
                           coef_mat.6 = mean(coef_mat[6,]),
                           coef_mat.7 = mean(coef_mat[7,]),
                           se_mat.1 = mean(se_mat[1,]),
                           se_mat.2 = mean(se_mat[2,]),
                           se_mat.3 = mean(se_mat[3,]),
                           se_mat.4 = mean(se_mat[4,]),
                           se_mat.5 = mean(se_mat[5,]),
                           se_mat.6 = mean(se_mat[6,]),
                           se_mat.7 = mean(se_mat[7,])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.sl3.tte.n1400 <- t(res.re.sl3[1:12, ])
se.re.sl3.tte.n1400 <- t(res.re.sl3[13:19, ])

#### stratified intercept ####
m <-1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.si <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+strata(trial)+(0+treat|trial),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.si[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.si),
                           coef_mat.1 = mean(coef_mat[1, ]),
                           coef_mat.2 = mean(coef_mat[2, ]),
                           coef_mat.3 = mean(coef_mat[3, ]),
                           coef_mat.4 = mean(coef_mat[4, ]),
                           coef_mat.5 = mean(coef_mat[5, ]),
                           coef_mat.6 = mean(coef_mat[6, ]),
                           coef_mat.7 = mean(coef_mat[7, ]),
                           se_mat.1 = mean(se_mat[1, ]),
                           se_mat.2 = mean(se_mat[2, ]),
                           se_mat.3 = mean(se_mat[3, ]),
                           se_mat.4 = mean(se_mat[4, ]),
                           se_mat.5 = mean(se_mat[5, ]),
                           se_mat.6 = mean(se_mat[6, ]),
                           se_mat.7 = mean(se_mat[7, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.sl3.tte.n1400 <- t(res.si.sl3[1:12, ])
se.si.sl3.tte.n1400 <- t(res.si.sl3[13:19, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       source("coxvc.R")
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.r1 <- c()
                       coef_mat <- matrix(NA, nrow = 8, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           Ft <- cbind(rep(1,nrow(train)),bs(train$time))
                           mod <- coxvc(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+trial, Ft, rank = 1, data = train)
                           coef <- mod$b # maybe t(mod$b) ?
                           coef_mat[,i] <- coef
                           
                           #recalibration of event rates adapted to train
                           lp <- unlist(as.matrix(train[2:5]) %*% coef[1:4] + train[5] * (as.matrix(train[2:4]) %*% coef[6:8]))
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[6:8])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[6:8])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben <- cstat[1]
                           c.ben.se <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.r1[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.r1),
                         coef_mat.1 = mean(coef_mat[1, ]),
                         coef_mat.2 = mean(coef_mat[2, ]),
                         coef_mat.3 = mean(coef_mat[3, ]),
                         coef_mat.4 = mean(coef_mat[4, ]),
                         coef_mat.5 = mean(coef_mat[6, ]),
                         coef_mat.6 = mean(coef_mat[7, ]),
                         coef_mat.7 = mean(coef_mat[8, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.r1.sl3.tte.n1400 <- t(res.r1.sl3[1:12, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.fs <- c()
                       coef_mat <- matrix(NA, nrow = 27, ncol = k)
                       se_mat <- matrix(NA, nrow = 27, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictors
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test0$trial+
                                                                    mean(coef[c(6,13:17)])*test0$age+mean(coef[c(7,18:22)])*test0$sex+
                                                                    mean(coef[c(8,23:27)])*test0$sbp+coef[9]*test0$treat+
                                                                    (coef[10]*test0$age+coef[11]*test0$sex+coef[12]*test0$sbp)*test0$treat))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test1$trial+
                                                                    mean(coef[c(6,13:17)])*test1$age+mean(coef[c(7,18:22)])*test1$sex+
                                                                    mean(coef[c(8,23:27)])*test1$sbp+coef[9]*test1$treat+
                                                                    (coef[10]*test1$age+coef[11]*test1$sex+coef[12]*test1$sbp)*test1$treat))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.fs[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.fs),
                           coef_mat.1 = mean(c(coef_mat[6,]+coef_mat[13:17,])),
                           coef_mat.2 = mean(c(coef_mat[7,]+coef_mat[18:22,])),
                           coef_mat.3 = mean(c(coef_mat[8,]+coef_mat[23:27,])),
                           coef_mat.4 = mean(coef_mat[9,]),
                           coef_mat.5 = mean(coef_mat[10,]),
                           coef_mat.6 = mean(coef_mat[11,]),
                           coef_mat.7 = mean(coef_mat[12,]),
                           se_mat.1 = mean(c(se_mat[6,]+se_mat[13:17,])),
                           se_mat.2 = mean(c(se_mat[7,]+se_mat[18:22,])),
                           se_mat.3 = mean(c(se_mat[8,]+se_mat[23:27,])),
                           se_mat.4 = mean(se_mat[9,]),
                           se_mat.5 = mean(se_mat[10,]),
                           se_mat.6 = mean(se_mat[11,]),
                           se_mat.7 = mean(se_mat[12,])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.fs.sl3.tte.n1400 <- t(res.fs.sl3[1:12, ])
se.fs.sl3.tte.n1400 <- t(res.fs.sl3[13:19, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train0)
                           mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.na[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.na),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.tl3.tte.n1400 <- t(res.na.tl3[1:12, ])
se.na.tl3.tte.n1400 <- t(res.na.tl3[13:18, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxme(Surv(time,status)~age+factor(sex)+sbp+(1|trial),
                                         data = train0)
                           mod1 <- coxme(Surv(time,status)~age+factor(sex)+sbp+(1|trial),
                                         data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictor
                           lp1 <- mod1$linear.predictor
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.re[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.re),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.tl3.tte.n1400 <- t(res.re.tl3[1:12, ])
se.re.tl3.tte.n1400 <- t(res.re.tl3[13:18, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.si <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+strata(trial),data = train0)
                           mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+strata(trial),data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.si[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.si),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.tl3.tte.n1400 <- t(res.si.tl3[1:12, ])
se.si.tl3.tte.n1400 <- t(res.si.tl3[13:18, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       source("coxvc.R")
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.r1 <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           Ft0<-cbind(rep(1,nrow(train0)),bs(train0$time))
                           mod0 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+trial, Ft0, rank = 1, data = train0)
                           Ft1<-cbind(rep(1,nrow(train0)),bs(train1$time))
                           mod1 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+trial, Ft1, rank = 1, data = train1)
                           
                           coef0 <- mod0$b
                           coef1 <- mod1$b
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           
                           #recalibration of event rates adapted to train
                           lp0 <- unlist(as.matrix(train0[2:4]) %*% coef[1:3])
                           lp1 <- unlist(as.matrix(train1[2:4]) %*% coef[1:3])
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.r1[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.r1),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                       
                     }
stopCluster(cl)
res.r1.tl3.tte.n1400 <- t(res.r1.tl3[1:12, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.fs <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 23, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 23, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 23, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 23, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(trial), data = train0)
                           mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(trial), data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(mean(coef0[4:8])*test$trial+
                                                                    mean(coef0[c(1,9:13)])*test$age+mean(coef0[c(2,14:18)])*test$sex+
                                                                    mean(coef0[c(3,19:23)])*test$sbp))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(mean(coef1[4:8])*test$trial+
                                                                    mean(coef1[c(1,9:13)])*test$age+mean(coef1[c(2,14:18)])*test$sex+
                                                                    mean(coef1[c(3,19:23)])*test$sbp))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.fs[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.fs),
                           coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[9:13, ]),
                           coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[14:18, ]),
                           coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[19:23, ]),
                           coef_mat.4 = mean(int1)-mean(int0),
                           coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[9:13, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[9:13, ])),
                           coef_mat.6 = (mean(coef_mat1[2, ]) + mean(coef_mat1[14:18, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[14:18, ])),
                           coef_mat.7 = (mean(coef_mat1[3, ]) + mean(coef_mat1[19:23, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[19:23, ])),
                           se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[9:13, ]),
                           se_mat.2 = mean(se_mat0[2, ]) + mean(se_mat0[14:18, ]),
                           se_mat.3 = mean(se_mat0[3, ]) + mean(se_mat0[19:23, ]),
                           se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[9:13, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[9:13, ])),
                           se_mat.6 = (mean(se_mat1[2, ]) + mean(se_mat1[14:18, ])) - (mean(se_mat0[2, ]) + mean(se_mat0[14:18, ])),
                           se_mat.7 = (mean(se_mat1[3, ]) + mean(se_mat1[19:23, ])) - (mean(se_mat0[3, ]) + mean(se_mat0[19:23, ]))),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.fs.tl3.tte.n1400 <- t(res.fs.tl3[1:12, ])
se.fs.tl3.tte.n1400 <- t(res.fs.tl3[13:18, ])

save(res.na.sl3.tte.n1400,se.na.sl3.tte.n1400,
     res.re.sl3.tte.n1400,se.re.sl3.tte.n1400,
     res.si.sl3.tte.n1400,se.si.sl3.tte.n1400,
     res.fs.sl3.tte.n1400,se.fs.sl3.tte.n1400,
     res.na.sl3.tte.n1400,se.na.sl3.tte.n1400,
     res.re.tl3.tte.n1400,se.re.tl3.tte.n1400,
     res.si.tl3.tte.n1400,se.si.tl3.tte.n1400,
     res.fs.tl3.tte.n1400,se.fs.tl3.tte.n1400,
     file = "res_scenario14.Rdata") #res.r1.sl3.tte.n1400, res.r1.tl3.tte.n1400,


#### scenario 15: 3 covariates & tte outcome & total sample size = 700 ####

#### s-learner ####

#### naive model ####
m <-1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.na[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.na),
                           coef_mat.1 = mean(coef_mat[1, ]),
                           coef_mat.2 = mean(coef_mat[2, ]),
                           coef_mat.3 = mean(coef_mat[3, ]),
                           coef_mat.4 = mean(coef_mat[4, ]),
                           coef_mat.5 = mean(coef_mat[5, ]),
                           coef_mat.6 = mean(coef_mat[6, ]),
                           coef_mat.7 = mean(coef_mat[7, ]),
                           se_mat.1 = mean(se_mat[1, ]),
                           se_mat.2 = mean(se_mat[2, ]),
                           se_mat.3 = mean(se_mat[3, ]),
                           se_mat.4 = mean(se_mat[4, ]),
                           se_mat.5 = mean(se_mat[5, ]),
                           se_mat.6 = mean(se_mat[6, ]),
                           se_mat.7 = mean(se_mat[7, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.sl3.tte.n700 <- t(res.na.sl3[1:12, ])
se.na.sl3.tte.n700 <- t(res.na.sl3[13:19, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying the model to train
                           mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+(1+treat|trial),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.re[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.re),
                           coef_mat.1 = mean(coef_mat[1,]),
                           coef_mat.2 = mean(coef_mat[2,]),
                           coef_mat.3 = mean(coef_mat[3,]),
                           coef_mat.4 = mean(coef_mat[4,]),
                           coef_mat.5 = mean(coef_mat[5,]),
                           coef_mat.6 = mean(coef_mat[6,]),
                           coef_mat.7 = mean(coef_mat[7,]),
                           se_mat.1 = mean(se_mat[1,]),
                           se_mat.2 = mean(se_mat[2,]),
                           se_mat.3 = mean(se_mat[3,]),
                           se_mat.4 = mean(se_mat[4,]),
                           se_mat.5 = mean(se_mat[5,]),
                           se_mat.6 = mean(se_mat[6,]),
                           se_mat.7 = mean(se_mat[7,])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.sl3.tte.n700 <- t(res.re.sl3[1:12, ])
se.re.sl3.tte.n700 <- t(res.re.sl3[13:19, ])

#### stratified intercept ####
m <-1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.si <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+strata(trial)+(0+treat|trial),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictor
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.si[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.si),
                           coef_mat.1 = mean(coef_mat[1, ]),
                           coef_mat.2 = mean(coef_mat[2, ]),
                           coef_mat.3 = mean(coef_mat[3, ]),
                           coef_mat.4 = mean(coef_mat[4, ]),
                           coef_mat.5 = mean(coef_mat[5, ]),
                           coef_mat.6 = mean(coef_mat[6, ]),
                           coef_mat.7 = mean(coef_mat[7, ]),
                           se_mat.1 = mean(se_mat[1, ]),
                           se_mat.2 = mean(se_mat[2, ]),
                           se_mat.3 = mean(se_mat[3, ]),
                           se_mat.4 = mean(se_mat[4, ]),
                           se_mat.5 = mean(se_mat[5, ]),
                           se_mat.6 = mean(se_mat[6, ]),
                           se_mat.7 = mean(se_mat[7, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.sl3.tte.n700 <- t(res.si.sl3[1:12, ])
se.si.sl3.tte.n700 <- t(res.si.sl3[13:19, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       source("coxvc.R")
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.r1 <- c()
                       coef_mat <- matrix(NA, nrow = 8, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           Ft <- cbind(rep(1,nrow(train)),bs(train$time))
                           mod <- coxvc(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat)+trial, Ft, rank = 1, data = train)
                           coef <- mod$b # maybe t(mod$b) ?
                           coef_mat[,i] <- coef
                           
                           #recalibration of event rates adapted to train
                           lp <- unlist(as.matrix(train[2:5]) %*% coef[1:4] + train[5] * (as.matrix(train[2:4]) %*% coef[6:8]))
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[6:8])))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[6:8])))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben <- cstat[1]
                           c.ben.se <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.r1[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.r1),
                         coef_mat.1 = mean(coef_mat[1, ]),
                         coef_mat.2 = mean(coef_mat[2, ]),
                         coef_mat.3 = mean(coef_mat[3, ]),
                         coef_mat.4 = mean(coef_mat[4, ]),
                         coef_mat.5 = mean(coef_mat[6, ]),
                         coef_mat.6 = mean(coef_mat[7, ]),
                         coef_mat.7 = mean(coef_mat[8, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.r1.sl3.tte.n700 <- t(res.r1.sl3[1:12, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.fs <- c()
                       coef_mat <- matrix(NA, nrow = 27, ncol = k)
                       se_mat <- matrix(NA, nrow = 27, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train <- df3[df3$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial),
                                        data = train)
                           coef <- mod$coefficients
                           coef_mat[,i] <- coef
                           se_mat[,i] <- sqrt(diag(vcov(mod)))
                           
                           #recalibration of event rates adapted to train
                           lp <- mod$linear.predictors
                           temp <- coxph(Surv(train$time,train$status)~offset(lp))
                           bh <- basehaz(temp, centered=F)
                           
                           #ite estimation
                           pred0 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test0$trial+
                                                                    mean(coef[c(6,13:17)])*test0$age+mean(coef[c(7,18:22)])*test0$sex+
                                                                    mean(coef[c(8,23:27)])*test0$sbp+coef[9]*test0$treat+
                                                                    (coef[10]*test0$age+coef[11]*test0$sex+coef[12]*test0$sbp)*test0$treat))
                           pred1 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test1$trial+
                                                                    mean(coef[c(6,13:17)])*test1$age+mean(coef[c(7,18:22)])*test1$sex+
                                                                    mean(coef[c(8,23:27)])*test1$sbp+coef[9]*test1$treat+
                                                                    (coef[10]*test1$age+coef[11]*test1$sex+coef[12]*test1$sbp)*test1$treat))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.fs[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.fs),
                           coef_mat.1 = mean(c(coef_mat[6,]+coef_mat[13:17,])),
                           coef_mat.2 = mean(c(coef_mat[7,]+coef_mat[18:22,])),
                           coef_mat.3 = mean(c(coef_mat[8,]+coef_mat[23:27,])),
                           coef_mat.4 = mean(coef_mat[9,]),
                           coef_mat.5 = mean(coef_mat[10,]),
                           coef_mat.6 = mean(coef_mat[11,]),
                           coef_mat.7 = mean(coef_mat[12,]),
                           se_mat.1 = mean(c(se_mat[6,]+se_mat[13:17,])),
                           se_mat.2 = mean(c(se_mat[7,]+se_mat[18:22,])),
                           se_mat.3 = mean(c(se_mat[8,]+se_mat[23:27,])),
                           se_mat.4 = mean(se_mat[9,]),
                           se_mat.5 = mean(se_mat[10,]),
                           se_mat.6 = mean(se_mat[11,]),
                           se_mat.7 = mean(se_mat[12,])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.fs.sl3.tte.n700 <- t(res.fs.sl3[1:12, ])
se.fs.sl3.tte.n700 <- t(res.fs.sl3[13:19, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train0)
                           mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.na[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.na),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.tl3.tte.n700 <- t(res.na.tl3[1:12, ])
se.na.tl3.tte.n700 <- t(res.na.tl3[13:18, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxme(Surv(time,status)~age+factor(sex)+sbp+(1|trial),
                                         data = train0)
                           mod1 <- coxme(Surv(time,status)~age+factor(sex)+sbp+(1|trial),
                                         data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictor
                           lp1 <- mod1$linear.predictor
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.re[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.re),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.tl3.tte.n700 <- t(res.re.tl3[1:12, ])
se.re.tl3.tte.n700 <- t(res.re.tl3[13:18, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.si <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+strata(trial),data = train0)
                           mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+strata(trial),data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.si[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.si),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                         se_mat.1 = mean(se_mat0[1, ]),
                         se_mat.2 = mean(se_mat0[2, ]),
                         se_mat.3 = mean(se_mat0[3, ]),
                         se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                         se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                         se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.tl3.tte.n700 <- t(res.si.tl3[1:12, ])
se.si.tl3.tte.n700 <- t(res.si.tl3[13:18, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       source("coxvc.R")
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.r1 <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           Ft0<-cbind(rep(1,nrow(train0)),bs(train0$time))
                           mod0 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+trial, Ft0, rank = 1, data = train0)
                           Ft1<-cbind(rep(1,nrow(train0)),bs(train1$time))
                           mod1 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+trial, Ft1, rank = 1, data = train1)
                           
                           coef0 <- mod0$b
                           coef1 <- mod1$b
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           
                           #recalibration of event rates adapted to train
                           lp0 <- unlist(as.matrix(train0[2:4]) %*% coef[1:3])
                           lp1 <- unlist(as.matrix(train1[2:4]) %*% coef[1:3])
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:4]) %*% coef0[1:3]))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:4]) %*% coef1[1:3]))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.r1[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.r1),
                         coef_mat.1 = mean(coef_mat0[1, ]),
                         coef_mat.2 = mean(coef_mat0[2, ]),
                         coef_mat.3 = mean(coef_mat0[3, ]),
                         coef_mat.4 = mean(int1)-mean(int0),
                         coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                         coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                         coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ])),
                         seed = j,
                         success = TRUE)
                       
                     }
stopCluster(cl)
res.r1.tl3.tte.n700 <- t(res.r1.tl3[1:12, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       
                       trial <- rep(1:k, each = floor(n/k))
                       df3 <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df3$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df3$age <- df3$age-70
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df3$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df3$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df3$sbp <- df3$sbp-185
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.1
                       b4 <- -0.3
                       b5 <- 0.01
                       b6 <- 0.04
                       b7 <- -0.2
                       
                       df3$treat <- rep(c(0,1), times = n/(2*k))
                       
                       trial_eff <- rep(0,7) 
                       
                       log_hr <- with(df3, b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
                       log_hr <- log_hr - mean(log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df3$time <- pmin(t,c)  # observed time is min of censored and true
                       df3$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df3[df3$treat==0,]
                       df1 <- df3[df3$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.fs <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 23, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 23, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 23, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 23, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df3[df3$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(trial), data = train0)
                           mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(trial), data = train1)
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef_mat0[,i] <- coef0
                           coef_mat1[,i] <- coef1
                           se_mat0[,i] <- sqrt(diag(vcov(mod0)))
                           se_mat1[,i] <- sqrt(diag(vcov(mod1)))
                           
                           #recalibration of event rates adapted to train
                           lp0 <- mod0$linear.predictors
                           lp1 <- mod1$linear.predictors
                           temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                           temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                           bh0 <- basehaz(temp0, centered=F)
                           bh1 <- basehaz(temp1, centered=F)
                           int0[i] <- bh0[nrow(bh0),]$hazard
                           int1[i] <- bh1[nrow(bh1),]$hazard
                           
                           #ite estimation
                           pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(mean(coef0[4:8])*test$trial+
                                                                    mean(coef0[c(1,9:13)])*test$age+mean(coef0[c(2,14:18)])*test$sex+
                                                                    mean(coef0[c(3,19:23)])*test$sbp))
                           pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(mean(coef1[4:8])*test$trial+
                                                                    mean(coef1[c(1,9:13)])*test$age+mean(coef1[c(2,14:18)])*test$sex+
                                                                    mean(coef1[c(3,19:23)])*test$sbp))
                           ite <- unlist(pred1 - pred0)
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration plot
                           predq5 <- quantcut(1-ite, q=5)
                           pben <-  tapply(1-ite, predq5, mean)
                           oben <- sapply(levels(predq5), function(x){1-diff(summary(survfit(Surv(time,status)~treat, data=test, subset=predq5==x), times=5, extend=TRUE)$surv)})
                           lm_ <- lm(oben~pben)
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.fs[i] <- mse(oben,pben)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(
                         res = c(
                           c.ben = mean(c.ben),
                           c.ben.se = mean(c.ben.se),
                           a = mean(a),
                           b = mean(b),
                           mse = mean(mse.fs),
                           coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[9:13, ]),
                           coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[14:18, ]),
                           coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[19:23, ]),
                           coef_mat.4 = mean(int1)-mean(int0),
                           coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[9:13, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[9:13, ])),
                           coef_mat.6 = (mean(coef_mat1[2, ]) + mean(coef_mat1[14:18, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[14:18, ])),
                           coef_mat.7 = (mean(coef_mat1[3, ]) + mean(coef_mat1[19:23, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[19:23, ])),
                           se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[9:13, ]),
                           se_mat.2 = mean(se_mat0[2, ]) + mean(se_mat0[14:18, ]),
                           se_mat.3 = mean(se_mat0[3, ]) + mean(se_mat0[19:23, ]),
                           se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[9:13, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[9:13, ])),
                           se_mat.6 = (mean(se_mat1[2, ]) + mean(se_mat1[14:18, ])) - (mean(se_mat0[2, ]) + mean(se_mat0[14:18, ])),
                           se_mat.7 = (mean(se_mat1[3, ]) + mean(se_mat1[19:23, ])) - (mean(se_mat0[3, ]) + mean(se_mat0[19:23, ]))),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.fs.tl3.tte.n700 <- t(res.fs.tl3[1:12, ])
se.fs.tl3.tte.n700 <- t(res.fs.tl3[13:18, ])

save(res.na.sl3.tte.n700,se.na.sl3.tte.n700,
     res.re.sl3.tte.n700,se.re.sl3.tte.n700,
     res.si.sl3.tte.n700,se.si.sl3.tte.n700,
     res.fs.sl3.tte.n700,se.fs.sl3.tte.n700,
     res.na.tl3.tte.n700,se.na.tl3.tte.n700,
     res.re.tl3.tte.n700,se.re.tl3.tte.n700,
     res.si.tl3.tte.n700,se.si.tl3.tte.n700,
     res.fs.tl3.tte.n700,se.fs.tl3.tte.n700,
     file = "res_scenario15.Rdata") #res.r1.sl3.tte.n700, res.r1.tl3.tte.n700,
