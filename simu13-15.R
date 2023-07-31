library(tidyverse) # data.frame manipulation
library(Hmisc) # cstat computation
library(magrittr) # cstat computation
library(doParallel) # parallel computing
library(gtools) # calibration computation
library(survival) # TTE outcome
library(coxme) # cox regression with mixed effects
library(mvmeta) # meta-analysis
library(msm) # delta method

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
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                       b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                                  expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 7, ncol = k)
                       se_mat <- matrix(NA, nrow = 7, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                  
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
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                               (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                           
                           true_val <- unlist(ite_true[i])
                           pred_int <- data.frame(lwr = ite - qt(.975,length(mod$residuals)) * se,
                                                  upr = ite + qt(.975,length(mod$residuals)) * se)
                           
                           in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration for benefit
                           lm_ <- lm((ite~unlist(ite_true[i])))
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- 1+lm_$coefficients[[2]]
                           mse.na[i] <- mse(unlist(ite_true[i]),ite)
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
                           in_int = mean(in_int),
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
res.na.sl.sim13 <- t(res.na.sl13[1:13, ])
se.na.sl.sim13 <- t(res.na.sl13[14:20, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                         b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                         expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 7, ncol = k)
  se_mat <- matrix(NA, nrow = 7, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                        (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]][1:4]),
                           upr = ite + qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]][1:4]))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- 1+lm_$coefficients[[2]]
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
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
    in_int = mean(in_int),
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
res.re.sl.sim13 <- t(res.re.sl13[1:13, ])
se.re.sl.sim13 <- t(res.re.sl13[14:20, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                         b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                         expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 7, ncol = k)
  se_mat <- matrix(NA, nrow = 7, ncol = k)
  
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                        (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]]),
                           upr = ite + qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]]))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- 1+lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
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
    in_int = mean(in_int),
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
res.si.sl.sim13 <- t(res.si.sl13[1:13, ])
se.si.sl.sim13 <- t(res.si.sl13[14:20, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl13 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 2800
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.r1 <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 8, ncol = k)
                        se_mat <- matrix(NA, nrow = 8, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat),
                                          data = train)
                            lin.pred <- mod1$linear.predictors
                            mod2 <- coxme(Surv(time,status)~(1+lin.pred|trial),
                                          data = train)
                            coef <- c(mod1$coefficients,mod2$frail[[1]][7:12])
                            
                            K.r1 <- array(0, dim=c(8, length(coef), 6))
                            ind <- c(rep(T,7),rep(F,length(coef)-7))
                            for(p in 1:6){
                              diag(K.r1[1:7,ind,p]) <- 1
                              K.r1[8,7+p,p] <- 1
                            }
                            allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
                            allcoef <- allcoef[,c(8,1:7)]
                            X <- model.matrix(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(treat),data = train)
                            
                            allvar <- apply(K.r1, 3, function(x){x %*% cov(X) %*% t(x)}, simplify=F)
                            
                            r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
                            coef = r1_meta$coefficients
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(r1_meta$vcov)
                            
                            #recalibration of event rates adapted to train
                            lp <- unlist(as.matrix(train[1:5]) %*% coef[1:5] + train[5] * (as.matrix(train[2:4]) %*% coef[6:8]))
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[1:5]) %*% coef[1:5] + test0[5] * (as.matrix(test0[2:4]) %*% coef[6:8])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[1:5]) %*% coef[1:5] + test1[5] * (as.matrix(test1[2:4]) %*% coef[6:8])))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                                                (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, r1_meta$vcov)
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod1$residuals)) * se * sqrt(1 + diag(r1_meta$Psi)),
                                                   upr = ite + qt(.975,length(mod1$residuals)) * se * sqrt(1 + diag(r1_meta$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.r1[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.r1),
                          in_int = mean(in_int),
                          coef_mat.1 = mean(coef_mat[1, ]),
                          coef_mat.2 = mean(coef_mat[2, ]),
                          coef_mat.3 = mean(coef_mat[3, ]),
                          coef_mat.4 = mean(coef_mat[4, ]),
                          coef_mat.5 = mean(coef_mat[5, ]),
                          coef_mat.6 = mean(coef_mat[6, ]),
                          coef_mat.7 = mean(coef_mat[7, ]),
                          coef_mat.8 = mean(coef_mat[8, ]),
                          se_mat.1 = mean(se_mat[1, ]),
                          se_mat.2 = mean(se_mat[2, ]),
                          se_mat.3 = mean(se_mat[3, ]),
                          se_mat.4 = mean(se_mat[4, ]),
                          se_mat.5 = mean(se_mat[5, ]),
                          se_mat.6 = mean(se_mat[6, ]),
                          se_mat.7 = mean(se_mat[7, ]),
                          se_mat.8 = mean(se_mat[8, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.r1.sl.sim13 <- t(res.r1.sl13[1:14, ])
se.r1.sl.sim13 <- t(res.r1.sl13[15:22, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                         b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                         expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 7, ncol = k)
  se_mat <- matrix(NA, nrow = 7, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp)*treat*factor(trial),
                 data = train, control = coxph.control(iter.max = 10000))
    coef <- mod$coefficients
    
    K.fs <- array(0, dim=c(8, length(coef), 6))
    ind <- c(rep(F,5),rep(T,7),rep(F,length(coef)-12))
    if(i==1){vec <- 3:7} else {vec <- 2:7}
    for(p in 1:6){
      if(p==1){K.fs[1,p,p] <- 1} else {K.fs[1,p-1,p] <- 1}
      diag(K.fs[2:8,ind,p]) <- 1
      if(p == 1) next
      K.fs[2,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
      K.fs[3,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
      K.fs[4,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
      K.fs[5,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":treat")),p] <- 1
      K.fs[6,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age:treat")),p] <- 1
      K.fs[7,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1:treat")),p] <- 1
      K.fs[8,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp:treat")),p] <- 1
    }
    allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
    allvar <- apply(K.fs, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
    
    fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = fs_meta$coefficients[c(-1)]
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(fs_meta$vcov[c(-1),c(-1)])
    
    #recalibration of event rates adapted to train
    lp <- mod$linear.predictors
    temp <- coxph(Surv(train$time,train$status)~offset(lp))
    bh <- basehaz(temp, centered=F)
    
    #ite estimation
    pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
    pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
    ite <- unlist(pred1 - pred0)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                        (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, fs_meta$vcov[-c(1), -c(1)])
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,length(mod$residuals)) * se * sqrt(1 + diag(fs_meta$Psi)),
                           upr = ite + qt(.975,length(mod$residuals)) * se * sqrt(1 + diag(fs_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- 1+lm_$coefficients[[2]]
    mse.fs[i] <- mse(unlist(ite_true[i]),ite)
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
    in_int = mean(in_int),
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
res.fs.sl.sim13 <- t(res.fs.sl13[1:13, ])
se.fs.sl.sim13 <- t(res.fs.sl13[14:20, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                                   (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                      (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                                 expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df[df$treat==0,]
                       df1 <- df[df$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       in_int <- c()
                       int0 <- c()
                       int1 <- c()
                       coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
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
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                           se <- max(se0,se1)
                           true_val <- unlist(ite_true[i])
                           pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$residuals)) * se,
                                                  upr = ite + qt(.975,length(mod0$residuals)) * se)
                           
                           in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$status),
                                             treatment = test$treat == 1,
                                             score = 1-ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration for benefit
                           lm_ <- lm(ite~unlist(ite_true[i]))
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- 1+lm_$coefficients[[2]]
                           mse.na[i] <- mse(unlist(ite_true[i]),ite)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.na),
                         in_int = mean(in_int),
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
res.na.tl.sim13 <- t(res.na.tl13[1:13, ])
se.na.tl.sim13 <- t(res.na.tl13[14:20, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                         (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                         (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                         expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df[df$treat==0,]
                       df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
  se_mat0 <- matrix(NA, nrow = 3, ncol = k)
  se_mat1 <- matrix(NA, nrow = 3, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + VarCorr(mod0)[[1]]),
                           upr = ite + qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + VarCorr(mod0)[[1]]))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- 1+lm_$coefficients[[2]]
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    in_int = mean(in_int),
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
res.re.tl.sim13 <- t(res.re.tl13[1:13, ])
se.re.tl.sim13<- t(res.re.tl13[14:19, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                         (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                         (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                         expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df[df$treat==0,]
                       df1 <- df[df$treat==1,]
                       
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
  se_mat0 <- matrix(NA, nrow = 3, ncol = k)
  se_mat1 <- matrix(NA, nrow = 3, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se,
                           upr = ite + qt(.975,length(mod0$linear.predictor)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- 1+lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    in_int = mean(in_int),
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
res.si.tl.sim13 <- t(res.si.tl13[1:13, ])
se.si.tl.sim13 <- t(res.si.tl13[13:19, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl13 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 2800
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.r1 <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 4, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 4, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod1.0 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train0)
                            mod1.1 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train1)
                            lin.pred0 <- mod1.0$linear.predictors
                            lin.pred1 <- mod1.1$linear.predictors
                            
                            mod2.0 <- coxme(Surv(time,status)~(1+lin.pred0|trial), data = train0)
                            mod2.1 <- coxme(Surv(time,status)~(1+lin.pred1|trial), data = train0)
                            coef0 <- c(mod1.0$coefficients,mod2.0$frail[[1]][7:12])
                            coef1 <- c(mod1.1$coefficients,mod2.1$frail[[1]][7:12])
                            
                            K.r1.0 <- array(0, dim=c(4, length(coef0), 6))
                            ind <- c(rep(T,3),rep(F,length(coef0)-3))
                            for(p in 1:6){
                              diag(K.r1.0[1:3,ind,p]) <- 1
                              K.r1.0[4,3+p,p] <- 1
                            }
                            allcoef0 <- t(apply(K.r1.0, 3, function(x){x %*% coef0}))
                            allcoef0 <- allcoef0[,c(4,1:3)]
                            X0 <- model.matrix(Surv(time,status)~-1+factor(trial)+age+factor(sex)+sbp,data = train0)
                            allvar0 <- apply(K.r1.0, 3, function(x){x %*% cov(X0) %*% t(x)}, simplify=F)
                            
                            r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                            coef0 = r1_meta0$coefficients
                            coef_mat0[,i] <- coef0
                            se_mat0[,i] <- diag(r1_meta0$vcov)
                            
                            K.r1.1 <- array(0, dim=c(4, length(coef1), 6))
                            ind <- c(rep(T,3),rep(F,length(coef1)-3))
                            for(p in 1:6){
                              diag(K.r1.1[1:3,ind,p]) <- 1
                              K.r1.1[4,3+p,p] <- 1
                            }
                            allcoef1 <- t(apply(K.r1.1, 3, function(x){x %*% coef1}))
                            allcoef1 <- allcoef1[,c(4,1:3)]
                            X1 <- model.matrix(Surv(time,status)~-1+factor(trial)+age+factor(sex)+sbp,data = train1)
                            allvar1 <- apply(K.r1.1, 3, function(x){x %*% cov(X1) %*% t(x)}, simplify=F)
                            
                            r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                            coef1 = r1_meta1$coefficients
                            coef_mat1[,i] <- coef1
                            se_mat1[,i] <- diag(r1_meta1$vcov)
                            
                            #recalibration of event rates adapted to train
                            lp0 <- unlist(as.matrix(train0[1:4]) %*% coef0[1:4])
                            lp1 <- unlist(as.matrix(train1[1:4]) %*% coef1[1:4])
                            temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                            temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                            bh0 <- basehaz(temp0, centered=F)
                            bh1 <- basehaz(temp1, centered=F)
                            int0[i] <- bh0[nrow(bh0),]$hazard
                            int1[i] <- bh1[nrow(bh1),]$hazard
                            
                            #ite estimation
                            pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[1:4]) %*% coef0[1:4]))
                            pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[1:4]) %*% coef1[1:4]))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, r1_meta0$vcov)
                            se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, r1_meta1$vcov)
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod1.0$residuals)) * se * sqrt(1 + diag(r1_meta0$Psi)),
                                                   upr = ite + qt(.975,length(mod1.0$residuals)) * se * sqrt(1 + diag(r1_meta0$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.r1[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.r1),
                          in_int = mean(in_int),
                          coef_mat.1 = mean(coef_mat0[1, ]),
                          coef_mat.2 = mean(coef_mat0[2, ]),
                          coef_mat.3 = mean(coef_mat0[3, ]),
                          coef_mat.4 = mean(coef_mat0[4, ]),
                          coef_mat.5 = mean(int1)-mean(int0),
                          coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                          coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                          coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                          se_mat.1 = mean(se_mat0[1, ]),
                          se_mat.2 = mean(se_mat0[2, ]),
                          se_mat.3 = mean(se_mat0[3, ]),
                          se_mat.4 = mean(se_mat0[4, ]),
                          se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                          se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                          se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ])),
                          seed = j,
                          success = TRUE)
                        
                      }
stopCluster(cl)
res.r1.tl.sim13 <- t(res.r1.tl13[1:14, ])
se.r1.tl.sim13 <- t(res.r1.tl13[15:21, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl13 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82)
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age - mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197)
                       sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp - mean(df$sbp)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- 50
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                         (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                       df$log_hr <- df$log_hr - mean(df$log_hr)
                       df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                         (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                         expit(b1*df$age + b2*df$sex + b3*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       sp <- 1.15
                       t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                       c <- runif(n, 0.5, 10)
                       c[c>5] <- 5 # censoring
                       df$time <- pmin(t,c)  # observed time is min of censored and true
                       df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                       
                       df0 <- df[df$treat==0,]
                       df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  in_int <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
  se_mat0 <- matrix(NA, nrow = 3, ncol = k)
  se_mat1 <- matrix(NA, nrow = 3, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- coxph(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(trial), data = train0)
    mod1 <- coxph(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(trial), data = train1)
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    
    K.fs0 <- array(0, dim=c(4, length(coef0), 6))
    ind <- c(rep(F,5),rep(T,3),rep(F,length(coef0)-8))
    if(i==1){vec <- 3:7} else {vec <- 2:7}
    for(p in 1:6){
      if(p==1){K.fs0[1,p,p] <- 1} else {K.fs0[1,p-1,p] <- 1}
      diag(K.fs0[2:4,ind,p]) <- 1
      if(p == 1) next
      K.fs0[2,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
      K.fs0[3,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
      K.fs0[4,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
    }
    allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.fs0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
    
    fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = fs_meta0$coefficients[c(-1)]
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(fs_meta0$vcov[c(-1),c(-1)])
    
    K.fs1 <- array(0, dim=c(4, length(coef1), 6))
    ind <- c(rep(F,5),rep(T,3),rep(F,length(coef1)-8))
    if(i==1){vec <- 3:7} else {vec <- 2:7}
    for(p in 1:6){
      if(p==1){K.fs1[1,p,p] <- 1} else {K.fs1[1,p-1,p] <- 1}
      diag(K.fs1[2:4,ind,p]) <- 1
      if(p == 1) next
      K.fs1[2,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
      K.fs1[3,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
      K.fs1[4,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
    }
    allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.fs1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
    
    fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = fs_meta1$coefficients[c(-1)]
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(fs_meta1$vcov[c(-1),c(-1)])
    
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
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, fs_meta0$vcov[c(-1),c(-1)])
    se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, fs_meta1$vcov[c(-1),c(-1)])
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + diag(fs_meta0$Psi)),
                           upr = ite + qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + diag(fs_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$status),
                      treatment = test$treat == 1,
                      score = 1-ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- 1+lm_$coefficients[[2]]
    mse.fs[i] <- mse(unlist(ite_true[i]),ite)
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
    in_int = mean(in_int),
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
res.fs.tl.sim13 <- t(res.fs.tl13[1:13, ])
se.fs.tl.sim13 <- t(res.fs.tl13[14:19, ])

save(res.na.sl.sim13,se.na.sl.sim13,
     res.re.sl.sim13,se.re.sl.sim13,
     res.si.sl.sim13,se.si.sl.sim13,
     res.r1.sl.sim13,se.r1.sl.sim13,
     res.fs.sl.sim13,se.fs.sl.sim13,
     res.na.tl.sim13,se.na.tl.sim13,
     res.re.tl.sim13,se.re.tl.sim13,
     res.si.tl.sim13,se.si.tl.sim13,
     res.r1.tl.sim13,se.r1.tl.sim13,
     res.fs.tl.sim13,se.fs.tl.sim13,
     file = "res_scenario13.Rdata")

#### scenario 14: 3 covariates & tte outcome & total sample size = 1400 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
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
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$residuals)) * se,
                                                   upr = ite + qt(.975,length(mod$residuals)) * se)
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm((ite~unlist(ite_true[i])))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.na[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.na.sl.sim14 <- t(res.na.sl14[1:13, ])
se.na.sl.sim14 <- t(res.na.sl14[14:20, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.re <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
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
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]][1:4]),
                                                   upr = ite + qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]][1:4]))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.re[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.re.sl.sim14 <- t(res.re.sl14[1:13, ])
se.re.sl.sim14 <- t(res.re.sl14[14:20, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.si <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
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
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]]),
                                                   upr = ite + qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]]))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.si[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.si.sl.sim14 <- t(res.si.sl14[1:13, ])
se.si.sl.sim14 <- t(res.si.sl14[14:20, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.r1 <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 8, ncol = k)
                        se_mat <- matrix(NA, nrow = 8, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat),
                                          data = train)
                            lin.pred <- mod1$linear.predictors
                            mod2 <- coxme(Surv(time,status)~(1+lin.pred|trial),
                                          data = train)
                            coef <- c(mod1$coefficients,mod2$frail[[1]][7:12])
                            
                            K.r1 <- array(0, dim=c(8, length(coef), 6))
                            ind <- c(rep(T,7),rep(F,length(coef)-7))
                            for(p in 1:6){
                              diag(K.r1[1:7,ind,p]) <- 1
                              K.r1[8,7+p,p] <- 1
                            }
                            allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
                            allcoef <- allcoef[,c(8,1:7)]
                            X <- model.matrix(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(treat),data = train)
                            
                            allvar <- apply(K.r1, 3, function(x){x %*% cov(X) %*% t(x)}, simplify=F)
                            
                            r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
                            coef = r1_meta$coefficients
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(r1_meta$vcov)
                            
                            #recalibration of event rates adapted to train
                            lp <- unlist(as.matrix(train[1:5]) %*% coef[1:5] + train[5] * (as.matrix(train[2:4]) %*% coef[6:8]))
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[1:5]) %*% coef[1:5] + test0[5] * (as.matrix(test0[2:4]) %*% coef[6:8])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[1:5]) %*% coef[1:5] + test1[5] * (as.matrix(test1[2:4]) %*% coef[6:8])))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                                                (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, r1_meta$vcov)
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod1$residuals)) * se * sqrt(1 + diag(r1_meta$Psi)),
                                                   upr = ite + qt(.975,length(mod1$residuals)) * se * sqrt(1 + diag(r1_meta$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.r1[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.r1),
                          in_int = mean(in_int),
                          coef_mat.1 = mean(coef_mat[1, ]),
                          coef_mat.2 = mean(coef_mat[2, ]),
                          coef_mat.3 = mean(coef_mat[3, ]),
                          coef_mat.4 = mean(coef_mat[4, ]),
                          coef_mat.5 = mean(coef_mat[5, ]),
                          coef_mat.6 = mean(coef_mat[6, ]),
                          coef_mat.7 = mean(coef_mat[7, ]),
                          coef_mat.8 = mean(coef_mat[8, ]),
                          se_mat.1 = mean(se_mat[1, ]),
                          se_mat.2 = mean(se_mat[2, ]),
                          se_mat.3 = mean(se_mat[3, ]),
                          se_mat.4 = mean(se_mat[4, ]),
                          se_mat.5 = mean(se_mat[5, ]),
                          se_mat.6 = mean(se_mat[6, ]),
                          se_mat.7 = mean(se_mat[7, ]),
                          se_mat.8 = mean(se_mat[8, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.r1.sl.sim14 <- t(res.r1.sl14[1:14, ])
se.r1.sl.sim14 <- t(res.r1.sl14[15:22, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp)*treat*factor(trial),
                                         data = train, control = coxph.control(iter.max = 10000))
                            coef <- mod$coefficients
                            
                            K.fs <- array(0, dim=c(8, length(coef), 6))
                            ind <- c(rep(F,5),rep(T,7),rep(F,length(coef)-12))
                            if(i==1){vec <- 3:7} else {vec <- 2:7}
                            for(p in 1:6){
                              if(p==1){K.fs[1,p,p] <- 1} else {K.fs[1,p-1,p] <- 1}
                              diag(K.fs[2:8,ind,p]) <- 1
                              if(p == 1) next
                              K.fs[2,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
                              K.fs[3,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
                              K.fs[4,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
                              K.fs[5,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":treat")),p] <- 1
                              K.fs[6,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age:treat")),p] <- 1
                              K.fs[7,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1:treat")),p] <- 1
                              K.fs[8,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp:treat")),p] <- 1
                            }
                            allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
                            allvar <- apply(K.fs, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
                            
                            fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
                            coef = fs_meta$coefficients[c(-1)]
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(fs_meta$vcov[c(-1),c(-1)])
                            
                            #recalibration of event rates adapted to train
                            lp <- mod$linear.predictors
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, fs_meta$vcov[-c(1), -c(1)])
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$residuals)) * se * sqrt(1 + diag(fs_meta$Psi)),
                                                   upr = ite + qt(.975,length(mod$residuals)) * se * sqrt(1 + diag(fs_meta$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.fs[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.fs.sl.sim14 <- t(res.fs.sl14[1:13, ])
se.fs.sl.sim14 <- t(res.fs.sl14[14:20, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$residuals)) * se,
                                                   upr = ite + qt(.975,length(mod0$residuals)) * se)
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.na[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.na),
                          in_int = mean(in_int),
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
res.na.tl.sim14 <- t(res.na.tl14[1:13, ])
se.na.tl.sim14 <- t(res.na.tl14[14:20, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.re <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + VarCorr(mod0)[[1]]),
                                                   upr = ite + qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + VarCorr(mod0)[[1]]))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.re[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.re),
                          in_int = mean(in_int),
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
res.re.tl.sim14 <- t(res.re.tl14[1:13, ])
se.re.tl.sim14<- t(res.re.tl14[14:19, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.si <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se,
                                                   upr = ite + qt(.975,length(mod0$linear.predictor)) * se)
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.si[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.si),
                          in_int = mean(in_int),
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
res.si.tl.sim14 <- t(res.si.tl14[1:13, ])
se.si.tl.sim14 <- t(res.si.tl14[13:19, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.r1 <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 4, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 4, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod1.0 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train0)
                            mod1.1 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train1)
                            lin.pred0 <- mod1.0$linear.predictors
                            lin.pred1 <- mod1.1$linear.predictors
                            
                            mod2.0 <- coxme(Surv(time,status)~(1+lin.pred0|trial), data = train0)
                            mod2.1 <- coxme(Surv(time,status)~(1+lin.pred1|trial), data = train0)
                            coef0 <- c(mod1.0$coefficients,mod2.0$frail[[1]][7:12])
                            coef1 <- c(mod1.1$coefficients,mod2.1$frail[[1]][7:12])
                            
                            K.r1.0 <- array(0, dim=c(4, length(coef0), 6))
                            ind <- c(rep(T,3),rep(F,length(coef0)-3))
                            for(p in 1:6){
                              diag(K.r1.0[1:3,ind,p]) <- 1
                              K.r1.0[4,3+p,p] <- 1
                            }
                            allcoef0 <- t(apply(K.r1.0, 3, function(x){x %*% coef0}))
                            allcoef0 <- allcoef0[,c(4,1:3)]
                            X0 <- model.matrix(Surv(time,status)~-1+factor(trial)+age+factor(sex)+sbp,data = train0)
                            allvar0 <- apply(K.r1.0, 3, function(x){x %*% cov(X0) %*% t(x)}, simplify=F)
                            
                            r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                            coef0 = r1_meta0$coefficients
                            coef_mat0[,i] <- coef0
                            se_mat0[,i] <- diag(r1_meta0$vcov)
                            
                            K.r1.1 <- array(0, dim=c(4, length(coef1), 6))
                            ind <- c(rep(T,3),rep(F,length(coef1)-3))
                            for(p in 1:6){
                              diag(K.r1.1[1:3,ind,p]) <- 1
                              K.r1.1[4,3+p,p] <- 1
                            }
                            allcoef1 <- t(apply(K.r1.1, 3, function(x){x %*% coef1}))
                            allcoef1 <- allcoef1[,c(4,1:3)]
                            X1 <- model.matrix(Surv(time,status)~-1+factor(trial)+age+factor(sex)+sbp,data = train1)
                            allvar1 <- apply(K.r1.1, 3, function(x){x %*% cov(X1) %*% t(x)}, simplify=F)
                            
                            r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                            coef1 = r1_meta1$coefficients
                            coef_mat1[,i] <- coef1
                            se_mat1[,i] <- diag(r1_meta1$vcov)
                            
                            #recalibration of event rates adapted to train
                            lp0 <- unlist(as.matrix(train0[1:4]) %*% coef0[1:4])
                            lp1 <- unlist(as.matrix(train1[1:4]) %*% coef1[1:4])
                            temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                            temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                            bh0 <- basehaz(temp0, centered=F)
                            bh1 <- basehaz(temp1, centered=F)
                            int0[i] <- bh0[nrow(bh0),]$hazard
                            int1[i] <- bh1[nrow(bh1),]$hazard
                            
                            #ite estimation
                            pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[1:4]) %*% coef0[1:4]))
                            pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[1:4]) %*% coef1[1:4]))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, r1_meta0$vcov)
                            se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, r1_meta1$vcov)
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod1.0$residuals)) * se * sqrt(1 + diag(r1_meta0$Psi)),
                                                   upr = ite + qt(.975,length(mod1.0$residuals)) * se * sqrt(1 + diag(r1_meta0$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.r1[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.r1),
                          in_int = mean(in_int),
                          coef_mat.1 = mean(coef_mat0[1, ]),
                          coef_mat.2 = mean(coef_mat0[2, ]),
                          coef_mat.3 = mean(coef_mat0[3, ]),
                          coef_mat.4 = mean(coef_mat0[4, ]),
                          coef_mat.5 = mean(int1)-mean(int0),
                          coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                          coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                          coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                          se_mat.1 = mean(se_mat0[1, ]),
                          se_mat.2 = mean(se_mat0[2, ]),
                          se_mat.3 = mean(se_mat0[3, ]),
                          se_mat.4 = mean(se_mat0[4, ]),
                          se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                          se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                          se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ])),
                          seed = j,
                          success = TRUE)
                        
                      }
stopCluster(cl)
res.r1.tl.sim14 <- t(res.r1.tl14[1:14, ])
se.r1.tl.sim14 <- t(res.r1.tl14[15:21, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl14 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- coxph(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(trial), data = train0)
                            mod1 <- coxph(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(trial), data = train1)
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            
                            K.fs0 <- array(0, dim=c(4, length(coef0), 6))
                            ind <- c(rep(F,5),rep(T,3),rep(F,length(coef0)-8))
                            if(i==1){vec <- 3:7} else {vec <- 2:7}
                            for(p in 1:6){
                              if(p==1){K.fs0[1,p,p] <- 1} else {K.fs0[1,p-1,p] <- 1}
                              diag(K.fs0[2:4,ind,p]) <- 1
                              if(p == 1) next
                              K.fs0[2,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
                              K.fs0[3,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
                              K.fs0[4,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
                            }
                            allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
                            allvar0 <- apply(K.fs0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
                            
                            fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                            coef0 = fs_meta0$coefficients[c(-1)]
                            
                            coef_mat0[,i] <- coef0
                            se_mat0[,i] <- diag(fs_meta0$vcov[c(-1),c(-1)])
                            
                            K.fs1 <- array(0, dim=c(4, length(coef1), 6))
                            ind <- c(rep(F,5),rep(T,3),rep(F,length(coef1)-8))
                            if(i==1){vec <- 3:7} else {vec <- 2:7}
                            for(p in 1:6){
                              if(p==1){K.fs1[1,p,p] <- 1} else {K.fs1[1,p-1,p] <- 1}
                              diag(K.fs1[2:4,ind,p]) <- 1
                              if(p == 1) next
                              K.fs1[2,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
                              K.fs1[3,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
                              K.fs1[4,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
                            }
                            allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
                            allvar1 <- apply(K.fs1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
                            
                            fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                            coef1 = fs_meta1$coefficients[c(-1)]
                            
                            coef_mat1[,i] <- coef1
                            se_mat1[,i] <- diag(fs_meta1$vcov[c(-1),c(-1)])
                            
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, fs_meta0$vcov[c(-1),c(-1)])
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, fs_meta1$vcov[c(-1),c(-1)])
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + diag(fs_meta0$Psi)),
                                                   upr = ite + qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + diag(fs_meta0$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.fs[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.fs.tl.sim14 <- t(res.fs.tl14[1:13, ])
se.fs.tl.sim14 <- t(res.fs.tl14[14:19, ])

save(res.na.sl.sim14,se.na.sl.sim14,
     res.re.sl.sim14,se.re.sl.sim14,
     res.si.sl.sim14,se.si.sl.sim14,
     res.r1.sl.sim14,se.r1.sl.sim14,
     res.fs.sl.sim14,se.fs.sl.sim14,
     res.na.tl.sim14,se.na.tl.sim14,
     res.re.tl.sim14,se.re.tl.sim14,
     res.r1.tl.sim14,se.r1.tl.sim14,
     res.si.tl.sim14,se.si.tl.sim14,
     res.fs.tl.sim14,se.fs.tl.sim14,
     file = "res_scenario14.Rdata")


#### scenario 15: 3 covariates & tte outcome & total sample size = 700 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
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
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$residuals)) * se,
                                                   upr = ite + qt(.975,length(mod$residuals)) * se)
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm((ite~unlist(ite_true[i])))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.na[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.na.sl.sim15 <- t(res.na.sl15[1:13, ])
se.na.sl.sim15 <- t(res.na.sl15[14:20, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.re <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
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
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]][1:4]),
                                                   upr = ite + qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]][1:4]))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.re[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.re.sl.sim15 <- t(res.re.sl15[1:13, ])
se.re.sl.sim15 <- t(res.re.sl15[14:20, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.si <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
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
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, vcov(mod))
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]]),
                                                   upr = ite + qt(.975,length(mod$linear.predictor)) * se * sqrt(1 + VarCorr(mod)[[1]]))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.si[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.si.sl.sim15 <- t(res.si.sl15[1:13, ])
se.si.sl.sim15 <- t(res.si.sl15[14:20, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.r1 <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 8, ncol = k)
                        se_mat <- matrix(NA, nrow = 8, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp)*factor(treat),
                                          data = train)
                            lin.pred <- mod1$linear.predictors
                            mod2 <- coxme(Surv(time,status)~(1+lin.pred|trial),
                                          data = train)
                            coef <- c(mod1$coefficients,mod2$frail[[1]][7:12])
                            
                            K.r1 <- array(0, dim=c(8, length(coef), 6))
                            ind <- c(rep(T,7),rep(F,length(coef)-7))
                            for(p in 1:6){
                              diag(K.r1[1:7,ind,p]) <- 1
                              K.r1[8,7+p,p] <- 1
                            }
                            allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
                            allcoef <- allcoef[,c(8,1:7)]
                            X <- model.matrix(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(treat),data = train)
                            
                            allvar <- apply(K.r1, 3, function(x){x %*% cov(X) %*% t(x)}, simplify=F)
                            
                            r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
                            coef = r1_meta$coefficients
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(r1_meta$vcov)
                            
                            #recalibration of event rates adapted to train
                            lp <- unlist(as.matrix(train[1:5]) %*% coef[1:5] + train[5] * (as.matrix(train[2:4]) %*% coef[6:8]))
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[1:5]) %*% coef[1:5] + test0[5] * (as.matrix(test0[2:4]) %*% coef[6:8])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[1:5]) %*% coef[1:5] + test1[5] * (as.matrix(test1[2:4]) %*% coef[6:8])))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                                                (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, r1_meta$vcov)
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod1$residuals)) * se * sqrt(1 + diag(r1_meta$Psi)),
                                                   upr = ite + qt(.975,length(mod1$residuals)) * se * sqrt(1 + diag(r1_meta$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.r1[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.r1),
                          in_int = mean(in_int),
                          coef_mat.1 = mean(coef_mat[1, ]),
                          coef_mat.2 = mean(coef_mat[2, ]),
                          coef_mat.3 = mean(coef_mat[3, ]),
                          coef_mat.4 = mean(coef_mat[4, ]),
                          coef_mat.5 = mean(coef_mat[5, ]),
                          coef_mat.6 = mean(coef_mat[6, ]),
                          coef_mat.7 = mean(coef_mat[7, ]),
                          coef_mat.8 = mean(coef_mat[8, ]),
                          se_mat.1 = mean(se_mat[1, ]),
                          se_mat.2 = mean(se_mat[2, ]),
                          se_mat.3 = mean(se_mat[3, ]),
                          se_mat.4 = mean(se_mat[4, ]),
                          se_mat.5 = mean(se_mat[5, ]),
                          se_mat.6 = mean(se_mat[6, ]),
                          se_mat.7 = mean(se_mat[7, ]),
                          se_mat.8 = mean(se_mat[8, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.r1.sl.sim15 <- t(res.r1.sl15[1:14, ])
se.r1.sl.sim15 <- t(res.r1.sl15[15:22, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        nt <- n/k
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(df$log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 7, ncol = k)
                        se_mat <- matrix(NA, nrow = 7, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp)*treat*factor(trial),
                                         data = train, control = coxph.control(iter.max = 10000))
                            coef <- mod$coefficients
                            
                            K.fs <- array(0, dim=c(8, length(coef), 6))
                            ind <- c(rep(F,5),rep(T,7),rep(F,length(coef)-12))
                            if(i==1){vec <- 3:7} else {vec <- 2:7}
                            for(p in 1:6){
                              if(p==1){K.fs[1,p,p] <- 1} else {K.fs[1,p-1,p] <- 1}
                              diag(K.fs[2:8,ind,p]) <- 1
                              if(p == 1) next
                              K.fs[2,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
                              K.fs[3,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
                              K.fs[4,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
                              K.fs[5,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":treat")),p] <- 1
                              K.fs[6,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age:treat")),p] <- 1
                              K.fs[7,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1:treat")),p] <- 1
                              K.fs[8,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp:treat")),p] <- 1
                            }
                            allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
                            allvar <- apply(K.fs, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
                            
                            fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
                            coef = fs_meta$coefficients[c(-1)]
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(fs_meta$vcov[c(-1),c(-1)])
                            
                            #recalibration of event rates adapted to train
                            lp <- mod$linear.predictors
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:5]) %*% coef[1:4] + test0[5] * (as.matrix(test0[2:4]) %*% coef[5:7])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:5]) %*% coef[1:4] + test1[5] * (as.matrix(test1[2:4]) %*% coef[5:7])))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7)/(1+exp(x1+x2+x3+x4+x5+x6+x7))) - 
                                                (exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef, fs_meta$vcov[-c(1), -c(1)])
                            
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod$residuals)) * se * sqrt(1 + diag(fs_meta$Psi)),
                                                   upr = ite + qt(.975,length(mod$residuals)) * se * sqrt(1 + diag(fs_meta$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.fs[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.fs.sl.sim15 <- t(res.fs.sl15[1:13, ])
se.fs.sl.sim15 <- t(res.fs.sl15[14:20, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$residuals)) * se,
                                                   upr = ite + qt(.975,length(mod0$residuals)) * se)
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.na[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.na),
                          in_int = mean(in_int),
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
res.na.tl.sim15 <- t(res.na.tl15[1:13, ])
se.na.tl.sim15 <- t(res.na.tl15[14:20, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.re <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + VarCorr(mod0)[[1]]),
                                                   upr = ite + qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + VarCorr(mod0)[[1]]))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.re[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.re),
                          in_int = mean(in_int),
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
res.re.tl.sim15 <- t(res.re.tl15[1:13, ])
se.re.tl.sim15<- t(res.re.tl15[14:19, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.si <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, vcov(mod0))
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, vcov(mod1))
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se,
                                                   upr = ite + qt(.975,length(mod0$linear.predictor)) * se)
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.si[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.si),
                          in_int = mean(in_int),
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
res.si.tl.sim15 <- t(res.si.tl15[1:13, ])
se.si.tl.sim15 <- t(res.si.tl15[13:19, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.r1 <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 4, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 4, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod1.0 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train0)
                            mod1.1 <- coxph(Surv(time,status)~age+factor(sex)+sbp,data = train1)
                            lin.pred0 <- mod1.0$linear.predictors
                            lin.pred1 <- mod1.1$linear.predictors
                            
                            mod2.0 <- coxme(Surv(time,status)~(1+lin.pred0|trial), data = train0)
                            mod2.1 <- coxme(Surv(time,status)~(1+lin.pred1|trial), data = train0)
                            coef0 <- c(mod1.0$coefficients,mod2.0$frail[[1]][7:12])
                            coef1 <- c(mod1.1$coefficients,mod2.1$frail[[1]][7:12])
                            
                            K.r1.0 <- array(0, dim=c(4, length(coef0), 6))
                            ind <- c(rep(T,3),rep(F,length(coef0)-3))
                            for(p in 1:6){
                              diag(K.r1.0[1:3,ind,p]) <- 1
                              K.r1.0[4,3+p,p] <- 1
                            }
                            allcoef0 <- t(apply(K.r1.0, 3, function(x){x %*% coef0}))
                            allcoef0 <- allcoef0[,c(4,1:3)]
                            X0 <- model.matrix(Surv(time,status)~-1+factor(trial)+age+factor(sex)+sbp,data = train0)
                            allvar0 <- apply(K.r1.0, 3, function(x){x %*% cov(X0) %*% t(x)}, simplify=F)
                            
                            r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                            coef0 = r1_meta0$coefficients
                            coef_mat0[,i] <- coef0
                            se_mat0[,i] <- diag(r1_meta0$vcov)
                            
                            K.r1.1 <- array(0, dim=c(4, length(coef1), 6))
                            ind <- c(rep(T,3),rep(F,length(coef1)-3))
                            for(p in 1:6){
                              diag(K.r1.1[1:3,ind,p]) <- 1
                              K.r1.1[4,3+p,p] <- 1
                            }
                            allcoef1 <- t(apply(K.r1.1, 3, function(x){x %*% coef1}))
                            allcoef1 <- allcoef1[,c(4,1:3)]
                            X1 <- model.matrix(Surv(time,status)~-1+factor(trial)+age+factor(sex)+sbp,data = train1)
                            allvar1 <- apply(K.r1.1, 3, function(x){x %*% cov(X1) %*% t(x)}, simplify=F)
                            
                            r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                            coef1 = r1_meta1$coefficients
                            coef_mat1[,i] <- coef1
                            se_mat1[,i] <- diag(r1_meta1$vcov)
                            
                            #recalibration of event rates adapted to train
                            lp0 <- unlist(as.matrix(train0[1:4]) %*% coef0[1:4])
                            lp1 <- unlist(as.matrix(train1[1:4]) %*% coef1[1:4])
                            temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
                            temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
                            bh0 <- basehaz(temp0, centered=F)
                            bh1 <- basehaz(temp1, centered=F)
                            int0[i] <- bh0[nrow(bh0),]$hazard
                            int1[i] <- bh1[nrow(bh1),]$hazard
                            
                            #ite estimation
                            pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[1:4]) %*% coef0[1:4]))
                            pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[1:4]) %*% coef1[1:4]))
                            ite <- unlist(pred1 - pred0)
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, r1_meta0$vcov)
                            se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, r1_meta1$vcov)
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod1.0$residuals)) * se * sqrt(1 + diag(r1_meta0$Psi)),
                                                   upr = ite + qt(.975,length(mod1.0$residuals)) * se * sqrt(1 + diag(r1_meta0$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.r1[i] <- mse(unlist(ite_true[i]),ite)
                          }
                        })
                        if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.r1),
                          in_int = mean(in_int),
                          coef_mat.1 = mean(coef_mat0[1, ]),
                          coef_mat.2 = mean(coef_mat0[2, ]),
                          coef_mat.3 = mean(coef_mat0[3, ]),
                          coef_mat.4 = mean(coef_mat0[4, ]),
                          coef_mat.5 = mean(int1)-mean(int0),
                          coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                          coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                          coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                          se_mat.1 = mean(se_mat0[1, ]),
                          se_mat.2 = mean(se_mat0[2, ]),
                          se_mat.3 = mean(se_mat0[3, ]),
                          se_mat.4 = mean(se_mat0[4, ]),
                          se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                          se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                          se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ])),
                          seed = j,
                          success = TRUE)
                        
                      }
stopCluster(cl)
res.r1.tl.sim15 <- t(res.r1.tl15[1:14, ])
se.r1.tl.sim15 <- t(res.r1.tl15[15:21, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl15 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival","msm","mvmeta")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = nt)
                        df <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82)
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age - mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5)
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197)
                        sd_sbp <- c(13,16.5,9.4,12,9,10,21)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp - mean(df$sbp)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- 50
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- -0.3
                        b5 <- 0.015
                        b6 <- 0.1
                        b7 <- -0.008
                        
                        df$log_hr <- b1*df$age + b2*df$sex + b3*df$sbp + b4*df$treat +
                          (b5*df$age +b6*(df$sex-0.5) + b7*df$sbp)*df$treat
                        df$log_hr <- df$log_hr - mean(df$log_hr)
                        df$ite <- expit(b1*df$age + b2*df$sex + b3*df$sbp + b4 +
                                          (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)) -
                          expit(b1*df$age + b2*df$sex + b3*df$sbp)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0.5, 10)
                        c[c>5] <- 5 # censoring
                        df$time <- pmin(t,c)  # observed time is min of censored and true
                        df$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df[df$treat==0,]
                        df1 <- df[df$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        in_int <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 3, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 3, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- coxph(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(trial), data = train0)
                            mod1 <- coxph(Surv(time,status)~-1+factor(trial)+(age+factor(sex)+sbp)*factor(trial), data = train1)
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            
                            K.fs0 <- array(0, dim=c(4, length(coef0), 6))
                            ind <- c(rep(F,5),rep(T,3),rep(F,length(coef0)-8))
                            if(i==1){vec <- 3:7} else {vec <- 2:7}
                            for(p in 1:6){
                              if(p==1){K.fs0[1,p,p] <- 1} else {K.fs0[1,p-1,p] <- 1}
                              diag(K.fs0[2:4,ind,p]) <- 1
                              if(p == 1) next
                              K.fs0[2,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
                              K.fs0[3,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
                              K.fs0[4,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
                            }
                            allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
                            allvar0 <- apply(K.fs0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
                            
                            fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                            coef0 = fs_meta0$coefficients[c(-1)]
                            
                            coef_mat0[,i] <- coef0
                            se_mat0[,i] <- diag(fs_meta0$vcov[c(-1),c(-1)])
                            
                            K.fs1 <- array(0, dim=c(4, length(coef1), 6))
                            ind <- c(rep(F,5),rep(T,3),rep(F,length(coef1)-8))
                            if(i==1){vec <- 3:7} else {vec <- 2:7}
                            for(p in 1:6){
                              if(p==1){K.fs1[1,p,p] <- 1} else {K.fs1[1,p-1,p] <- 1}
                              diag(K.fs1[2:4,ind,p]) <- 1
                              if(p == 1) next
                              K.fs1[2,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":age")),p] <- 1
                              K.fs1[3,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1")),p] <- 1
                              K.fs1[4,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":sbp")),p] <- 1
                            }
                            allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
                            allvar1 <- apply(K.fs1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
                            
                            fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                            coef1 = fs_meta1$coefficients[c(-1)]
                            
                            coef_mat1[,i] <- coef1
                            se_mat1[,i] <- diag(fs_meta1$vcov[c(-1),c(-1)])
                            
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
                            
                            #ite prediction interval
                            se0 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef0, fs_meta0$vcov[c(-1),c(-1)])
                            se1 <- deltamethod(~(exp(x1+x2+x3)/(1+exp(x1+x2+x3))), coef1, fs_meta1$vcov[c(-1),c(-1)])
                            se <- max(se0,se1)
                            true_val <- unlist(ite_true[i])
                            pred_int <- data.frame(lwr = ite - qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + diag(fs_meta0$Psi)),
                                                   upr = ite + qt(.975,length(mod0$linear.predictor)) * se * sqrt(1 + diag(fs_meta0$Psi)))
                            
                            in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                            
                            #c-statistic for benefit
                            cstat = cstat4ben(outcome = as.numeric(test$status),
                                              treatment = test$treat == 1,
                                              score = 1-ite)
                            
                            c.ben[i] <- cstat[1]
                            c.ben.se[i] <- cstat[2]
                            
                            #calibration for benefit
                            lm_ <- lm(ite~unlist(ite_true[i]))
                            
                            a[i] <- lm_$coefficients[[1]]
                            b[i] <- 1+lm_$coefficients[[2]]
                            mse.fs[i] <- mse(unlist(ite_true[i]),ite)
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
                            in_int = mean(in_int),
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
res.fs.tl.sim15 <- t(res.fs.tl15[1:13, ])
se.fs.tl.sim15 <- t(res.fs.tl15[14:19, ])

save(res.na.sl.sim15,se.na.sl.sim15,
     res.re.sl.sim15,se.re.sl.sim15,
     res.si.sl.sim15,se.si.sl.sim15,
     res.r1.sl.sim15,se.r1.sl.sim15,
     res.fs.sl.sim15,se.fs.sl.sim15,
     res.na.tl.sim15,se.na.tl.sim15,
     res.re.tl.sim15,se.re.tl.sim15,
     res.si.tl.sim15,se.si.tl.sim15,
     res.r1.tl.sim15,se.r1.tl.sim15,
     res.fs.tl.sim15,se.fs.tl.sim15,
     file = "res_scenario15.Rdata")
