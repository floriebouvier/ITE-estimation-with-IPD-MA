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

#### scenario 22: 10 covariates & tte outcome & total sample size = 2800 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        set.seed(j)
                        n <- 2800
                        k <- 7
                        
                        trial <- rep(1:k, each = floor(n/k))
                        df10 <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df10$age <- df10$age-70
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df10$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df10$sbp <- df10$sbp-185
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                        df10$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df10$height <- df10$height-168
                        
                        mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
                        sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
                        df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
                        df10$chol <- df10$chol-6.3
                        
                        b0 <- 50-((mean_age-60)/100)
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- -0.3-((mean_age-60)/100)
                        b5 <- 0.82
                        b6 <- 0.8
                        b7 <- 0.7
                        b8 <- 0.1
                        b9 <- 0.33
                        b10 <- -0.02
                        b11 <- 0.06
                        b12 <- 0.01+((mean_age-70)/1000)
                        b13 <- 0.04+0.005*pman
                        b14 <- -0.2-((mean_sbp-170)/100)
                        b15 <- -0.01+0.003*psmoke
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        log_hr <- log_hr - mean(log_hr)
                        
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0, 5)
                        df10$time <- pmin(t,c)  #observed time is min of censored and true
                        df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat <- matrix(NA, nrow = 21, ncol = k)
                        se_mat <- matrix(NA, nrow = 21, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat),
                                         data = train)
                            coef <- mod$coefficients
                            coef[is.na(coef)] <- 0
                            coef_mat[,i] <- coef
                            se_mat[,i] <- sqrt(diag(vcov(mod)))
                            
                            #recalibration of event rates adapted to train
                            lp <- mod$linear.predictors
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
                            coef_mat.8 = mean(coef_mat[8, ]),
                            coef_mat.9 = mean(coef_mat[9, ]),
                            coef_mat.10 = mean(coef_mat[10, ]),
                            coef_mat.11 = mean(coef_mat[11, ]),
                            coef_mat.12 = mean(coef_mat[12, ]),
                            coef_mat.13 = mean(coef_mat[13, ]),
                            coef_mat.14 = mean(coef_mat[14, ]),
                            coef_mat.15 = mean(coef_mat[15, ]),
                            coef_mat.16 = mean(coef_mat[16, ]),
                            coef_mat.17 = mean(coef_mat[17, ]),
                            coef_mat.18 = mean(coef_mat[18, ]),
                            coef_mat.19 = mean(coef_mat[19, ]),
                            coef_mat.20 = mean(coef_mat[20, ]),
                            coef_mat.21 = mean(coef_mat[21, ]),
                            se_mat.1 = mean(se_mat[1, ]),
                            se_mat.2 = mean(se_mat[2, ]),
                            se_mat.3 = mean(se_mat[3, ]),
                            se_mat.4 = mean(se_mat[4, ]),
                            se_mat.5 = mean(se_mat[5, ]),
                            se_mat.6 = mean(se_mat[6, ]),
                            se_mat.7 = mean(se_mat[7, ]),
                            se_mat.8 = mean(se_mat[8, ]),
                            se_mat.9 = mean(se_mat[9, ]),
                            se_mat.10 = mean(se_mat[10, ]),
                            se_mat.11 = mean(se_mat[11, ]),
                            se_mat.12 = mean(se_mat[12, ]),
                            se_mat.13 = mean(se_mat[13, ]),
                            se_mat.14 = mean(se_mat[14, ]),
                            se_mat.15 = mean(se_mat[15, ]),
                            se_mat.16 = mean(se_mat[16, ]),
                            se_mat.17 = mean(se_mat[17, ]),
                            se_mat.18 = mean(se_mat[18, ]),
                            se_mat.19 = mean(se_mat[19, ]),
                            se_mat.20 = mean(se_mat[20, ]),
                            se_mat.21 = mean(se_mat[21, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.sl10.tte.hte.n2800 <- t(res.na.sl10[1:26, ])
se.na.sl10.tte.hte.n2800 <- t(res.na.sl10[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.re <- c()
    coef_mat <- matrix(NA, nrow = 21, ncol = k)
    se_mat <- matrix(NA, nrow = 21, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying the model to train
        mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+
                                          factor(diab)+lvhn1+height+chol)*factor(treat)+(1+treat|trial),
                     data = train)
        coef <- mod$coefficients
        coef_mat[,i] <- coef
        se_mat[,i] <- sqrt(diag(vcov(mod)))
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictor
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
        coef_mat.1 = mean(coef_mat[1, ]),
        coef_mat.2 = mean(coef_mat[2, ]),
        coef_mat.3 = mean(coef_mat[3, ]),
        coef_mat.4 = mean(coef_mat[4, ]),
        coef_mat.5 = mean(coef_mat[5, ]),
        coef_mat.6 = mean(coef_mat[6, ]),
        coef_mat.7 = mean(coef_mat[7, ]),
        coef_mat.8 = mean(coef_mat[8, ]),
        coef_mat.9 = mean(coef_mat[9, ]),
        coef_mat.10 = mean(coef_mat[10, ]),
        coef_mat.11 = mean(coef_mat[11, ]),
        coef_mat.12 = mean(coef_mat[12, ]),
        coef_mat.13 = mean(coef_mat[13, ]),
        coef_mat.14 = mean(coef_mat[14, ]),
        coef_mat.15 = mean(coef_mat[15, ]),
        coef_mat.16 = mean(coef_mat[16, ]),
        coef_mat.17 = mean(coef_mat[17, ]),
        coef_mat.18 = mean(coef_mat[18, ]),
        coef_mat.19 = mean(coef_mat[19, ]),
        coef_mat.20 = mean(coef_mat[20, ]),
        coef_mat.21 = mean(coef_mat[21, ]),
        se_mat.1 = mean(se_mat[1, ]),
        se_mat.2 = mean(se_mat[2, ]),
        se_mat.3 = mean(se_mat[3, ]),
        se_mat.4 = mean(se_mat[4, ]),
        se_mat.5 = mean(se_mat[5, ]),
        se_mat.6 = mean(se_mat[6, ]),
        se_mat.7 = mean(se_mat[7, ]),
        se_mat.8 = mean(se_mat[8, ]),
        se_mat.9 = mean(se_mat[9, ]),
        se_mat.10 = mean(se_mat[10, ]),
        se_mat.11 = mean(se_mat[11, ]),
        se_mat.12 = mean(se_mat[12, ]),
        se_mat.13 = mean(se_mat[13, ]),
        se_mat.14 = mean(se_mat[14, ]),
        se_mat.15 = mean(se_mat[15, ]),
        se_mat.16 = mean(se_mat[16, ]),
        se_mat.17 = mean(se_mat[17, ]),
        se_mat.18 = mean(se_mat[18, ]),
        se_mat.19 = mean(se_mat[19, ]),
        se_mat.20 = mean(se_mat[20, ]),
        se_mat.21 = mean(se_mat[21, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.re.sl10.tte.hte.n2800 <- t(res.re.sl10[1:26, ])
se.re.sl10.tte.hte.n2800 <- t(res.re.sl10[27:47 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.si <- c()
    coef_mat <- matrix(NA, nrow = 21, ncol = k)
    se_mat <- matrix(NA, nrow = 21, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+strata(trial),
                     data = train) #+(0+treat|trial)
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        se_mat[,i] <- sqrt(diag(vcov(mod)))
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictors
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
        coef_mat.8 = mean(coef_mat[8, ]),
        coef_mat.9 = mean(coef_mat[9, ]),
        coef_mat.10 = mean(coef_mat[10, ]),
        coef_mat.11 = mean(coef_mat[11, ]),
        coef_mat.12 = mean(coef_mat[12, ]),
        coef_mat.13 = mean(coef_mat[13, ]),
        coef_mat.14 = mean(coef_mat[14, ]),
        coef_mat.15 = mean(coef_mat[15, ]),
        coef_mat.16 = mean(coef_mat[16, ]),
        coef_mat.17 = mean(coef_mat[17, ]),
        coef_mat.18 = mean(coef_mat[18, ]),
        coef_mat.19 = mean(coef_mat[19, ]),
        coef_mat.20 = mean(coef_mat[20, ]),
        coef_mat.21 = mean(coef_mat[21, ]),
        se_mat.1 = mean(se_mat[1, ]),
        se_mat.2 = mean(se_mat[2, ]),
        se_mat.3 = mean(se_mat[3, ]),
        se_mat.4 = mean(se_mat[4, ]),
        se_mat.5 = mean(se_mat[5, ]),
        se_mat.6 = mean(se_mat[6, ]),
        se_mat.7 = mean(se_mat[7, ]),
        se_mat.8 = mean(se_mat[8, ]),
        se_mat.9 = mean(se_mat[9, ]),
        se_mat.10 = mean(se_mat[10, ]),
        se_mat.11 = mean(se_mat[11, ]),
        se_mat.12 = mean(se_mat[12, ]),
        se_mat.13 = mean(se_mat[13, ]),
        se_mat.14 = mean(se_mat[14, ]),
        se_mat.15 = mean(se_mat[15, ]),
        se_mat.16 = mean(se_mat[16, ]),
        se_mat.17 = mean(se_mat[17, ]),
        se_mat.18 = mean(se_mat[18, ]),
        se_mat.19 = mean(se_mat[19, ]),
        se_mat.20 = mean(se_mat[20, ]),
        se_mat.21 = mean(se_mat[21, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.si.sl10.tte.hte.n2800 <- t(res.si.sl10[1:26, ])
se.si.sl10.tte.hte.n2800 <- t(res.si.sl10[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        
    source("coxvc.R")
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.r1 <- c()
    coef_mat <- matrix(NA, nrow = 22, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        Ft<-cbind(rep(1,nrow(train)),bs(train$time))
        mod <- coxvc(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                          factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial,
                     Ft, rank = 1, data = train)
        coef <- mod$b
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        lp <- unlist(as.matrix(train[2:12]) %*% coef[1:11] + train[12] * (as.matrix(train[2:11]) %*% coef[13:22]))
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[13:22])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[13:22])))
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
      coef_mat.5 = mean(coef_mat[5, ]),
      coef_mat.6 = mean(coef_mat[6, ]),
      coef_mat.7 = mean(coef_mat[7, ]),
      coef_mat.8 = mean(coef_mat[8, ]),
      coef_mat.9 = mean(coef_mat[9, ]),
      coef_mat.10 = mean(coef_mat[10, ]),
      coef_mat.11 = mean(coef_mat[11, ]),
      coef_mat.12 = mean(coef_mat[12, ]),
      coef_mat.13 = mean(coef_mat[14, ]),
      coef_mat.14 = mean(coef_mat[15, ]),
      coef_mat.15 = mean(coef_mat[16, ]),
      coef_mat.16 = mean(coef_mat[17, ]),
      coef_mat.17 = mean(coef_mat[18, ]),
      coef_mat.18 = mean(coef_mat[19, ]),
      coef_mat.19 = mean(coef_mat[20, ]),
      coef_mat.20 = mean(coef_mat[21, ]),
      coef_mat.21 = mean(coef_mat[22, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.r1.sl10.tte.hte.n2800 <- t(res.r1.sl10[1:26, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat <- matrix(NA, nrow = 76, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+
                       (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol):factor(trial),
                     data = train)
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictors
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test0$trial+
                                                 mean(coef[c(6,27:31)])*test0$age+mean(coef[c(7,32:36)])*test0$sex+
                                                 mean(coef[c(8,37:41)])*test0$sbp+mean(coef[c(9,42:46)])*test0$mi+
                                                 mean(coef[c(10,47:51)])*test0$stroke+mean(coef[c(11,52:56)])*test0$smoke+
                                                 mean(coef[c(12,57:61)])*test0$diab+mean(coef[c(13,62:66)])*test0$lvhn1+
                                                 mean(coef[c(14,67:71)])*test0$height+mean(coef[c(15,72:76)])*test0$chol+
                                                 coef[16]*test0$treat+(coef[17]*test0$age+coef[18]*test0$sex+
                                                                         coef[19]*test0$sbp+coef[20]*test0$mi+coef[21]*test0$stroke+
                                                                         coef[22]*test0$smoke+coef[23]*test0$diab+coef[24]*test0$lvhn1+
                                                                         coef[25]*test0$height+coef[26]*test0$chol)*test0$treat))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test1$trial+
                                                 mean(coef[c(6,27:31)])*test1$age+mean(coef[c(7,32:36)])*test1$sex+
                                                 mean(coef[c(8,37:41)])*test1$sbp+mean(coef[c(9,42:46)])*test1$mi+
                                                 mean(coef[c(10,47:51)])*test1$stroke+mean(coef[c(11,52:56)])*test1$smoke+
                                                 mean(coef[c(12,57:61)])*test1$diab+mean(coef[c(13,62:66)])*test1$lvhn1+
                                                 mean(coef[c(14,67:71)])*test1$height+mean(coef[c(15,72:76)])*test1$chol+
                                                 coef[16]*test1$treat+(coef[17]*test1$age+coef[18]*test1$sex+
                                                                         coef[19]*test1$sbp+coef[20]*test1$mi+coef[21]*test1$stroke+
                                                                         coef[22]*test1$smoke+coef[23]*test1$diab+coef[24]*test1$lvhn1+
                                                                         coef[25]*test1$height+coef[26]*test1$chol)*test1$treat))
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
    
    list(res = c(
      c.ben = mean(c.ben),
      c.ben.se = mean(c.ben.se),
      a = mean(a),
      b = mean(b),
      mse = mean(mse.fs),
      coef_mat.1 = mean(coef_mat[6, ])+coef_mat[27:31, ]),
      coef_mat.2 = mean(c(coef_mat[7, ]+coef_mat[32:36, ])),
      coef_mat.3 = mean(c(coef_mat[8, ]+coef_mat[37:41, ])),
      coef_mat.4 = mean(c(coef_mat[9, ]+coef_mat[42:46, ])),
      coef_mat.5 = mean(c(coef_mat[10, ]+coef_mat[47:51, ])),
      coef_mat.6 = mean(c(coef_mat[11, ]+coef_mat[52:56, ])),
      coef_mat.7 = mean(c(coef_mat[12, ]+coef_mat[57:61, ])),
      coef_mat.8 = mean(c(coef_mat[13, ]+coef_mat[62:66, ])),
      coef_mat.9 = mean(c(coef_mat[14, ]+coef_mat[67:71, ])),
      coef_mat.10 = mean(c(coef_mat[15, ]+coef_mat[72:76, ])),
      coef_mat.11 = mean(coef_mat[16, ]),
      coef_mat.12 = mean(coef_mat[17, ]),
      coef_mat.13 = mean(coef_mat[18, ]),
      coef_mat.14 = mean(coef_mat[19, ]),
      coef_mat.15 = mean(coef_mat[20, ]),
      coef_mat.16 = mean(coef_mat[21, ]),
      coef_mat.17 = mean(coef_mat[22, ]),
      coef_mat.18 = mean(coef_mat[23, ]),
      coef_mat.19 = mean(coef_mat[24, ]),
      coef_mat.20 = mean(coef_mat[25, ]),
      coef_mat.21 = mean(coef_mat[26, ]),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.fs.sl10.tte.hte.n2800 <- t(res.fs.sl10[1:26, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        set.seed(j)
                        n <- 2800
                        k <- 7
                        
                        trial <- rep(1:k, each = floor(n/k))
                        df10 <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df10$age <- df10$age-70
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df10$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df10$sbp <- df10$sbp-185
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                        df10$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df10$height <- df10$height-168
                        
                        mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
                        sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
                        df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
                        df10$chol <- df10$chol-6.3
                        
                        b0 <- 50-((mean_age-60)/100)
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- -0.3-((mean_age-60)/100)
                        b5 <- 0.82
                        b6 <- 0.8
                        b7 <- 0.7
                        b8 <- 0.1
                        b9 <- 0.33
                        b10 <- -0.02
                        b11 <- 0.06
                        b12 <- 0.01+((mean_age-70)/1000)
                        b13 <- 0.04+0.005*pman
                        b14 <- -0.2-((mean_sbp-170)/100)
                        b15 <- -0.01+0.003*psmoke
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        log_hr <- log_hr - mean(log_hr)
                        
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0, 5)
                        df10$time <- pmin(t,c)  #observed time is min of censored and true
                        df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df10[df10$treat==0,]
                        df1 <- df10[df10$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                            factor(smoke)+factor(diab)+lvhn1+height+chol,
                                          data = train0)
                            mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                            factor(smoke)+factor(diab)+lvhn1+height+chol,
                                          data = train1)
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            coef0[is.na(coef0)] <- 0
                            coef1[is.na(coef1)] <- 0
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
                            pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
                            pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
                            coef_mat.1 = mean(coef_mat0[1, ]),
                            coef_mat.2 = mean(coef_mat0[2, ]),
                            coef_mat.3 = mean(coef_mat0[3, ]),
                            coef_mat.4 = mean(coef_mat0[4, ]),
                            coef_mat.5 = mean(coef_mat0[5, ]),
                            coef_mat.6 = mean(coef_mat0[6, ]),
                            coef_mat.7 = mean(coef_mat0[7, ]),
                            coef_mat.8 = mean(coef_mat0[8, ]),
                            coef_mat.9 = mean(coef_mat0[9, ]),
                            coef_mat.10 = mean(coef_mat0[10, ]),
                            coef_mat.11 = mean(int1)-mean(int0),
                            coef_mat.12 = mean(coef_mat1[1, ])- mean(coef_mat0[1, ]),
                            coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                            coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                            coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                            coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                            coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                            coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                            coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                            coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                            coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
                            se_mat.1 = mean(se_mat0[1, ]),
                            se_mat.2 = mean(se_mat0[2, ]),
                            se_mat.3 = mean(se_mat0[3, ]),
                            se_mat.4 = mean(se_mat0[4, ]),
                            se_mat.5 = mean(se_mat0[5, ]),
                            se_mat.6 = mean(se_mat0[6, ]),
                            se_mat.7 = mean(se_mat0[7, ]),
                            se_mat.8 = mean(se_mat0[8, ]),
                            se_mat.9 = mean(se_mat0[9, ]),
                            se_mat.10 = mean(se_mat0[10, ]),
                            se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                            se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                            se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                            se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                            se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                            se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                            se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                            se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                            se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                            se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.tl10.tte.hte.n2800 <- t(res.na.tl10[1:26, ])
se.na.tl10.tte.hte.n2800 <- t(res.na.tl10[27:46, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.re <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxme(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),data = train0)
        mod1 <- coxme(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),data = train1)
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        coef_mat.1 = mean(coef_mat0[1, ]),
        coef_mat.2 = mean(coef_mat0[2, ]),
        coef_mat.3 = mean(coef_mat0[3, ]),
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
        coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
        se_mat.1 = mean(se_mat0[1, ]),
        se_mat.2 = mean(se_mat0[2, ]),
        se_mat.3 = mean(se_mat0[3, ]),
        se_mat.4 = mean(se_mat0[4, ]),
        se_mat.5 = mean(se_mat0[5, ]),
        se_mat.6 = mean(se_mat0[6, ]),
        se_mat.7 = mean(se_mat0[7, ]),
        se_mat.8 = mean(se_mat0[8, ]),
        se_mat.9 = mean(se_mat0[9, ]),
        se_mat.10 = mean(se_mat0[10, ]),
        se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.re.tl10.tte.hte.n2800 <- t(res.re.tl10[1:26, ])
se.re.tl10.tte.hte.n2800 <- t(res.re.tl10[27:46 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.si <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+strata(trial),
                      data = train0)
        mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+strata(trial),
                      data = train1)
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        coef_mat.1 = mean(coef_mat0[1, ]),
        coef_mat.2 = mean(coef_mat0[2, ]),
        coef_mat.3 = mean(coef_mat0[3, ]),
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = mean(coef_mat1[1, ])- mean(coef_mat0[1, ]),
        coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
        se_mat.1 = mean(se_mat0[1, ]),
        se_mat.2 = mean(se_mat0[2, ]),
        se_mat.3 = mean(se_mat0[3, ]),
        se_mat.4 = mean(se_mat0[4, ]),
        se_mat.5 = mean(se_mat0[5, ]),
        se_mat.6 = mean(se_mat0[6, ]),
        se_mat.7 = mean(se_mat0[7, ]),
        se_mat.8 = mean(se_mat0[8, ]),
        se_mat.9 = mean(se_mat0[9, ]),
        se_mat.10 = mean(se_mat0[10, ]),
        se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.si.tl10.tte.hte.n2800 <- t(res.si.tl10[1:26, ])
se.si.tl10.tte.hte.n2800 <- t(res.si.tl10[27:46, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    source("coxvc.R")
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.r1 <- c()
    coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        Ft0<-cbind(rep(1,nrow(train0)),bs(train0$time))
        mod0 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+trial, Ft0, rank = 1,
                      data = train0)
        Ft1<-cbind(rep(1,nrow(train0)),bs(train1$time))
        mod1 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+trial, Ft1, rank = 1,
                      data = train1)
        
        coef0 <- mod0$b
        coef1 <- mod1$b
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
        #recalibration of event rates adapted to train
        lp0 <- unlist(as.matrix(train0[2:11]) %*% coef[1:10])
        lp1 <- unlist(as.matrix(train1[2:11]) %*% coef[1:10])
        temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
        temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
        bh0 <- basehaz(temp0, centered=F)
        bh1 <- basehaz(temp1, centered=F)
        int0[i] <- bh0[nrow(bh0),]$hazard
        int1[i] <- bh1[nrow(bh1),]$hazard
        
        #ite estimation
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        oben <- unlist(ifelse(sapply(oben, length) == 0, 0, oben))
        lm_ <- lm(oben~pben)
        
        a[i] <- lm_$coefficients[[1]]
        b[i] <- lm_$coefficients[[2]]
        mse.r1[i] <- mse(oben,pben)
      }
    })
    if(any(class(testerror) == "try-error")) return(list(success = FALSE))
    list(
      res = c(c.ben = mean(c.ben),
              c.ben.se = mean(c.ben.se),
              a = mean(a),
              b = mean(b),
              mse = mean(mse.r1),
              coef_mat.1 = mean(coef_mat0[1, ]),
              coef_mat.2 = mean(coef_mat0[2, ]),
              coef_mat.3 = mean(coef_mat0[3, ]),
              coef_mat.4 = mean(coef_mat0[4, ]),
              coef_mat.5 = mean(coef_mat0[5, ]),
              coef_mat.6 = mean(coef_mat0[6, ]),
              coef_mat.7 = mean(coef_mat0[7, ]),
              coef_mat.8 = mean(coef_mat0[8, ]),
              coef_mat.9 = mean(coef_mat0[9, ]),
              coef_mat.10 = mean(coef_mat0[10, ]),
              coef_mat.11 = mean(int1)-mean(int0),
              coef_mat.12 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
              coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
              coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
              coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
              coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
              coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
              coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
              coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
              coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
              coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.r1.tl10.tte.hte.n2800 <- t(res.r1.tl10[1:26, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl10 = foreach(
  j = 1:m,
  .final = function(l) do.call("cbind",
                               lapply(
                                 l[sapply(l, `[[`, "success")],
                                 function(subl) c(subl$res, subl$seed))),
  .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar%{
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 65, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 65, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                           factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(trial),
                      data = train0)
        mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                           factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(trial),
                      data = train1)
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(mean(coef0[11:15])*test$trial+
                                                 mean(coef0[c(1,16:20)])*test$age+mean(coef0[c(2,21:25)])*test$sex+
                                                 mean(coef0[c(3,26:30)])*test$sbp+mean(coef0[c(4,31:35)])*test$mi+
                                                 mean(coef0[c(5,36:40)])*test$stroke+mean(coef0[c(6,41:45)])*test$smoke+
                                                 mean(coef0[c(7,46:50)])*test$diab+mean(coef0[c(8,51:55)])*test$lvhn1+
                                                 mean(coef0[c(9,56:60)])*test$height+mean(coef0[c(10,61:65)])*test$chol))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(mean(coef1[11:15])*test$trial+
                                                 mean(coef1[c(1,16:20)])*test$age+mean(coef1[c(2,21:25)])*test$sex+
                                                 mean(coef1[c(3,26:30)])*test$sbp+mean(coef1[c(4,31:35)])*test$mi+
                                                 mean(coef1[c(5,36:40)])*test$stroke+mean(coef1[c(6,41:45)])*test$smoke+
                                                 mean(coef1[c(7,46:50)])*test$diab+mean(coef1[c(8,51:55)])*test$lvhn1+
                                                 mean(coef1[c(9,56:60)])*test$height+mean(coef1[c(10,61:65)])*test$chol))
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
    if (any(class(testerror) == "try-error")) return(list(success = FALSE))
    
    list(
      res = c(
        c.ben = mean(c.ben),
        c.ben.se = mean(c.ben.se),
        a = mean(a),
        b = mean(b),
        mse = mean(mse.fs),
        coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[16:20, ]),
        coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[21:25, ]),
        coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[26:30, ]),
        coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[31:35, ]),
        coef_mat.5 = mean(coef_mat0[5, ]) + mean(coef_mat0[36:40, ]),
        coef_mat.6 = mean(coef_mat0[6, ]) + mean(coef_mat0[41:45, ]),
        coef_mat.7 = mean(coef_mat0[7, ]) + mean(coef_mat0[46:50, ]),
        coef_mat.8 = mean(coef_mat0[8, ]) + mean(coef_mat0[51:55, ]),
        coef_mat.9 = mean(coef_mat0[9, ]) + mean(coef_mat0[56:60, ]),
        coef_mat.10 = mean(coef_mat0[10, ]) + mean(coef_mat0[61:65, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = (mean(coef_mat1[1]) + mean(coef_mat1[16:20, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[16:20, ])),
        coef_mat.13 = (mean(coef_mat1[2]) + mean(coef_mat1[17:21, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[21:25, ])),
        coef_mat.14 = (mean(coef_mat1[3]) + mean(coef_mat1[22:26, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[26:30, ])),
        coef_mat.15 = (mean(coef_mat1[4]) + mean(coef_mat1[27:31, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[31:35, ])),
        coef_mat.16 = (mean(coef_mat1[5]) + mean(coef_mat1[32:36, ])) - (mean(coef_mat0[5, ]) + mean(coef_mat0[36:40, ])),
        coef_mat.17 = (mean(coef_mat1[6]) + mean(coef_mat1[37:41, ])) - (mean(coef_mat0[6, ]) + mean(coef_mat0[41:45, ])),
        coef_mat.18 = (mean(coef_mat1[7]) + mean(coef_mat1[42:46, ])) - (mean(coef_mat0[7, ]) + mean(coef_mat0[46:50, ])),
        coef_mat.19 = (mean(coef_mat1[8]) + mean(coef_mat1[47:51, ])) - (mean(coef_mat0[8, ]) + mean(coef_mat0[51:55, ])),
        coef_mat.20 = (mean(coef_mat1[9]) + mean(coef_mat1[52:56, ])) - (mean(coef_mat0[9, ]) + mean(coef_mat0[56:60, ])),
        coef_mat.21 = (mean(coef_mat1[10]) + mean(coef_mat1[57:61, ])) - (mean(coef_mat0[10, ]) + mean(coef_mat0[61:65, ]))),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl10.tte.hte.n2800 <- t(res.fs.tl10[1:26, ])

save(res.na.sl10.tte.hte.n2800,se.na.sl10.tte.hte.n2800,
     res.re.sl10.tte.hte.n2800,se.re.sl10.tte.hte.n2800,
     res.si.sl10.tte.hte.n2800,se.si.sl10.tte.hte.n2800,
     res.fs.sl10.tte.hte.n2800,
     res.na.tl10.tte.hte.n2800,se.na.tl10.tte.hte.n2800,
     res.re.tl10.tte.hte.n2800,se.re.tl10.tte.hte.n2800,
     res.si.tl10.tte.hte.n2800,se.si.tl10.tte.hte.n2800,
     res.fs.tl10.tte.hte.n2800,
     file = "res_scenario22.Rdata") #res.r1.sl10.tte.hte.n2800,res.r1.tl10.tte.hte.n2800

#### scenario 23: 10 covariates & tte outcome & total sample size = 1400 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = floor(n/k))
                        df10 <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df10$age <- df10$age-70
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df10$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df10$sbp <- df10$sbp-185
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                        df10$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df10$height <- df10$height-168
                        
                        mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
                        sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
                        df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
                        df10$chol <- df10$chol-6.3
                        
                        b0 <- 50-((mean_age-60)/100)
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- -0.3-((mean_age-60)/100)
                        b5 <- 0.82
                        b6 <- 0.8
                        b7 <- 0.7
                        b8 <- 0.1
                        b9 <- 0.33
                        b10 <- -0.02
                        b11 <- 0.06
                        b12 <- 0.01+((mean_age-70)/1000)
                        b13 <- 0.04+0.005*pman
                        b14 <- -0.2-((mean_sbp-170)/100)
                        b15 <- -0.01+0.003*psmoke
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        log_hr <- log_hr - mean(log_hr)
                        
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0, 5)
                        df10$time <- pmin(t,c)  #observed time is min of censored and true
                        df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat <- matrix(NA, nrow = 21, ncol = k)
                        se_mat <- matrix(NA, nrow = 21, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat),
                                         data = train)
                            coef <- mod$coefficients
                            coef[is.na(coef)] <- 0
                            coef_mat[,i] <- coef
                            se_mat[,i] <- sqrt(diag(vcov(mod)))
                            
                            #recalibration of event rates adapted to train
                            lp <- mod$linear.predictors
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
                            coef_mat.8 = mean(coef_mat[8, ]),
                            coef_mat.9 = mean(coef_mat[9, ]),
                            coef_mat.10 = mean(coef_mat[10, ]),
                            coef_mat.11 = mean(coef_mat[11, ]),
                            coef_mat.12 = mean(coef_mat[12, ]),
                            coef_mat.13 = mean(coef_mat[13, ]),
                            coef_mat.14 = mean(coef_mat[14, ]),
                            coef_mat.15 = mean(coef_mat[15, ]),
                            coef_mat.16 = mean(coef_mat[16, ]),
                            coef_mat.17 = mean(coef_mat[17, ]),
                            coef_mat.18 = mean(coef_mat[18, ]),
                            coef_mat.19 = mean(coef_mat[19, ]),
                            coef_mat.20 = mean(coef_mat[20, ]),
                            coef_mat.21 = mean(coef_mat[21, ]),
                            se_mat.1 = mean(se_mat[1, ]),
                            se_mat.2 = mean(se_mat[2, ]),
                            se_mat.3 = mean(se_mat[3, ]),
                            se_mat.4 = mean(se_mat[4, ]),
                            se_mat.5 = mean(se_mat[5, ]),
                            se_mat.6 = mean(se_mat[6, ]),
                            se_mat.7 = mean(se_mat[7, ]),
                            se_mat.8 = mean(se_mat[8, ]),
                            se_mat.9 = mean(se_mat[9, ]),
                            se_mat.10 = mean(se_mat[10, ]),
                            se_mat.11 = mean(se_mat[11, ]),
                            se_mat.12 = mean(se_mat[12, ]),
                            se_mat.13 = mean(se_mat[13, ]),
                            se_mat.14 = mean(se_mat[14, ]),
                            se_mat.15 = mean(se_mat[15, ]),
                            se_mat.16 = mean(se_mat[16, ]),
                            se_mat.17 = mean(se_mat[17, ]),
                            se_mat.18 = mean(se_mat[18, ]),
                            se_mat.19 = mean(se_mat[19, ]),
                            se_mat.20 = mean(se_mat[20, ]),
                            se_mat.21 = mean(se_mat[21, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.sl10.tte.hte.n1400 <- t(res.na.sl10[1:26, ])
se.na.sl10.tte.hte.n1400 <- t(res.na.sl10[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.re <- c()
    coef_mat <- matrix(NA, nrow = 21, ncol = k)
    se_mat <- matrix(NA, nrow = 21, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying the model to train
        mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+
                                          factor(diab)+lvhn1+height+chol)*factor(treat)+(1+treat|trial),
                     data = train)
        coef <- mod$coefficients
        coef_mat[,i] <- coef
        se_mat[,i] <- sqrt(diag(vcov(mod)))
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictor
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
        coef_mat.1 = mean(coef_mat[1, ]),
        coef_mat.2 = mean(coef_mat[2, ]),
        coef_mat.3 = mean(coef_mat[3, ]),
        coef_mat.4 = mean(coef_mat[4, ]),
        coef_mat.5 = mean(coef_mat[5, ]),
        coef_mat.6 = mean(coef_mat[6, ]),
        coef_mat.7 = mean(coef_mat[7, ]),
        coef_mat.8 = mean(coef_mat[8, ]),
        coef_mat.9 = mean(coef_mat[9, ]),
        coef_mat.10 = mean(coef_mat[10, ]),
        coef_mat.11 = mean(coef_mat[11, ]),
        coef_mat.12 = mean(coef_mat[12, ]),
        coef_mat.13 = mean(coef_mat[13, ]),
        coef_mat.14 = mean(coef_mat[14, ]),
        coef_mat.15 = mean(coef_mat[15, ]),
        coef_mat.16 = mean(coef_mat[16, ]),
        coef_mat.17 = mean(coef_mat[17, ]),
        coef_mat.18 = mean(coef_mat[18, ]),
        coef_mat.19 = mean(coef_mat[19, ]),
        coef_mat.20 = mean(coef_mat[20, ]),
        coef_mat.21 = mean(coef_mat[21, ]),
        se_mat.1 = mean(se_mat[1, ]),
        se_mat.2 = mean(se_mat[2, ]),
        se_mat.3 = mean(se_mat[3, ]),
        se_mat.4 = mean(se_mat[4, ]),
        se_mat.5 = mean(se_mat[5, ]),
        se_mat.6 = mean(se_mat[6, ]),
        se_mat.7 = mean(se_mat[7, ]),
        se_mat.8 = mean(se_mat[8, ]),
        se_mat.9 = mean(se_mat[9, ]),
        se_mat.10 = mean(se_mat[10, ]),
        se_mat.11 = mean(se_mat[11, ]),
        se_mat.12 = mean(se_mat[12, ]),
        se_mat.13 = mean(se_mat[13, ]),
        se_mat.14 = mean(se_mat[14, ]),
        se_mat.15 = mean(se_mat[15, ]),
        se_mat.16 = mean(se_mat[16, ]),
        se_mat.17 = mean(se_mat[17, ]),
        se_mat.18 = mean(se_mat[18, ]),
        se_mat.19 = mean(se_mat[19, ]),
        se_mat.20 = mean(se_mat[20, ]),
        se_mat.21 = mean(se_mat[21, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.re.sl10.tte.hte.n1400 <- t(res.re.sl10[1:26, ])
se.re.sl10.tte.hte.n1400 <- t(res.re.sl10[27:47 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.si <- c()
    coef_mat <- matrix(NA, nrow = 21, ncol = k)
    se_mat <- matrix(NA, nrow = 21, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+strata(trial),
                     data = train) #+(0+treat|trial)
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        se_mat[,i] <- sqrt(diag(vcov(mod)))
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictors
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
        coef_mat.8 = mean(coef_mat[8, ]),
        coef_mat.9 = mean(coef_mat[9, ]),
        coef_mat.10 = mean(coef_mat[10, ]),
        coef_mat.11 = mean(coef_mat[11, ]),
        coef_mat.12 = mean(coef_mat[12, ]),
        coef_mat.13 = mean(coef_mat[13, ]),
        coef_mat.14 = mean(coef_mat[14, ]),
        coef_mat.15 = mean(coef_mat[15, ]),
        coef_mat.16 = mean(coef_mat[16, ]),
        coef_mat.17 = mean(coef_mat[17, ]),
        coef_mat.18 = mean(coef_mat[18, ]),
        coef_mat.19 = mean(coef_mat[19, ]),
        coef_mat.20 = mean(coef_mat[20, ]),
        coef_mat.21 = mean(coef_mat[21, ]),
        se_mat.1 = mean(se_mat[1, ]),
        se_mat.2 = mean(se_mat[2, ]),
        se_mat.3 = mean(se_mat[3, ]),
        se_mat.4 = mean(se_mat[4, ]),
        se_mat.5 = mean(se_mat[5, ]),
        se_mat.6 = mean(se_mat[6, ]),
        se_mat.7 = mean(se_mat[7, ]),
        se_mat.8 = mean(se_mat[8, ]),
        se_mat.9 = mean(se_mat[9, ]),
        se_mat.10 = mean(se_mat[10, ]),
        se_mat.11 = mean(se_mat[11, ]),
        se_mat.12 = mean(se_mat[12, ]),
        se_mat.13 = mean(se_mat[13, ]),
        se_mat.14 = mean(se_mat[14, ]),
        se_mat.15 = mean(se_mat[15, ]),
        se_mat.16 = mean(se_mat[16, ]),
        se_mat.17 = mean(se_mat[17, ]),
        se_mat.18 = mean(se_mat[18, ]),
        se_mat.19 = mean(se_mat[19, ]),
        se_mat.20 = mean(se_mat[20, ]),
        se_mat.21 = mean(se_mat[21, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.si.sl10.tte.hte.n1400 <- t(res.si.sl10[1:26, ])
se.si.sl10.tte.hte.n1400 <- t(res.si.sl10[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    
    source("coxvc.R")
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.r1 <- c()
    coef_mat <- matrix(NA, nrow = 22, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        Ft<-cbind(rep(1,nrow(train)),bs(train$time))
        mod <- coxvc(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                          factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial,
                     Ft, rank = 1, data = train)
        coef <- mod$b
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        lp <- unlist(as.matrix(train[2:12]) %*% coef[1:11] + train[12] * (as.matrix(train[2:11]) %*% coef[13:22]))
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[13:22])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[13:22])))
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
      coef_mat.5 = mean(coef_mat[5, ]),
      coef_mat.6 = mean(coef_mat[6, ]),
      coef_mat.7 = mean(coef_mat[7, ]),
      coef_mat.8 = mean(coef_mat[8, ]),
      coef_mat.9 = mean(coef_mat[9, ]),
      coef_mat.10 = mean(coef_mat[10, ]),
      coef_mat.11 = mean(coef_mat[11, ]),
      coef_mat.12 = mean(coef_mat[12, ]),
      coef_mat.13 = mean(coef_mat[14, ]),
      coef_mat.14 = mean(coef_mat[15, ]),
      coef_mat.15 = mean(coef_mat[16, ]),
      coef_mat.16 = mean(coef_mat[17, ]),
      coef_mat.17 = mean(coef_mat[18, ]),
      coef_mat.18 = mean(coef_mat[19, ]),
      coef_mat.19 = mean(coef_mat[20, ]),
      coef_mat.20 = mean(coef_mat[21, ]),
      coef_mat.21 = mean(coef_mat[22, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.r1.sl10.tte.hte.n1400 <- t(res.r1.sl10[1:26, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat <- matrix(NA, nrow = 76, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+
                       (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol):factor(trial),
                     data = train)
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictors
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test0$trial+
                                                 mean(coef[c(6,27:31)])*test0$age+mean(coef[c(7,32:36)])*test0$sex+
                                                 mean(coef[c(8,37:41)])*test0$sbp+mean(coef[c(9,42:46)])*test0$mi+
                                                 mean(coef[c(10,47:51)])*test0$stroke+mean(coef[c(11,52:56)])*test0$smoke+
                                                 mean(coef[c(12,57:61)])*test0$diab+mean(coef[c(13,62:66)])*test0$lvhn1+
                                                 mean(coef[c(14,67:71)])*test0$height+mean(coef[c(15,72:76)])*test0$chol+
                                                 coef[16]*test0$treat+(coef[17]*test0$age+coef[18]*test0$sex+
                                                                         coef[19]*test0$sbp+coef[20]*test0$mi+coef[21]*test0$stroke+
                                                                         coef[22]*test0$smoke+coef[23]*test0$diab+coef[24]*test0$lvhn1+
                                                                         coef[25]*test0$height+coef[26]*test0$chol)*test0$treat))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test1$trial+
                                                 mean(coef[c(6,27:31)])*test1$age+mean(coef[c(7,32:36)])*test1$sex+
                                                 mean(coef[c(8,37:41)])*test1$sbp+mean(coef[c(9,42:46)])*test1$mi+
                                                 mean(coef[c(10,47:51)])*test1$stroke+mean(coef[c(11,52:56)])*test1$smoke+
                                                 mean(coef[c(12,57:61)])*test1$diab+mean(coef[c(13,62:66)])*test1$lvhn1+
                                                 mean(coef[c(14,67:71)])*test1$height+mean(coef[c(15,72:76)])*test1$chol+
                                                 coef[16]*test1$treat+(coef[17]*test1$age+coef[18]*test1$sex+
                                                                         coef[19]*test1$sbp+coef[20]*test1$mi+coef[21]*test1$stroke+
                                                                         coef[22]*test1$smoke+coef[23]*test1$diab+coef[24]*test1$lvhn1+
                                                                         coef[25]*test1$height+coef[26]*test1$chol)*test1$treat))
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
    
    list(res = c(
      c.ben = mean(c.ben),
      c.ben.se = mean(c.ben.se),
      a = mean(a),
      b = mean(b),
      mse = mean(mse.fs),
      coef_mat.1 = mean(coef_mat[6, ])+coef_mat[27:31, ]),
      coef_mat.2 = mean(c(coef_mat[7, ]+coef_mat[32:36, ])),
      coef_mat.3 = mean(c(coef_mat[8, ]+coef_mat[37:41, ])),
      coef_mat.4 = mean(c(coef_mat[9, ]+coef_mat[42:46, ])),
      coef_mat.5 = mean(c(coef_mat[10, ]+coef_mat[47:51, ])),
      coef_mat.6 = mean(c(coef_mat[11, ]+coef_mat[52:56, ])),
      coef_mat.7 = mean(c(coef_mat[12, ]+coef_mat[57:61, ])),
      coef_mat.8 = mean(c(coef_mat[13, ]+coef_mat[62:66, ])),
      coef_mat.9 = mean(c(coef_mat[14, ]+coef_mat[67:71, ])),
      coef_mat.10 = mean(c(coef_mat[15, ]+coef_mat[72:76, ])),
      coef_mat.11 = mean(coef_mat[16, ]),
      coef_mat.12 = mean(coef_mat[17, ]),
      coef_mat.13 = mean(coef_mat[18, ]),
      coef_mat.14 = mean(coef_mat[19, ]),
      coef_mat.15 = mean(coef_mat[20, ]),
      coef_mat.16 = mean(coef_mat[21, ]),
      coef_mat.17 = mean(coef_mat[22, ]),
      coef_mat.18 = mean(coef_mat[23, ]),
      coef_mat.19 = mean(coef_mat[24, ]),
      coef_mat.20 = mean(coef_mat[25, ]),
      coef_mat.21 = mean(coef_mat[26, ]),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.fs.sl10.tte.hte.n1400 <- t(res.fs.sl10[1:26, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        set.seed(j)
                        n <- 1400
                        k <- 7
                        
                        trial <- rep(1:k, each = floor(n/k))
                        df10 <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df10$age <- df10$age-70
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df10$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df10$sbp <- df10$sbp-185
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                        df10$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df10$height <- df10$height-168
                        
                        mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
                        sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
                        df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
                        df10$chol <- df10$chol-6.3
                        
                        b0 <- 50-((mean_age-60)/100)
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- -0.3-((mean_age-60)/100)
                        b5 <- 0.82
                        b6 <- 0.8
                        b7 <- 0.7
                        b8 <- 0.1
                        b9 <- 0.33
                        b10 <- -0.02
                        b11 <- 0.06
                        b12 <- 0.01+((mean_age-70)/1000)
                        b13 <- 0.04+0.005*pman
                        b14 <- -0.2-((mean_sbp-170)/100)
                        b15 <- -0.01+0.003*psmoke
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        log_hr <- log_hr - mean(log_hr)
                        
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0, 5)
                        df10$time <- pmin(t,c)  #observed time is min of censored and true
                        df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df10[df10$treat==0,]
                        df1 <- df10[df10$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                            factor(smoke)+factor(diab)+lvhn1+height+chol,
                                          data = train0)
                            mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                            factor(smoke)+factor(diab)+lvhn1+height+chol,
                                          data = train1)
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            coef0[is.na(coef0)] <- 0
                            coef1[is.na(coef1)] <- 0
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
                            pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
                            pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
                            coef_mat.1 = mean(coef_mat0[1, ]),
                            coef_mat.2 = mean(coef_mat0[2, ]),
                            coef_mat.3 = mean(coef_mat0[3, ]),
                            coef_mat.4 = mean(coef_mat0[4, ]),
                            coef_mat.5 = mean(coef_mat0[5, ]),
                            coef_mat.6 = mean(coef_mat0[6, ]),
                            coef_mat.7 = mean(coef_mat0[7, ]),
                            coef_mat.8 = mean(coef_mat0[8, ]),
                            coef_mat.9 = mean(coef_mat0[9, ]),
                            coef_mat.10 = mean(coef_mat0[10, ]),
                            coef_mat.11 = mean(int1)-mean(int0),
                            coef_mat.12 = mean(coef_mat1[1, ])- mean(coef_mat0[1, ]),
                            coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                            coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                            coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                            coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                            coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                            coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                            coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                            coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                            coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
                            se_mat.1 = mean(se_mat0[1, ]),
                            se_mat.2 = mean(se_mat0[2, ]),
                            se_mat.3 = mean(se_mat0[3, ]),
                            se_mat.4 = mean(se_mat0[4, ]),
                            se_mat.5 = mean(se_mat0[5, ]),
                            se_mat.6 = mean(se_mat0[6, ]),
                            se_mat.7 = mean(se_mat0[7, ]),
                            se_mat.8 = mean(se_mat0[8, ]),
                            se_mat.9 = mean(se_mat0[9, ]),
                            se_mat.10 = mean(se_mat0[10, ]),
                            se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                            se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                            se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                            se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                            se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                            se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                            se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                            se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                            se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                            se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.tl10.tte.hte.n1400 <- t(res.na.tl10[1:26, ])
se.na.tl10.tte.hte.n1400 <- t(res.na.tl10[27:46, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.re <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxme(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),data = train0)
        mod1 <- coxme(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),data = train1)
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        coef_mat.1 = mean(coef_mat0[1, ]),
        coef_mat.2 = mean(coef_mat0[2, ]),
        coef_mat.3 = mean(coef_mat0[3, ]),
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
        coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
        se_mat.1 = mean(se_mat0[1, ]),
        se_mat.2 = mean(se_mat0[2, ]),
        se_mat.3 = mean(se_mat0[3, ]),
        se_mat.4 = mean(se_mat0[4, ]),
        se_mat.5 = mean(se_mat0[5, ]),
        se_mat.6 = mean(se_mat0[6, ]),
        se_mat.7 = mean(se_mat0[7, ]),
        se_mat.8 = mean(se_mat0[8, ]),
        se_mat.9 = mean(se_mat0[9, ]),
        se_mat.10 = mean(se_mat0[10, ]),
        se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.re.tl10.tte.hte.n1400 <- t(res.re.tl10[1:26, ])
se.re.tl10.tte.hte.n1400 <- t(res.re.tl10[27:46 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.si <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+strata(trial),
                      data = train0)
        mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+strata(trial),
                      data = train1)
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        coef_mat.1 = mean(coef_mat0[1, ]),
        coef_mat.2 = mean(coef_mat0[2, ]),
        coef_mat.3 = mean(coef_mat0[3, ]),
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = mean(coef_mat1[1, ])- mean(coef_mat0[1, ]),
        coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
        se_mat.1 = mean(se_mat0[1, ]),
        se_mat.2 = mean(se_mat0[2, ]),
        se_mat.3 = mean(se_mat0[3, ]),
        se_mat.4 = mean(se_mat0[4, ]),
        se_mat.5 = mean(se_mat0[5, ]),
        se_mat.6 = mean(se_mat0[6, ]),
        se_mat.7 = mean(se_mat0[7, ]),
        se_mat.8 = mean(se_mat0[8, ]),
        se_mat.9 = mean(se_mat0[9, ]),
        se_mat.10 = mean(se_mat0[10, ]),
        se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.si.tl10.tte.hte.n1400 <- t(res.si.tl10[1:26, ])
se.si.tl10.tte.hte.n1400 <- t(res.si.tl10[27:46, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    source("coxvc.R")
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.r1 <- c()
    coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        Ft0<-cbind(rep(1,nrow(train0)),bs(train0$time))
        mod0 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+trial, Ft0, rank = 1,
                      data = train0)
        Ft1<-cbind(rep(1,nrow(train0)),bs(train1$time))
        mod1 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+trial, Ft1, rank = 1,
                      data = train1)
        
        coef0 <- mod0$b
        coef1 <- mod1$b
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
        #recalibration of event rates adapted to train
        lp0 <- unlist(as.matrix(train0[2:11]) %*% coef[1:10])
        lp1 <- unlist(as.matrix(train1[2:11]) %*% coef[1:10])
        temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
        temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
        bh0 <- basehaz(temp0, centered=F)
        bh1 <- basehaz(temp1, centered=F)
        int0[i] <- bh0[nrow(bh0),]$hazard
        int1[i] <- bh1[nrow(bh1),]$hazard
        
        #ite estimation
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        oben <- unlist(ifelse(sapply(oben, length) == 0, 0, oben))
        lm_ <- lm(oben~pben)
        
        a[i] <- lm_$coefficients[[1]]
        b[i] <- lm_$coefficients[[2]]
        mse.r1[i] <- mse(oben,pben)
      }
    })
    if(any(class(testerror) == "try-error")) return(list(success = FALSE))
    list(
      res = c(c.ben = mean(c.ben),
              c.ben.se = mean(c.ben.se),
              a = mean(a),
              b = mean(b),
              mse = mean(mse.r1),
              coef_mat.1 = mean(coef_mat0[1, ]),
              coef_mat.2 = mean(coef_mat0[2, ]),
              coef_mat.3 = mean(coef_mat0[3, ]),
              coef_mat.4 = mean(coef_mat0[4, ]),
              coef_mat.5 = mean(coef_mat0[5, ]),
              coef_mat.6 = mean(coef_mat0[6, ]),
              coef_mat.7 = mean(coef_mat0[7, ]),
              coef_mat.8 = mean(coef_mat0[8, ]),
              coef_mat.9 = mean(coef_mat0[9, ]),
              coef_mat.10 = mean(coef_mat0[10, ]),
              coef_mat.11 = mean(int1)-mean(int0),
              coef_mat.12 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
              coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
              coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
              coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
              coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
              coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
              coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
              coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
              coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
              coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.r1.tl10.tte.hte.n1400 <- t(res.r1.tl10[1:26, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl10 = foreach(
  j = 1:m,
  .final = function(l) do.call("cbind",
                               lapply(
                                 l[sapply(l, `[[`, "success")],
                                 function(subl) c(subl$res, subl$seed))),
  .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar%
  {
  set.seed(j)
  n <- 1400
  k <- 7
  
  trial <- rep(1:k, each = floor(n/k))
  df10 <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df10$age <- df10$age-70
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df10$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df10$sbp <- df10$sbp-185
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df10$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df10$height <- df10$height-168
  
  mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
  sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
  df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
  df10$chol <- df10$chol-6.3
  
  b0 <- 50-((mean_age-60)/100)
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3-((mean_age-60)/100)
  b5 <- 0.82
  b6 <- 0.8
  b7 <- 0.7
  b8 <- 0.1
  b9 <- 0.33
  b10 <- -0.02
  b11 <- 0.06
  b12 <- 0.01+((mean_age-70)/1000)
  b13 <- 0.04+0.005*pman
  b14 <- -0.2-((mean_sbp-170)/100)
  b15 <- -0.01+0.003*psmoke
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  log_hr <- log_hr - mean(log_hr)
  
  sp <- 1.15
  t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
  c <- runif(n, 0, 5)
  df10$time <- pmin(t,c)  #observed time is min of censored and true
  df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  int0 <- c()
  int1 <- c()
  coef_mat0 <- matrix(NA, nrow = 65, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 65, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                         factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(trial),
                    data = train0)
      mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                         factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(trial),
                    data = train1)
      coef0 <- mod0$coefficients
      coef1 <- mod1$coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      coef_mat0[,i] <- coef0
      coef_mat1[,i] <- coef1
      
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
      pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(mean(coef0[11:15])*test$trial+
                                               mean(coef0[c(1,16:20)])*test$age+mean(coef0[c(2,21:25)])*test$sex+
                                               mean(coef0[c(3,26:30)])*test$sbp+mean(coef0[c(4,31:35)])*test$mi+
                                               mean(coef0[c(5,36:40)])*test$stroke+mean(coef0[c(6,41:45)])*test$smoke+
                                               mean(coef0[c(7,46:50)])*test$diab+mean(coef0[c(8,51:55)])*test$lvhn1+
                                               mean(coef0[c(9,56:60)])*test$height+mean(coef0[c(10,61:65)])*test$chol))
      pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(mean(coef1[11:15])*test$trial+
                                               mean(coef1[c(1,16:20)])*test$age+mean(coef1[c(2,21:25)])*test$sex+
                                               mean(coef1[c(3,26:30)])*test$sbp+mean(coef1[c(4,31:35)])*test$mi+
                                               mean(coef1[c(5,36:40)])*test$stroke+mean(coef1[c(6,41:45)])*test$smoke+
                                               mean(coef1[c(7,46:50)])*test$diab+mean(coef1[c(8,51:55)])*test$lvhn1+
                                               mean(coef1[c(9,56:60)])*test$height+mean(coef1[c(10,61:65)])*test$chol))
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
  if (any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
      c.ben = mean(c.ben),
      c.ben.se = mean(c.ben.se),
      a = mean(a),
      b = mean(b),
      mse = mean(mse.fs),
      coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[16:20, ]),
      coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[21:25, ]),
      coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[26:30, ]),
      coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[31:35, ]),
      coef_mat.5 = mean(coef_mat0[5, ]) + mean(coef_mat0[36:40, ]),
      coef_mat.6 = mean(coef_mat0[6, ]) + mean(coef_mat0[41:45, ]),
      coef_mat.7 = mean(coef_mat0[7, ]) + mean(coef_mat0[46:50, ]),
      coef_mat.8 = mean(coef_mat0[8, ]) + mean(coef_mat0[51:55, ]),
      coef_mat.9 = mean(coef_mat0[9, ]) + mean(coef_mat0[56:60, ]),
      coef_mat.10 = mean(coef_mat0[10, ]) + mean(coef_mat0[61:65, ]),
      coef_mat.11 = mean(int1)-mean(int0),
      coef_mat.12 = (mean(coef_mat1[1]) + mean(coef_mat1[16:20, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[16:20, ])),
      coef_mat.13 = (mean(coef_mat1[2]) + mean(coef_mat1[17:21, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[21:25, ])),
      coef_mat.14 = (mean(coef_mat1[3]) + mean(coef_mat1[22:26, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[26:30, ])),
      coef_mat.15 = (mean(coef_mat1[4]) + mean(coef_mat1[27:31, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[31:35, ])),
      coef_mat.16 = (mean(coef_mat1[5]) + mean(coef_mat1[32:36, ])) - (mean(coef_mat0[5, ]) + mean(coef_mat0[36:40, ])),
      coef_mat.17 = (mean(coef_mat1[6]) + mean(coef_mat1[37:41, ])) - (mean(coef_mat0[6, ]) + mean(coef_mat0[41:45, ])),
      coef_mat.18 = (mean(coef_mat1[7]) + mean(coef_mat1[42:46, ])) - (mean(coef_mat0[7, ]) + mean(coef_mat0[46:50, ])),
      coef_mat.19 = (mean(coef_mat1[8]) + mean(coef_mat1[47:51, ])) - (mean(coef_mat0[8, ]) + mean(coef_mat0[51:55, ])),
      coef_mat.20 = (mean(coef_mat1[9]) + mean(coef_mat1[52:56, ])) - (mean(coef_mat0[9, ]) + mean(coef_mat0[56:60, ])),
      coef_mat.21 = (mean(coef_mat1[10]) + mean(coef_mat1[57:61, ])) - (mean(coef_mat0[10, ]) + mean(coef_mat0[61:65, ]))),
    seed = j,
    success = TRUE
  )
}
stopCluster(cl)
res.fs.tl10.tte.hte.n1400 <- t(res.fs.tl10[1:26, ])

save(res.na.sl10.tte.hte.n1400,se.na.sl10.tte.hte.n1400,
     res.re.sl10.tte.hte.n1400,se.re.sl10.tte.hte.n1400,
     res.si.sl10.tte.hte.n1400,se.si.sl10.tte.hte.n1400,
     res.fs.sl10.tte.hte.n1400,
     res.na.tl10.tte.hte.n1400,se.na.tl10.tte.hte.n1400,
     res.re.tl10.tte.hte.n1400,se.re.tl10.tte.hte.n1400,
     res.si.tl10.tte.hte.n1400,se.si.tl10.tte.hte.n1400,
     res.fs.tl10.tte.hte.n1400,
     file = "res_scenario23.Rdata") #res.r1.sl10.tte.hte.n1400,res.r1.tl10.tte.hte.n1400


#### scenario 24: 10 covariates & tte outcome & total sample size = 700 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = floor(n/k))
                        df10 <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df10$age <- df10$age-70
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df10$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df10$sbp <- df10$sbp-185
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                        df10$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df10$height <- df10$height-168
                        
                        mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
                        sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
                        df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
                        df10$chol <- df10$chol-6.3
                        
                        b0 <- 50-((mean_age-60)/100)
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- -0.3-((mean_age-60)/100)
                        b5 <- 0.82
                        b6 <- 0.8
                        b7 <- 0.7
                        b8 <- 0.1
                        b9 <- 0.33
                        b10 <- -0.02
                        b11 <- 0.06
                        b12 <- 0.01+((mean_age-70)/1000)
                        b13 <- 0.04+0.005*pman
                        b14 <- -0.2-((mean_sbp-170)/100)
                        b15 <- -0.01+0.003*psmoke
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        log_hr <- log_hr - mean(log_hr)
                        
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0, 5)
                        df10$time <- pmin(t,c)  #observed time is min of censored and true
                        df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat <- matrix(NA, nrow = 21, ncol = k)
                        se_mat <- matrix(NA, nrow = 21, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat),
                                         data = train)
                            coef <- mod$coefficients
                            coef[is.na(coef)] <- 0
                            coef_mat[,i] <- coef
                            se_mat[,i] <- sqrt(diag(vcov(mod)))
                            
                            #recalibration of event rates adapted to train
                            lp <- mod$linear.predictors
                            temp <- coxph(Surv(train$time,train$status)~offset(lp))
                            bh <- basehaz(temp, centered=F)
                            
                            #ite estimation
                            pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
                            pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
                            coef_mat.8 = mean(coef_mat[8, ]),
                            coef_mat.9 = mean(coef_mat[9, ]),
                            coef_mat.10 = mean(coef_mat[10, ]),
                            coef_mat.11 = mean(coef_mat[11, ]),
                            coef_mat.12 = mean(coef_mat[12, ]),
                            coef_mat.13 = mean(coef_mat[13, ]),
                            coef_mat.14 = mean(coef_mat[14, ]),
                            coef_mat.15 = mean(coef_mat[15, ]),
                            coef_mat.16 = mean(coef_mat[16, ]),
                            coef_mat.17 = mean(coef_mat[17, ]),
                            coef_mat.18 = mean(coef_mat[18, ]),
                            coef_mat.19 = mean(coef_mat[19, ]),
                            coef_mat.20 = mean(coef_mat[20, ]),
                            coef_mat.21 = mean(coef_mat[21, ]),
                            se_mat.1 = mean(se_mat[1, ]),
                            se_mat.2 = mean(se_mat[2, ]),
                            se_mat.3 = mean(se_mat[3, ]),
                            se_mat.4 = mean(se_mat[4, ]),
                            se_mat.5 = mean(se_mat[5, ]),
                            se_mat.6 = mean(se_mat[6, ]),
                            se_mat.7 = mean(se_mat[7, ]),
                            se_mat.8 = mean(se_mat[8, ]),
                            se_mat.9 = mean(se_mat[9, ]),
                            se_mat.10 = mean(se_mat[10, ]),
                            se_mat.11 = mean(se_mat[11, ]),
                            se_mat.12 = mean(se_mat[12, ]),
                            se_mat.13 = mean(se_mat[13, ]),
                            se_mat.14 = mean(se_mat[14, ]),
                            se_mat.15 = mean(se_mat[15, ]),
                            se_mat.16 = mean(se_mat[16, ]),
                            se_mat.17 = mean(se_mat[17, ]),
                            se_mat.18 = mean(se_mat[18, ]),
                            se_mat.19 = mean(se_mat[19, ]),
                            se_mat.20 = mean(se_mat[20, ]),
                            se_mat.21 = mean(se_mat[21, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.sl10.tte.hte.n700 <- t(res.na.sl10[1:26, ])
se.na.sl10.tte.hte.n700 <- t(res.na.sl10[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.re <- c()
    coef_mat <- matrix(NA, nrow = 21, ncol = k)
    se_mat <- matrix(NA, nrow = 21, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying the model to train
        mod <- coxme(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+
                                          factor(diab)+lvhn1+height+chol)*factor(treat)+(1+treat|trial),
                     data = train)
        coef <- mod$coefficients
        coef_mat[,i] <- coef
        se_mat[,i] <- sqrt(diag(vcov(mod)))
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictor
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
        coef_mat.1 = mean(coef_mat[1, ]),
        coef_mat.2 = mean(coef_mat[2, ]),
        coef_mat.3 = mean(coef_mat[3, ]),
        coef_mat.4 = mean(coef_mat[4, ]),
        coef_mat.5 = mean(coef_mat[5, ]),
        coef_mat.6 = mean(coef_mat[6, ]),
        coef_mat.7 = mean(coef_mat[7, ]),
        coef_mat.8 = mean(coef_mat[8, ]),
        coef_mat.9 = mean(coef_mat[9, ]),
        coef_mat.10 = mean(coef_mat[10, ]),
        coef_mat.11 = mean(coef_mat[11, ]),
        coef_mat.12 = mean(coef_mat[12, ]),
        coef_mat.13 = mean(coef_mat[13, ]),
        coef_mat.14 = mean(coef_mat[14, ]),
        coef_mat.15 = mean(coef_mat[15, ]),
        coef_mat.16 = mean(coef_mat[16, ]),
        coef_mat.17 = mean(coef_mat[17, ]),
        coef_mat.18 = mean(coef_mat[18, ]),
        coef_mat.19 = mean(coef_mat[19, ]),
        coef_mat.20 = mean(coef_mat[20, ]),
        coef_mat.21 = mean(coef_mat[21, ]),
        se_mat.1 = mean(se_mat[1, ]),
        se_mat.2 = mean(se_mat[2, ]),
        se_mat.3 = mean(se_mat[3, ]),
        se_mat.4 = mean(se_mat[4, ]),
        se_mat.5 = mean(se_mat[5, ]),
        se_mat.6 = mean(se_mat[6, ]),
        se_mat.7 = mean(se_mat[7, ]),
        se_mat.8 = mean(se_mat[8, ]),
        se_mat.9 = mean(se_mat[9, ]),
        se_mat.10 = mean(se_mat[10, ]),
        se_mat.11 = mean(se_mat[11, ]),
        se_mat.12 = mean(se_mat[12, ]),
        se_mat.13 = mean(se_mat[13, ]),
        se_mat.14 = mean(se_mat[14, ]),
        se_mat.15 = mean(se_mat[15, ]),
        se_mat.16 = mean(se_mat[16, ]),
        se_mat.17 = mean(se_mat[17, ]),
        se_mat.18 = mean(se_mat[18, ]),
        se_mat.19 = mean(se_mat[19, ]),
        se_mat.20 = mean(se_mat[20, ]),
        se_mat.21 = mean(se_mat[21, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.re.sl10.tte.hte.n700 <- t(res.re.sl10[1:26, ])
se.re.sl10.tte.hte.n700 <- t(res.re.sl10[27:47 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.si <- c()
    coef_mat <- matrix(NA, nrow = 21, ncol = k)
    se_mat <- matrix(NA, nrow = 21, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+strata(trial),
                     data = train) #+(0+treat|trial)
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        se_mat[,i] <- sqrt(diag(vcov(mod)))
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictors
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[12:21])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[12:21])))
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
        coef_mat.8 = mean(coef_mat[8, ]),
        coef_mat.9 = mean(coef_mat[9, ]),
        coef_mat.10 = mean(coef_mat[10, ]),
        coef_mat.11 = mean(coef_mat[11, ]),
        coef_mat.12 = mean(coef_mat[12, ]),
        coef_mat.13 = mean(coef_mat[13, ]),
        coef_mat.14 = mean(coef_mat[14, ]),
        coef_mat.15 = mean(coef_mat[15, ]),
        coef_mat.16 = mean(coef_mat[16, ]),
        coef_mat.17 = mean(coef_mat[17, ]),
        coef_mat.18 = mean(coef_mat[18, ]),
        coef_mat.19 = mean(coef_mat[19, ]),
        coef_mat.20 = mean(coef_mat[20, ]),
        coef_mat.21 = mean(coef_mat[21, ]),
        se_mat.1 = mean(se_mat[1, ]),
        se_mat.2 = mean(se_mat[2, ]),
        se_mat.3 = mean(se_mat[3, ]),
        se_mat.4 = mean(se_mat[4, ]),
        se_mat.5 = mean(se_mat[5, ]),
        se_mat.6 = mean(se_mat[6, ]),
        se_mat.7 = mean(se_mat[7, ]),
        se_mat.8 = mean(se_mat[8, ]),
        se_mat.9 = mean(se_mat[9, ]),
        se_mat.10 = mean(se_mat[10, ]),
        se_mat.11 = mean(se_mat[11, ]),
        se_mat.12 = mean(se_mat[12, ]),
        se_mat.13 = mean(se_mat[13, ]),
        se_mat.14 = mean(se_mat[14, ]),
        se_mat.15 = mean(se_mat[15, ]),
        se_mat.16 = mean(se_mat[16, ]),
        se_mat.17 = mean(se_mat[17, ]),
        se_mat.18 = mean(se_mat[18, ]),
        se_mat.19 = mean(se_mat[19, ]),
        se_mat.20 = mean(se_mat[20, ]),
        se_mat.21 = mean(se_mat[21, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.si.sl10.tte.hte.n700 <- t(res.si.sl10[1:26, ])
se.si.sl10.tte.hte.n700 <- t(res.si.sl10[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        
    source("coxvc.R")
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.r1 <- c()
    coef_mat <- matrix(NA, nrow = 22, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        Ft<-cbind(rep(1,nrow(train)),bs(train$time))
        mod <- coxvc(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                          factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial,
                     Ft, rank = 1, data = train)
        coef <- mod$b
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        lp <- unlist(as.matrix(train[2:12]) %*% coef[1:11] + train[12] * (as.matrix(train[2:11]) %*% coef[13:22]))
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test0[2:12]) %*% coef[1:11] + test0[12] * (as.matrix(test0[2:11]) %*% coef[13:22])))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(as.matrix(test1[2:12]) %*% coef[1:11] + test1[12] * (as.matrix(test1[2:11]) %*% coef[13:22])))
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
      coef_mat.5 = mean(coef_mat[5, ]),
      coef_mat.6 = mean(coef_mat[6, ]),
      coef_mat.7 = mean(coef_mat[7, ]),
      coef_mat.8 = mean(coef_mat[8, ]),
      coef_mat.9 = mean(coef_mat[9, ]),
      coef_mat.10 = mean(coef_mat[10, ]),
      coef_mat.11 = mean(coef_mat[11, ]),
      coef_mat.12 = mean(coef_mat[12, ]),
      coef_mat.13 = mean(coef_mat[14, ]),
      coef_mat.14 = mean(coef_mat[15, ]),
      coef_mat.15 = mean(coef_mat[16, ]),
      coef_mat.16 = mean(coef_mat[17, ]),
      coef_mat.17 = mean(coef_mat[18, ]),
      coef_mat.18 = mean(coef_mat[19, ]),
      coef_mat.19 = mean(coef_mat[20, ]),
      coef_mat.20 = mean(coef_mat[21, ]),
      coef_mat.21 = mean(coef_mat[22, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.r1.sl10.tte.hte.n700 <- t(res.r1.sl10[1:26, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat <- matrix(NA, nrow = 76, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- coxph(Surv(time,status)~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+
                       (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol):factor(trial),
                     data = train)
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        lp <- mod$linear.predictors
        temp <- coxph(Surv(train$time,train$status)~offset(lp))
        bh <- basehaz(temp, centered=F)
        
        #ite estimation
        pred0 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test0$trial+
                                                 mean(coef[c(6,27:31)])*test0$age+mean(coef[c(7,32:36)])*test0$sex+
                                                 mean(coef[c(8,37:41)])*test0$sbp+mean(coef[c(9,42:46)])*test0$mi+
                                                 mean(coef[c(10,47:51)])*test0$stroke+mean(coef[c(11,52:56)])*test0$smoke+
                                                 mean(coef[c(12,57:61)])*test0$diab+mean(coef[c(13,62:66)])*test0$lvhn1+
                                                 mean(coef[c(14,67:71)])*test0$height+mean(coef[c(15,72:76)])*test0$chol+
                                                 coef[16]*test0$treat+(coef[17]*test0$age+coef[18]*test0$sex+
                                                                         coef[19]*test0$sbp+coef[20]*test0$mi+coef[21]*test0$stroke+
                                                                         coef[22]*test0$smoke+coef[23]*test0$diab+coef[24]*test0$lvhn1+
                                                                         coef[25]*test0$height+coef[26]*test0$chol)*test0$treat))
        pred1 <- exp(-bh[nrow(bh),]$hazard*exp(mean(coef[1:5])*test1$trial+
                                                 mean(coef[c(6,27:31)])*test1$age+mean(coef[c(7,32:36)])*test1$sex+
                                                 mean(coef[c(8,37:41)])*test1$sbp+mean(coef[c(9,42:46)])*test1$mi+
                                                 mean(coef[c(10,47:51)])*test1$stroke+mean(coef[c(11,52:56)])*test1$smoke+
                                                 mean(coef[c(12,57:61)])*test1$diab+mean(coef[c(13,62:66)])*test1$lvhn1+
                                                 mean(coef[c(14,67:71)])*test1$height+mean(coef[c(15,72:76)])*test1$chol+
                                                 coef[16]*test1$treat+(coef[17]*test1$age+coef[18]*test1$sex+
                                                                         coef[19]*test1$sbp+coef[20]*test1$mi+coef[21]*test1$stroke+
                                                                         coef[22]*test1$smoke+coef[23]*test1$diab+coef[24]*test1$lvhn1+
                                                                         coef[25]*test1$height+coef[26]*test1$chol)*test1$treat))
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
    
    list(res = c(
      c.ben = mean(c.ben),
      c.ben.se = mean(c.ben.se),
      a = mean(a),
      b = mean(b),
      mse = mean(mse.fs),
      coef_mat.1 = mean(coef_mat[6, ])+coef_mat[27:31, ]),
      coef_mat.2 = mean(c(coef_mat[7, ]+coef_mat[32:36, ])),
      coef_mat.3 = mean(c(coef_mat[8, ]+coef_mat[37:41, ])),
      coef_mat.4 = mean(c(coef_mat[9, ]+coef_mat[42:46, ])),
      coef_mat.5 = mean(c(coef_mat[10, ]+coef_mat[47:51, ])),
      coef_mat.6 = mean(c(coef_mat[11, ]+coef_mat[52:56, ])),
      coef_mat.7 = mean(c(coef_mat[12, ]+coef_mat[57:61, ])),
      coef_mat.8 = mean(c(coef_mat[13, ]+coef_mat[62:66, ])),
      coef_mat.9 = mean(c(coef_mat[14, ]+coef_mat[67:71, ])),
      coef_mat.10 = mean(c(coef_mat[15, ]+coef_mat[72:76, ])),
      coef_mat.11 = mean(coef_mat[16, ]),
      coef_mat.12 = mean(coef_mat[17, ]),
      coef_mat.13 = mean(coef_mat[18, ]),
      coef_mat.14 = mean(coef_mat[19, ]),
      coef_mat.15 = mean(coef_mat[20, ]),
      coef_mat.16 = mean(coef_mat[21, ]),
      coef_mat.17 = mean(coef_mat[22, ]),
      coef_mat.18 = mean(coef_mat[23, ]),
      coef_mat.19 = mean(coef_mat[24, ]),
      coef_mat.20 = mean(coef_mat[25, ]),
      coef_mat.21 = mean(coef_mat[26, ]),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.fs.sl10.tte.hte.n700 <- t(res.fs.sl10[1:26, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        set.seed(j)
                        n <- 700
                        k <- 7
                        
                        trial <- rep(1:k, each = floor(n/k))
                        df10 <- as.data.frame(trial)
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df10$age <- df10$age-70
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df10$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df10$sbp <- df10$sbp-185
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                        df10$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df10$height <- df10$height-168
                        
                        mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
                        sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
                        df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
                        df10$chol <- df10$chol-6.3
                        
                        b0 <- 50-((mean_age-60)/100)
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- -0.3-((mean_age-60)/100)
                        b5 <- 0.82
                        b6 <- 0.8
                        b7 <- 0.7
                        b8 <- 0.1
                        b9 <- 0.33
                        b10 <- -0.02
                        b11 <- 0.06
                        b12 <- 0.01+((mean_age-70)/1000)
                        b13 <- 0.04+0.005*pman
                        b14 <- -0.2-((mean_sbp-170)/100)
                        b15 <- -0.01+0.003*psmoke
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        log_hr <- log_hr - mean(log_hr)
                        
                        sp <- 1.15
                        t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
                        c <- runif(n, 0, 5)
                        df10$time <- pmin(t,c)  #observed time is min of censored and true
                        df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
                        
                        df0 <- df10[df10$treat==0,]
                        df1 <- df10[df10$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        int0 <- c()
                        int1 <- c()
                        coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                            factor(smoke)+factor(diab)+lvhn1+height+chol,
                                          data = train0)
                            mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                            factor(smoke)+factor(diab)+lvhn1+height+chol,
                                          data = train1)
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            coef0[is.na(coef0)] <- 0
                            coef1[is.na(coef1)] <- 0
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
                            pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
                            pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
                            coef_mat.1 = mean(coef_mat0[1, ]),
                            coef_mat.2 = mean(coef_mat0[2, ]),
                            coef_mat.3 = mean(coef_mat0[3, ]),
                            coef_mat.4 = mean(coef_mat0[4, ]),
                            coef_mat.5 = mean(coef_mat0[5, ]),
                            coef_mat.6 = mean(coef_mat0[6, ]),
                            coef_mat.7 = mean(coef_mat0[7, ]),
                            coef_mat.8 = mean(coef_mat0[8, ]),
                            coef_mat.9 = mean(coef_mat0[9, ]),
                            coef_mat.10 = mean(coef_mat0[10, ]),
                            coef_mat.11 = mean(int1)-mean(int0),
                            coef_mat.12 = mean(coef_mat1[1, ])- mean(coef_mat0[1, ]),
                            coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                            coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                            coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                            coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                            coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                            coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                            coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                            coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                            coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
                            se_mat.1 = mean(se_mat0[1, ]),
                            se_mat.2 = mean(se_mat0[2, ]),
                            se_mat.3 = mean(se_mat0[3, ]),
                            se_mat.4 = mean(se_mat0[4, ]),
                            se_mat.5 = mean(se_mat0[5, ]),
                            se_mat.6 = mean(se_mat0[6, ]),
                            se_mat.7 = mean(se_mat0[7, ]),
                            se_mat.8 = mean(se_mat0[8, ]),
                            se_mat.9 = mean(se_mat0[9, ]),
                            se_mat.10 = mean(se_mat0[10, ]),
                            se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                            se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                            se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                            se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                            se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                            se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                            se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                            se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                            se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                            se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.tl10.tte.hte.n700 <- t(res.na.tl10[1:26, ])
se.na.tl10.tte.hte.n700 <- t(res.na.tl10[27:46, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.re <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxme(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),data = train0)
        mod1 <- coxme(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),data = train1)
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        coef_mat.1 = mean(coef_mat0[1, ]),
        coef_mat.2 = mean(coef_mat0[2, ]),
        coef_mat.3 = mean(coef_mat0[3, ]),
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
        coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
        se_mat.1 = mean(se_mat0[1, ]),
        se_mat.2 = mean(se_mat0[2, ]),
        se_mat.3 = mean(se_mat0[3, ]),
        se_mat.4 = mean(se_mat0[4, ]),
        se_mat.5 = mean(se_mat0[5, ]),
        se_mat.6 = mean(se_mat0[6, ]),
        se_mat.7 = mean(se_mat0[7, ]),
        se_mat.8 = mean(se_mat0[8, ]),
        se_mat.9 = mean(se_mat0[9, ]),
        se_mat.10 = mean(se_mat0[10, ]),
        se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.re.tl10.tte.hte.n700 <- t(res.re.tl10[1:26, ])
se.re.tl10.tte.hte.n700 <- t(res.re.tl10[27:46 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.si <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+strata(trial),
                      data = train0)
        mod1 <- coxph(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+strata(trial),
                      data = train1)
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        coef_mat.1 = mean(coef_mat0[1, ]),
        coef_mat.2 = mean(coef_mat0[2, ]),
        coef_mat.3 = mean(coef_mat0[3, ]),
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = mean(coef_mat1[1, ])- mean(coef_mat0[1, ]),
        coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
        se_mat.1 = mean(se_mat0[1, ]),
        se_mat.2 = mean(se_mat0[2, ]),
        se_mat.3 = mean(se_mat0[3, ]),
        se_mat.4 = mean(se_mat0[4, ]),
        se_mat.5 = mean(se_mat0[5, ]),
        se_mat.6 = mean(se_mat0[6, ]),
        se_mat.7 = mean(se_mat0[7, ]),
        se_mat.8 = mean(se_mat0[8, ]),
        se_mat.9 = mean(se_mat0[9, ]),
        se_mat.10 = mean(se_mat0[10, ]),
        se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.si.tl10.tte.hte.n700 <- t(res.si.tl10[1:26, ])
se.si.tl10.tte.hte.n700 <- t(res.si.tl10[27:46, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar% {
                        source("coxvc.R")
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.r1 <- c()
    coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        Ft0<-cbind(rep(1,nrow(train0)),bs(train0$time))
        mod0 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+trial, Ft0, rank = 1,
                      data = train0)
        Ft1<-cbind(rep(1,nrow(train0)),bs(train1$time))
        mod1 <- coxvc(Surv(time,status)~age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                        factor(smoke)+factor(diab)+lvhn1+height+chol+trial, Ft1, rank = 1,
                      data = train1)
        
        coef0 <- mod0$b
        coef1 <- mod1$b
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
        #recalibration of event rates adapted to train
        lp0 <- unlist(as.matrix(train0[2:11]) %*% coef[1:10])
        lp1 <- unlist(as.matrix(train1[2:11]) %*% coef[1:10])
        temp0 <- coxph(Surv(train0$time,train0$status)~offset(lp0))
        temp1 <- coxph(Surv(train1$time,train1$status)~offset(lp1))
        bh0 <- basehaz(temp0, centered=F)
        bh1 <- basehaz(temp1, centered=F)
        int0[i] <- bh0[nrow(bh0),]$hazard
        int1[i] <- bh1[nrow(bh1),]$hazard
        
        #ite estimation
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(as.matrix(test[2:11]) %*% coef0[1:10]))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(as.matrix(test[2:11]) %*% coef1[1:10]))
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
        oben <- unlist(ifelse(sapply(oben, length) == 0, 0, oben))
        lm_ <- lm(oben~pben)
        
        a[i] <- lm_$coefficients[[1]]
        b[i] <- lm_$coefficients[[2]]
        mse.r1[i] <- mse(oben,pben)
      }
    })
    if(any(class(testerror) == "try-error")) return(list(success = FALSE))
    list(
      res = c(c.ben = mean(c.ben),
              c.ben.se = mean(c.ben.se),
              a = mean(a),
              b = mean(b),
              mse = mean(mse.r1),
              coef_mat.1 = mean(coef_mat0[1, ]),
              coef_mat.2 = mean(coef_mat0[2, ]),
              coef_mat.3 = mean(coef_mat0[3, ]),
              coef_mat.4 = mean(coef_mat0[4, ]),
              coef_mat.5 = mean(coef_mat0[5, ]),
              coef_mat.6 = mean(coef_mat0[6, ]),
              coef_mat.7 = mean(coef_mat0[7, ]),
              coef_mat.8 = mean(coef_mat0[8, ]),
              coef_mat.9 = mean(coef_mat0[9, ]),
              coef_mat.10 = mean(coef_mat0[10, ]),
              coef_mat.11 = mean(int1)-mean(int0),
              coef_mat.12 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
              coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
              coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
              coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
              coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
              coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
              coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
              coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
              coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
              coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.r1.tl10.tte.hte.n700 <- t(res.r1.tl10[1:26, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl10 = foreach(
  j = 1:m,
  .final = function(l) do.call("cbind",
                               lapply(
                                 l[sapply(l, `[[`, "success")],
                                 function(subl) c(subl$res, subl$seed))),
  .packages = c("tidyverse","Hmisc","coxme","magrittr","gtools","survival")) %dopar%{
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n/k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df10$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df10$age <- df10$age-70
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df10$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df10$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df10$sbp <- df10$sbp-185
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
    df10$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob=pdiab[trial]) 
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df10$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df10$height <- df10$height-168
    
    mean_chol <- c(6.6,6.5,6.4,6.1,6.4,6,6.4) #Cholesterol
    sd_chol <- c(0.01,0.015,0.011,0.012,0.012,0.012,0.01)
    df10$chol <- round(rnorm(n,mean=mean_chol[trial], sd=sd_chol[trial]),0)
    df10$chol <- df10$chol-6.3
    
    b0 <- 50-((mean_age-60)/100)
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- -0.3-((mean_age-60)/100)
    b5 <- 0.82
    b6 <- 0.8
    b7 <- 0.7
    b8 <- 0.1
    b9 <- 0.33
    b10 <- -0.02
    b11 <- 0.06
    b12 <- 0.01+((mean_age-70)/1000)
    b13 <- 0.04+0.005*pman
    b14 <- -0.2-((mean_sbp-170)/100)
    b15 <- -0.01+0.003*psmoke
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_hr <- with(df10, b1*age+b2*(sex==1)+b3*sbp+b4*(treat==1)+b5*(mi==1)+b6*(stroke==1)+b7*(smoke==1)+b8*(diab==1)+b9*(lvhn1==1)+b10*height+b11*chol+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    log_hr <- log_hr - mean(log_hr)
    
    sp <- 1.15
    t <- rweibull(n, shape = 1.15, scale = b0/exp(log_hr)^sp)
    c <- runif(n, 0, 5)
    df10$time <- pmin(t,c)  #observed time is min of censored and true
    df10$status <- as.numeric(t<=c) # set to 1 if event is observed 0 if censored
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    int0 <- c()
    int1 <- c()
    coef_mat0 <- matrix(NA, nrow = 65, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 65, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                           factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(trial),
                      data = train0)
        mod1 <- coxph(Surv(time,status)~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+
                                           factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(trial),
                      data = train1)
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
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
        pred0 <- exp(-bh0[nrow(bh0),]$hazard*exp(mean(coef0[11:15])*test$trial+
                                                 mean(coef0[c(1,16:20)])*test$age+mean(coef0[c(2,21:25)])*test$sex+
                                                 mean(coef0[c(3,26:30)])*test$sbp+mean(coef0[c(4,31:35)])*test$mi+
                                                 mean(coef0[c(5,36:40)])*test$stroke+mean(coef0[c(6,41:45)])*test$smoke+
                                                 mean(coef0[c(7,46:50)])*test$diab+mean(coef0[c(8,51:55)])*test$lvhn1+
                                                 mean(coef0[c(9,56:60)])*test$height+mean(coef0[c(10,61:65)])*test$chol))
        pred1 <- exp(-bh1[nrow(bh1),]$hazard*exp(mean(coef1[11:15])*test$trial+
                                                 mean(coef1[c(1,16:20)])*test$age+mean(coef1[c(2,21:25)])*test$sex+
                                                 mean(coef1[c(3,26:30)])*test$sbp+mean(coef1[c(4,31:35)])*test$mi+
                                                 mean(coef1[c(5,36:40)])*test$stroke+mean(coef1[c(6,41:45)])*test$smoke+
                                                 mean(coef1[c(7,46:50)])*test$diab+mean(coef1[c(8,51:55)])*test$lvhn1+
                                                 mean(coef1[c(9,56:60)])*test$height+mean(coef1[c(10,61:65)])*test$chol))
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
    if (any(class(testerror) == "try-error")) return(list(success = FALSE))
    
    list(
      res = c(
        c.ben = mean(c.ben),
        c.ben.se = mean(c.ben.se),
        a = mean(a),
        b = mean(b),
        mse = mean(mse.fs),
        coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[16:20, ]),
        coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[21:25, ]),
        coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[26:30, ]),
        coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[31:35, ]),
        coef_mat.5 = mean(coef_mat0[5, ]) + mean(coef_mat0[36:40, ]),
        coef_mat.6 = mean(coef_mat0[6, ]) + mean(coef_mat0[41:45, ]),
        coef_mat.7 = mean(coef_mat0[7, ]) + mean(coef_mat0[46:50, ]),
        coef_mat.8 = mean(coef_mat0[8, ]) + mean(coef_mat0[51:55, ]),
        coef_mat.9 = mean(coef_mat0[9, ]) + mean(coef_mat0[56:60, ]),
        coef_mat.10 = mean(coef_mat0[10, ]) + mean(coef_mat0[61:65, ]),
        coef_mat.11 = mean(int1)-mean(int0),
        coef_mat.12 = (mean(coef_mat1[1]) + mean(coef_mat1[16:20, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[16:20, ])),
        coef_mat.13 = (mean(coef_mat1[2]) + mean(coef_mat1[17:21, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[21:25, ])),
        coef_mat.14 = (mean(coef_mat1[3]) + mean(coef_mat1[22:26, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[26:30, ])),
        coef_mat.15 = (mean(coef_mat1[4]) + mean(coef_mat1[27:31, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[31:35, ])),
        coef_mat.16 = (mean(coef_mat1[5]) + mean(coef_mat1[32:36, ])) - (mean(coef_mat0[5, ]) + mean(coef_mat0[36:40, ])),
        coef_mat.17 = (mean(coef_mat1[6]) + mean(coef_mat1[37:41, ])) - (mean(coef_mat0[6, ]) + mean(coef_mat0[41:45, ])),
        coef_mat.18 = (mean(coef_mat1[7]) + mean(coef_mat1[42:46, ])) - (mean(coef_mat0[7, ]) + mean(coef_mat0[46:50, ])),
        coef_mat.19 = (mean(coef_mat1[8]) + mean(coef_mat1[47:51, ])) - (mean(coef_mat0[8, ]) + mean(coef_mat0[51:55, ])),
        coef_mat.20 = (mean(coef_mat1[9]) + mean(coef_mat1[52:56, ])) - (mean(coef_mat0[9, ]) + mean(coef_mat0[56:60, ])),
        coef_mat.21 = (mean(coef_mat1[10]) + mean(coef_mat1[57:61, ])) - (mean(coef_mat0[10, ]) + mean(coef_mat0[61:65, ]))),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl10.tte.hte.n700 <- t(res.fs.tl10[1:26, ])

save(res.na.sl10.tte.hte.n700,se.na.sl10.tte.hte.n700,
     res.re.sl10.tte.hte.n700,se.re.sl10.tte.hte.n700,
     res.si.sl10.tte.hte.n700,se.si.sl10.tte.hte.n700,
     res.fs.sl10.tte.hte.n700,
     res.na.tl10.tte.hte.n700,se.na.tl10.tte.hte.n700,
     res.re.tl10.tte.hte.n700,se.re.tl10.tte.hte.n700,
     res.si.tl10.tte.hte.n700,se.si.tl10.tte.hte.n700,
     res.fs.tl10.tte.hte.n700,
     file = "res_scenario24.Rdata") #res.r1.sl10.tte.hte.n700,res.r1.tl10.tte.hte.n700
