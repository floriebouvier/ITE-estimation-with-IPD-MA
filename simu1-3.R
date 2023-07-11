library(tidyverse) # data.frame manipulation
library(Hmisc) # cstat computation
library(magrittr) # cstat computation
library(lme4) # glmer
library(doParallel) # parallel computing
library(gtools) # calibration computation
library(VGAM) # rrvglm

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

#### scenario 1: 3 covariates (1 binary and 2 continuous) & total sample size = 2800 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glm(death~(age+factor(sex)+sbp)*factor(treat),
               data = train,family = binomial,control = glm.control(maxit = 100000))
    coef <- mod$coefficients
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- X0 %*% coef[1:5]
    lp1 <- X1 %*% coef[1:5] + X1[,2:4] %*% coef[6:8] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(oben,pben)
  }
  c(
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
    se_mat.1 = mean(se_mat[1, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[6, ]),
    se_mat.7 = mean(se_mat[7, ]),
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.na.sl3.n2800 <- t(res.na.sl3[1:13, ])
se.na.sl3.n2800 <- t(res.na.sl3[14:21, ])
colnames(res.na.sl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.sl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- glmer(death~(age+factor(sex)+sbp)*factor(treat)+(1|trial)+(0+treat|trial),
                 data = train, family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- X0 %*% coef[1:5]
    lp1 <- X1 %*% coef[1:5] + X1[,2:4] %*% coef[6:8] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  c(
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
    coef_mat.8 = mean(coef_mat[8,]),
    se_mat.1 = mean(se_mat[1,]),
    se_mat.2 = mean(se_mat[2,]),
    se_mat.3 = mean(se_mat[3,]),
    se_mat.4 = mean(se_mat[4,]),
    se_mat.5 = mean(se_mat[5,]),
    se_mat.6 = mean(se_mat[6,]),
    se_mat.7 = mean(se_mat[7,]),
    se_mat.8 = mean(se_mat[8,]))
}
stopCluster(cl)
res.re.sl3.n2800 <- t(res.re.sl3[1:13, ])
se.re.sl3.n2800 <- t(res.re.sl3[14:21 , ])
colnames(res.re.sl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.sl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 9, ncol = k)
  se_mat <- matrix(NA, nrow = 9, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glmer(death~(age+factor(sex)+sbp)*factor(treat)+trial+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:6]
    lp1 <- X1 %*% coef[1:6] + X1[,2:4] %*% coef[7:9] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[6, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ]),
    coef_mat.8 = mean(coef_mat[9, ]),
    se_mat.1 = mean(se_mat[1, ]) + mean(se_mat[6, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[7, ]),
    se_mat.7 = mean(se_mat[8, ]),
    se_mat.8 = mean(se_mat[9, ]))
}
stopCluster(cl)
res.si.sl3.n2800 <- t(res.si.sl3[1:13, ])
se.si.sl3.n2800 <- t(res.si.sl3[14:21, ])
colnames(res.si.sl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.sl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 9, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- rrvglm(death ~ (age+factor(sex)+sbp)*factor(treat)+trial, family = binomialff, data = train,Rank = 1)
    coef <- mod@coefficients
    coef_mat[,i] <- coef
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:6]
    lp1 <- X1 %*% coef[1:6] + X1[,2:4] %*% coef[7:9] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[6, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ]),
    coef_mat.8 = mean(coef_mat[9, ]))
}
stopCluster(cl)
res.r1.sl3.n2800 <- t(res.r1.sl3[1:13, ])
colnames(res.r1.sl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat <- matrix(NA, nrow = 28, ncol = k)
  se_mat <- matrix(NA, nrow = 28, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glm(death~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial),
               data = train,family = binomial,control = glm.control(maxit = 100000))
    coef <- mod$coefficients
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial)"),
                      data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~trial+age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~trial+age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- coef[1]+mean(coef[2:6])*X0[,2]+mean(c(coef[7],coef[14:18]))*X0[,3]+mean(c(coef[8],coef[19:23]))*X0[,4]+
      mean(c(coef[9],coef[24:28]))*X0[,5]+coef[10]*X0[,6]
    lp1 <- coef[1]+mean(coef[2:6])*X1[,2]+mean(c(coef[7],coef[14:18]))*X1[,3]+mean(c(coef[8],coef[19:23]))*X1[,4]+
      mean(c(coef[9],coef[24:28]))*X1[,5]+coef[10]*X1[,6]+coef[11]*X1[,3]+coef[12]*X1[,4]+coef[13]*X1[,5]
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(oben,pben)
  } 
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
    coef_mat.1 = mean(coef_mat[1:6,]),
    coef_mat.2 = mean(c(coef_mat[7,]+coef_mat[14:18,])),
    coef_mat.3 = mean(c(coef_mat[8,]+coef_mat[19:23,])),
    coef_mat.4 = mean(c(coef_mat[9,]+coef_mat[24:28,])),
    coef_mat.5 = mean(coef_mat[10,]),
    coef_mat.6 = mean(coef_mat[11,]),
    coef_mat.7 = mean(coef_mat[12,]),
    coef_mat.8 = mean(coef_mat[13,]),
    se_mat.1 = mean(se_mat[1:6,]),
    se_mat.2 = mean(c(se_mat[7,]+se_mat[14:18,])),
    se_mat.3 = mean(c(se_mat[8,]+se_mat[19:23,])),
    se_mat.4 = mean(c(se_mat[9,]+se_mat[24:28,])),
    se_mat.5 = mean(se_mat[10,]),
    se_mat.6 = mean(se_mat[11,]),
    se_mat.7 = mean(se_mat[12,]),
    se_mat.8 = mean(se_mat[13,]))
}
stopCluster(cl)
res.fs.sl3.n2800 <- t(res.fs.sl3[1:13, ])
se.fs.sl3.n2800 <- t(res.fs.sl3[14:21 , ])
colnames(res.fs.sl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.sl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.na),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.na.tl3.n2800 <- t(res.na.tl3[1:13, ])
se.na.tl3.n2800 <- t(res.na.tl3[14:21 , ])
colnames(res.na.tl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.tl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")


#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glmer(death~age+factor(sex)+sbp+(1|trial),
                  data = train0,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    mod1 <- glmer(death~age+factor(sex)+sbp+(1|trial),
                  data = train1,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef0 <- mod0@beta
    coef1 <- mod1@beta
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.re.tl3.n2800 <- t(res.re.tl3[1:13, ])
se.re.tl3.n2800 <- t(res.re.tl3[14:21, ])
colnames(res.re.tl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.tl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat0 <- matrix(NA, nrow = 5, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 5, ncol = k)
  se_mat0 <- matrix(NA, nrow = 5, ncol = k)
  se_mat1 <- matrix(NA, nrow = 5, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp+trial,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp+trial,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5, ])),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[5, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[5, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[5, ])),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.si.tl3.n2800 <- t(res.si.tl3[1:13, ])
se.si.tl3.n2800 <- t(res.si.tl3[14:21 , ])
colnames(res.si.tl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.tl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat0 <- matrix(NA, nrow = 5, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 5, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- rrvglm(death ~ age+factor(sex)+sbp+trial, family = binomialff, data = train0,Rank = 1)
    mod1 <- rrvglm(death ~ age+factor(sex)+sbp+trial, family = binomialff, data = train1,Rank = 1)
    coef0 <- mod0@coefficients
    coef1 <- mod1@coefficients
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5, ])),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]))
}
stopCluster(cl)
res.r1.tl3.n2800 <- t(res.r1.tl3[1:13, ])
colnames(res.r1.tl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat0 <- matrix(NA, nrow = 24, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 24, ncol = k)
  se_mat0 <- matrix(NA, nrow = 24, ncol = k)
  se_mat1 <- matrix(NA, nrow = 24, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~(age+factor(sex)+sbp)*factor(trial),
                data = train0,family = binomial,control = glm.control(maxit = 100000))
    mod1 <- glm(death~(age+factor(sex)+sbp)*factor(trial),
                data = train1,family = binomial,control = glm.control(maxit = 100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(trial)"), data=train0)
    X1 <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(trial)"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- coef0[1]+mean(coef0[5:8])*X[,5]+mean(c(coef0[2],coef0[9:12]))*X[,2]+
      mean(c(coef0[3],coef0[13:16]))*X[,3]+mean(c(coef0[4],coef0[17:20]))*X[,4]
    lp1 <- coef1[1]+mean(coef1[5:8])*X[,5]+mean(c(coef1[2],coef1[9:12]))*X[,2]+
      mean(c(coef1[3],coef1[13:16]))*X[,3]+mean(c(coef1[4],coef1[17:20]))*X[,4]
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5:9, ]),
    coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[10:14, ]),
    coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[15:19, ]),
    coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[20:24, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5:9, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5:9, ])),
    coef_mat.6 = (mean(coef_mat1[2, ]) + mean(coef_mat1[10:14, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[10:14, ])),
    coef_mat.7 = (mean(coef_mat1[3, ]) + mean(coef_mat1[15:19, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[15:19, ])),
    coef_mat.8 = (mean(coef_mat1[4, ]) + mean(coef_mat1[20:24, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[20:24, ])),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[5:9, ]),
    se_mat.2 = mean(se_mat0[2, ]) + mean(se_mat0[10:14, ]),
    se_mat.3 = mean(se_mat0[3, ]) + mean(se_mat0[15:19, ]),
    se_mat.4 = mean(se_mat0[4, ]) + mean(se_mat0[20:24, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[5:9, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[5:9, ])),
    se_mat.6 = (mean(se_mat1[2, ]) + mean(se_mat1[10:14, ])) - (mean(se_mat0[2, ]) + mean(se_mat0[10:14, ])),
    se_mat.7 = (mean(se_mat1[3, ]) + mean(se_mat1[15:19, ])) - (mean(se_mat0[3, ]) + mean(se_mat0[15:19, ])),
    se_mat.8 = (mean(se_mat1[4, ]) + mean(se_mat1[20:24, ])) - (mean(se_mat0[4, ]) + mean(se_mat0[20:24, ])))
}
stopCluster(cl)
res.fs.tl3.n2800 <- t(res.fs.tl3[1:13, ])
se.fs.tl3.n2800 <- t(res.fs.tl3[14:21, ])
colnames(res.fs.tl3.n2800) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.tl3.n2800) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

save(res.na.sl3.n2800,se.na.sl3.n2800,
     res.re.sl3.n2800,se.re.sl3.n2800,
     res.si.sl3.n2800,se.si.sl3.n2800,
     res.r1.sl3.n2800,
     res.fs.sl3.n2800,se.fs.sl3.n2800,
     res.na.tl3.n2800,se.na.tl3.n2800,
     res.re.tl3.n2800,se.re.tl3.n2800,
     res.si.tl3.n2800,se.si.tl3.n2800,
     res.r1.tl3.n2800,
     res.fs.tl3.n2800,se.fs.tl3.n2800,
     file = "res_scenario1.Rdata")

#### scenario 2: 3 covariates (1 binary and 2 continuous) & total sample size = 1400 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glm(death~(age+factor(sex)+sbp)*factor(treat),
               data = train,family = binomial,control = glm.control(maxit = 100000))
    coef <- mod$coefficients
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- X0 %*% coef[1:5]
    lp1 <- X1 %*% coef[1:5] + X1[,2:4] %*% coef[6:8] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(oben,pben)
  }
  c(
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
    se_mat.1 = mean(se_mat[1, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[6, ]),
    se_mat.7 = mean(se_mat[7, ]),
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.na.sl3.n1400 <- t(res.na.sl3[1:13, ])
se.na.sl3.n1400 <- t(res.na.sl3[14:21, ])
colnames(res.na.sl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.sl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- glmer(death~(age+factor(sex)+sbp)*factor(treat)+(1|trial)+(0+treat|trial),
                 data = train, family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- X0 %*% coef[1:5]
    lp1 <- X1 %*% coef[1:5] + X1[,2:4] %*% coef[6:8] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  c(
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
    coef_mat.8 = mean(coef_mat[8,]),
    se_mat.1 = mean(se_mat[1,]),
    se_mat.2 = mean(se_mat[2,]),
    se_mat.3 = mean(se_mat[3,]),
    se_mat.4 = mean(se_mat[4,]),
    se_mat.5 = mean(se_mat[5,]),
    se_mat.6 = mean(se_mat[6,]),
    se_mat.7 = mean(se_mat[7,]),
    se_mat.8 = mean(se_mat[8,]))
}
stopCluster(cl)
res.re.sl3.n1400 <- t(res.re.sl3[1:13, ])
se.re.sl3.n1400 <- t(res.re.sl3[14:21 , ])
colnames(res.re.sl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.sl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 9, ncol = k)
  se_mat <- matrix(NA, nrow = 9, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glmer(death~(age+factor(sex)+sbp)*factor(treat)+trial+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:6]
    lp1 <- X1 %*% coef[1:6] + X1[,2:4] %*% coef[7:9] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[6, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ]),
    coef_mat.8 = mean(coef_mat[9, ]),
    se_mat.1 = mean(se_mat[1, ]) + mean(se_mat[6, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[7, ]),
    se_mat.7 = mean(se_mat[8, ]),
    se_mat.8 = mean(se_mat[9, ]))
}
stopCluster(cl)
res.si.sl3.n1400 <- t(res.si.sl3[1:13, ])
se.si.sl3.n1400 <- t(res.si.sl3[14:21, ])
colnames(res.si.sl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.sl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 9, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- rrvglm(death ~ (age+factor(sex)+sbp)*factor(treat)+trial, family = binomialff, data = train,Rank = 1)
    coef <- mod@coefficients
    coef_mat[,i] <- coef
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:6]
    lp1 <- X1 %*% coef[1:6] + X1[,2:4] %*% coef[7:9] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[6, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ]),
    coef_mat.8 = mean(coef_mat[9, ]))
}
stopCluster(cl)
res.r1.sl3.n1400 <- t(res.r1.sl3[1:13, ])
colnames(res.r1.sl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat <- matrix(NA, nrow = 28, ncol = k)
  se_mat <- matrix(NA, nrow = 28, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glm(death~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial),
               data = train,family = binomial,control = glm.control(maxit = 100000))
    coef <- mod$coefficients
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial)"),
                      data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~trial+age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~trial+age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- coef[1]+mean(coef[2:6])*X0[,2]+mean(c(coef[7],coef[14:18]))*X0[,3]+mean(c(coef[8],coef[19:23]))*X0[,4]+
      mean(c(coef[9],coef[24:28]))*X0[,5]+coef[10]*X0[,6]
    lp1 <- coef[1]+mean(coef[2:6])*X1[,2]+mean(c(coef[7],coef[14:18]))*X1[,3]+mean(c(coef[8],coef[19:23]))*X1[,4]+
      mean(c(coef[9],coef[24:28]))*X1[,5]+coef[10]*X1[,6]+coef[11]*X1[,3]+coef[12]*X1[,4]+coef[13]*X1[,5]
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
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
    coef_mat.1 = mean(coef_mat[1:6,]),
    coef_mat.2 = mean(c(coef_mat[7,]+coef_mat[14:18,])),
    coef_mat.3 = mean(c(coef_mat[8,]+coef_mat[19:23,])),
    coef_mat.4 = mean(c(coef_mat[9,]+coef_mat[24:28,])),
    coef_mat.5 = mean(coef_mat[10,]),
    coef_mat.6 = mean(coef_mat[11,]),
    coef_mat.7 = mean(coef_mat[12,]),
    coef_mat.8 = mean(coef_mat[13,]),
    se_mat.1 = mean(se_mat[1:6,]),
    se_mat.2 = mean(c(se_mat[7,]+se_mat[14:18,])),
    se_mat.3 = mean(c(se_mat[8,]+se_mat[19:23,])),
    se_mat.4 = mean(c(se_mat[9,]+se_mat[24:28,])),
    se_mat.5 = mean(se_mat[10,]),
    se_mat.6 = mean(se_mat[11,]),
    se_mat.7 = mean(se_mat[12,]),
    se_mat.8 = mean(se_mat[13,])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.fs.sl3.n1400 <- t(res.fs.sl3[1:13, ])
se.fs.sl3.n1400 <- t(res.fs.sl3[14:21 , ])
colnames(res.fs.sl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.sl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.na),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.na.tl3.n1400 <- t(res.na.tl3[1:13, ])
se.na.tl3.n1400 <- t(res.na.tl3[14:21 , ])
colnames(res.na.tl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.tl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glmer(death~age+factor(sex)+sbp+(1|trial),
                  data = train0,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    mod1 <- glmer(death~age+factor(sex)+sbp+(1|trial),
                  data = train1,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef0 <- mod0@beta
    coef1 <- mod1@beta
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.re.tl3.n1400 <- t(res.re.tl3[1:13, ])
se.re.tl3.n1400 <- t(res.re.tl3[14:21, ])
colnames(res.re.tl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.tl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat0 <- matrix(NA, nrow = 5, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 5, ncol = k)
  se_mat0 <- matrix(NA, nrow = 5, ncol = k)
  se_mat1 <- matrix(NA, nrow = 5, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp+trial,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp+trial,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5, ])),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[5, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[5, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[5, ])),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.si.tl3.n1400 <- t(res.si.tl3[1:13, ])
se.si.tl3.n1400 <- t(res.si.tl3[14:21 , ])
colnames(res.si.tl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.tl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat0 <- matrix(NA, nrow = 5, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 5, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- rrvglm(death ~ age+factor(sex)+sbp+trial, family = binomialff, data = train0,Rank = 1)
    mod1 <- rrvglm(death ~ age+factor(sex)+sbp+trial, family = binomialff, data = train1,Rank = 1)
    coef0 <- mod0@coefficients
    coef1 <- mod1@coefficients
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5, ])),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]))
}
stopCluster(cl)
res.r1.tl3.n1400 <- t(res.r1.tl3[1:13, ])
colnames(res.r1.tl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat0 <- matrix(NA, nrow = 24, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 24, ncol = k)
  se_mat0 <- matrix(NA, nrow = 24, ncol = k)
  se_mat1 <- matrix(NA, nrow = 24, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~(age+factor(sex)+sbp)*factor(trial),
                data = train0,family = binomial,control = glm.control(maxit = 100000))
    mod1 <- glm(death~(age+factor(sex)+sbp)*factor(trial),
                data = train1,family = binomial,control = glm.control(maxit = 100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(trial)"), data=train0)
    X1 <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(trial)"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- coef0[1]+mean(coef0[5:8])*X[,5]+mean(c(coef0[2],coef0[9:12]))*X[,2]+
      mean(c(coef0[3],coef0[13:16]))*X[,3]+mean(c(coef0[4],coef0[17:20]))*X[,4]
    lp1 <- coef1[1]+mean(coef1[5:8])*X[,5]+mean(c(coef1[2],coef1[9:12]))*X[,2]+
      mean(c(coef1[3],coef1[13:16]))*X[,3]+mean(c(coef1[4],coef1[17:20]))*X[,4]
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
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
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5:9, ]),
    coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[10:14, ]),
    coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[15:19, ]),
    coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[20:24, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5:9, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5:9, ])),
    coef_mat.6 = (mean(coef_mat1[2, ]) + mean(coef_mat1[10:14, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[10:14, ])),
    coef_mat.7 = (mean(coef_mat1[3, ]) + mean(coef_mat1[15:19, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[15:19, ])),
    coef_mat.8 = (mean(coef_mat1[4, ]) + mean(coef_mat1[20:24, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[20:24, ])),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[5:9, ]),
    se_mat.2 = mean(se_mat0[2, ]) + mean(se_mat0[10:14, ]),
    se_mat.3 = mean(se_mat0[3, ]) + mean(se_mat0[15:19, ]),
    se_mat.4 = mean(se_mat0[4, ]) + mean(se_mat0[20:24, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[5:9, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[5:9, ])),
    se_mat.6 = (mean(se_mat1[2, ]) + mean(se_mat1[10:14, ])) - (mean(se_mat0[2, ]) + mean(se_mat0[10:14, ])),
    se_mat.7 = (mean(se_mat1[3, ]) + mean(se_mat1[15:19, ])) - (mean(se_mat0[3, ]) + mean(se_mat0[15:19, ])),
    se_mat.8 = (mean(se_mat1[4, ]) + mean(se_mat1[20:24, ])) - (mean(se_mat0[4, ]) + mean(se_mat0[20:24, ]))),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.fs.tl3.n1400 <- t(res.fs.tl3[1:13, ])
se.fs.tl3.n1400 <- t(res.fs.tl3[14:21, ])
colnames(res.fs.tl3.n1400) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.tl3.n1400) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

save(res.na.sl3.n1400,se.na.sl3.n1400,
     res.re.sl3.n1400,se.re.sl3.n1400,
     res.si.sl3.n1400,se.si.sl3.n1400,
     res.r1.sl3.n1400,
     res.fs.sl3.n1400,se.fs.sl3.n1400,
     res.na.tl3.n1400,se.na.tl3.n1400,
     res.re.tl3.n1400,se.re.tl3.n1400,
     res.si.tl3.n1400,se.si.tl3.n1400,
     res.r1.tl3.n1400,
     res.fs.tl3.n1400,se.fs.tl3.n1400,
     file = "res_scenario2.Rdata")


#### scenario 3: 3 covariates (1 binary and 2 continuous) & total sample size = 700 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glm(death~(age+factor(sex)+sbp)*factor(treat),
               data = train,family = binomial,control = glm.control(maxit = 100000))
    coef <- mod$coefficients
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- X0 %*% coef[1:5]
    lp1 <- X1 %*% coef[1:5] + X1[,2:4] %*% coef[6:8] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(oben,pben)
  }
  c(
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
    se_mat.1 = mean(se_mat[1, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[6, ]),
    se_mat.7 = mean(se_mat[7, ]),
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.na.sl3.n700 <- t(res.na.sl3[1:13, ])
se.na.sl3.n700 <- t(res.na.sl3[14:21, ])
colnames(res.na.sl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.sl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- glmer(death~(age+factor(sex)+sbp)*factor(treat)+(1|trial)+(0+treat|trial),
                 data = train, family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- X0 %*% coef[1:5]
    lp1 <- X1 %*% coef[1:5] + X1[,2:4] %*% coef[6:8] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  c(
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
    coef_mat.8 = mean(coef_mat[8,]),
    se_mat.1 = mean(se_mat[1,]),
    se_mat.2 = mean(se_mat[2,]),
    se_mat.3 = mean(se_mat[3,]),
    se_mat.4 = mean(se_mat[4,]),
    se_mat.5 = mean(se_mat[5,]),
    se_mat.6 = mean(se_mat[6,]),
    se_mat.7 = mean(se_mat[7,]),
    se_mat.8 = mean(se_mat[8,]))
}
stopCluster(cl)
res.re.sl3.n700 <- t(res.re.sl3[1:13, ])
se.re.sl3.n700 <- t(res.re.sl3[14:21 , ])
colnames(res.re.sl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.sl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 9, ncol = k)
  se_mat <- matrix(NA, nrow = 9, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glmer(death~(age+factor(sex)+sbp)*factor(treat)+trial+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:6]
    lp1 <- X1 %*% coef[1:6] + X1[,2:4] %*% coef[7:9] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[6, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ]),
    coef_mat.8 = mean(coef_mat[9, ]),
    se_mat.1 = mean(se_mat[1, ]) + mean(se_mat[6, ]),
    se_mat.2 = mean(se_mat[2, ]),
    se_mat.3 = mean(se_mat[3, ]),
    se_mat.4 = mean(se_mat[4, ]),
    se_mat.5 = mean(se_mat[5, ]),
    se_mat.6 = mean(se_mat[7, ]),
    se_mat.7 = mean(se_mat[8, ]),
    se_mat.8 = mean(se_mat[9, ]))
}
stopCluster(cl)
res.si.sl3.n700 <- t(res.si.sl3[1:13, ])
se.si.sl3.n700 <- t(res.si.sl3[14:21, ])
colnames(res.si.sl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.sl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 9, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- rrvglm(death ~ (age+factor(sex)+sbp)*factor(treat)+trial, family = binomialff, data = train,Rank = 1)
    coef <- mod@coefficients
    coef_mat[,i] <- coef
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:6]
    lp1 <- X1 %*% coef[1:6] + X1[,2:4] %*% coef[7:9] 
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[6, ]),
    coef_mat.2 = mean(coef_mat[2, ]),
    coef_mat.3 = mean(coef_mat[3, ]),
    coef_mat.4 = mean(coef_mat[4, ]),
    coef_mat.5 = mean(coef_mat[5, ]),
    coef_mat.6 = mean(coef_mat[7, ]),
    coef_mat.7 = mean(coef_mat[8, ]),
    coef_mat.8 = mean(coef_mat[9, ]))
}
stopCluster(cl)
res.r1.sl3.n700 <- t(res.r1.sl3[1:13, ])
colnames(res.r1.sl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat <- matrix(NA, nrow = 28, ncol = k)
  se_mat <- matrix(NA, nrow = 28, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train <- df3[df3$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glm(death~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial),
               data = train,family = binomial,control = glm.control(maxit = 100000))
    coef <- mod$coefficients
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~factor(trial)+(age+factor(sex)+sbp)*factor(treat)+(age+factor(sex)+sbp):factor(trial)"),
                      data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~trial+age+factor(sex)+sbp+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~trial+age+factor(sex)+sbp+treat"), data=test1)
    lp0 <- coef[1]+mean(coef[2:6])*X0[,2]+mean(c(coef[7],coef[14:18]))*X0[,3]+mean(c(coef[8],coef[19:23]))*X0[,4]+
      mean(c(coef[9],coef[24:28]))*X0[,5]+coef[10]*X0[,6]
    lp1 <- coef[1]+mean(coef[2:6])*X1[,2]+mean(c(coef[7],coef[14:18]))*X1[,3]+mean(c(coef[8],coef[19:23]))*X1[,4]+
      mean(c(coef[9],coef[24:28]))*X1[,5]+coef[10]*X1[,6]+coef[11]*X1[,3]+coef[12]*X1[,4]+coef[13]*X1[,5]
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
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
    coef_mat.1 = mean(coef_mat[1:6,]),
    coef_mat.2 = mean(c(coef_mat[7,]+coef_mat[14:18,])),
    coef_mat.3 = mean(c(coef_mat[8,]+coef_mat[19:23,])),
    coef_mat.4 = mean(c(coef_mat[9,]+coef_mat[24:28,])),
    coef_mat.5 = mean(coef_mat[10,]),
    coef_mat.6 = mean(coef_mat[11,]),
    coef_mat.7 = mean(coef_mat[12,]),
    coef_mat.8 = mean(coef_mat[13,]),
    se_mat.1 = mean(se_mat[1:6,]),
    se_mat.2 = mean(c(se_mat[7,]+se_mat[14:18,])),
    se_mat.3 = mean(c(se_mat[8,]+se_mat[19:23,])),
    se_mat.4 = mean(c(se_mat[9,]+se_mat[24:28,])),
    se_mat.5 = mean(se_mat[10,]),
    se_mat.6 = mean(se_mat[11,]),
    se_mat.7 = mean(se_mat[12,]),
    se_mat.8 = mean(se_mat[13,])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.fs.sl3.n700 <- t(res.fs.sl3[1:13, ])
se.fs.sl3.n700 <- t(res.fs.sl3[14:21 , ])
colnames(res.fs.sl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.sl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.na),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.na.tl3.n700 <- t(res.na.tl3[1:13, ])
se.na.tl3.n700 <- t(res.na.tl3[14:21 , ])
colnames(res.na.tl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.tl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glmer(death~age+factor(sex)+sbp+(1|trial),
                  data = train0,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    mod1 <- glmer(death~age+factor(sex)+sbp+(1|trial),
                  data = train1,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef0 <- mod0@beta
    coef1 <- mod1@beta
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(oben,pben)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    coef_mat.1 = mean(coef_mat0[1, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.re.tl3.n700 <- t(res.re.tl3[1:13, ])
se.re.tl3.n700 <- t(res.re.tl3[14:21, ])
colnames(res.re.tl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.tl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat0 <- matrix(NA, nrow = 5, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 5, ncol = k)
  se_mat0 <- matrix(NA, nrow = 5, ncol = k)
  se_mat1 <- matrix(NA, nrow = 5, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp+trial,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp+trial,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5, ])),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[5, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[5, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[5, ])),
    se_mat.6 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.7 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]))
}
stopCluster(cl)
res.si.tl3.n700 <- t(res.si.tl3[1:13, ])
se.si.tl3.n700 <- t(res.si.tl3[14:21 , ])
colnames(res.si.tl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.tl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl3 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat0 <- matrix(NA, nrow = 5, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 5, ncol = k)
  
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- rrvglm(death ~ age+factor(sex)+sbp+trial, family = binomialff, data = train0,Rank = 1)
    mod1 <- rrvglm(death ~ age+factor(sex)+sbp+trial, family = binomialff, data = train1,Rank = 1)
    coef0 <- mod0@coefficients
    coef1 <- mod1@coefficients
    coef_mat0[,i] <- coef0
    coef_mat1[,i] <- coef1
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- X %*% coef0
    lp1 <- X %*% coef1
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
    lm_ <- lm(oben~pben)
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(oben,pben)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5, ])),
    coef_mat.6 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.7 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.8 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]))
}
stopCluster(cl)
res.r1.tl3.n700 <- t(res.r1.tl3[1:13, ])
colnames(res.r1.tl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl3 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- -0.3
  b5 <- 0.01
  b6 <- 0.04
  b7 <- -0.2
  
  df3$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7) 
  
  log_odds <- with(df3, b0 + b1*age + b2*(sex==1)+ b3*sbp + b4*(treat==1) + trial_eff[trial] + (b5*age + b6*(sex==1) + b7*sbp)*(treat==1))
  p <- plogis(log_odds)
  df3$death <- rbinom(n, 1, p)
  
  df0 <- df3[df3$treat==0,]
  df1 <- df3[df3$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  coef_mat0 <- matrix(NA, nrow = 24, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 24, ncol = k)
  se_mat0 <- matrix(NA, nrow = 24, ncol = k)
  se_mat1 <- matrix(NA, nrow = 24, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df3[df3$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~(age+factor(sex)+sbp)*factor(trial),
                data = train0,family = binomial,control = glm.control(maxit = 100000))
    mod1 <- glm(death~(age+factor(sex)+sbp)*factor(trial),
                data = train1,family = binomial,control = glm.control(maxit = 100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(trial)"), data=train0)
    X1 <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(trial)"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+factor(sex)+sbp+trial"), data=test)
    lp0 <- coef0[1]+mean(coef0[5:8])*X[,5]+mean(c(coef0[2],coef0[9:12]))*X[,2]+
      mean(c(coef0[3],coef0[13:16]))*X[,3]+mean(c(coef0[4],coef0[17:20]))*X[,4]
    lp1 <- coef1[1]+mean(coef1[5:8])*X[,5]+mean(c(coef1[2],coef1[9:12]))*X[,2]+
      mean(c(coef1[3],coef1[13:16]))*X[,3]+mean(c(coef1[4],coef1[17:20]))*X[,4]
    ite <- expit(lp1) - expit(lp0)
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration plot
    predq5 <- quantcut(ite, q=5)
    pben <-  tapply(ite, predq5, mean)
    oben <- sapply(levels(predq5), function(x){lm(death ~ treat, data=test, subset=predq5==x)$coef[2]})
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
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[5:9, ]),
    coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[10:14, ]),
    coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[15:19, ]),
    coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[20:24, ]),
    coef_mat.5 = (mean(coef_mat1[1, ]) + mean(coef_mat1[5:9, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[5:9, ])),
    coef_mat.6 = (mean(coef_mat1[2, ]) + mean(coef_mat1[10:14, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[10:14, ])),
    coef_mat.7 = (mean(coef_mat1[3, ]) + mean(coef_mat1[15:19, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[15:19, ])),
    coef_mat.8 = (mean(coef_mat1[4, ]) + mean(coef_mat1[20:24, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[20:24, ])),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[5:9, ]),
    se_mat.2 = mean(se_mat0[2, ]) + mean(se_mat0[10:14, ]),
    se_mat.3 = mean(se_mat0[3, ]) + mean(se_mat0[15:19, ]),
    se_mat.4 = mean(se_mat0[4, ]) + mean(se_mat0[20:24, ]),
    se_mat.5 = (mean(se_mat1[1, ]) + mean(se_mat1[5:9, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[5:9, ])),
    se_mat.6 = (mean(se_mat1[2, ]) + mean(se_mat1[10:14, ])) - (mean(se_mat0[2, ]) + mean(se_mat0[10:14, ])),
    se_mat.7 = (mean(se_mat1[3, ]) + mean(se_mat1[15:19, ])) - (mean(se_mat0[3, ]) + mean(se_mat0[15:19, ])),
    se_mat.8 = (mean(se_mat1[4, ]) + mean(se_mat1[20:24, ])) - (mean(se_mat0[4, ]) + mean(se_mat0[20:24, ]))),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.fs.tl3.n700 <- t(res.fs.tl3[1:13, ])
se.fs.tl3.n700 <- t(res.fs.tl3[14:21, ])
colnames(res.fs.tl3.n700) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.tl3.n700) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

save(res.na.sl3.n700,se.na.sl3.n700,
     res.re.sl3.n700,se.re.sl3.n700,
     res.si.sl3.n700,se.si.sl3.n700,
     res.r1.sl3.n700,
     res.fs.sl3.n700,se.fs.sl3.n700,
     res.na.tl3.n700,se.na.tl3.n700,
     res.re.tl3.n700,se.re.tl3.n700,
     res.si.tl3.n700,se.si.tl3.n700,
     res.r1.tl3.n700,
     res.fs.tl3.n700,se.fs.tl3.n700,
     file = "res_scenario3.Rdata")
