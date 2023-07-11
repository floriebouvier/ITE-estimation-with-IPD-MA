library(tidyverse) # data.frame manipulation
library(Hmisc) # cstat computation
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

#### scenario 4: 10 covariates (6 binary and 4 continuous) & total sample size = 2800 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat <- matrix(NA, nrow = 22, ncol = k)
                        se_mat <- matrix(NA, nrow = 22, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                            se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)"), data=train)
                            lp <- X %*% coef
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
                            X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
                            lp0 <- X0 %*% coef[1:12]
                            lp1 <- X1 %*% coef[1:12] + X1[,2:11] %*% coef[13:22] 
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
                            coef_mat.22 = mean(coef_mat[22, ]),
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
                            se_mat.21 = mean(se_mat[21, ]),
                            se_mat.22 = mean(se_mat[22, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.sl10.n2800 <- t(res.na.sl10[1:27, ])
se.na.sl10.n2800 <- t(res.na.sl10[28:49, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 22, ncol = k)
  se_mat <- matrix(NA, nrow = 22, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train <- df10[df10$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+(1|trial)+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
    lp0 <- X0 %*% coef[1:12]
    lp1 <- X1 %*% coef[1:12] + X1[,2:11] %*% coef[13:22] 
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
    coef_mat.22 = mean(coef_mat[22, ]),
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
    se_mat.21 = mean(se_mat[21, ]),
    se_mat.22 = mean(se_mat[22, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.sl10.n2800 <- t(res.re.sl10[1:27, ])
se.re.sl10.n2800 <- t(res.re.sl10[28:49 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 23, ncol = k)
  se_mat <- matrix(NA, nrow = 23, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train <- df10[df10$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:13]
    lp1 <- X1 %*% coef[1:13] + X1[,2:11] %*% coef[14:23] 
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[13, ]),
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
    coef_mat.21 = mean(coef_mat[22, ]),
    coef_mat.22 = mean(coef_mat[23, ]),
    se_mat.1 = mean(se_mat[1, ]) + mean(se_mat[13, ]),
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
    se_mat.13 = mean(se_mat[14, ]),
    se_mat.14 = mean(se_mat[15, ]),
    se_mat.15 = mean(se_mat[16, ]),
    se_mat.16 = mean(se_mat[17, ]),
    se_mat.17 = mean(se_mat[18, ]),
    se_mat.18 = mean(se_mat[19, ]),
    se_mat.19 = mean(se_mat[20, ]),
    se_mat.20 = mean(se_mat[21, ]),
    se_mat.21 = mean(se_mat[22, ]),
    se_mat.22 = mean(se_mat[23, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.sl10.n2800 <- t(res.si.sl10[1:27, ])
se.si.sl10.n2800 <- t(res.si.sl10[28:49, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")
) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 23, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train <- df10[df10$trial!=i,]
      
      test0 <- test
      test0$treat <- 0
      test1 <- test
      test1$treat <- 1
      
      #applying model to train
      mod <- rrvglm(death ~ (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+
                               chol)*factor(treat)+trial, family = binomialff, data = train,Rank = 1)
      coef <- mod@coefficients
      coef_mat[,i] <- coef
      
      #recalibration of event rates adapted to train
      X <- model.matrix(as.formula("death ~ (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial"), data=train)
      lp <- X %*% coef
      coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
      
      #ite estimation
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test1)
      lp0 <- X0 %*% coef[1:13]
      lp1 <- X1 %*% coef[1:13] + X1[,2:11] %*% coef[14:23] 
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[13, ]),
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
    coef_mat.21 = mean(coef_mat[22, ]),
    coef_mat.22 = mean(coef_mat[23, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl10.n2800 <- t(res.r1.sl10[1:27, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar% {
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
    
    b0 <- -1.4
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- 0.82
    b5 <- 0.8
    b6 <- 0.7
    b7 <- 0.1
    b8 <- 0.33
    b9 <- -0.02
    b10 <- 0.06
    b11 <- -0.3
    b12 <- 0.01
    b13 <- 0.04
    b14 <- -0.2
    b15 <- -0.01
    
    df10$treat <- rep(c(0,1), times = n/(2*k))
    
    trial_eff <- rep(0,7)
    
    log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
    p <- plogis(log_odds)
    df10$death <- rbinom(n, 1, p)
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat <- matrix(NA, nrow = 77, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train <- df10[df10$trial!=i,]
        
        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1
        
        #applying model to train
        mod <- glm(death~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+
                                          height+chol)*factor(treat)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+
                                                                        lvhn1+height+chol):factor(trial),
                   data = train,family = binomial,control = glm.control(maxit = 100000))
        coef <- mod$coefficients
        coef[is.na(coef)] <- 0
        coef_mat[,i] <- coef
        
        #recalibration of event rates adapted to train
        X <- model.matrix(as.formula("death~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol):factor(trial)"),
                          data=train)
        lp <- X %*% coef
        coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
        
        #ite estimation
        X0 <- model.matrix(as.formula("death~trial+age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
        X1 <- model.matrix(as.formula("death~trial+age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
        lp0 <- coef[1]+mean(coef[2:6])*X0[,2]+mean(c(coef[7],coef[28:32]))*X0[,3]+mean(c(coef[8],coef[33:37]))*X0[,4]+
          mean(c(coef[9],coef[38:42]))*X0[,5]+mean(c(coef[9],coef[38:42]))*X0[,5]+mean(c(coef[9],coef[38:42]))*X0[,5]+
          mean(c(coef[10],coef[43:47]))*X0[,6]+mean(c(coef[11],coef[48:52]))*X0[,7]+mean(c(coef[12],coef[53:57]))*X0[,8]+
          mean(c(coef[13],coef[58:62]))*X0[,9]+mean(c(coef[14],coef[63:67]))*X0[,10]+mean(c(coef[15],coef[68:72]))*X0[,11]+
          mean(c(coef[16],coef[73:77]))*X0[,12]+coef[17]*X0[,13]
        lp1 <- coef[1]+mean(coef[2:6])*X1[,2]+mean(c(coef[7],coef[28:32]))*X1[,3]+mean(c(coef[8],coef[33:37]))*X1[,4]+
          mean(c(coef[9],coef[38:42]))*X1[,5]+mean(c(coef[9],coef[38:42]))*X1[,5]+mean(c(coef[9],coef[38:42]))*X1[,5]+
          mean(c(coef[10],coef[43:47]))*X1[,6]+mean(c(coef[11],coef[48:52]))*X1[,7]+mean(c(coef[12],coef[53:57]))*X1[,8]+
          mean(c(coef[13],coef[58:62]))*X1[,9]+mean(c(coef[14],coef[63:67]))*X1[,10]+mean(c(coef[15],coef[68:72]))*X1[,11]+
          mean(c(coef[16],coef[73:77]))*X1[,12]+coef[17]*X1[,13]+coef[18]*X1[,3]+coef[19]*X1[,4]+coef[20]*X1[,5]+
          coef[21]*X1[,6]+coef[22]*X1[,7]+coef[23]*X1[,8]+coef[24]*X1[,9]+coef[25]*X1[,10]+coef[26]*X1[,11]+coef[27]*X1[,12]
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
    
    list(res = c(
      c.ben = mean(c.ben),
      c.ben.se = mean(c.ben.se),
      a = mean(a),
      b = mean(b),
      mse = mean(mse.fs),
      coef_mat.1 = mean(coef_mat[1:6, ]),
      coef_mat.2 = mean(c(coef_mat[7, ]+coef_mat[28:32, ])),
      coef_mat.3 = mean(c(coef_mat[8, ]+coef_mat[33:37, ])),
      coef_mat.4 = mean(c(coef_mat[9, ]+coef_mat[38:42, ])),
      coef_mat.5 = mean(c(coef_mat[10, ]+coef_mat[43:47, ])),
      coef_mat.6 = mean(c(coef_mat[11, ]+coef_mat[48:52, ])),
      coef_mat.7 = mean(c(coef_mat[12, ]+coef_mat[53:57, ])),
      coef_mat.8 = mean(c(coef_mat[13, ]+coef_mat[58:62, ])),
      coef_mat.9 = mean(c(coef_mat[14, ]+coef_mat[63:67, ])),
      coef_mat.10 = mean(c(coef_mat[15, ]+coef_mat[68:72, ])),
      coef_mat.11 = mean(c(coef_mat[16, ]+coef_mat[73:77, ])),
      coef_mat.12 = mean(c(coef_mat[17, ])),
      coef_mat.13 = mean(coef_mat[18, ]),
      coef_mat.14 = mean(coef_mat[19, ]),
      coef_mat.15 = mean(coef_mat[20, ]),
      coef_mat.16 = mean(coef_mat[21, ]),
      coef_mat.17 = mean(coef_mat[22, ]),
      coef_mat.18 = mean(c(coef_mat[23, ])),
      coef_mat.19 = mean(c(coef_mat[24, ])),
      coef_mat.20 = mean(c(coef_mat[25, ])),
      coef_mat.21 = mean(coef_mat[26, ]),
      coef_mat.22 = mean(coef_mat[27, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.fs.sl10.n2800 <- t(res.fs.sl10[1:27, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        df0 <- df10[df10$treat==0,]
                        df1 <- df10[df10$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 11, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 11, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol,
                                        data = train0,family = binomial,
                                        control = glm.control(maxit=100000))
                            mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol,
                                        data = train1,family = binomial,
                                        control = glm.control(maxit=100000))
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            coef0[is.na(coef0)] <- 0
                            coef1[is.na(coef1)] <- 0
                            coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
                            coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
                            se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
                            se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
                            
                            #recalibration of event rates adapted to train
                            X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train0)
                            X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train1)
                            lp0 <- X0 %*% coef0
                            lp1 <- X1 %*% coef1
                            
                            coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                            coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                            
                            #ite estimaiton
                            X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=test)
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
                            coef_mat.11 = mean(coef_mat0[11, ]),
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
                            coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
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
                            se_mat.11 = mean(se_mat0[11, ]),
                            se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                            se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                            se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                            se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                            se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                            se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                            se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                            se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                            se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                            se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
                            se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.tl10.n2800 <- t(res.na.tl10[1:27, ])
se.na.tl10.n2800 <- t(res.na.tl10[28:49, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
  se_mat0 <- matrix(NA, nrow = 11, ncol = k)
  se_mat1 <- matrix(NA, nrow = 11, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),
                  data = train0,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    mod1 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),
                  data = train1,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef0 <- mod0@beta
    coef1 <- mod1@beta
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train1)
    lp0 <- X0[ ,1:length(coef0)] %*% coef0
    lp1 <- X1[ ,1:length(coef1)] %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=test)
    lp0 <- X[ ,1:length(coef0)] %*% coef0
    lp1 <- X[ ,1:length(coef1)] %*% coef1
    
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
    coef_mat.11 = mean(coef_mat0[11, ]),
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
    coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
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
    se_mat.11 = mean(se_mat0[11, ]),
    se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
    se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
    se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
    se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
    se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
    se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
    se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
    se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.tl10.n2800 <- t(res.re.tl10[1:27, ])
se.re.tl10.n2800 <- t(res.re.tl10[28:49 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat0 <- matrix(NA, nrow = 12, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 12, ncol = k)
  se_mat0 <- matrix(NA, nrow = 12, ncol = k)
  se_mat1 <- matrix(NA, nrow = 12, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+trial,
                data = train0,family = binomial,
                control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+trial,
                data = train1,family = binomial,
                control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef0[is.na(coef0)] <- 0
    coef1[is.na(coef1)] <- 0
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=test)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat0[5, ]),
    coef_mat.6 = mean(coef_mat0[6, ]),
    coef_mat.7 = mean(coef_mat0[7, ]),
    coef_mat.8 = mean(coef_mat0[8, ]),
    coef_mat.9 = mean(coef_mat0[9, ]),
    coef_mat.10 = mean(coef_mat0[10, ]),
    coef_mat.11 = mean(coef_mat0[11, ]),
    coef_mat.12 = (mean(coef_mat1[1, ]) + mean(coef_mat1[12, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12, ])),
    coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
    coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
    coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
    coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
    coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
    coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
    coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[12, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat0[5, ]),
    se_mat.6 = mean(se_mat0[6, ]),
    se_mat.7 = mean(se_mat0[7, ]),
    se_mat.8 = mean(se_mat0[8, ]),
    se_mat.9 = mean(se_mat0[9, ]),
    se_mat.10 = mean(se_mat0[10, ]),
    se_mat.11 = mean(se_mat0[11, ]),
    se_mat.12 = (mean(se_mat1[1, ]) + mean(se_mat1[12, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[12, ])),
    se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
    se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
    se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
    se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
    se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
    se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
    se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
    se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.tl10.n2800 <- t(res.si.tl10[1:27, ])
se.si.tl10.n2800 <- t(res.si.tl10[28:49, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")
) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat0 <- matrix(NA, nrow = 12, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 12, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- rrvglm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height+chol+trial, family = binomialff, data = train0,Rank = 1)
      mod1 <- rrvglm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height+chol+trial, family = binomialff, data = train1,Rank = 1)
      
      coef0 <- mod0@coefficients
      coef1 <- mod1@coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      coef_mat0[,i] <- coef0
      coef_mat1[,i] <- coef1
      
      #recalibration of event rates adapted to train
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train1)
      lp0 <- X0 %*% coef0
      lp1 <- X1 %*% coef1
      
      coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
      coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
      
      #ite estimation
      X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=test)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(
    res = c(c.ben = mean(c.ben),
            c.ben.se = mean(c.ben.se),
            a = mean(a),
            b = mean(b),
            mse = mean(mse.r1),
            coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12, ]),
            coef_mat.2 = mean(coef_mat0[2, ]),
            coef_mat.3 = mean(coef_mat0[3, ]),
            coef_mat.4 = mean(coef_mat0[4, ]),
            coef_mat.5 = mean(coef_mat0[5, ]),
            coef_mat.6 = mean(coef_mat0[6, ]),
            coef_mat.7 = mean(coef_mat0[7, ]),
            coef_mat.8 = mean(coef_mat0[8, ]),
            coef_mat.9 = mean(coef_mat0[9, ]),
            coef_mat.10 = mean(coef_mat0[10, ]),
            coef_mat.11 = mean(coef_mat0[11, ]),
            coef_mat.12 = (mean(coef_mat1[1, ]) + mean(coef_mat1[12, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12, ])),
            coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
            coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
            coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
            coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
            coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
            coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
            coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
            coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
            coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
            coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.tl10.n2800 <- t(res.r1.tl10[1:27, ])

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
  .packages = c("tidyverse", "Hmisc", "lme4", "magrittr", "gtools", "VGAM")) %dopar%
  {
    set.seed(j)
    n <- 2800
    k <- 7
    
    trial <- rep(1:k, each = floor(n / k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52, 56, 64, 70, 77, 78, 82) #Age
    sd_age <- c(4, 2, 1, 3, 4, 6, 2)
    df10$age <- round(rnorm(n, mean = mean_age[trial], sd = sd_age[trial]), 0)
    df10$age <- df10$age - 70
    
    pman <- c(0.8, 0.4, 0.5, 0.6, 0.5, 0.7, 0.5) #Sex
    df10$sex <- rbinom(n, 1, prob = pman[trial])
    
    mean_sbp <- c(186, 182, 170, 185, 190, 188, 197) #SBP
    sd_sbp <- c(9, 11, 5, 12, 9, 10, 16)
    df10$sbp <- round(rnorm(n, mean = mean_sbp[trial], sd = sd_sbp[trial]), 0)
    df10$sbp <- df10$sbp - 185
    
    pmi <- c(0.1, 0.005, 0.01, 0.02, 0.05, 0.01, 0.04) #Previous myocardial infarction
    df10$mi <- rbinom(n, 1, prob = pmi[trial])
    
    pstroke <- c(0.002, 0.06, 0.02, 0.02, 0.001, 0.008, 0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob = pstroke[trial])
    
    psmoke <- c(0.5, 0.2, 0.3, 0.4, 0.3, 0.25, 0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob = psmoke[trial])
    
    pdiab <- c(0.03, 0.001, 0.002, 0.07, 0.003, 0.01, 0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob = pdiab[trial])
    
    plvhn1 <- c(0.13, 0.11, 0.05, 0.25, 0.05, 0.06, 0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob = plvhn1[trial])
    
    mean_height <- c(176, 162, 167, 169, 168, 170, 167) #Height
    sd_height <- c(6, 9, 10, 10, 10, 9, 9)
    df10$height <- round(rnorm(n, mean = mean_height[trial], sd = sd_height[trial]), 0)
    df10$height <- df10$height - 168
    
    mean_chol <- c(6.6, 6.5, 6.4, 6.1, 6.4, 6, 6.4) #Cholesterol
    sd_chol <- c(0.01, 0.015, 0.011, 0.012, 0.012, 0.012, 0.01)
    df10$chol <- round(rnorm(n, mean = mean_chol[trial], sd = sd_chol[trial]), 0)
    df10$chol <- df10$chol - 6.3
    
    b0 <- -1.4
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- 0.82
    b5 <- 0.8
    b6 <- 0.7
    b7 <- 0.1
    b8 <- 0.33
    b9 <- -0.02
    b10 <- 0.06
    b11 <- -0.3
    b12 <- 0.01
    b13 <- 0.04
    b14 <- -0.2
    b15 <- -0.01
    
    df10$treat <- rep(c(0, 1), times = n / (2 * k))
    
    trial_eff <- rep(0, 7)
    
    log_odds <- with(df10,b0 + b1 * age + b2 * (sex == 1) + b3 * sbp + b4 * (mi == 1) + 
                       b5 * (stroke == 1) + b6 * (smoke == 1) + b7 * (diab == 1) + b8 * (lvhn1 == 1) +
                       b9 * height + b10 * chol + b11 * (treat == 1) + trial_eff[trial] +
                       (b12 * age + b13 * (sex == 1) + b14 * sbp + b15 * (smoke == 1)) * (treat == 1))
    p <- plogis(log_odds)
    df10$death <- rbinom(n, 1, p)
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat0 <- matrix(NA, nrow = 66, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 66, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- glm(death ~ (age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial),
                    data = train0,family = binomial,control = glm.control(maxit = 100000))
        mod1 <- glm(death ~ (age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial),
                    data = train1,family = binomial,control = glm.control(maxit = 100000))
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
        #recalibration of event rates adapted to train
        X0 <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial)"),
                           data = train0)
        X1 <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial)"),
                           data = train1)
        lp0 <- X0 %*% coef0
        lp1 <- X1 %*% coef1
        
        coef0[1] <- coef0[1] + glm(death ~ offset(lp0), data = train0, family = "binomial")$coef[1]
        coef1[1] <- coef1[1] + glm(death ~ offset(lp1), data = train1, family = "binomial")$coef[1]
        
        #ite estimation
        X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"),
                          data = test)
        lp0 <- coef0[1] + mean(coef0[12:15])*X[, 12] + mean(c(coef0[2],coef0[16:19]))*X[, 2] +
          mean(c(coef0[3],coef0[20:23]))*X[, 3] + mean(c(coef0[4],coef0[24:27]))*X[, 4] +
          mean(c(coef0[5],coef0[28:31]))*X[, 5] + mean(c(coef0[6],coef0[32:35]))*X[, 6] +
          mean(c(coef0[7],coef0[36:39]))*X[, 7] + mean(c(coef0[8],coef0[40:43]))*X[, 8] +
          mean(c(coef0[9],coef0[44:47]))*X[, 9] + mean(c(coef0[10],coef0[48:51]))*X[, 10] +
          mean(c(coef0[11],coef0[52:55]))*X[, 11]
        lp1 <- coef1[1] + mean(coef1[12:15])*X[, 12] + mean(c(coef1[2],coef1[16:19]))*X[, 2] +
          mean(c(coef1[3],coef1[20:23]))*X[, 3] + mean(c(coef1[4],coef1[24:27]))*X[, 4] +
          mean(c(coef1[5],coef1[28:31]))*X[, 5] + mean(c(coef1[6],coef1[32:35]))*X[, 6] +
          mean(c(coef1[7],coef1[36:39]))*X[, 7] + mean(c(coef1[8],coef1[40:43]))*X[, 8] +
          mean(c(coef1[9],coef1[44:47]))*X[, 9] + mean(c(coef1[10],coef1[48:51]))*X[, 10] +
          mean(c(coef1[11],coef1[52:55]))*X[, 11]
        ite <- expit(lp1) - expit(lp0)
        
        #c-statistic for benefit
        cstat = cstat4ben(outcome = as.numeric(test$death),
                          treatment = test$treat == 1,
                          score = ite)
        
        c.ben[i] <- cstat[1]
        c.ben.se[i] <- cstat[2]
        
        #calibration plot
        predq5 <- quantcut(ite, q = 5)
        pben <- tapply(ite, predq5, mean)
        oben <- sapply(levels(predq5), function(x) {
          lm(death ~ treat, data = test, subset = predq5 == x)$coef[2]
        })
        lm_ <- lm(oben ~ pben)
        
        a[i] <- lm_$coefficients[[1]]
        b[i] <- lm_$coefficients[[2]]
        mse.fs[i] <- mse(oben, pben)
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
        coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12:16, ]),
        coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[17:21, ]),
        coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[22:26, ]),
        coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[27:31, ]),
        coef_mat.5 = mean(coef_mat0[5, ]) + mean(coef_mat0[32:36, ]),
        coef_mat.6 = mean(coef_mat0[6, ]) + mean(coef_mat0[37:41, ]),
        coef_mat.7 = mean(coef_mat0[7, ]) + mean(coef_mat0[42:46, ]),
        coef_mat.8 = mean(coef_mat0[8, ]) + mean(coef_mat0[47:51, ]),
        coef_mat.9 = mean(coef_mat0[9, ]) + mean(coef_mat0[52:56, ]),
        coef_mat.10 = mean(coef_mat0[10, ]) + mean(coef_mat0[57:61, ]),
        coef_mat.11 = mean(coef_mat0[11, ]) + mean(coef_mat0[62:66, ]),
        coef_mat.12 = (mean(coef_mat1[1]) + mean(coef_mat1[12:16, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12:16, ])),
        coef_mat.13 = (mean(coef_mat1[2]) + mean(coef_mat1[17:21, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[17:21, ])),
        coef_mat.14 = (mean(coef_mat1[3]) + mean(coef_mat1[22:26, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[22:26, ])),
        coef_mat.15 = (mean(coef_mat1[4]) + mean(coef_mat1[27:31, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[27:31, ])),
        coef_mat.16 = (mean(coef_mat1[5]) + mean(coef_mat1[32:36, ])) - (mean(coef_mat0[5, ]) + mean(coef_mat0[32:36, ])),
        coef_mat.17 = (mean(coef_mat1[6]) + mean(coef_mat1[37:41, ])) - (mean(coef_mat0[6, ]) + mean(coef_mat0[37:41, ])),
        coef_mat.18 = (mean(coef_mat1[7]) + mean(coef_mat1[42:46, ])) - (mean(coef_mat0[7, ]) + mean(coef_mat0[42:46, ])),
        coef_mat.19 = (mean(coef_mat1[8]) + mean(coef_mat1[47:51, ])) - (mean(coef_mat0[8, ]) + mean(coef_mat0[47:51, ])),
        coef_mat.20 = (mean(coef_mat1[9]) + mean(coef_mat1[52:56, ])) - (mean(coef_mat0[9, ]) + mean(coef_mat0[52:56, ])),
        coef_mat.21 = (mean(coef_mat1[10]) + mean(coef_mat1[57:61, ])) - (mean(coef_mat0[10, ]) + mean(coef_mat0[57:61, ])),
        coef_mat.22 = (mean(coef_mat1[11]) + mean(coef_mat1[62:66, ])) - (mean(coef_mat0[11, ]) + mean(coef_mat0[62:66, ]))),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl10.n2800 <- t(res.fs.tl10[1:27, ])

save(res.na.sl10.n2800,se.na.sl10.n2800,
     res.re.sl10.n2800,se.re.sl10.n2800,
     res.si.sl10.n2800,se.si.sl10.n2800,
     res.r1.sl10.n2800,
     res.fs.sl10.n2800,
     res.na.tl10.n2800,se.na.tl10.n2800,
     res.re.tl10.n2800,se.re.tl10.n2800,
     res.si.tl10.n2800,se.si.tl10.n2800,
     res.r1.tl10.n2800,
     res.fs.tl10.n2800,
     file = "res_scenario4.Rdata")

#### scenario 5: 10 covariates (6 binary and 4 continuous) & total sample size = 1400 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat <- matrix(NA, nrow = 22, ncol = k)
                        se_mat <- matrix(NA, nrow = 22, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                            se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)"), data=train)
                            lp <- X %*% coef
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
                            X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
                            lp0 <- X0 %*% coef[1:12]
                            lp1 <- X1 %*% coef[1:12] + X1[,2:11] %*% coef[13:22] 
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
                            coef_mat.22 = mean(coef_mat[22, ]),
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
                            se_mat.21 = mean(se_mat[21, ]),
                            se_mat.22 = mean(se_mat[22, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.sl10.n1400 <- t(res.na.sl10[1:27, ])
se.na.sl10.n1400 <- t(res.na.sl10[28:49, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 22, ncol = k)
  se_mat <- matrix(NA, nrow = 22, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train <- df10[df10$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+(1|trial)+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
    lp0 <- X0 %*% coef[1:12]
    lp1 <- X1 %*% coef[1:12] + X1[,2:11] %*% coef[13:22] 
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
    coef_mat.22 = mean(coef_mat[22, ]),
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
    se_mat.21 = mean(se_mat[21, ]),
    se_mat.22 = mean(se_mat[22, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.sl10.n1400 <- t(res.re.sl10[1:27, ])
se.re.sl10.n1400 <- t(res.re.sl10[28:49 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 23, ncol = k)
  se_mat <- matrix(NA, nrow = 23, ncol = k)
  
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train <- df10[df10$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:13]
    lp1 <- X1 %*% coef[1:13] + X1[,2:11] %*% coef[14:23] 
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[13, ]),
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
    coef_mat.21 = mean(coef_mat[22, ]),
    coef_mat.22 = mean(coef_mat[23, ]),
    se_mat.1 = mean(se_mat[1, ]) + mean(se_mat[13, ]),
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
    se_mat.13 = mean(se_mat[14, ]),
    se_mat.14 = mean(se_mat[15, ]),
    se_mat.15 = mean(se_mat[16, ]),
    se_mat.16 = mean(se_mat[17, ]),
    se_mat.17 = mean(se_mat[18, ]),
    se_mat.18 = mean(se_mat[19, ]),
    se_mat.19 = mean(se_mat[20, ]),
    se_mat.20 = mean(se_mat[21, ]),
    se_mat.21 = mean(se_mat[22, ]),
    se_mat.22 = mean(se_mat[23, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.sl10.n1400 <- t(res.si.sl10[1:27, ])
se.si.sl10.n1400 <- t(res.si.sl10[28:49, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")
) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 23, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train <- df10[df10$trial!=i,]
      
      test0 <- test
      test0$treat <- 0
      test1 <- test
      test1$treat <- 1
      
      #applying model to train
      mod <- rrvglm(death ~ (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+
                               chol)*factor(treat)+trial, family = binomialff, data = train,Rank = 1)
      coef <- mod@coefficients
      coef_mat[,i] <- coef
      
      #recalibration of event rates adapted to train
      X <- model.matrix(as.formula("death ~ (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial"), data=train)
      lp <- X %*% coef
      coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
      
      #ite estimation
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test1)
      lp0 <- X0 %*% coef[1:13]
      lp1 <- X1 %*% coef[1:13] + X1[,2:11] %*% coef[14:23] 
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[13, ]),
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
    coef_mat.21 = mean(coef_mat[22, ]),
    coef_mat.22 = mean(coef_mat[23, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl10.n1400 <- t(res.r1.sl10[1:27, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        coef_mat <- matrix(NA, nrow = 77, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+
                                                              height+chol)*factor(treat)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+
                                                                                            lvhn1+height+chol):factor(trial),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef[is.na(coef)] <- 0
                            coef_mat[,i] <- coef
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol):factor(trial)"),
                                              data=train)
                            lp <- X %*% coef
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            X0 <- model.matrix(as.formula("death~trial+age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
                            X1 <- model.matrix(as.formula("death~trial+age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
                            lp0 <- coef[1]+mean(coef[2:6])*X0[,2]+mean(c(coef[7],coef[28:32]))*X0[,3]+mean(c(coef[8],coef[33:37]))*X0[,4]+
                              mean(c(coef[9],coef[38:42]))*X0[,5]+mean(c(coef[9],coef[38:42]))*X0[,5]+mean(c(coef[9],coef[38:42]))*X0[,5]+
                              mean(c(coef[10],coef[43:47]))*X0[,6]+mean(c(coef[11],coef[48:52]))*X0[,7]+mean(c(coef[12],coef[53:57]))*X0[,8]+
                              mean(c(coef[13],coef[58:62]))*X0[,9]+mean(c(coef[14],coef[63:67]))*X0[,10]+mean(c(coef[15],coef[68:72]))*X0[,11]+
                              mean(c(coef[16],coef[73:77]))*X0[,12]+coef[17]*X0[,13]
                            lp1 <- coef[1]+mean(coef[2:6])*X1[,2]+mean(c(coef[7],coef[28:32]))*X1[,3]+mean(c(coef[8],coef[33:37]))*X1[,4]+
                              mean(c(coef[9],coef[38:42]))*X1[,5]+mean(c(coef[9],coef[38:42]))*X1[,5]+mean(c(coef[9],coef[38:42]))*X1[,5]+
                              mean(c(coef[10],coef[43:47]))*X1[,6]+mean(c(coef[11],coef[48:52]))*X1[,7]+mean(c(coef[12],coef[53:57]))*X1[,8]+
                              mean(c(coef[13],coef[58:62]))*X1[,9]+mean(c(coef[14],coef[63:67]))*X1[,10]+mean(c(coef[15],coef[68:72]))*X1[,11]+
                              mean(c(coef[16],coef[73:77]))*X1[,12]+coef[17]*X1[,13]+coef[18]*X1[,3]+coef[19]*X1[,4]+coef[20]*X1[,5]+
                              coef[21]*X1[,6]+coef[22]*X1[,7]+coef[23]*X1[,8]+coef[24]*X1[,9]+coef[25]*X1[,10]+coef[26]*X1[,11]+coef[27]*X1[,12]
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
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.fs),
                          coef_mat.1 = mean(coef_mat[1:6, ]),
                          coef_mat.2 = mean(c(coef_mat[7, ]+coef_mat[28:32, ])),
                          coef_mat.3 = mean(c(coef_mat[8, ]+coef_mat[33:37, ])),
                          coef_mat.4 = mean(c(coef_mat[9, ]+coef_mat[38:42, ])),
                          coef_mat.5 = mean(c(coef_mat[10, ]+coef_mat[43:47, ])),
                          coef_mat.6 = mean(c(coef_mat[11, ]+coef_mat[48:52, ])),
                          coef_mat.7 = mean(c(coef_mat[12, ]+coef_mat[53:57, ])),
                          coef_mat.8 = mean(c(coef_mat[13, ]+coef_mat[58:62, ])),
                          coef_mat.9 = mean(c(coef_mat[14, ]+coef_mat[63:67, ])),
                          coef_mat.10 = mean(c(coef_mat[15, ]+coef_mat[68:72, ])),
                          coef_mat.11 = mean(c(coef_mat[16, ]+coef_mat[73:77, ])),
                          coef_mat.12 = mean(c(coef_mat[17, ])),
                          coef_mat.13 = mean(coef_mat[18, ]),
                          coef_mat.14 = mean(coef_mat[19, ]),
                          coef_mat.15 = mean(coef_mat[20, ]),
                          coef_mat.16 = mean(coef_mat[21, ]),
                          coef_mat.17 = mean(coef_mat[22, ]),
                          coef_mat.18 = mean(c(coef_mat[23, ])),
                          coef_mat.19 = mean(c(coef_mat[24, ])),
                          coef_mat.20 = mean(c(coef_mat[25, ])),
                          coef_mat.21 = mean(coef_mat[26, ]),
                          coef_mat.22 = mean(coef_mat[27, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.fs.sl10.n1400 <- t(res.fs.sl10[1:27, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        df0 <- df10[df10$treat==0,]
                        df1 <- df10[df10$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 11, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 11, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol,
                                        data = train0,family = binomial,
                                        control = glm.control(maxit=100000))
                            mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol,
                                        data = train1,family = binomial,
                                        control = glm.control(maxit=100000))
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            coef0[is.na(coef0)] <- 0
                            coef1[is.na(coef1)] <- 0
                            coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
                            coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
                            se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
                            se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
                            
                            #recalibration of event rates adapted to train
                            X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train0)
                            X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train1)
                            lp0 <- X0 %*% coef0
                            lp1 <- X1 %*% coef1
                            
                            coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                            coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                            
                            #ite estimaiton
                            X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=test)
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
                            coef_mat.11 = mean(coef_mat0[11, ]),
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
                            coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
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
                            se_mat.11 = mean(se_mat0[11, ]),
                            se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                            se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                            se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                            se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                            se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                            se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                            se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                            se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                            se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                            se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
                            se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.tl10.n1400 <- t(res.na.tl10[1:27, ])
se.na.tl10.n1400 <- t(res.na.tl10[28:49, ])


#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
  se_mat0 <- matrix(NA, nrow = 11, ncol = k)
  se_mat1 <- matrix(NA, nrow = 11, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),
                  data = train0,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    mod1 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),
                  data = train1,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef0 <- mod0@beta
    coef1 <- mod1@beta
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train1)
    lp0 <- X0[ ,1:length(coef0)] %*% coef0
    lp1 <- X1[ ,1:length(coef1)] %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=test)
    lp0 <- X[ ,1:length(coef0)] %*% coef0
    lp1 <- X[ ,1:length(coef1)] %*% coef1
    
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
    coef_mat.11 = mean(coef_mat0[11, ]),
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
    coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
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
    se_mat.11 = mean(se_mat0[11, ]),
    se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
    se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
    se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
    se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
    se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
    se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
    se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
    se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.tl10.n1400 <- t(res.re.tl10[1:27, ])
se.re.tl10.n1400 <- t(res.re.tl10[28:49 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat0 <- matrix(NA, nrow = 12, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 12, ncol = k)
  se_mat0 <- matrix(NA, nrow = 12, ncol = k)
  se_mat1 <- matrix(NA, nrow = 12, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+trial,
                data = train0,family = binomial,
                control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+trial,
                data = train1,family = binomial,
                control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef0[is.na(coef0)] <- 0
    coef1[is.na(coef1)] <- 0
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=test)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat0[5, ]),
    coef_mat.6 = mean(coef_mat0[6, ]),
    coef_mat.7 = mean(coef_mat0[7, ]),
    coef_mat.8 = mean(coef_mat0[8, ]),
    coef_mat.9 = mean(coef_mat0[9, ]),
    coef_mat.10 = mean(coef_mat0[10, ]),
    coef_mat.11 = mean(coef_mat0[11, ]),
    coef_mat.12 = (mean(coef_mat1[1, ]) + mean(coef_mat1[12, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12, ])),
    coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
    coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
    coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
    coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
    coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
    coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
    coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[12, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat0[5, ]),
    se_mat.6 = mean(se_mat0[6, ]),
    se_mat.7 = mean(se_mat0[7, ]),
    se_mat.8 = mean(se_mat0[8, ]),
    se_mat.9 = mean(se_mat0[9, ]),
    se_mat.10 = mean(se_mat0[10, ]),
    se_mat.11 = mean(se_mat0[11, ]),
    se_mat.12 = (mean(se_mat1[1, ]) + mean(se_mat1[12, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[12, ])),
    se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
    se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
    se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
    se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
    se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
    se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
    se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
    se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.tl10.n1400 <- t(res.si.tl10[1:27, ])
se.si.tl10.n1400 <- t(res.si.tl10[28:49, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")
) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat0 <- matrix(NA, nrow = 12, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 12, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- rrvglm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height+chol+trial, family = binomialff, data = train0,Rank = 1)
      mod1 <- rrvglm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height+chol+trial, family = binomialff, data = train1,Rank = 1)
      
      coef0 <- mod0@coefficients
      coef1 <- mod1@coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      coef_mat0[,i] <- coef0
      coef_mat1[,i] <- coef1
      
      #recalibration of event rates adapted to train
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train1)
      lp0 <- X0 %*% coef0
      lp1 <- X1 %*% coef1
      
      coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
      coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
      
      #ite estimation
      X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=test)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(
    res = c(c.ben = mean(c.ben),
            c.ben.se = mean(c.ben.se),
            a = mean(a),
            b = mean(b),
            mse = mean(mse.r1),
            coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12, ]),
            coef_mat.2 = mean(coef_mat0[2, ]),
            coef_mat.3 = mean(coef_mat0[3, ]),
            coef_mat.4 = mean(coef_mat0[4, ]),
            coef_mat.5 = mean(coef_mat0[5, ]),
            coef_mat.6 = mean(coef_mat0[6, ]),
            coef_mat.7 = mean(coef_mat0[7, ]),
            coef_mat.8 = mean(coef_mat0[8, ]),
            coef_mat.9 = mean(coef_mat0[9, ]),
            coef_mat.10 = mean(coef_mat0[10, ]),
            coef_mat.11 = mean(coef_mat0[11, ]),
            coef_mat.12 = (mean(coef_mat1[1, ]) + mean(coef_mat1[12, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12, ])),
            coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
            coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
            coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
            coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
            coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
            coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
            coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
            coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
            coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
            coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.tl10.n1400 <- t(res.r1.tl10[1:27, ])

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
  .packages = c("tidyverse", "Hmisc", "lme4", "magrittr", "gtools", "VGAM")) %dopar%
  {
    set.seed(j)
    n <- 1400
    k <- 7
    
    trial <- rep(1:k, each = floor(n / k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52, 56, 64, 70, 77, 78, 82) #Age
    sd_age <- c(4, 2, 1, 3, 4, 6, 2)
    df10$age <- round(rnorm(n, mean = mean_age[trial], sd = sd_age[trial]), 0)
    df10$age <- df10$age - 70
    
    pman <- c(0.8, 0.4, 0.5, 0.6, 0.5, 0.7, 0.5) #Sex
    df10$sex <- rbinom(n, 1, prob = pman[trial])
    
    mean_sbp <- c(186, 182, 170, 185, 190, 188, 197) #SBP
    sd_sbp <- c(9, 11, 5, 12, 9, 10, 16)
    df10$sbp <- round(rnorm(n, mean = mean_sbp[trial], sd = sd_sbp[trial]), 0)
    df10$sbp <- df10$sbp - 185
    
    pmi <- c(0.1, 0.005, 0.01, 0.02, 0.05, 0.01, 0.04) #Previous myocardial infarction
    df10$mi <- rbinom(n, 1, prob = pmi[trial])
    
    pstroke <- c(0.002, 0.06, 0.02, 0.02, 0.001, 0.008, 0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob = pstroke[trial])
    
    psmoke <- c(0.5, 0.2, 0.3, 0.4, 0.3, 0.25, 0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob = psmoke[trial])
    
    pdiab <- c(0.03, 0.001, 0.002, 0.07, 0.003, 0.01, 0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob = pdiab[trial])
    
    plvhn1 <- c(0.13, 0.11, 0.05, 0.25, 0.05, 0.06, 0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob = plvhn1[trial])
    
    mean_height <- c(176, 162, 167, 169, 168, 170, 167) #Height
    sd_height <- c(6, 9, 10, 10, 10, 9, 9)
    df10$height <- round(rnorm(n, mean = mean_height[trial], sd = sd_height[trial]), 0)
    df10$height <- df10$height - 168
    
    mean_chol <- c(6.6, 6.5, 6.4, 6.1, 6.4, 6, 6.4) #Cholesterol
    sd_chol <- c(0.01, 0.015, 0.011, 0.012, 0.012, 0.012, 0.01)
    df10$chol <- round(rnorm(n, mean = mean_chol[trial], sd = sd_chol[trial]), 0)
    df10$chol <- df10$chol - 6.3
    
    b0 <- -1.4
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- 0.82
    b5 <- 0.8
    b6 <- 0.7
    b7 <- 0.1
    b8 <- 0.33
    b9 <- -0.02
    b10 <- 0.06
    b11 <- -0.3
    b12 <- 0.01
    b13 <- 0.04
    b14 <- -0.2
    b15 <- -0.01
    
    df10$treat <- rep(c(0, 1), times = n / (2 * k))
    
    trial_eff <- rep(0, 7)
    
    log_odds <- with(df10,b0 + b1 * age + b2 * (sex == 1) + b3 * sbp + b4 * (mi == 1) + 
                       b5 * (stroke == 1) + b6 * (smoke == 1) + b7 * (diab == 1) + b8 * (lvhn1 == 1) +
                       b9 * height + b10 * chol + b11 * (treat == 1) + trial_eff[trial] +
                       (b12 * age + b13 * (sex == 1) + b14 * sbp + b15 * (smoke == 1)) * (treat == 1))
    p <- plogis(log_odds)
    df10$death <- rbinom(n, 1, p)
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat0 <- matrix(NA, nrow = 66, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 66, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- glm(death ~ (age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial),
                    data = train0,family = binomial,control = glm.control(maxit = 100000))
        mod1 <- glm(death ~ (age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial),
                    data = train1,family = binomial,control = glm.control(maxit = 100000))
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
        #recalibration of event rates adapted to train
        X0 <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial)"),
                           data = train0)
        X1 <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial)"),
                           data = train1)
        lp0 <- X0 %*% coef0
        lp1 <- X1 %*% coef1
        
        coef0[1] <- coef0[1] + glm(death ~ offset(lp0), data = train0, family = "binomial")$coef[1]
        coef1[1] <- coef1[1] + glm(death ~ offset(lp1), data = train1, family = "binomial")$coef[1]
        
        #ite estimation
        X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"),
                          data = test)
        lp0 <- coef0[1] + mean(coef0[12:15])*X[, 12] + mean(c(coef0[2],coef0[16:19]))*X[, 2] +
          mean(c(coef0[3],coef0[20:23]))*X[, 3] + mean(c(coef0[4],coef0[24:27]))*X[, 4] +
          mean(c(coef0[5],coef0[28:31]))*X[, 5] + mean(c(coef0[6],coef0[32:35]))*X[, 6] +
          mean(c(coef0[7],coef0[36:39]))*X[, 7] + mean(c(coef0[8],coef0[40:43]))*X[, 8] +
          mean(c(coef0[9],coef0[44:47]))*X[, 9] + mean(c(coef0[10],coef0[48:51]))*X[, 10] +
          mean(c(coef0[11],coef0[52:55]))*X[, 11]
        lp1 <- coef1[1] + mean(coef1[12:15])*X[, 12] + mean(c(coef1[2],coef1[16:19]))*X[, 2] +
          mean(c(coef1[3],coef1[20:23]))*X[, 3] + mean(c(coef1[4],coef1[24:27]))*X[, 4] +
          mean(c(coef1[5],coef1[28:31]))*X[, 5] + mean(c(coef1[6],coef1[32:35]))*X[, 6] +
          mean(c(coef1[7],coef1[36:39]))*X[, 7] + mean(c(coef1[8],coef1[40:43]))*X[, 8] +
          mean(c(coef1[9],coef1[44:47]))*X[, 9] + mean(c(coef1[10],coef1[48:51]))*X[, 10] +
          mean(c(coef1[11],coef1[52:55]))*X[, 11]
        ite <- expit(lp1) - expit(lp0)
        
        #c-statistic for benefit
        cstat = cstat4ben(outcome = as.numeric(test$death),
                          treatment = test$treat == 1,
                          score = ite)
        
        c.ben[i] <- cstat[1]
        c.ben.se[i] <- cstat[2]
        
        #calibration plot
        predq5 <- quantcut(ite, q = 5)
        pben <- tapply(ite, predq5, mean)
        oben <- sapply(levels(predq5), function(x) {
          lm(death ~ treat, data = test, subset = predq5 == x)$coef[2]
        })
        lm_ <- lm(oben ~ pben)
        
        a[i] <- lm_$coefficients[[1]]
        b[i] <- lm_$coefficients[[2]]
        mse.fs[i] <- mse(oben, pben)
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
        coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12:16, ]),
        coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[17:21, ]),
        coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[22:26, ]),
        coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[27:31, ]),
        coef_mat.5 = mean(coef_mat0[5, ]) + mean(coef_mat0[32:36, ]),
        coef_mat.6 = mean(coef_mat0[6, ]) + mean(coef_mat0[37:41, ]),
        coef_mat.7 = mean(coef_mat0[7, ]) + mean(coef_mat0[42:46, ]),
        coef_mat.8 = mean(coef_mat0[8, ]) + mean(coef_mat0[47:51, ]),
        coef_mat.9 = mean(coef_mat0[9, ]) + mean(coef_mat0[52:56, ]),
        coef_mat.10 = mean(coef_mat0[10, ]) + mean(coef_mat0[57:61, ]),
        coef_mat.11 = mean(coef_mat0[11, ]) + mean(coef_mat0[62:66, ]),
        coef_mat.12 = (mean(coef_mat1[1]) + mean(coef_mat1[12:16, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12:16, ])),
        coef_mat.13 = (mean(coef_mat1[2]) + mean(coef_mat1[17:21, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[17:21, ])),
        coef_mat.14 = (mean(coef_mat1[3]) + mean(coef_mat1[22:26, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[22:26, ])),
        coef_mat.15 = (mean(coef_mat1[4]) + mean(coef_mat1[27:31, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[27:31, ])),
        coef_mat.16 = (mean(coef_mat1[5]) + mean(coef_mat1[32:36, ])) - (mean(coef_mat0[5, ]) + mean(coef_mat0[32:36, ])),
        coef_mat.17 = (mean(coef_mat1[6]) + mean(coef_mat1[37:41, ])) - (mean(coef_mat0[6, ]) + mean(coef_mat0[37:41, ])),
        coef_mat.18 = (mean(coef_mat1[7]) + mean(coef_mat1[42:46, ])) - (mean(coef_mat0[7, ]) + mean(coef_mat0[42:46, ])),
        coef_mat.19 = (mean(coef_mat1[8]) + mean(coef_mat1[47:51, ])) - (mean(coef_mat0[8, ]) + mean(coef_mat0[47:51, ])),
        coef_mat.20 = (mean(coef_mat1[9]) + mean(coef_mat1[52:56, ])) - (mean(coef_mat0[9, ]) + mean(coef_mat0[52:56, ])),
        coef_mat.21 = (mean(coef_mat1[10]) + mean(coef_mat1[57:61, ])) - (mean(coef_mat0[10, ]) + mean(coef_mat0[57:61, ])),
        coef_mat.22 = (mean(coef_mat1[11]) + mean(coef_mat1[62:66, ])) - (mean(coef_mat0[11, ]) + mean(coef_mat0[62:66, ]))),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl10.n1400 <- t(res.fs.tl10[1:27, ])

save(res.na.sl10.n1400,se.na.sl10.n1400,
     res.re.sl10.n1400,se.re.sl10.n1400,
     res.si.sl10.n1400,se.si.sl10.n1400,
     res.r1.sl10.n1400,
     res.fs.sl10.n1400,
     res.na.tl10.n1400,se.na.tl10.n1400,
     res.re.tl10.n1400,se.re.tl10.n1400,
     res.si.tl10.n1400,se.si.tl10.n1400,
     res.r1.tl10.n1400,
     res.fs.tl10.n1400,
     file = "res_scenario5.Rdata")


#### scenario 6: 10 covariates (6 binary and 4 continuous) & total sample size = 700 ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat <- matrix(NA, nrow = 22, ncol = k)
                        se_mat <- matrix(NA, nrow = 22, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                            se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)"), data=train)
                            lp <- X %*% coef
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
                            X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
                            lp0 <- X0 %*% coef[1:12]
                            lp1 <- X1 %*% coef[1:12] + X1[,2:11] %*% coef[13:22] 
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
                            coef_mat.22 = mean(coef_mat[22, ]),
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
                            se_mat.21 = mean(se_mat[21, ]),
                            se_mat.22 = mean(se_mat[22, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.sl10.n700 <- t(res.na.sl10[1:27, ])
se.na.sl10.n700 <- t(res.na.sl10[28:49, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat <- matrix(NA, nrow = 22, ncol = k)
  se_mat <- matrix(NA, nrow = 22, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train <- df10[df10$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying the model to train
    mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+(1|trial)+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
    lp0 <- X0 %*% coef[1:12]
    lp1 <- X1 %*% coef[1:12] + X1[,2:11] %*% coef[13:22] 
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
    coef_mat.22 = mean(coef_mat[22, ]),
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
    se_mat.21 = mean(se_mat[21, ]),
    se_mat.22 = mean(se_mat[22, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.sl10.n700 <- t(res.re.sl10[1:27, ])
se.re.sl10.n700 <- t(res.re.sl10[28:49 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat <- matrix(NA, nrow = 23, ncol = k)
  se_mat <- matrix(NA, nrow = 23, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train <- df10[df10$trial!=i,]
    
    test0 <- test
    test0$treat <- 0
    test1 <- test
    test1$treat <- 1
    
    #applying model to train
    mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial+(0+treat|trial),
                 data = train,family = binomial,
                 control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef <- mod@beta
    coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
    se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial"), data=train)
    lp <- X %*% coef
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test1)
    lp0 <- X0 %*% coef[1:13]
    lp1 <- X1 %*% coef[1:13] + X1[,2:11] %*% coef[14:23] 
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[13, ]),
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
    coef_mat.21 = mean(coef_mat[22, ]),
    coef_mat.22 = mean(coef_mat[23, ]),
    se_mat.1 = mean(se_mat[1, ]) + mean(se_mat[13, ]),
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
    se_mat.13 = mean(se_mat[14, ]),
    se_mat.14 = mean(se_mat[15, ]),
    se_mat.15 = mean(se_mat[16, ]),
    se_mat.16 = mean(se_mat[17, ]),
    se_mat.17 = mean(se_mat[18, ]),
    se_mat.18 = mean(se_mat[19, ]),
    se_mat.19 = mean(se_mat[20, ]),
    se_mat.20 = mean(se_mat[21, ]),
    se_mat.21 = mean(se_mat[22, ]),
    se_mat.22 = mean(se_mat[23, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.sl10.n700 <- t(res.si.sl10[1:27, ])
se.si.sl10.n700 <- t(res.si.sl10[28:49, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")
) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat <- matrix(NA, nrow = 23, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train <- df10[df10$trial!=i,]
      
      test0 <- test
      test0$treat <- 0
      test1 <- test
      test1$treat <- 1
      
      #applying model to train
      mod <- rrvglm(death ~ (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+
                               chol)*factor(treat)+trial, family = binomialff, data = train,Rank = 1)
      coef <- mod@coefficients
      coef_mat[,i] <- coef
      
      #recalibration of event rates adapted to train
      X <- model.matrix(as.formula("death ~ (age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+trial"), data=train)
      lp <- X %*% coef
      coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
      
      #ite estimation
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat+trial"), data=test1)
      lp0 <- X0 %*% coef[1:13]
      lp1 <- X1 %*% coef[1:13] + X1[,2:11] %*% coef[14:23] 
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.r1),
    coef_mat.1 = mean(coef_mat[1, ]) + mean(coef_mat[13, ]),
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
    coef_mat.21 = mean(coef_mat[22, ]),
    coef_mat.22 = mean(coef_mat[23, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl10.n700 <- t(res.r1.sl10[1:27, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        coef_mat <- matrix(NA, nrow = 77, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train <- df10[df10$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+
                                                              height+chol)*factor(treat)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+
                                                                                            lvhn1+height+chol):factor(trial),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef[is.na(coef)] <- 0
                            coef_mat[,i] <- coef
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol)*factor(treat)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol):factor(trial)"),
                                              data=train)
                            lp <- X %*% coef
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            X0 <- model.matrix(as.formula("death~trial+age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test0)
                            X1 <- model.matrix(as.formula("death~trial+age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+treat"), data=test1)
                            lp0 <- coef[1]+mean(coef[2:6])*X0[,2]+mean(c(coef[7],coef[28:32]))*X0[,3]+mean(c(coef[8],coef[33:37]))*X0[,4]+
                              mean(c(coef[9],coef[38:42]))*X0[,5]+mean(c(coef[9],coef[38:42]))*X0[,5]+mean(c(coef[9],coef[38:42]))*X0[,5]+
                              mean(c(coef[10],coef[43:47]))*X0[,6]+mean(c(coef[11],coef[48:52]))*X0[,7]+mean(c(coef[12],coef[53:57]))*X0[,8]+
                              mean(c(coef[13],coef[58:62]))*X0[,9]+mean(c(coef[14],coef[63:67]))*X0[,10]+mean(c(coef[15],coef[68:72]))*X0[,11]+
                              mean(c(coef[16],coef[73:77]))*X0[,12]+coef[17]*X0[,13]
                            lp1 <- coef[1]+mean(coef[2:6])*X1[,2]+mean(c(coef[7],coef[28:32]))*X1[,3]+mean(c(coef[8],coef[33:37]))*X1[,4]+
                              mean(c(coef[9],coef[38:42]))*X1[,5]+mean(c(coef[9],coef[38:42]))*X1[,5]+mean(c(coef[9],coef[38:42]))*X1[,5]+
                              mean(c(coef[10],coef[43:47]))*X1[,6]+mean(c(coef[11],coef[48:52]))*X1[,7]+mean(c(coef[12],coef[53:57]))*X1[,8]+
                              mean(c(coef[13],coef[58:62]))*X1[,9]+mean(c(coef[14],coef[63:67]))*X1[,10]+mean(c(coef[15],coef[68:72]))*X1[,11]+
                              mean(c(coef[16],coef[73:77]))*X1[,12]+coef[17]*X1[,13]+coef[18]*X1[,3]+coef[19]*X1[,4]+coef[20]*X1[,5]+
                              coef[21]*X1[,6]+coef[22]*X1[,7]+coef[23]*X1[,8]+coef[24]*X1[,9]+coef[25]*X1[,10]+coef[26]*X1[,11]+coef[27]*X1[,12]
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
                        
                        list(res = c(
                          c.ben = mean(c.ben),
                          c.ben.se = mean(c.ben.se),
                          a = mean(a),
                          b = mean(b),
                          mse = mean(mse.fs),
                          coef_mat.1 = mean(coef_mat[1:6, ]),
                          coef_mat.2 = mean(c(coef_mat[7, ]+coef_mat[28:32, ])),
                          coef_mat.3 = mean(c(coef_mat[8, ]+coef_mat[33:37, ])),
                          coef_mat.4 = mean(c(coef_mat[9, ]+coef_mat[38:42, ])),
                          coef_mat.5 = mean(c(coef_mat[10, ]+coef_mat[43:47, ])),
                          coef_mat.6 = mean(c(coef_mat[11, ]+coef_mat[48:52, ])),
                          coef_mat.7 = mean(c(coef_mat[12, ]+coef_mat[53:57, ])),
                          coef_mat.8 = mean(c(coef_mat[13, ]+coef_mat[58:62, ])),
                          coef_mat.9 = mean(c(coef_mat[14, ]+coef_mat[63:67, ])),
                          coef_mat.10 = mean(c(coef_mat[15, ]+coef_mat[68:72, ])),
                          coef_mat.11 = mean(c(coef_mat[16, ]+coef_mat[73:77, ])),
                          coef_mat.12 = mean(c(coef_mat[17, ])),
                          coef_mat.13 = mean(coef_mat[18, ]),
                          coef_mat.14 = mean(coef_mat[19, ]),
                          coef_mat.15 = mean(coef_mat[20, ]),
                          coef_mat.16 = mean(coef_mat[21, ]),
                          coef_mat.17 = mean(coef_mat[22, ]),
                          coef_mat.18 = mean(c(coef_mat[23, ])),
                          coef_mat.19 = mean(c(coef_mat[24, ])),
                          coef_mat.20 = mean(c(coef_mat[25, ])),
                          coef_mat.21 = mean(coef_mat[26, ]),
                          coef_mat.22 = mean(coef_mat[27, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.fs.sl10.n700 <- t(res.fs.sl10[1:27, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","magrittr","gtools")) %dopar% {
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
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.1
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- 0.06
                        b11 <- -0.3
                        b12 <- 0.01
                        b13 <- 0.04
                        b14 <- -0.2
                        b15 <- -0.01
                        
                        df10$treat <- rep(c(0,1), times = n/(2*k))
                        
                        trial_eff <- rep(0,7)
                        
                        log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
                        p <- plogis(log_odds)
                        df10$death <- rbinom(n, 1, p)
                        
                        df0 <- df10[df10$treat==0,]
                        df1 <- df10[df10$treat==1,]
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.na <- c()
                        coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
                        coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
                        se_mat0 <- matrix(NA, nrow = 11, ncol = k)
                        se_mat1 <- matrix(NA, nrow = 11, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df10[df10$trial==i,]
                            train0 <- df0[df0$trial!=i,]
                            train1 <- df1[df1$trial!=i,]
                            
                            #applying the model to train
                            mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol,
                                        data = train0,family = binomial,
                                        control = glm.control(maxit=100000))
                            mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol,
                                        data = train1,family = binomial,
                                        control = glm.control(maxit=100000))
                            coef0 <- mod0$coefficients
                            coef1 <- mod1$coefficients
                            coef0[is.na(coef0)] <- 0
                            coef1[is.na(coef1)] <- 0
                            coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
                            coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
                            se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
                            se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
                            
                            #recalibration of event rates adapted to train
                            X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train0)
                            X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train1)
                            lp0 <- X0 %*% coef0
                            lp1 <- X1 %*% coef1
                            
                            coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                            coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                            
                            #ite estimaiton
                            X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=test)
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
                            coef_mat.11 = mean(coef_mat0[11, ]),
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
                            coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
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
                            se_mat.11 = mean(se_mat0[11, ]),
                            se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                            se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                            se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                            se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                            se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                            se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                            se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                            se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                            se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                            se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
                            se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.na.tl10.n700 <- t(res.na.tl10[1:27, ])
se.na.tl10.n700 <- t(res.na.tl10[28:49, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  coef_mat0 <- matrix(NA, nrow = 11, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 11, ncol = k)
  se_mat0 <- matrix(NA, nrow = 11, ncol = k)
  se_mat1 <- matrix(NA, nrow = 11, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),
                  data = train0,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    mod1 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+(1|trial),
                  data = train1,family = binomial,
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
    coef0 <- mod0@beta
    coef1 <- mod1@beta
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=train1)
    lp0 <- X0[ ,1:length(coef0)] %*% coef0
    lp1 <- X1[ ,1:length(coef1)] %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol"), data=test)
    lp0 <- X[ ,1:length(coef0)] %*% coef0
    lp1 <- X[ ,1:length(coef1)] %*% coef1
    
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
    coef_mat.11 = mean(coef_mat0[11, ]),
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
    coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
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
    se_mat.11 = mean(se_mat0[11, ]),
    se_mat.12 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
    se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
    se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
    se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
    se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
    se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
    se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
    se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
    se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.re.tl10.n700 <- t(res.re.tl10[1:27, ])
se.re.tl10.n700 <- t(res.re.tl10[28:49 , ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  coef_mat0 <- matrix(NA, nrow = 12, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 12, ncol = k)
  se_mat0 <- matrix(NA, nrow = 12, ncol = k)
  se_mat1 <- matrix(NA, nrow = 12, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df10[df10$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+trial,
                data = train0,family = binomial,
                control = glm.control(maxit=100000))
    mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+chol+trial,
                data = train1,family = binomial,
                control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    coef0[is.na(coef0)] <- 0
    coef1[is.na(coef1)] <- 0
    coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
    coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
    se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
    se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train0)
    X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train1)
    lp0 <- X0 %*% coef0
    lp1 <- X1 %*% coef1
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimaiton
    X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=test)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  
  list(
    res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12, ]),
    coef_mat.2 = mean(coef_mat0[2, ]),
    coef_mat.3 = mean(coef_mat0[3, ]),
    coef_mat.4 = mean(coef_mat0[4, ]),
    coef_mat.5 = mean(coef_mat0[5, ]),
    coef_mat.6 = mean(coef_mat0[6, ]),
    coef_mat.7 = mean(coef_mat0[7, ]),
    coef_mat.8 = mean(coef_mat0[8, ]),
    coef_mat.9 = mean(coef_mat0[9, ]),
    coef_mat.10 = mean(coef_mat0[10, ]),
    coef_mat.11 = mean(coef_mat0[11, ]),
    coef_mat.12 = (mean(coef_mat1[1, ]) + mean(coef_mat1[12, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12, ])),
    coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
    coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
    coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
    coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
    coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
    coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
    coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
    coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
    coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
    coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ]),
    se_mat.1 = mean(se_mat0[1, ]) + mean(se_mat0[12, ]),
    se_mat.2 = mean(se_mat0[2, ]),
    se_mat.3 = mean(se_mat0[3, ]),
    se_mat.4 = mean(se_mat0[4, ]),
    se_mat.5 = mean(se_mat0[5, ]),
    se_mat.6 = mean(se_mat0[6, ]),
    se_mat.7 = mean(se_mat0[7, ]),
    se_mat.8 = mean(se_mat0[8, ]),
    se_mat.9 = mean(se_mat0[9, ]),
    se_mat.10 = mean(se_mat0[10, ]),
    se_mat.11 = mean(se_mat0[11, ]),
    se_mat.12 = (mean(se_mat1[1, ]) + mean(se_mat1[12, ])) - (mean(se_mat0[1, ]) + mean(se_mat0[12, ])),
    se_mat.13 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
    se_mat.14 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
    se_mat.15 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
    se_mat.16 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
    se_mat.17 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
    se_mat.18 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
    se_mat.19 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
    se_mat.20 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
    se_mat.21 = mean(se_mat1[10, ]) - mean(se_mat0[10, ]),
    se_mat.22 = mean(se_mat1[11, ]) - mean(se_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.si.tl10.n700 <- t(res.si.tl10[1:27, ])
se.si.tl10.n700 <- t(res.si.tl10[28:49, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM")
) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.1
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- 0.06
  b11 <- -0.3
  b12 <- 0.01
  b13 <- 0.04
  b14 <- -0.2
  b15 <- -0.01
  
  df10$treat <- rep(c(0,1), times = n/(2*k))
  
  trial_eff <- rep(0,7)
  
  log_odds <- with(df10, b0+b1*age+b2*(sex==1)+b3*sbp+b4*(mi==1)+b5*(stroke==1)+b6*(smoke==1)+b7*(diab==1)+b8*(lvhn1==1)+b9*height+b10*chol+b11*(treat==1)+trial_eff[trial]+(b12*age+b13*(sex==1)+b14*sbp+b15*(smoke==1))*(treat==1))
  p <- plogis(log_odds)
  df10$death <- rbinom(n, 1, p)
  
  df0 <- df10[df10$treat==0,]
  df1 <- df10[df10$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  coef_mat0 <- matrix(NA, nrow = 12, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 12, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df10[df10$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- rrvglm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height+chol+trial, family = binomialff, data = train0,Rank = 1)
      mod1 <- rrvglm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height+chol+trial, family = binomialff, data = train1,Rank = 1)
      
      coef0 <- mod0@coefficients
      coef1 <- mod1@coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      coef_mat0[,i] <- coef0
      coef_mat1[,i] <- coef1
      
      #recalibration of event rates adapted to train
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=train1)
      lp0 <- X0 %*% coef0
      lp1 <- X1 %*% coef1
      
      coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
      coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
      
      #ite estimation
      X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"), data=test)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(
    res = c(c.ben = mean(c.ben),
            c.ben.se = mean(c.ben.se),
            a = mean(a),
            b = mean(b),
            mse = mean(mse.r1),
            coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12, ]),
            coef_mat.2 = mean(coef_mat0[2, ]),
            coef_mat.3 = mean(coef_mat0[3, ]),
            coef_mat.4 = mean(coef_mat0[4, ]),
            coef_mat.5 = mean(coef_mat0[5, ]),
            coef_mat.6 = mean(coef_mat0[6, ]),
            coef_mat.7 = mean(coef_mat0[7, ]),
            coef_mat.8 = mean(coef_mat0[8, ]),
            coef_mat.9 = mean(coef_mat0[9, ]),
            coef_mat.10 = mean(coef_mat0[10, ]),
            coef_mat.11 = mean(coef_mat0[11, ]),
            coef_mat.12 = (mean(coef_mat1[1, ]) + mean(coef_mat1[12, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12, ])),
            coef_mat.13 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
            coef_mat.14 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
            coef_mat.15 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
            coef_mat.16 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
            coef_mat.17 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
            coef_mat.18 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
            coef_mat.19 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
            coef_mat.20 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
            coef_mat.21 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
            coef_mat.22 = mean(coef_mat1[11, ]) - mean(coef_mat0[11, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.tl10.n700 <- t(res.r1.tl10[1:27, ])

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
  .packages = c("tidyverse", "Hmisc", "lme4", "magrittr", "gtools", "VGAM")) %dopar%
  {
    set.seed(j)
    n <- 700
    k <- 7
    
    trial <- rep(1:k, each = floor(n / k))
    df10 <- as.data.frame(trial)
    
    mean_age <- c(52, 56, 64, 70, 77, 78, 82) #Age
    sd_age <- c(4, 2, 1, 3, 4, 6, 2)
    df10$age <- round(rnorm(n, mean = mean_age[trial], sd = sd_age[trial]), 0)
    df10$age <- df10$age - 70
    
    pman <- c(0.8, 0.4, 0.5, 0.6, 0.5, 0.7, 0.5) #Sex
    df10$sex <- rbinom(n, 1, prob = pman[trial])
    
    mean_sbp <- c(186, 182, 170, 185, 190, 188, 197) #SBP
    sd_sbp <- c(9, 11, 5, 12, 9, 10, 16)
    df10$sbp <- round(rnorm(n, mean = mean_sbp[trial], sd = sd_sbp[trial]), 0)
    df10$sbp <- df10$sbp - 185
    
    pmi <- c(0.1, 0.005, 0.01, 0.02, 0.05, 0.01, 0.04) #Previous myocardial infarction
    df10$mi <- rbinom(n, 1, prob = pmi[trial])
    
    pstroke <- c(0.002, 0.06, 0.02, 0.02, 0.001, 0.008, 0.04) #Previous stroke
    df10$stroke <- rbinom(n, 1, prob = pstroke[trial])
    
    psmoke <- c(0.5, 0.2, 0.3, 0.4, 0.3, 0.25, 0.3) #Current smoker
    df10$smoke <- rbinom(n, 1, prob = psmoke[trial])
    
    pdiab <- c(0.03, 0.001, 0.002, 0.07, 0.003, 0.01, 0.002) #Diabetes
    df10$diab <- rbinom(n, 1, prob = pdiab[trial])
    
    plvhn1 <- c(0.13, 0.11, 0.05, 0.25, 0.05, 0.06, 0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df10$lvhn1 <- rbinom(n, 1, prob = plvhn1[trial])
    
    mean_height <- c(176, 162, 167, 169, 168, 170, 167) #Height
    sd_height <- c(6, 9, 10, 10, 10, 9, 9)
    df10$height <- round(rnorm(n, mean = mean_height[trial], sd = sd_height[trial]), 0)
    df10$height <- df10$height - 168
    
    mean_chol <- c(6.6, 6.5, 6.4, 6.1, 6.4, 6, 6.4) #Cholesterol
    sd_chol <- c(0.01, 0.015, 0.011, 0.012, 0.012, 0.012, 0.01)
    df10$chol <- round(rnorm(n, mean = mean_chol[trial], sd = sd_chol[trial]), 0)
    df10$chol <- df10$chol - 6.3
    
    b0 <- -1.4
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.1
    b4 <- 0.82
    b5 <- 0.8
    b6 <- 0.7
    b7 <- 0.1
    b8 <- 0.33
    b9 <- -0.02
    b10 <- 0.06
    b11 <- -0.3
    b12 <- 0.01
    b13 <- 0.04
    b14 <- -0.2
    b15 <- -0.01
    
    df10$treat <- rep(c(0, 1), times = n / (2 * k))
    
    trial_eff <- rep(0, 7)
    
    log_odds <- with(df10,b0 + b1 * age + b2 * (sex == 1) + b3 * sbp + b4 * (mi == 1) + 
                       b5 * (stroke == 1) + b6 * (smoke == 1) + b7 * (diab == 1) + b8 * (lvhn1 == 1) +
                       b9 * height + b10 * chol + b11 * (treat == 1) + trial_eff[trial] +
                       (b12 * age + b13 * (sex == 1) + b14 * sbp + b15 * (smoke == 1)) * (treat == 1))
    p <- plogis(log_odds)
    df10$death <- rbinom(n, 1, p)
    
    df0 <- df10[df10$treat==0,]
    df1 <- df10[df10$treat==1,]
    
    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    coef_mat0 <- matrix(NA, nrow = 66, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 66, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df10[df10$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- glm(death ~ (age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial),
                    data = train0,family = binomial,control = glm.control(maxit = 100000))
        mod1 <- glm(death ~ (age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial),
                    data = train1,family = binomial,control = glm.control(maxit = 100000))
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0[is.na(coef0)] <- 0
        coef1[is.na(coef1)] <- 0
        coef_mat0[,i] <- coef0
        coef_mat1[,i] <- coef1
        
        #recalibration of event rates adapted to train
        X0 <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial)"),
                           data = train0)
        X1 <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol)*factor(trial)"),
                           data = train1)
        lp0 <- X0 %*% coef0
        lp1 <- X1 %*% coef1
        
        coef0[1] <- coef0[1] + glm(death ~ offset(lp0), data = train0, family = "binomial")$coef[1]
        coef1[1] <- coef1[1] + glm(death ~ offset(lp1), data = train1, family = "binomial")$coef[1]
        
        #ite estimation
        X <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height+chol+trial"),
                          data = test)
        lp0 <- coef0[1] + mean(coef0[12:15])*X[, 12] + mean(c(coef0[2],coef0[16:19]))*X[, 2] +
          mean(c(coef0[3],coef0[20:23]))*X[, 3] + mean(c(coef0[4],coef0[24:27]))*X[, 4] +
          mean(c(coef0[5],coef0[28:31]))*X[, 5] + mean(c(coef0[6],coef0[32:35]))*X[, 6] +
          mean(c(coef0[7],coef0[36:39]))*X[, 7] + mean(c(coef0[8],coef0[40:43]))*X[, 8] +
          mean(c(coef0[9],coef0[44:47]))*X[, 9] + mean(c(coef0[10],coef0[48:51]))*X[, 10] +
          mean(c(coef0[11],coef0[52:55]))*X[, 11]
        lp1 <- coef1[1] + mean(coef1[12:15])*X[, 12] + mean(c(coef1[2],coef1[16:19]))*X[, 2] +
          mean(c(coef1[3],coef1[20:23]))*X[, 3] + mean(c(coef1[4],coef1[24:27]))*X[, 4] +
          mean(c(coef1[5],coef1[28:31]))*X[, 5] + mean(c(coef1[6],coef1[32:35]))*X[, 6] +
          mean(c(coef1[7],coef1[36:39]))*X[, 7] + mean(c(coef1[8],coef1[40:43]))*X[, 8] +
          mean(c(coef1[9],coef1[44:47]))*X[, 9] + mean(c(coef1[10],coef1[48:51]))*X[, 10] +
          mean(c(coef1[11],coef1[52:55]))*X[, 11]
        ite <- expit(lp1) - expit(lp0)
        
        #c-statistic for benefit
        cstat = cstat4ben(outcome = as.numeric(test$death),
                          treatment = test$treat == 1,
                          score = ite)
        
        c.ben[i] <- cstat[1]
        c.ben.se[i] <- cstat[2]
        
        #calibration plot
        predq5 <- quantcut(ite, q = 5)
        pben <- tapply(ite, predq5, mean)
        oben <- sapply(levels(predq5), function(x) {
          lm(death ~ treat, data = test, subset = predq5 == x)$coef[2]
        })
        lm_ <- lm(oben ~ pben)
        
        a[i] <- lm_$coefficients[[1]]
        b[i] <- lm_$coefficients[[2]]
        mse.fs[i] <- mse(oben, pben)
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
        coef_mat.1 = mean(coef_mat0[1, ]) + mean(coef_mat0[12:16, ]),
        coef_mat.2 = mean(coef_mat0[2, ]) + mean(coef_mat0[17:21, ]),
        coef_mat.3 = mean(coef_mat0[3, ]) + mean(coef_mat0[22:26, ]),
        coef_mat.4 = mean(coef_mat0[4, ]) + mean(coef_mat0[27:31, ]),
        coef_mat.5 = mean(coef_mat0[5, ]) + mean(coef_mat0[32:36, ]),
        coef_mat.6 = mean(coef_mat0[6, ]) + mean(coef_mat0[37:41, ]),
        coef_mat.7 = mean(coef_mat0[7, ]) + mean(coef_mat0[42:46, ]),
        coef_mat.8 = mean(coef_mat0[8, ]) + mean(coef_mat0[47:51, ]),
        coef_mat.9 = mean(coef_mat0[9, ]) + mean(coef_mat0[52:56, ]),
        coef_mat.10 = mean(coef_mat0[10, ]) + mean(coef_mat0[57:61, ]),
        coef_mat.11 = mean(coef_mat0[11, ]) + mean(coef_mat0[62:66, ]),
        coef_mat.12 = (mean(coef_mat1[1]) + mean(coef_mat1[12:16, ])) - (mean(coef_mat0[1, ]) + mean(coef_mat0[12:16, ])),
        coef_mat.13 = (mean(coef_mat1[2]) + mean(coef_mat1[17:21, ])) - (mean(coef_mat0[2, ]) + mean(coef_mat0[17:21, ])),
        coef_mat.14 = (mean(coef_mat1[3]) + mean(coef_mat1[22:26, ])) - (mean(coef_mat0[3, ]) + mean(coef_mat0[22:26, ])),
        coef_mat.15 = (mean(coef_mat1[4]) + mean(coef_mat1[27:31, ])) - (mean(coef_mat0[4, ]) + mean(coef_mat0[27:31, ])),
        coef_mat.16 = (mean(coef_mat1[5]) + mean(coef_mat1[32:36, ])) - (mean(coef_mat0[5, ]) + mean(coef_mat0[32:36, ])),
        coef_mat.17 = (mean(coef_mat1[6]) + mean(coef_mat1[37:41, ])) - (mean(coef_mat0[6, ]) + mean(coef_mat0[37:41, ])),
        coef_mat.18 = (mean(coef_mat1[7]) + mean(coef_mat1[42:46, ])) - (mean(coef_mat0[7, ]) + mean(coef_mat0[42:46, ])),
        coef_mat.19 = (mean(coef_mat1[8]) + mean(coef_mat1[47:51, ])) - (mean(coef_mat0[8, ]) + mean(coef_mat0[47:51, ])),
        coef_mat.20 = (mean(coef_mat1[9]) + mean(coef_mat1[52:56, ])) - (mean(coef_mat0[9, ]) + mean(coef_mat0[52:56, ])),
        coef_mat.21 = (mean(coef_mat1[10]) + mean(coef_mat1[57:61, ])) - (mean(coef_mat0[10, ]) + mean(coef_mat0[57:61, ])),
        coef_mat.22 = (mean(coef_mat1[11]) + mean(coef_mat1[62:66, ])) - (mean(coef_mat0[11, ]) + mean(coef_mat0[62:66, ]))),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl10.n700 <- t(res.fs.tl10[1:27, ])

save(res.na.sl10.n700,se.na.sl10.n700,
     res.re.sl10.n700,se.re.sl10.n700,
     res.si.sl10.n700,se.si.sl10.n700,
     res.r1.sl10.n700,
     res.fs.sl10.n700,
     res.na.tl10.n700,se.na.tl10.n700,
     res.re.tl10.n700,se.re.tl10.n700,
     res.si.tl10.n700,se.si.tl10.n700,
     res.r1.tl10.n700,
     res.fs.tl10.n700,
     file = "res_scenario6.Rdata")
