library(tidyverse) # data.frame manipulation
library(Hmisc) # cstat computation
library(magrittr) # cstat computation
library(lme4) # glmer
library(doParallel) # parallel computing
library(gtools) # calibration computation
library(VGAM) # rrvglm
library(mvmeta) # meta-analysis
library(msm) # delta method

setwd("...")

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

#### scenario 7: 3 covariates (1 binary and 2 continuous) & sample size = 2800 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) -
            expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    
    #ite estimation
    ite <- predict(mod, newdata=test1, type="response") -
      predict(mod, newdata=test0, type="response")
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se,
                           upr = ite + qt(.975,df.residual(mod)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.na.sl.sim7 <- t(res.na.sl7[1:14, ])
se.na.sl.sim7 <- t(res.na.sl7[15:22, ])
colnames(res.na.sl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.sl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + as.data.frame(VarCorr(mod))$vcov),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + as.data.frame(VarCorr(mod))$vcov))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.re.sl.sim7 <- t(res.re.sl7[1:14, ])
se.re.sl.sim7 <- t(res.re.sl7[15:22, ])
colnames(res.re.sl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.sl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- glm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat +
                 treat:factor(trial), data = train,family = "binomial")
    coef <- mod$coefficients
    
    K.si <- array(0, dim=c(8, 18, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
    for(p in 1:6){
      diag(K.si[2:8,ind,p]) <- 1
      K.si[1,p,p] <- 1
      if(p == 1) next
      K.si[5,12+p,p] <- 1
    }
    allcoef <- t(apply(K.si, 3, function(x){x %*% coef}))
    allvar <- apply(K.si, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
    
    si_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = si_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(si_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, si_meta$vcov)
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(si_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(si_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.si.sl.sim7 <- t(res.si.sl7[1:14, ])
se.si.sl.sim7 <- t(res.si.sl7[15:22, ])
colnames(res.si.sl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.sl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- rrvglm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat + treat:factor(trial),
                  family = binomialff, data = train,Rank = 1) 
    coef <- mod@coefficients
    
    K.r1 <- array(0, dim=c(8, 18, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
    for(p in 1:6){
      diag(K.r1[2:8,ind,p]) <- 1
      K.r1[1,p,p] <- 1
      if(p == 1) next
      K.r1[5,12+p,p] <- 1
    }
    allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
    allvar <- apply(K.r1, 3, function(x){x %*% vcovvlm(mod) %*% t(x)}, simplify=F)
    
    r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = r1_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(si_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, r1_meta$vcov)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(r1_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(r1_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(c.ben = mean(c.ben),
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
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.r1.sl.sim7 <- t(res.r1.sl7[1:14, ])
se.r1.sl.sim7 <- t(res.r1.sl7[15:22, ])
colnames(res.r1.sl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.r1.sl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- glm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat * factor(trial),
               data=train, family="binomial")
    coef <- mod$coefficients
    
    K.fs <- array(0, dim=c(8, 48, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_age=F,
             t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
             t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
             t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F,t2_age_treat=F,
             t3_age_treat=F, t4_age_treat=F, t5_age_treat=F, t6_age_treat=F,t2_sex_treat=F,
             t3_sex_treat=F, t4_sex_treat=F, t5_sex_treat=F, t6_sex_treat=F,t2_sbp_treat=F,
             t3_sbp_treat=F, t4_sbp_treat=F, t5_sbp_treat=F, t6_sbp_treat=F)
    for(p in 1:6){
      diag(K.fs[2:8,ind,p]) <- 1
      K.fs[1,p,p] <- 1
      if(p == 1) next
      K.fs[2,12+p,p] <- 1
      K.fs[3,17+p,p] <- 1
      K.fs[4,22+p,p] <- 1
      K.fs[5,27+p,p] <- 1
      K.fs[6,32+p,p] <- 1
      K.fs[7,37+p,p] <- 1
      K.fs[8,42+p,p] <- 1
    }
    allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
    allvar <- apply(K.fs, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
    
    fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = fs_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(fs_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"),
                      data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, fs_meta$vcov)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(fs_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(fs_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(unlist(ite_true[i]),ite)
  } 
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
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
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.fs.sl.sim7 <- t(res.fs.sl7[1:14, ])
se.fs.sl.sim7 <- t(res.fs.sl7[15:22, ])
colnames(res.fs.sl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.sl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    
    #ite estimation
    ite <- predict(mod1, newdata=test, type="response") -
      predict(mod0, newdata=test, type="response")
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se,
                           upr = ite + qt(.975,df.residual(mod0)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.na),
    in_int = mean(in_int),
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
res.na.tl.sim7 <- t(res.na.tl7[1:14, ])
se.na.tl.sim7 <- t(res.na.tl7[15:22, ])
colnames(res.na.tl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.tl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")


#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + as.data.frame(VarCorr(mod0))$vcov),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + as.data.frame(VarCorr(mod0))$vcov))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    in_int = mean(in_int),
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
res.re.tl.sim7 <- t(res.re.tl7[1:14, ])
se.re.tl.sim7 <- t(res.re.tl7[15:22, ])
colnames(res.re.tl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.tl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death ~ -1 + factor(trial) + age + factor(sex) + sbp,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death ~ -1 + factor(trial) + age + factor(sex) + sbp,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    
    K.si0 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.si0[2:4,ind,p]) <- 1
      K.si0[1,p,p] <- 1
    }
    allcoef0 <- t(apply(K.si0, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.si0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
    
    si_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = si_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(si_meta0$vcov)
    
    K.si1 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.si1[2:4,ind,p]) <- 1
      K.si1[1,p,p] <- 1
    }
    allcoef1 <- t(apply(K.si1, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.si1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
    
    si_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = si_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(si_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, si_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, si_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(si_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(si_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    in_int = mean(in_int),
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
res.si.tl.sim7 <- t(res.si.tl7[1:14, ])
se.si.tl.sim7 <- t(res.si.tl7[15:22, ])
colnames(res.si.tl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.tl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- rrvglm(death ~ -1 + factor(trial) + age + factor(sex) + sbp, family = binomialff, data = train0,Rank = 1)
    mod1 <- rrvglm(death ~ -1 + factor(trial) + age + factor(sex) + sbp, family = binomialff, data = train1,Rank = 1)
    coef0 <- mod0@coefficients
    coef1 <- mod1@coefficients
    
    K.r10 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.r10[2:4,ind,p]) <- 1
      K.r10[1,p,p] <- 1
    }
    allcoef0 <- t(apply(K.r10, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.r10, 3, function(x){x %*% vcovvlm(mod0) %*% t(x)}, simplify=F)
    
    r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = r1_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(r1_meta0$vcov)
    
    K.r11 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.r11[2:4,ind,p]) <- 1
      K.r11[1,p,p] <- 1
    }
    allcoef1 <- t(apply(K.r11, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.r11, 3, function(x){x %*% vcovvlm(mod1) %*% t(x)}, simplify=F)
    
    r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = r1_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(r1_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, r1_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, r1_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(r1_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(r1_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.r1.tl.sim7 <- t(res.r1.tl7[1:14, ])
se.r1.tl.sim7 <- t(res.r1.tl7[15:22, ])
colnames(res.r1.tl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.r1.tl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl7 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death ~ -1 + factor(trial) + (age+factor(sex)+sbp)*factor(trial),
                data = train0,family = binomial,control = glm.control(maxit = 100000))
    mod1 <- glm(death ~ -1 + factor(trial) + (age+factor(sex)+sbp)*factor(trial),
                data = train1,family = binomial,control = glm.control(maxit = 100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    
    K.fs0 <- array(0, dim=c(4, 24, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, t2_age=F,
             t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
             t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
             t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F)
    for(p in 1:6){
      diag(K.fs0[2:4,ind,p]) <- 1
      K.fs0[1,p,p] <- 1
      if(p == 1) next
      K.fs0[2,8+p,p] <- 1
      K.fs0[3,13+p,p] <- 1
      K.fs0[4,18+p,p] <- 1
    }
    allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.fs0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
    
    fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = fs_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(fs_meta0$vcov)
    
    K.fs1 <- array(0, dim=c(4, 24, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, t2_age=F,
             t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
             t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
             t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F)
    for(p in 1:6){
      diag(K.fs1[2:4,ind,p]) <- 1
      K.fs1[1,p,p] <- 1
      if(p == 1) next
      K.fs1[2,8+p,p] <- 1
      K.fs1[3,13+p,p] <- 1
      K.fs1[4,18+p,p] <- 1
    }
    allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.fs1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
    
    fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = fs_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(fs_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, fs_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, fs_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(fs_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(fs_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
    in_int = mean(in_int),
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
res.fs.tl.sim7 <- t(res.fs.tl7[1:14, ])
se.fs.tl.sim7 <- t(res.fs.tl7[15:22, ])
colnames(res.fs.tl.sim7) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.tl.sim7) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

save(res.na.sl.sim7,se.na.sl.sim7,
     res.re.sl.sim7,se.re.sl.sim7,
     res.si.sl.sim7,se.si.sl.sim7,
     res.r1.sl.sim7,se.r1.sl.sim7,
     res.fs.sl.sim7,se.fs.sl.sim7,
     res.na.tl.sim7,se.na.tl.sim7,
     res.re.tl.sim7,se.re.tl.sim7,
     res.si.tl.sim7,se.si.tl.sim7,
     res.r1.tl.sim7,se.r1.tl.sim7,
     res.fs.tl.sim7,se.fs.tl.sim7,
     file = "res_scenario7.Rdata")

#### scenario 8: 3 covariates (1 binary and 2 continuous) & sample size = 1400 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    
    #ite estimation
    ite <- predict(mod, newdata=test1, type="response") -
      predict(mod, newdata=test0, type="response")
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se,
                           upr = ite + qt(.975,df.residual(mod)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.na.sl.sim8 <- t(res.na.sl8[1:14, ])
se.na.sl.sim8 <- t(res.na.sl8[15:22, ])
colnames(res.na.sl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.sl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + as.data.frame(VarCorr(mod))$vcov),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + as.data.frame(VarCorr(mod))$vcov))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.re.sl.sim8 <- t(res.re.sl8[1:14, ])
se.re.sl.sim8 <- t(res.re.sl8[15:22, ])
colnames(res.re.sl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.sl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- glm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat +
                 treat:factor(trial), data = train,family = "binomial")
    coef <- mod$coefficients
    
    K.si <- array(0, dim=c(8, 18, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
    for(p in 1:6){
      diag(K.si[2:8,ind,p]) <- 1
      K.si[1,p,p] <- 1
      if(p == 1) next
      K.si[5,12+p,p] <- 1
    }
    allcoef <- t(apply(K.si, 3, function(x){x %*% coef}))
    allvar <- apply(K.si, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
    
    si_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = si_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(si_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, si_meta$vcov)
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(si_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(si_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.si.sl.sim8 <- t(res.si.sl8[1:14, ])
se.si.sl.sim8 <- t(res.si.sl8[15:22, ])
colnames(res.si.sl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.sl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- rrvglm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat + treat:factor(trial),
                  family = binomialff, data = train,Rank = 1) 
    coef <- mod@coefficients
    
    K.r1 <- array(0, dim=c(8, 18, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
    for(p in 1:6){
      diag(K.r1[2:8,ind,p]) <- 1
      K.r1[1,p,p] <- 1
      if(p == 1) next
      K.r1[5,12+p,p] <- 1
    }
    allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
    allvar <- apply(K.r1, 3, function(x){x %*% vcovvlm(mod) %*% t(x)}, simplify=F)
    
    r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = r1_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(si_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, r1_meta$vcov)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(r1_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(r1_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(c.ben = mean(c.ben),
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
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.r1.sl.sim8 <- t(res.r1.sl8[1:14, ])
se.r1.sl.sim8 <- t(res.r1.sl8[15:22, ])
colnames(res.r1.sl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.r1.sl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl8 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")) %dopar% {
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  testerror = try({
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- glm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat * factor(trial),
               data=train, family="binomial")
    coef <- mod$coefficients
    
    K.fs <- array(0, dim=c(8, 48, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_age=F,
             t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
             t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
             t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F,t2_age_treat=F,
             t3_age_treat=F, t4_age_treat=F, t5_age_treat=F, t6_age_treat=F,t2_sex_treat=F,
             t3_sex_treat=F, t4_sex_treat=F, t5_sex_treat=F, t6_sex_treat=F,t2_sbp_treat=F,
             t3_sbp_treat=F, t4_sbp_treat=F, t5_sbp_treat=F, t6_sbp_treat=F)
    for(p in 1:6){
      diag(K.fs[2:8,ind,p]) <- 1
      K.fs[1,p,p] <- 1
      if(p == 1) next
      K.fs[2,12+p,p] <- 1
      K.fs[3,17+p,p] <- 1
      K.fs[4,22+p,p] <- 1
      K.fs[5,27+p,p] <- 1
      K.fs[6,32+p,p] <- 1
      K.fs[7,37+p,p] <- 1
      K.fs[8,42+p,p] <- 1
    }
    allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
    allvar <- apply(K.fs, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
    
    fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = fs_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(fs_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"),
                      data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, fs_meta$vcov)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(fs_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(fs_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(unlist(ite_true[i]),ite)
  } 
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(res = c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
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
res.fs.sl.sim8 <- t(res.fs.sl8[1:14, ])
se.fs.sl.sim8 <- t(res.fs.sl8[15:22, ])
colnames(res.fs.sl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.sl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    
    #ite estimation
    ite <- predict(mod1, newdata=test, type="response") -
      predict(mod0, newdata=test, type="response")
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se,
                           upr = ite + qt(.975,df.residual(mod0)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.na),
    in_int = mean(in_int),
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
res.na.tl.sim8 <- t(res.na.tl8[1:14, ])
se.na.tl.sim8 <- t(res.na.tl8[15:22, ])
colnames(res.na.tl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.tl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")


#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + as.data.frame(VarCorr(mod0))$vcov),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + as.data.frame(VarCorr(mod0))$vcov))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    in_int = mean(in_int),
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
res.re.tl.sim8 <- t(res.re.tl8[1:14, ])
se.re.tl.sim8 <- t(res.re.tl8[15:22, ])
colnames(res.re.tl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.tl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death ~ -1 + factor(trial) + age + factor(sex) + sbp,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death ~ -1 + factor(trial) + age + factor(sex) + sbp,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    
    K.si0 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.si0[2:4,ind,p]) <- 1
      K.si0[1,p,p] <- 1
    }
    allcoef0 <- t(apply(K.si0, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.si0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
    
    si_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = si_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(si_meta0$vcov)
    
    K.si1 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.si1[2:4,ind,p]) <- 1
      K.si1[1,p,p] <- 1
    }
    allcoef1 <- t(apply(K.si1, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.si1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
    
    si_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = si_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(si_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, si_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, si_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(si_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(si_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    in_int = mean(in_int),
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
res.si.tl.sim8 <- t(res.si.tl8[1:14, ])
se.si.tl.sim8 <- t(res.si.tl8[15:22, ])
colnames(res.si.tl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.tl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- rrvglm(death ~ -1 + factor(trial) + age + factor(sex) + sbp, family = binomialff, data = train0,Rank = 1)
    mod1 <- rrvglm(death ~ -1 + factor(trial) + age + factor(sex) + sbp, family = binomialff, data = train1,Rank = 1)
    coef0 <- mod0@coefficients
    coef1 <- mod1@coefficients
    
    K.r10 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.r10[2:4,ind,p]) <- 1
      K.r10[1,p,p] <- 1
    }
    allcoef0 <- t(apply(K.r10, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.r10, 3, function(x){x %*% vcovvlm(mod0) %*% t(x)}, simplify=F)
    
    r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = r1_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(r1_meta0$vcov)
    
    K.r11 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.r11[2:4,ind,p]) <- 1
      K.r11[1,p,p] <- 1
    }
    allcoef1 <- t(apply(K.r11, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.r11, 3, function(x){x %*% vcovvlm(mod1) %*% t(x)}, simplify=F)
    
    r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = r1_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(r1_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, r1_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, r1_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(r1_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(r1_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.r1.tl.sim8 <- t(res.r1.tl8[1:14, ])
se.r1.tl.sim8 <- t(res.r1.tl8[15:22, ])
colnames(res.r1.tl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.r1.tl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl8 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.fs <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death ~ -1 + factor(trial) + (age+factor(sex)+sbp)*factor(trial),
                data = train0,family = binomial,control = glm.control(maxit = 100000))
    mod1 <- glm(death ~ -1 + factor(trial) + (age+factor(sex)+sbp)*factor(trial),
                data = train1,family = binomial,control = glm.control(maxit = 100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    
    K.fs0 <- array(0, dim=c(4, 24, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, t2_age=F,
             t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
             t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
             t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F)
    for(p in 1:6){
      diag(K.fs0[2:4,ind,p]) <- 1
      K.fs0[1,p,p] <- 1
      if(p == 1) next
      K.fs0[2,8+p,p] <- 1
      K.fs0[3,13+p,p] <- 1
      K.fs0[4,18+p,p] <- 1
    }
    allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.fs0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
    
    fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = fs_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(fs_meta0$vcov)
    
    K.fs1 <- array(0, dim=c(4, 24, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, t2_age=F,
             t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
             t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
             t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F)
    for(p in 1:6){
      diag(K.fs1[2:4,ind,p]) <- 1
      K.fs1[1,p,p] <- 1
      if(p == 1) next
      K.fs1[2,8+p,p] <- 1
      K.fs1[3,13+p,p] <- 1
      K.fs1[4,18+p,p] <- 1
    }
    allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.fs1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
    
    fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = fs_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(fs_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, fs_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, fs_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(fs_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(fs_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.fs[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.fs),
    in_int = mean(in_int),
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
res.fs.tl.sim8 <- t(res.fs.tl8[1:14, ])
se.fs.tl.sim8 <- t(res.fs.tl8[15:22, ])
colnames(res.fs.tl.sim8) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.tl.sim8) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

save(res.na.sl.sim8,se.na.sl.sim8,
     res.re.sl.sim8,se.re.sl.sim8,
     res.si.sl.sim8,se.si.sl.sim8,
     res.r1.sl.sim8,se.r1.sl.sim8,
     res.fs.sl.sim8,se.fs.sl.sim8,
     res.na.tl.sim8,se.na.tl.sim8,
     res.re.tl.sim8,se.re.tl.sim8,
     res.si.tl.sim8,se.si.tl.sim8,
     res.r1.tl.sim8,se.r1.tl.sim8,
     res.fs.tl.sim8,se.fs.tl.sim8,
     file = "res_scenario8.Rdata")


#### scenario 9: 3 covariates (1 binary and 2 continuous) & sample size = 700 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    
    #ite estimation
    ite <- predict(mod, newdata=test1, type="response") -
      predict(mod, newdata=test0, type="response")
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se,
                           upr = ite + qt(.975,df.residual(mod)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.na.sl.sim9 <- t(res.na.sl9[1:14, ])
se.na.sl.sim9 <- t(res.na.sl9[15:22, ])
colnames(res.na.sl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.sl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
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
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, vcov(mod))
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + as.data.frame(VarCorr(mod))$vcov),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + as.data.frame(VarCorr(mod))$vcov))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.re.sl.sim9 <- t(res.re.sl9[1:14, ])
se.re.sl.sim9 <- t(res.re.sl9[15:22, ])
colnames(res.re.sl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.sl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- glm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat +
                 treat:factor(trial), data = train,family = "binomial")
    coef <- mod$coefficients
    
    K.si <- array(0, dim=c(8, 18, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
    for(p in 1:6){
      diag(K.si[2:8,ind,p]) <- 1
      K.si[1,p,p] <- 1
      if(p == 1) next
      K.si[5,12+p,p] <- 1
    }
    allcoef <- t(apply(K.si, 3, function(x){x %*% coef}))
    allvar <- apply(K.si, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
    
    si_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = si_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(si_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, si_meta$vcov)
    
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(si_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(si_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.si.sl.sim9 <- t(res.si.sl9[1:14, ])
se.si.sl.sim9 <- t(res.si.sl9[15:22, ])
colnames(res.si.sl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.sl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 8, ncol = k)
  se_mat <- matrix(NA, nrow = 8, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train <- df[df$trial!=i,]
    
    #applying model to train
    mod <- rrvglm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat + treat:factor(trial),
                  family = binomialff, data = train,Rank = 1) 
    coef <- mod@coefficients
    
    K.r1 <- array(0, dim=c(8, 18, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
             age_treat=T, sex_treat=T, sbp_treat=T, t2_treat=F,
             t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
    for(p in 1:6){
      diag(K.r1[2:8,ind,p]) <- 1
      K.r1[1,p,p] <- 1
      if(p == 1) next
      K.r1[5,12+p,p] <- 1
    }
    allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
    allvar <- apply(K.r1, 3, function(x){x %*% vcovvlm(mod) %*% t(x)}, simplify=F)
    
    r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
    coef = r1_meta$coefficients
    
    coef_mat[,i] <- coef
    se_mat[,i] <- diag(si_meta$vcov)
    
    #recalibration of event rates adapted to train
    X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"), data=train)
    lp <- X %*% c(coef)
    coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                   coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
      expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
    
    #ite prediction interval
    se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                        (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, r1_meta$vcov)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(r1_meta$Psi)),
                           upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(r1_meta$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(c.ben = mean(c.ben),
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
    se_mat.8 = mean(se_mat[8, ]))
}
stopCluster(cl)
res.r1.sl.sim9 <- t(res.r1.sl9[1:14, ])
se.r1.sl.sim9 <- t(res.r1.sl9[15:22, ])
colnames(res.r1.sl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.r1.sl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl9 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")) %dopar% {
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
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
                         b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
                         trt_het[df$trial]*df$treat
                       df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                                         b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                                         trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.fs <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 8, ncol = k)
                       se_mat <- matrix(NA, nrow = 8, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           #applying model to train
                           mod <- glm(death ~ -1 + factor(trial) + (age + factor(sex) + sbp) * treat * factor(trial),
                                      data=train, family="binomial")
                           coef <- mod$coefficients
                           
                           K.fs <- array(0, dim=c(8, 48, 6))
                           ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, treat=T,
                                    age_treat=T, sex_treat=T, sbp_treat=T, t2_age=F,
                                    t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
                                    t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
                                    t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F, t2_treat=F,
                                    t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F,t2_age_treat=F,
                                    t3_age_treat=F, t4_age_treat=F, t5_age_treat=F, t6_age_treat=F,t2_sex_treat=F,
                                    t3_sex_treat=F, t4_sex_treat=F, t5_sex_treat=F, t6_sex_treat=F,t2_sbp_treat=F,
                                    t3_sbp_treat=F, t4_sbp_treat=F, t5_sbp_treat=F, t6_sbp_treat=F)
                           for(p in 1:6){
                             diag(K.fs[2:8,ind,p]) <- 1
                             K.fs[1,p,p] <- 1
                             if(p == 1) next
                             K.fs[2,12+p,p] <- 1
                             K.fs[3,17+p,p] <- 1
                             K.fs[4,22+p,p] <- 1
                             K.fs[5,27+p,p] <- 1
                             K.fs[6,32+p,p] <- 1
                             K.fs[7,37+p,p] <- 1
                             K.fs[8,42+p,p] <- 1
                           }
                           allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
                           allvar <- apply(K.fs, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
                           
                           fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml")
                           coef = fs_meta$coefficients
                           
                           coef_mat[,i] <- coef
                           se_mat[,i] <- diag(fs_meta$vcov)
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+factor(sex)+sbp)*factor(treat)"),
                                             data=train)
                           lp <- X %*% c(coef)
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp + 
                                          coef[5] + coef[6]*test$age + coef[7]*test$sex + coef[8]*test$sbp) -
                             expit(coef[1] + coef[2]*test$age + coef[3]*test$sex + coef[4]*test$sbp)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8))) - 
                                               (exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef, fs_meta$vcov)
                           true_val <- unlist(ite_true[i])
                           pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod)) * se * sqrt(1 + diag(fs_meta$Psi)),
                                                  upr = ite + qt(.975,df.residual(mod)) * se * sqrt(1 + diag(fs_meta$Psi)))
                           
                           in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$death),
                                             treatment = test$treat == 1,
                                             score = ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration for benefit
                           lm_ <- lm(ite~unlist(ite_true[i]))
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.fs[i] <- mse(unlist(ite_true[i]),ite)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.fs),
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
res.fs.sl.sim9 <- t(res.fs.sl9[1:14, ])
se.fs.sl.sim9 <- t(res.fs.sl9[15:22, ])
colnames(res.fs.sl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.sl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.na <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    
    #ite estimation
    ite <- predict(mod1, newdata=test, type="response") -
      predict(mod0, newdata=test, type="response")
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se,
                           upr = ite + qt(.975,df.residual(mod0)) * se)
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.na[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.na),
    in_int = mean(in_int),
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
res.na.tl.sim9 <- t(res.na.tl9[1:14, ])
se.na.tl.sim9 <- t(res.na.tl9[15:22, ])
colnames(res.na.tl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.na.tl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")


#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.re <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
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
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, vcov(mod0))
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, vcov(mod1))
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + as.data.frame(VarCorr(mod0))$vcov),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + as.data.frame(VarCorr(mod0))$vcov))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    
    mse.re[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.re),
    in_int = mean(in_int),
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
res.re.tl.sim9 <- t(res.re.tl9[1:14, ])
se.re.tl.sim9 <- t(res.re.tl9[15:22, ])
colnames(res.re.tl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.re.tl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.si <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- glm(death ~ -1 + factor(trial) + age + factor(sex) + sbp,
                data = train0,family = binomial,control = glm.control(maxit=100000))
    mod1 <- glm(death ~ -1 + factor(trial) + age + factor(sex) + sbp,
                data = train1,family = binomial,control = glm.control(maxit=100000))
    coef0 <- mod0$coefficients
    coef1 <- mod1$coefficients
    
    K.si0 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.si0[2:4,ind,p]) <- 1
      K.si0[1,p,p] <- 1
    }
    allcoef0 <- t(apply(K.si0, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.si0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
    
    si_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = si_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(si_meta0$vcov)
    
    K.si1 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.si1[2:4,ind,p]) <- 1
      K.si1[1,p,p] <- 1
    }
    allcoef1 <- t(apply(K.si1, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.si1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
    
    si_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = si_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(si_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, si_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, si_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(si_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(si_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.si[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
    c.ben = mean(c.ben),
    c.ben.se = mean(c.ben.se),
    a = mean(a),
    b = mean(b),
    mse = mean(mse.si),
    in_int = mean(in_int),
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
res.si.tl.sim9 <- t(res.si.tl9[1:14, ])
se.si.tl.sim9 <- t(res.si.tl9[15:22, ])
colnames(res.si.tl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.si.tl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl9 = foreach(j = 1:m, .combine = "cbind", .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","mvmeta","msm")) %dopar%{
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
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- -0.3
  b5 <- 0.015
  b6 <- 0.1
  b7 <- -0.008
  
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
    b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
    trt_het[df$trial]*df$treat
  df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                    b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                    trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  df0 <- df[df$treat==0,]
  df1 <- df[df$treat==1,]
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat0 <- matrix(NA, nrow = 4, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 4, ncol = k)
  se_mat0 <- matrix(NA, nrow = 4, ncol = k)
  se_mat1 <- matrix(NA, nrow = 4, ncol = k)
  
  for(i in 1:k){
    test <- df[df$trial==i,]
    train0 <- df0[df0$trial!=i,]
    train1 <- df1[df1$trial!=i,]
    
    #applying the model to train
    mod0 <- rrvglm(death ~ -1 + factor(trial) + age + factor(sex) + sbp, family = binomialff, data = train0,Rank = 1)
    mod1 <- rrvglm(death ~ -1 + factor(trial) + age + factor(sex) + sbp, family = binomialff, data = train1,Rank = 1)
    coef0 <- mod0@coefficients
    coef1 <- mod1@coefficients
    
    K.r10 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.r10[2:4,ind,p]) <- 1
      K.r10[1,p,p] <- 1
    }
    allcoef0 <- t(apply(K.r10, 3, function(x){x %*% coef0}))
    allvar0 <- apply(K.r10, 3, function(x){x %*% vcovvlm(mod0) %*% t(x)}, simplify=F)
    
    r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
    coef0 = r1_meta0$coefficients
    
    coef_mat0[,i] <- coef0
    se_mat0[,i] <- diag(r1_meta0$vcov)
    
    K.r11 <- array(0, dim=c(4, 9, 6))
    ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T)
    for(p in 1:6){
      diag(K.r11[2:4,ind,p]) <- 1
      K.r11[1,p,p] <- 1
    }
    allcoef1 <- t(apply(K.r11, 3, function(x){x %*% coef1}))
    allvar1 <- apply(K.r11, 3, function(x){x %*% vcovvlm(mod1) %*% t(x)}, simplify=F)
    
    r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
    coef1 = r1_meta1$coefficients
    
    coef_mat1[,i] <- coef1
    se_mat1[,i] <- diag(r1_meta1$vcov)
    
    #recalibration of event rates adapted to train
    X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
    X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
    lp0 <- X0 %*% c(coef0)
    lp1 <- X1 %*% c(coef1)
    
    coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
    coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
    
    #ite estimation
    ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
      expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
    
    #ite prediction interval
    se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, r1_meta0$vcov)
    se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, r1_meta1$vcov)
    se <- max(se0,se1)
    true_val <- unlist(ite_true[i])
    pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(r1_meta0$Psi)),
                           upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(r1_meta0$Psi)))
    
    in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
    
    #c-statistic for benefit
    cstat = cstat4ben(outcome = as.numeric(test$death),
                      treatment = test$treat == 1,
                      score = ite)
    
    c.ben[i] <- cstat[1]
    c.ben.se[i] <- cstat[2]
    
    #calibration for benefit
    lm_ <- lm(ite~unlist(ite_true[i]))
    
    a[i] <- lm_$coefficients[[1]]
    b[i] <- lm_$coefficients[[2]]
    mse.r1[i] <- mse(unlist(ite_true[i]),ite)
  }
  c(
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
res.r1.tl.sim9 <- t(res.r1.tl9[1:14, ])
se.r1.tl.sim9 <- t(res.r1.tl9[15:22, ])
colnames(res.r1.tl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.r1.tl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl9 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")) %dopar% {
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
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- -0.3
                       b5 <- 0.015
                       b6 <- 0.1
                       b7 <- -0.008
                       
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       df$logodds <- b0 + trial_eff[df$trial] + b1*df$age + b2*df$sex + (sbp_het[df$trial]+b3)*df$sbp + 
                         b4*df$treat + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp)*df$treat +
                         trt_het[df$trial]*df$treat
                       df$ite <- expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp + 
                                         b4 + (b5*df$age + b6*(df$sex-0.5) + b7*df$sbp) +
                                         trt_het[df$trial]) - expit(b0 + trial_eff[df$trial] + b1*df$age + b2*(df$sex==1) + (sbp_het[df$trial]+b3)*df$sbp)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       df0 <- df[df$treat==0,]
                       df1 <- df[df$treat==1,]
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.fs <- c()
                       in_int <- c()
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
                           mod0 <- glm(death ~ -1 + factor(trial) + (age+factor(sex)+sbp)*factor(trial),
                                       data = train0,family = binomial,control = glm.control(maxit = 100000))
                           mod1 <- glm(death ~ -1 + factor(trial) + (age+factor(sex)+sbp)*factor(trial),
                                       data = train1,family = binomial,control = glm.control(maxit = 100000))
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           
                           K.fs0 <- array(0, dim=c(4, 24, 6))
                           ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, t2_age=F,
                                    t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
                                    t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
                                    t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F)
                           for(p in 1:6){
                             diag(K.fs0[2:4,ind,p]) <- 1
                             K.fs0[1,p,p] <- 1
                             if(p == 1) next
                             K.fs0[2,8+p,p] <- 1
                             K.fs0[3,13+p,p] <- 1
                             K.fs0[4,18+p,p] <- 1
                           }
                           allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
                           allvar0 <- apply(K.fs0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
                           
                           fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                           coef0 = fs_meta0$coefficients
                           
                           coef_mat0[,i] <- coef0
                           se_mat0[,i] <- diag(fs_meta0$vcov)
                           
                           K.fs1 <- array(0, dim=c(4, 24, 6))
                           ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, t2_age=F,
                                    t3_age=F, t4_age=F, t5_age=F, t6_age=F, t2_sex=F,
                                    t3_sex=F, t4_sex=F, t5_sex=F, t6_sex=F, t2_sbp=F,
                                    t3_sbp=F, t4_sbp=F, t5_sbp=F, t6_sbp=F)
                           for(p in 1:6){
                             diag(K.fs1[2:4,ind,p]) <- 1
                             K.fs1[1,p,p] <- 1
                             if(p == 1) next
                             K.fs1[2,8+p,p] <- 1
                             K.fs1[3,13+p,p] <- 1
                             K.fs1[4,18+p,p] <- 1
                           }
                           allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
                           allvar1 <- apply(K.fs1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
                           
                           fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                           coef1 = fs_meta1$coefficients
                           
                           coef_mat1[,i] <- coef1
                           se_mat1[,i] <- diag(fs_meta1$vcov)
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+factor(sex)+sbp"), data=train1)
                           lp0 <- X0 %*% c(coef0)
                           lp1 <- X1 %*% c(coef1)
                           
                           coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                           coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef1[1] + coef1[2]*test$age + coef1[3]*test$sex + coef1[4]*test$sbp) -
                             expit(coef0[1] + coef0[2]*test$age + coef0[3]*test$sex + coef0[4]*test$sbp)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef0, fs_meta0$vcov)
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4)/(1+exp(x1+x2+x3+x4))), coef1, fs_meta1$vcov)
                           se <- max(se0,se1)
                           true_val <- unlist(ite_true[i])
                           pred_int <- data.frame(lwr = ite - qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(fs_meta0$Psi)),
                                                  upr = ite + qt(.975,df.residual(mod0)) * se * sqrt(1 + diag(fs_meta0$Psi)))
                           
                           in_int[i] <- mean(with(pred_int, lwr <= true_val & upr >= true_val))
                           
                           #c-statistic for benefit
                           cstat = cstat4ben(outcome = as.numeric(test$death),
                                             treatment = test$treat == 1,
                                             score = ite)
                           
                           c.ben[i] <- cstat[1]
                           c.ben.se[i] <- cstat[2]
                           
                           #calibration for benefit
                           lm_ <- lm(ite~unlist(ite_true[i]))
                           
                           a[i] <- lm_$coefficients[[1]]
                           b[i] <- lm_$coefficients[[2]]
                           mse.fs[i] <- mse(unlist(ite_true[i]),ite)
                         }
                       })
                       if(any(class(testerror) == "try-error")) return(list(success = FALSE))
                       list(res = c(
                         c.ben = mean(c.ben),
                         c.ben.se = mean(c.ben.se),
                         a = mean(a),
                         b = mean(b),
                         mse = mean(mse.fs),
                         in_int = mean(in_int),
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
                         se_mat.8 = mean(se_mat1[4, ]) - mean(se_mat0[4, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.fs.tl.sim9 <- t(res.fs.tl9[1:14, ])
se.fs.tl.sim9 <- t(res.fs.tl9[15:22, ])
colnames(res.fs.tl.sim9) <- c("C-statistic for benefit","C-statistic for benefit SE","Intercept","Slope","MSE","Prop in int","b0","b1","b2","b3","b4","b5","b6","b7")
colnames(se.fs.tl.sim9) <- c("b0","b1","b2","b3","b4","b5","b6","b7")

save(res.na.sl.sim9,se.na.sl.sim9,
     res.re.sl.sim9,se.re.sl.sim9,
     res.si.sl.sim9,se.si.sl.sim9,
     res.r1.sl.sim9,se.r1.sl.sim9,
     res.fs.sl.sim9,se.fs.sl.sim9,
     res.na.tl.sim9,se.na.tl.sim9,
     res.re.tl.sim9,se.re.tl.sim9,
     res.si.tl.sim9,se.si.tl.sim9,
     res.r1.tl.sim9,se.r1.tl.sim9,
     res.fs.tl.sim9,se.fs.tl.sim9,
     file = "res_scenario9.Rdata")
