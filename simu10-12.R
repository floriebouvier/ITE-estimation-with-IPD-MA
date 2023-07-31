library(tidyverse) # data.frame manipulation
library(Hmisc) # cstat computation
library(lme4) # glmer
library(doParallel) # parallel computing
library(gtools) # calibration computation
library(VGAM) # rrvglm
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

#### scenario 10: 10 covariates (6 binary and 4 continuous) & sample size = 2800 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- glm(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat,
                                      data = train,family = binomial,control = glm.control(maxit = 100000))
                           coef <- mod$coefficients
                           coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                           se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                           
                           #ite estimation
                           ite <- predict(mod, newdata=test1, type="response") -
                             predict(mod, newdata=test0, type="response")
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, vcov(mod))
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.sl.sim10 <- t(res.na.sl10[1:26, ])
se.na.sl.sim10 <- t(res.na.sl10[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying the model to train
                           mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*factor(treat)+(1|trial)+(0+treat|trial),
                                        data = train,family = binomial,
                                        control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           coef <- mod@beta
                           coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                           se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                           lp <- X %*% coef
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                          coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                          coef[10]*test$height+coef[11]+
                                          coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                          coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                          coef[20]*test$height) -
                             expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                     coef[10]*test$height)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, vcov(mod))
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.sl.sim10 <- t(res.re.sl10[1:26, ])
se.re.sl.sim10 <- t(res.re.sl10[27:47, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.si <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- glm(death ~ -1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat+treat:factor(trial),
                                      data = train,family = "binomial",control = glm.control(maxit = 100000))
                           coef <- mod$coefficients
                           
                           K.si <- array(0, dim=c(20, 30, 6))
                           ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, mi=T,
                                    stroke=T,smoke=T,diab=T,lvhn1=T,height=T,treat=T,
                                    age_treat=T, sex_treat=T, sbp_treat=T, mi_treat=T,stroke_treat=T,
                                    smoke_treat=T,diab_treat=T,lvhn1_treat=T,height_treat=T,
                                    t2_treat=F,t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
                           for(p in 1:6){
                             diag(K.si[2:20,ind,p]) <- 1
                             K.si[1,p,p] <- 1
                             if(p == 1) next
                             K.si[11,24+p,p] <- 1
                           }
                           allcoef <- t(apply(K.si, 3, function(x){x %*% coef}))
                           allvar <- apply(K.si, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
                           
                           si_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
                           coef = si_meta$coefficients
                           
                           coef_mat[,i] <- coef
                           se_mat[,i] <- diag(si_meta$vcov)
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                           lp <- X %*% c(coef)
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                          coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                          coef[10]*test$height+coef[11]+
                                          coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                          coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                          coef[20]*test$height) -
                             expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                     coef[10]*test$height)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, si_meta$vcov)
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.sl.sim10 <- t(res.si.sl10[1:26, ])
se.si.sl.sim10 <- t(res.si.sl10[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")
) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  nt <- n/k
  
  trial <- rep(1:k, each = nt)
  df <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df$age <- df$age-mean(df$age)
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df$sbp <- df$sbp-mean(df$sbp)
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df$height <- df$height-mean(df$height)
  
  df$treat <- rep(c(0,1), times = n/(2*k))
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- -0.3
  b11 <- 0.015
  b12 <- 0.04
  b13 <- 0.1
  b14 <- -0.008
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  smoke_het <- 0.2 * pman
  smoke_het <- smoke_het - mean(smoke_het)
  
  diab_het <- (mean_age-60)/60
  diab_het <- diab_het - mean(diab_het)
  
  df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
    (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
  df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                    b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
    expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
            (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 20, ncol = k)
  se_mat <- matrix(NA, nrow = 20, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df[df$trial==i,]
      train <- df[df$trial!=i,]
      
      test0 <- test
      test0$treat <- 0
      test1 <- test
      test1$treat <- 1
      
      #applying model to train
      mod <- rrvglm(death ~ -1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat+treat:factor(trial),
                    family = binomialff, data = train,Rank = 1)
      coef <- mod@coefficients
      
      K.r1 <- array(0, dim=c(20, 30, 6))
      ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, mi=T,
               stroke=T,smoke=T,diab=T,lvhn1=T,height=T,treat=T,
               age_treat=T, sex_treat=T, sbp_treat=T, mi_treat=T,stroke_treat=T,
               smoke_treat=T,diab_treat=T,lvhn1_treat=T,height_treat=T,
               t2_treat=F,t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
      for(p in 1:6){
        diag(K.r1[2:20,ind,p]) <- 1
        K.r1[1,p,p] <- 1
        if(p == 1) next
        K.r1[11,24+p,p] <- 1
      }
      allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
      allvar <- apply(K.r1, 3, function(x){x %*% vcovvlm(mod) %*% t(x)}, simplify=F)
      
      r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
      coef = r1_meta$coefficients
      
      coef_mat[,i] <- coef
      se_mat[,i] <- diag(r1_meta$vcov)
      
      #recalibration of event rates adapted to train
      X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
      lp <- X %*% c(coef)
      coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
      
      #ite estimation
      ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                     coef[10]*test$height+coef[11]+
                     coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                     coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                     coef[20]*test$height) -
        expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                coef[10]*test$height)
      
      #ite prediction interval
      se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                          (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, r1_meta$vcov)
      
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
    se_mat.20 = mean(se_mat[20, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl.sim10 <- t(res.r1.sl10[1:26, ])
se.r1.sl.sim10 <- t(res.r1.sl10[27:47, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl10 = foreach(j = 1:m,
                      .final = function(l) {do.call("cbind", lapply(
                        l[sapply(l,`[[`, "success")],
                        function(subl) c(subl$res, subl$seed)))},
                      .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")) %dopar% {
    set.seed(j)
n <- 2800
k <- 7
nt <- n/k

trial <- rep(1:k, each = nt)
df <- as.data.frame(trial)

mean_age <- c(52,56,64,70,77,78,82) #Age
sd_age <- c(4,2,1,3,4,6,2)
df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
df$age <- df$age-mean(df$age)

pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
df$sex <- rbinom(n, 1, prob=pman[trial])

mean_sbp <- c(186,182,170,185,190,188,197) #SBP
sd_sbp <- c(9,11,5,12,9,10,16)
df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
df$sbp <- df$sbp-mean(df$sbp)

pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction
df$mi <- rbinom(n, 1, prob=pmi[trial])

pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
df$stroke <- rbinom(n, 1, prob=pstroke[trial])

psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
df$smoke <- rbinom(n, 1, prob=psmoke[trial])

pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
df$diab <- rbinom(n, 1, prob=pdiab[trial])

plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])

mean_height <- c(176,162,167,169,168,170,167) #Height
sd_height <- c(6,9,10,10,10,9,9)
df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
df$height <- df$height-mean(df$height)

df$treat <- rep(c(0,1), times = n/(2*k))

b0 <- -1.4
b1 <- 0.03
b2 <- 0.7
b3 <- 0.02
b4 <- 0.82
b5 <- 0.8
b6 <- 0.7
b7 <- 0.1
b8 <- 0.33
b9 <- -0.02
b10 <- -0.3
b11 <- 0.015
b12 <- 0.04
b13 <- 0.1
b14 <- -0.008

trial_eff <- (mean_age-60)/20
trial_eff <- trial_eff - mean(trial_eff)

trt_het <- (mean_age-60)/100
trt_het <- trt_het - mean(trt_het)

sbp_het <- 0.2 * pman
sbp_het <- sbp_het - mean(sbp_het)

smoke_het <- 0.2 * pman
smoke_het <- smoke_het - mean(smoke_het)

diab_het <- (mean_age-60)/60
diab_het <- diab_het - mean(diab_het)

df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
  (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
  (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                  (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                  b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
  expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
          (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
pout <- plogis(df$logodds)
df$death <- rbinom(n, 1, pout)

    c.ben <- c()
    c.ben.se <- c()
    a <- c()
    b <- c()
    mse.fs <- c()
    in_int <- c()
    coef_mat <- matrix(NA, nrow = 20, ncol = k)
    se_mat <- matrix(NA, nrow = 20, ncol = k)

    testerror = try({
      for(i in 1:k){
        test <- df[df$trial==i,]
        train <- df[df$trial!=i,]

        test0 <- test
        test0$treat <- 0
        test1 <- test
        test1$treat <- 1

        #applying model to train
        mod <- glm(death~-1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+
                                          height)*treat*factor(trial),
                   data = train,family = binomial,control = glm.control(maxit = 100000))
        coef <- mod$coefficients
        coef <- coef[!is.na(coef)]

        K.fs <- array(0, dim=c(20, length(coef), 6))
        ind <- c(rep(F,6),rep(T,19),rep(F,length(coef)-25))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs[2:20,ind,p]) <- 1
          K.fs[1,p,p] <- 1
          if(p == 1) next
          K.fs[2,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs[3,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs[4,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs[5,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs[6,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs[7,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs[8,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs[9,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs[10,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
          K.fs[11,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":treat"))-2,p] <- 1
          K.fs[12,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"age:treat"))-2,p] <- 1
          K.fs[13,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(sex)1:treat"))-2,p] <- 1
          K.fs[14,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp:treat"))-2,p] <- 1
          K.fs[15,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1:treat"))-2,p] <- 1
          K.fs[16,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1:treat"))-2,p] <- 1
          K.fs[17,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(smoke)1:treat"))-2,p] <- 1
          K.fs[18,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(diab)1:treat"))-2,p] <- 1
          K.fs[19,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"lvhn1:treat"))-2,p] <- 1
          K.fs[20,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"height:treat"))-2,p] <- 1
        }
        allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
        
        X <- model.matrix(death~-1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+
                                                    factor(diab)+lvhn1+height)*treat*factor(trial),data = train)
        idx <- which(colnames(X) %in% names(coef))
        
        X2 <- X[,idx]
        allvar <- apply(K.fs, 3, function(x){x %*% cov(X2) %*% t(x)}, simplify=F)
        
        fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
        coef = fs_meta$coefficients
        
        coef_mat[,i] <- coef
        se_mat[,i] <- diag(fs_meta$vcov)

        #recalibration of event rates adapted to train
        X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
        lp <- X %*% c(coef)
        coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]

        #ite estimation
        ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                       coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                       coef[10]*test$height+coef[11]+
                       coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                       coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                       coef[20]*test$height) -
          expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                  coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                  coef[10]*test$height)

        #ite prediction interval
        se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) -
                            (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, fs_meta$vcov)

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
      se_mat.20 = mean(se_mat[20, ])),
      seed = j,
      success = TRUE)
  }
stopCluster(cl)
res.fs.sl.sim10 <- t(res.fs.sl10[1:26, ])
se.fs.sl.sim10 <- t(res.fs.sl10[27:47, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train0,family = binomial,
                                       control = glm.control(maxit=100000))
                           mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
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
                           
                           #ite estimation
                           ite <- predict(mod1, newdata=test, type="response") -
                             predict(mod0, newdata=test, type="response")
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, vcov(mod1))
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.tl.sim10 <- t(res.na.tl10[1:26, ])
se.na.tl.sim10 <- t(res.na.tl10[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+(1|trial),
                                         data = train0,family = binomial,
                                         control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           mod1 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+(1|trial),
                                         data = train1,family = binomial,
                                         control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           coef0 <- mod0@beta
                           coef1 <- mod1@beta
                           coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
                           coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
                           se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
                           se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
                           lp0 <- X0 %*% coef0
                           lp1 <- X1 %*% coef1
                           
                           coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                           coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                                          coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                                          coef1[10]*test$height) -
                             expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                                     coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                                     coef0[10]*test$height)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, vcov(mod1))
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.tl.sim10 <- t(res.re.tl10[1:26, ])
se.re.tl.sim10 <- t(res.re.tl10[27:47, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 2800
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train0,family = binomial,
                                       control = glm.control(maxit=100000))
                           mod1 <- glm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train1,family = binomial,
                                       control = glm.control(maxit=100000))
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef0[is.na(coef0)] <- 0
                           coef1[is.na(coef1)] <- 0
                           
                           K.si0 <- array(0, dim=c(10, 15, 6))
                           ind <- c(rep(F,6),rep(T,9))
                           for(p in 1:6){
                             diag(K.si0[2:10,ind,p]) <- 1
                             K.si0[1,p,p] <- 1
                           }
                           allcoef0 <- t(apply(K.si0, 3, function(x){x %*% coef0}))
                           allvar0 <- apply(K.si0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
                           
                           si_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                           coef0 = si_meta0$coefficients
                           
                           coef_mat0[,i] <- coef0
                           se_mat0[,i] <- diag(si_meta0$vcov)
                           
                           K.si1 <- array(0, dim=c(10, 15, 6))
                           ind <- c(rep(F,6),rep(T,9))
                           for(p in 1:6){
                             diag(K.si1[2:10,ind,p]) <- 1
                             K.si1[1,p,p] <- 1
                           }
                           allcoef1 <- t(apply(K.si1, 3, function(x){x %*% coef1}))
                           allvar1 <- apply(K.si1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
                           
                           si_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                           coef1 = si_meta1$coefficients
                           
                           coef_mat1[,i] <- coef1
                           se_mat1[,i] <- diag(si_meta1$vcov)
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
                           lp0 <- X0 %*% c(coef0)
                           lp1 <- X1 %*% c(coef1)
                           
                           #ite estimation
                           ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                                          coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                                          coef1[10]*test$height) -
                             expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                                     coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                                     coef0[10]*test$height)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, si_meta0$vcov)
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, si_meta1$vcov)
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.tl.sim10 <- t(res.si.tl10[1:26, ])
se.si.tl.sim10 <- t(res.si.tl10[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl10 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")
) %dopar% {
  set.seed(j)
  n <- 2800
  k <- 7
  nt <- n/k
  
  trial <- rep(1:k, each = nt)
  df <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df$age <- df$age-mean(df$age)
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df$sbp <- df$sbp-mean(df$sbp)
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df$height <- df$height-mean(df$height)
  
  df$treat <- rep(c(0,1), times = n/(2*k))
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- -0.3
  b11 <- 0.015
  b12 <- 0.04
  b13 <- 0.1
  b14 <- -0.008
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  smoke_het <- 0.2 * pman
  smoke_het <- smoke_het - mean(smoke_het)
  
  diab_het <- (mean_age-60)/60
  diab_het <- diab_het - mean(diab_het)
  
  df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
    (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
  df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                    b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
    expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
            (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
  coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
  se_mat0 <- matrix(NA, nrow = 10, ncol = k)
  se_mat1 <- matrix(NA, nrow = 10, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df[df$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- rrvglm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height, family = binomialff, data = train0,Rank = 1)
      mod1 <- rrvglm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height, family = binomialff, data = train1,Rank = 1)
      
      coef0 <- mod0@coefficients
      coef1 <- mod1@coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      
      K.r10 <- array(0, dim=c(10, 15, 6))
      ind <- c(rep(F,6),rep(T,9))
      for(p in 1:6){
        diag(K.r10[2:10,ind,p]) <- 1
        K.r10[1,p,p] <- 1
      }
      allcoef0 <- t(apply(K.r10, 3, function(x){x %*% coef0}))
      allvar0 <- apply(K.r10, 3, function(x){x %*% vcovvlm(mod0) %*% t(x)}, simplify=F)
      
      r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
      coef0 = r1_meta0$coefficients
      
      coef_mat0[,i] <- coef0
      se_mat0[,i] <- diag(r1_meta0$vcov)
      
      K.r11 <- array(0, dim=c(10, 15, 6))
      ind <- c(rep(F,6),rep(T,9))
      for(p in 1:6){
        diag(K.r11[2:10,ind,p]) <- 1
        K.r11[1,p,p] <- 1
      }
      allcoef1 <- t(apply(K.r11, 3, function(x){x %*% coef1}))
      allvar1 <- apply(K.r11, 3, function(x){x %*% vcovvlm(mod1) %*% t(x)}, simplify=F)
      
      r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
      coef1 = r1_meta1$coefficients
      
      coef_mat1[,i] <- coef1
      se_mat1[,i] <- diag(r1_meta1$vcov)
      
      #recalibration of event rates adapted to train
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
      lp0 <- X0 %*% c(coef0)
      lp1 <- X1 %*% c(coef1)
      
      #ite estimation
      ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                     coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                     coef1[10]*test$height) -
        expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                coef0[10]*test$height)
      
      #ite prediction interval
      se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, r1_meta0$vcov)
      se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, r1_meta1$vcov)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(
    res = c(c.ben = mean(c.ben),
            c.ben.se = mean(c.ben.se),
            a = mean(a),
            b = mean(b),
            mse = mean(mse.r1),
            in_int = mean(in_int),
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
            coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
            coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
            coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
            coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
            coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
            coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
            coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
            coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
            coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
            coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
            se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
            se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
            se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
            se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
            se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
            se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
            se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
            se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
            se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
            se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.tl.sim10 <- t(res.r1.tl10[1:26, ])
se.r1.tl.sim10 <- t(res.r1.tl10[27:47, ])

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
  .packages = c("tidyverse", "Hmisc", "lme4", "magrittr", "gtools", "VGAM","msm","mvmeta")) %dopar%
  {
set.seed(j)
n <- 2800
k <- 7
nt <- n/k

trial <- rep(1:k, each = nt)
df <- as.data.frame(trial)

mean_age <- c(52,56,64,70,77,78,82) #Age
sd_age <- c(4,2,1,3,4,6,2)
df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
df$age <- df$age-mean(df$age)

pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
df$sex <- rbinom(n, 1, prob=pman[trial])

mean_sbp <- c(186,182,170,185,190,188,197) #SBP
sd_sbp <- c(9,11,5,12,9,10,16)
df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
df$sbp <- df$sbp-mean(df$sbp)

pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction
df$mi <- rbinom(n, 1, prob=pmi[trial])

pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
df$stroke <- rbinom(n, 1, prob=pstroke[trial])

psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
df$smoke <- rbinom(n, 1, prob=psmoke[trial])

pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
df$diab <- rbinom(n, 1, prob=pdiab[trial])

plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])

mean_height <- c(176,162,167,169,168,170,167) #Height
sd_height <- c(6,9,10,10,10,9,9)
df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
df$height <- df$height-mean(df$height)

df$treat <- rep(c(0,1), times = n/(2*k))

b0 <- -1.4
b1 <- 0.03
b2 <- 0.7
b3 <- 0.02
b4 <- 0.82
b5 <- 0.8
b6 <- 0.7
b7 <- 0.1
b8 <- 0.33
b9 <- -0.02
b10 <- -0.3
b11 <- 0.015
b12 <- 0.04
b13 <- 0.1
b14 <- -0.008

trial_eff <- (mean_age-60)/20
trial_eff <- trial_eff - mean(trial_eff)

trt_het <- (mean_age-60)/100
trt_het <- trt_het - mean(trt_het)

sbp_het <- 0.2 * pman
sbp_het <- sbp_het - mean(sbp_het)

smoke_het <- 0.2 * pman
smoke_het <- smoke_het - mean(smoke_het)

diab_het <- (mean_age-60)/60
diab_het <- diab_het - mean(diab_het)

df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
  (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
  (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                  (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                  b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
  expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
          (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)

    testerror = try({
      for(i in 1:k){
        test <- df[df$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]

        #applying the model to train
        mod0 <- glm(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),
                    data = train0,family = binomial,control = glm.control(maxit = 100000))
        mod1 <- glm(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),
                    data = train1,family = binomial,control = glm.control(maxit = 100000))
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0 <- coef0[!is.na(coef0)]
        coef1 <- coef1[!is.na(coef1)]
        
        K.fs0 <- array(0, dim=c(10, length(coef0), 6))
        ind <- c(rep(F,6),rep(T,9),rep(F,length(coef0)-15))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs0[2:10,ind,p]) <- 1
          K.fs0[1,p,p] <- 1
          if(p == 1) next
          K.fs0[2,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs0[3,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs0[4,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs0[5,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs0[6,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs0[7,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs0[8,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs0[9,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs0[10,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
        }
        allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
        
        X0 <- model.matrix(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),data = train0)
        idx0 <- which(colnames(X0) %in% names(coef0))
        X0.2 <- X0[,idx0]
        allvar0 <- apply(K.fs0, 3, function(x){x %*% cov(X0.2) %*% t(x)}, simplify=F)
        
        fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
        coef0 = fs_meta0$coefficients
        
        coef_mat0[,i] <- coef0
        se_mat0[,i] <- diag(fs_meta0$vcov)
        
        K.fs1 <- array(0, dim=c(10, length(coef1), 6))
        ind <- c(rep(F,6),rep(T,9),rep(F,length(coef1)-15))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs1[2:10,ind,p]) <- 1
          K.fs1[1,p,p] <- 1
          if(p == 1) next
          K.fs1[2,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs1[3,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs1[4,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs1[5,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs1[6,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs1[7,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs1[8,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs1[9,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs1[10,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
        }
        allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
        
        X1 <- model.matrix(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),data = train1)
        idx1 <- which(colnames(X1) %in% names(coef1))
        X1.2 <- X1[,idx1]
        allvar1 <- apply(K.fs1, 3, function(x){x %*% cov(X1.2) %*% t(x)}, simplify=F)
        
        fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
        coef1 = fs_meta1$coefficients
        
        coef_mat1[,i] <- coef1
        se_mat1[,i] <- diag(fs_meta1$vcov)
        
        #recalibration of event rates adapted to train
        X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
        X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
        lp0 <- X0 %*% c(coef0)
        lp1 <- X1 %*% c(coef1)
        
        coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
        coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]

        #ite estimation
        ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                       coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                       coef1[10]*test$height) -
          expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                  coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                  coef0[10]*test$height)

        #ite prediction interval
        se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, fs_meta0$vcov)
        se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, fs_meta1$vcov)
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
    if (any(class(testerror) == "try-error")) return(list(success = FALSE))

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
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
        coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
        se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl.sim10 <- t(res.fs.tl10[1:26, ])
se.fs.tl.sim10 <- t(res.fs.tl10[27:47, ])

save(res.na.sl.sim10,se.na.sl.sim10,
     res.re.sl.sim10,se.re.sl.sim10,
     res.si.sl.sim10,se.si.sl.sim10,
     res.r1.sl.sim10,se.r1.sl.sim10,
     res.fs.sl.sim10,se.fs.sl.sim10,
     res.na.tl.sim10,se.na.tl.sim10,
     res.re.tl.sim10,se.re.tl.sim10,
     res.si.tl.sim10,se.si.tl.sim10,
     res.r1.tl.sim10,se.r1.tl.sim10,
     res.fs.tl.sim10,se.fs.tl.sim10,
     file = "res_scenario10.Rdata")

#### scenario 11: 10 covariates (6 binary and 4 continuous) & sample size = 1400 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- glm(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat,
                                      data = train,family = binomial,control = glm.control(maxit = 100000))
                           coef <- mod$coefficients
                           coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                           se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                           
                           #ite estimation
                           ite <- predict(mod, newdata=test1, type="response") -
                             predict(mod, newdata=test0, type="response")
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, vcov(mod))
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.sl.sim11 <- t(res.na.sl11[1:26, ])
se.na.sl.sim11 <- t(res.na.sl11[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying the model to train
                           mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*factor(treat)+(1|trial)+(0+treat|trial),
                                        data = train,family = binomial,
                                        control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           coef <- mod@beta
                           coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                           se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                           lp <- X %*% coef
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                          coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                          coef[10]*test$height+coef[11]+
                                          coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                          coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                          coef[20]*test$height) -
                             expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                     coef[10]*test$height)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, vcov(mod))
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.sl.sim11 <- t(res.re.sl11[1:26, ])
se.re.sl.sim11 <- t(res.re.sl11[27:47, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.si <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- glm(death ~ -1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat+treat:factor(trial),
                                      data = train,family = "binomial",control = glm.control(maxit = 100000))
                           coef <- mod$coefficients
                           
                           K.si <- array(0, dim=c(20, 30, 6))
                           ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, mi=T,
                                    stroke=T,smoke=T,diab=T,lvhn1=T,height=T,treat=T,
                                    age_treat=T, sex_treat=T, sbp_treat=T, mi_treat=T,stroke_treat=T,
                                    smoke_treat=T,diab_treat=T,lvhn1_treat=T,height_treat=T,
                                    t2_treat=F,t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
                           for(p in 1:6){
                             diag(K.si[2:20,ind,p]) <- 1
                             K.si[1,p,p] <- 1
                             if(p == 1) next
                             K.si[11,24+p,p] <- 1
                           }
                           allcoef <- t(apply(K.si, 3, function(x){x %*% coef}))
                           allvar <- apply(K.si, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
                           
                           si_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
                           coef = si_meta$coefficients
                           
                           coef_mat[,i] <- coef
                           se_mat[,i] <- diag(si_meta$vcov)
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                           lp <- X %*% c(coef)
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                          coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                          coef[10]*test$height+coef[11]+
                                          coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                          coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                          coef[20]*test$height) -
                             expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                     coef[10]*test$height)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, si_meta$vcov)
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.sl.sim11 <- t(res.si.sl11[1:26, ])
se.si.sl.sim11 <- t(res.si.sl11[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")
) %dopar% {
  set.seed(j)
  n <- 1400
  k <- 7
  nt <- n/k
  
  trial <- rep(1:k, each = nt)
  df <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df$age <- df$age-mean(df$age)
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df$sbp <- df$sbp-mean(df$sbp)
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df$height <- df$height-mean(df$height)
  
  df$treat <- rep(c(0,1), times = n/(2*k))
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- -0.3
  b11 <- 0.015
  b12 <- 0.04
  b13 <- 0.1
  b14 <- -0.008
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  smoke_het <- 0.2 * pman
  smoke_het <- smoke_het - mean(smoke_het)
  
  diab_het <- (mean_age-60)/60
  diab_het <- diab_het - mean(diab_het)
  
  df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
    (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
  df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                    b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
    expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
            (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 20, ncol = k)
  se_mat <- matrix(NA, nrow = 20, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df[df$trial==i,]
      train <- df[df$trial!=i,]
      
      test0 <- test
      test0$treat <- 0
      test1 <- test
      test1$treat <- 1
      
      #applying model to train
      mod <- rrvglm(death ~ -1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat+treat:factor(trial),
                    family = binomialff, data = train,Rank = 1)
      coef <- mod@coefficients
      
      K.r1 <- array(0, dim=c(20, 30, 6))
      ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, mi=T,
               stroke=T,smoke=T,diab=T,lvhn1=T,height=T,treat=T,
               age_treat=T, sex_treat=T, sbp_treat=T, mi_treat=T,stroke_treat=T,
               smoke_treat=T,diab_treat=T,lvhn1_treat=T,height_treat=T,
               t2_treat=F,t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
      for(p in 1:6){
        diag(K.r1[2:20,ind,p]) <- 1
        K.r1[1,p,p] <- 1
        if(p == 1) next
        K.r1[11,24+p,p] <- 1
      }
      allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
      allvar <- apply(K.r1, 3, function(x){x %*% vcovvlm(mod) %*% t(x)}, simplify=F)
      
      r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
      coef = r1_meta$coefficients
      
      coef_mat[,i] <- coef
      se_mat[,i] <- diag(r1_meta$vcov)
      
      #recalibration of event rates adapted to train
      X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
      lp <- X %*% c(coef)
      coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
      
      #ite estimation
      ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                     coef[10]*test$height+coef[11]+
                     coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                     coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                     coef[20]*test$height) -
        expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                coef[10]*test$height)
      
      #ite prediction interval
      se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                          (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, r1_meta$vcov)
      
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
    se_mat.20 = mean(se_mat[20, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl.sim11 <- t(res.r1.sl11[1:26, ])
se.r1.sl.sim11 <- t(res.r1.sl11[27:47, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl11 = foreach(j = 1:m,
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
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age-mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp-mean(df$sbp)
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction
                        df$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df$diab <- rbinom(n, 1, prob=pdiab[trial])
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df$height <- df$height-mean(df$height)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- -0.3
                        b11 <- 0.015
                        b12 <- 0.04
                        b13 <- 0.1
                        b14 <- -0.008
                        
                        trial_eff <- (mean_age-60)/20
                        trial_eff <- trial_eff - mean(trial_eff)
                        
                        trt_het <- (mean_age-60)/100
                        trt_het <- trt_het - mean(trt_het)
                        
                        sbp_het <- 0.2 * pman
                        sbp_het <- sbp_het - mean(sbp_het)
                        
                        smoke_het <- 0.2 * pman
                        smoke_het <- smoke_het - mean(smoke_het)
                        
                        diab_het <- (mean_age-60)/60
                        diab_het <- diab_het - mean(diab_het)
                        
                        df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                          (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                          (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                        df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                          (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                          b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                          expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                  (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        pout <- plogis(df$logodds)
                        df$death <- rbinom(n, 1, pout)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 20, ncol = k)
                        se_mat <- matrix(NA, nrow = 20, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~-1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+
                                                                 height)*treat*factor(trial),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef <- coef[!is.na(coef)]
                            
                            K.fs <- array(0, dim=c(20, length(coef), 6))
                            ind <- c(rep(F,6),rep(T,19),rep(F,length(coef)-25))
                            if(i==1){vec <- 3:6} else {vec <- 2:7}
                            for(p in 1:6){
                              diag(K.fs[2:20,ind,p]) <- 1
                              K.fs[1,p,p] <- 1
                              if(p == 1) next
                              K.fs[2,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
                              K.fs[3,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
                              K.fs[4,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
                              K.fs[5,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
                              K.fs[6,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
                              K.fs[7,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
                              K.fs[8,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
                              K.fs[9,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
                              K.fs[10,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
                              K.fs[11,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":treat"))-2,p] <- 1
                              K.fs[12,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"age:treat"))-2,p] <- 1
                              K.fs[13,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(sex)1:treat"))-2,p] <- 1
                              K.fs[14,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp:treat"))-2,p] <- 1
                              K.fs[15,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1:treat"))-2,p] <- 1
                              K.fs[16,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1:treat"))-2,p] <- 1
                              K.fs[17,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(smoke)1:treat"))-2,p] <- 1
                              K.fs[18,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(diab)1:treat"))-2,p] <- 1
                              K.fs[19,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"lvhn1:treat"))-2,p] <- 1
                              K.fs[20,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"height:treat"))-2,p] <- 1
                            }
                            allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
                            
                            X <- model.matrix(death~-1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+
                                                                        factor(diab)+lvhn1+height)*treat*factor(trial),data = train)
                            idx <- which(colnames(X) %in% names(coef))
                            
                            X2 <- X[,idx]
                            allvar <- apply(K.fs, 3, function(x){x %*% cov(X2) %*% t(x)}, simplify=F)
                            
                            fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
                            coef = fs_meta$coefficients
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(fs_meta$vcov)
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                            lp <- X %*% c(coef)
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                           coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                           coef[10]*test$height+coef[11]+
                                           coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                           coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                           coef[20]*test$height) -
                              expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                      coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                      coef[10]*test$height)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) -
                                                (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, fs_meta$vcov)
                            
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
                          se_mat.20 = mean(se_mat[20, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.fs.sl.sim11 <- t(res.fs.sl11[1:26, ])
se.fs.sl.sim11 <- t(res.fs.sl11[27:47, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train0,family = binomial,
                                       control = glm.control(maxit=100000))
                           mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
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
                           
                           #ite estimation
                           ite <- predict(mod1, newdata=test, type="response") -
                             predict(mod0, newdata=test, type="response")
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, vcov(mod1))
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.tl.sim11 <- t(res.na.tl11[1:26, ])
se.na.tl.sim11 <- t(res.na.tl11[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+(1|trial),
                                         data = train0,family = binomial,
                                         control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           mod1 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+(1|trial),
                                         data = train1,family = binomial,
                                         control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           coef0 <- mod0@beta
                           coef1 <- mod1@beta
                           coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
                           coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
                           se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
                           se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
                           lp0 <- X0 %*% coef0
                           lp1 <- X1 %*% coef1
                           
                           coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                           coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                                          coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                                          coef1[10]*test$height) -
                             expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                                     coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                                     coef0[10]*test$height)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, vcov(mod1))
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.tl.sim11 <- t(res.re.tl11[1:26, ])
se.re.tl.sim11 <- t(res.re.tl11[27:47, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 1400
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train0,family = binomial,
                                       control = glm.control(maxit=100000))
                           mod1 <- glm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train1,family = binomial,
                                       control = glm.control(maxit=100000))
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef0[is.na(coef0)] <- 0
                           coef1[is.na(coef1)] <- 0
                           
                           K.si0 <- array(0, dim=c(10, 15, 6))
                           ind <- c(rep(F,6),rep(T,9))
                           for(p in 1:6){
                             diag(K.si0[2:10,ind,p]) <- 1
                             K.si0[1,p,p] <- 1
                           }
                           allcoef0 <- t(apply(K.si0, 3, function(x){x %*% coef0}))
                           allvar0 <- apply(K.si0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
                           
                           si_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                           coef0 = si_meta0$coefficients
                           
                           coef_mat0[,i] <- coef0
                           se_mat0[,i] <- diag(si_meta0$vcov)
                           
                           K.si1 <- array(0, dim=c(10, 15, 6))
                           ind <- c(rep(F,6),rep(T,9))
                           for(p in 1:6){
                             diag(K.si1[2:10,ind,p]) <- 1
                             K.si1[1,p,p] <- 1
                           }
                           allcoef1 <- t(apply(K.si1, 3, function(x){x %*% coef1}))
                           allvar1 <- apply(K.si1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
                           
                           si_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                           coef1 = si_meta1$coefficients
                           
                           coef_mat1[,i] <- coef1
                           se_mat1[,i] <- diag(si_meta1$vcov)
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
                           lp0 <- X0 %*% c(coef0)
                           lp1 <- X1 %*% c(coef1)
                           
                           #ite estimation
                           ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                                          coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                                          coef1[10]*test$height) -
                             expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                                     coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                                     coef0[10]*test$height)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, si_meta0$vcov)
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, si_meta1$vcov)
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.tl.sim11 <- t(res.si.tl11[1:26, ])
se.si.tl.sim11 <- t(res.si.tl11[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl11 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")
) %dopar% {
  set.seed(j)
  n <- 1400
  k <- 7
  nt <- n/k
  
  trial <- rep(1:k, each = nt)
  df <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df$age <- df$age-mean(df$age)
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df$sbp <- df$sbp-mean(df$sbp)
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df$height <- df$height-mean(df$height)
  
  df$treat <- rep(c(0,1), times = n/(2*k))
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- -0.3
  b11 <- 0.015
  b12 <- 0.04
  b13 <- 0.1
  b14 <- -0.008
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  smoke_het <- 0.2 * pman
  smoke_het <- smoke_het - mean(smoke_het)
  
  diab_het <- (mean_age-60)/60
  diab_het <- diab_het - mean(diab_het)
  
  df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
    (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
  df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                    b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
    expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
            (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
  coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
  se_mat0 <- matrix(NA, nrow = 10, ncol = k)
  se_mat1 <- matrix(NA, nrow = 10, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df[df$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- rrvglm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height, family = binomialff, data = train0,Rank = 1)
      mod1 <- rrvglm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height, family = binomialff, data = train1,Rank = 1)
      
      coef0 <- mod0@coefficients
      coef1 <- mod1@coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      
      K.r10 <- array(0, dim=c(10, 15, 6))
      ind <- c(rep(F,6),rep(T,9))
      for(p in 1:6){
        diag(K.r10[2:10,ind,p]) <- 1
        K.r10[1,p,p] <- 1
      }
      allcoef0 <- t(apply(K.r10, 3, function(x){x %*% coef0}))
      allvar0 <- apply(K.r10, 3, function(x){x %*% vcovvlm(mod0) %*% t(x)}, simplify=F)
      
      r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
      coef0 = r1_meta0$coefficients
      
      coef_mat0[,i] <- coef0
      se_mat0[,i] <- diag(r1_meta0$vcov)
      
      K.r11 <- array(0, dim=c(10, 15, 6))
      ind <- c(rep(F,6),rep(T,9))
      for(p in 1:6){
        diag(K.r11[2:10,ind,p]) <- 1
        K.r11[1,p,p] <- 1
      }
      allcoef1 <- t(apply(K.r11, 3, function(x){x %*% coef1}))
      allvar1 <- apply(K.r11, 3, function(x){x %*% vcovvlm(mod1) %*% t(x)}, simplify=F)
      
      r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
      coef1 = r1_meta1$coefficients
      
      coef_mat1[,i] <- coef1
      se_mat1[,i] <- diag(r1_meta1$vcov)
      
      #recalibration of event rates adapted to train
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
      lp0 <- X0 %*% c(coef0)
      lp1 <- X1 %*% c(coef1)
      
      #ite estimation
      ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                     coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                     coef1[10]*test$height) -
        expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                coef0[10]*test$height)
      
      #ite prediction interval
      se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, r1_meta0$vcov)
      se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, r1_meta1$vcov)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(
    res = c(c.ben = mean(c.ben),
            c.ben.se = mean(c.ben.se),
            a = mean(a),
            b = mean(b),
            mse = mean(mse.r1),
            in_int = mean(in_int),
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
            coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
            coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
            coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
            coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
            coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
            coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
            coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
            coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
            coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
            coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
            se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
            se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
            se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
            se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
            se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
            se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
            se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
            se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
            se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
            se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.tl.sim11 <- t(res.r1.tl11[1:26, ])
se.r1.tl.sim11 <- t(res.r1.tl11[27:47, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl11 = foreach(
  j = 1:m,
  .final = function(l) do.call("cbind",
                               lapply(
                                 l[sapply(l, `[[`, "success")],
                                 function(subl) c(subl$res, subl$seed))),
  .packages = c("tidyverse", "Hmisc", "lme4", "magrittr", "gtools", "VGAM","msm","mvmeta")) %dopar%
  {
    set.seed(j)
    n <- 1400
    k <- 7
    nt <- n/k
    
    trial <- rep(1:k, each = nt)
    df <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df$age <- df$age-mean(df$age)
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df$sbp <- df$sbp-mean(df$sbp)
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction
    df$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df$diab <- rbinom(n, 1, prob=pdiab[trial])
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df$height <- df$height-mean(df$height)
    
    df$treat <- rep(c(0,1), times = n/(2*k))
    
    b0 <- -1.4
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.02
    b4 <- 0.82
    b5 <- 0.8
    b6 <- 0.7
    b7 <- 0.1
    b8 <- 0.33
    b9 <- -0.02
    b10 <- -0.3
    b11 <- 0.015
    b12 <- 0.04
    b13 <- 0.1
    b14 <- -0.008
    
    trial_eff <- (mean_age-60)/20
    trial_eff <- trial_eff - mean(trial_eff)
    
    trt_het <- (mean_age-60)/100
    trt_het <- trt_het - mean(trt_het)
    
    sbp_het <- 0.2 * pman
    sbp_het <- sbp_het - mean(sbp_het)
    
    smoke_het <- 0.2 * pman
    smoke_het <- smoke_het - mean(smoke_het)
    
    diab_het <- (mean_age-60)/60
    diab_het <- diab_het - mean(diab_het)
    
    df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
      (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
      (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
    df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                      (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                      b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
      expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
              (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df[df$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- glm(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),
                    data = train0,family = binomial,control = glm.control(maxit = 100000))
        mod1 <- glm(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),
                    data = train1,family = binomial,control = glm.control(maxit = 100000))
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0 <- coef0[!is.na(coef0)]
        coef1 <- coef1[!is.na(coef1)]
        
        K.fs0 <- array(0, dim=c(10, length(coef0), 6))
        ind <- c(rep(F,6),rep(T,9),rep(F,length(coef0)-15))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs0[2:10,ind,p]) <- 1
          K.fs0[1,p,p] <- 1
          if(p == 1) next
          K.fs0[2,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs0[3,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs0[4,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs0[5,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs0[6,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs0[7,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs0[8,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs0[9,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs0[10,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
        }
        allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
        
        X0 <- model.matrix(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),data = train0)
        idx0 <- which(colnames(X0) %in% names(coef0))
        X0.2 <- X0[,idx0]
        allvar0 <- apply(K.fs0, 3, function(x){x %*% cov(X0.2) %*% t(x)}, simplify=F)
        
        fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
        coef0 = fs_meta0$coefficients
        
        coef_mat0[,i] <- coef0
        se_mat0[,i] <- diag(fs_meta0$vcov)
        
        K.fs1 <- array(0, dim=c(10, length(coef1), 6))
        ind <- c(rep(F,6),rep(T,9),rep(F,length(coef1)-15))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs1[2:10,ind,p]) <- 1
          K.fs1[1,p,p] <- 1
          if(p == 1) next
          K.fs1[2,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs1[3,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs1[4,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs1[5,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs1[6,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs1[7,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs1[8,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs1[9,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs1[10,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
        }
        allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
        
        X1 <- model.matrix(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),data = train1)
        idx1 <- which(colnames(X1) %in% names(coef1))
        X1.2 <- X1[,idx1]
        allvar1 <- apply(K.fs1, 3, function(x){x %*% cov(X1.2) %*% t(x)}, simplify=F)
        
        fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
        coef1 = fs_meta1$coefficients
        
        coef_mat1[,i] <- coef1
        se_mat1[,i] <- diag(fs_meta1$vcov)
        
        #recalibration of event rates adapted to train
        X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
        X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
        lp0 <- X0 %*% c(coef0)
        lp1 <- X1 %*% c(coef1)
        
        coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
        coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
        
        #ite estimation
        ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                       coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                       coef1[10]*test$height) -
          expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                  coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                  coef0[10]*test$height)
        
        #ite prediction interval
        se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, fs_meta0$vcov)
        se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, fs_meta1$vcov)
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
    if (any(class(testerror) == "try-error")) return(list(success = FALSE))
    
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
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
        coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
        se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl.sim11 <- t(res.fs.tl11[1:26, ])
se.fs.tl.sim11 <- t(res.fs.tl11[27:47, ])

save(res.na.sl.sim11,se.na.sl.sim11,
     res.re.sl.sim11,se.re.sl.sim11,
     res.si.sl.sim11,se.si.sl.sim11,
     res.r1.sl.sim11,se.r1.sl.sim11,
     res.fs.sl.sim11,se.fs.sl.sim11,
     res.na.tl.sim11,se.na.tl.sim11,
     res.re.tl.sim11,se.re.tl.sim11,
     res.si.tl.sim11,se.si.tl.sim11,
     res.r1.tl.sim11,se.r1.tl.sim11,
     res.fs.tl.sim11,se.fs.tl.sim11,
     file = "res_scenario11.Rdata")


#### scenario 12: 10 covariates (6 binary and 4 continuous) & sample size = 700 & variation ####

#### s-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.sl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.na <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- glm(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat,
                                      data = train,family = binomial,control = glm.control(maxit = 100000))
                           coef <- mod$coefficients
                           coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                           se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                           
                           #ite estimation
                           ite <- predict(mod, newdata=test1, type="response") -
                             predict(mod, newdata=test0, type="response")
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, vcov(mod))
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.sl.sim12 <- t(res.na.sl12[1:26, ])
se.na.sl.sim12 <- t(res.na.sl12[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.sl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       
                       c.ben <- c()
                       c.ben.se <- c()
                       a <- c()
                       b <- c()
                       mse.re <- c()
                       in_int <- c()
                       coef_mat <- matrix(NA, nrow = 20, ncol = k)
                       se_mat <- matrix(NA, nrow = 20, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying the model to train
                           mod <- glmer(death~(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*factor(treat)+(1|trial)+(0+treat|trial),
                                        data = train,family = binomial,
                                        control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           coef <- mod@beta
                           coef_mat[,i] <- summary(mod)$coef[ ,"Estimate"]
                           se_mat[,i] <- summary(mod)$coef[ ,"Std. Error"]
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                           lp <- X %*% coef
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                          coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                          coef[10]*test$height+coef[11]+
                                          coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                          coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                          coef[20]*test$height) -
                             expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                     coef[10]*test$height)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, vcov(mod))
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.sl.sim12 <- t(res.re.sl12[1:26, ])
se.re.sl.sim12 <- t(res.re.sl12[27:47, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.sl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                       ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                       pout <- plogis(df$logodds)
                       df$death <- rbinom(n, 1, pout)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train <- df[df$trial!=i,]
                           
                           test0 <- test
                           test0$treat <- 0
                           test1 <- test
                           test1$treat <- 1
                           
                           #applying model to train
                           mod <- glm(death ~ -1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat+treat:factor(trial),
                                      data = train,family = "binomial",control = glm.control(maxit = 100000))
                           coef <- mod$coefficients
                           
                           K.si <- array(0, dim=c(20, 30, 6))
                           ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, mi=T,
                                    stroke=T,smoke=T,diab=T,lvhn1=T,height=T,treat=T,
                                    age_treat=T, sex_treat=T, sbp_treat=T, mi_treat=T,stroke_treat=T,
                                    smoke_treat=T,diab_treat=T,lvhn1_treat=T,height_treat=T,
                                    t2_treat=F,t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
                           for(p in 1:6){
                             diag(K.si[2:20,ind,p]) <- 1
                             K.si[1,p,p] <- 1
                             if(p == 1) next
                             K.si[11,24+p,p] <- 1
                           }
                           allcoef <- t(apply(K.si, 3, function(x){x %*% coef}))
                           allvar <- apply(K.si, 3, function(x){x %*% vcov(mod) %*% t(x)}, simplify=F)
                           
                           si_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
                           coef = si_meta$coefficients
                           
                           coef_mat[,i] <- coef
                           se_mat[,i] <- diag(si_meta$vcov)
                           
                           #recalibration of event rates adapted to train
                           X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                           lp <- X %*% c(coef)
                           coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                          coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                          coef[10]*test$height+coef[11]+
                                          coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                          coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                          coef[20]*test$height) -
                             expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                     coef[10]*test$height)
                           
                           #ite prediction interval
                           se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                                               (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, si_meta$vcov)
                           
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
                           se_mat.20 = mean(se_mat[20, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.sl.sim12 <- t(res.si.sl12[1:26, ])
se.si.sl.sim12 <- t(res.si.sl12[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.sl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")
) %dopar% {
  set.seed(j)
  n <- 700
  k <- 7
  nt <- n/k
  
  trial <- rep(1:k, each = nt)
  df <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df$age <- df$age-mean(df$age)
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df$sbp <- df$sbp-mean(df$sbp)
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df$height <- df$height-mean(df$height)
  
  df$treat <- rep(c(0,1), times = n/(2*k))
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- -0.3
  b11 <- 0.015
  b12 <- 0.04
  b13 <- 0.1
  b14 <- -0.008
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  smoke_het <- 0.2 * pman
  smoke_het <- smoke_het - mean(smoke_het)
  
  diab_het <- (mean_age-60)/60
  diab_het <- diab_het - mean(diab_het)
  
  df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
    (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
  df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                    b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
    expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
            (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
  ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
  pout <- plogis(df$logodds)
  df$death <- rbinom(n, 1, pout)
  
  c.ben <- c()
  c.ben.se <- c()
  a <- c()
  b <- c()
  mse.r1 <- c()
  in_int <- c()
  coef_mat <- matrix(NA, nrow = 20, ncol = k)
  se_mat <- matrix(NA, nrow = 20, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df[df$trial==i,]
      train <- df[df$trial!=i,]
      
      test0 <- test
      test0$treat <- 0
      test1 <- test
      test1$treat <- 1
      
      #applying model to train
      mod <- rrvglm(death ~ -1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height)*treat+treat:factor(trial),
                    family = binomialff, data = train,Rank = 1)
      coef <- mod@coefficients
      
      K.r1 <- array(0, dim=c(20, 30, 6))
      ind <- c(t1=F, t2=F, t3=F, t4=F, t5=F, t6=F, age=T, sex=T,sbp=T, mi=T,
               stroke=T,smoke=T,diab=T,lvhn1=T,height=T,treat=T,
               age_treat=T, sex_treat=T, sbp_treat=T, mi_treat=T,stroke_treat=T,
               smoke_treat=T,diab_treat=T,lvhn1_treat=T,height_treat=T,
               t2_treat=F,t3_treat=F, t4_treat=F, t5_treat=F, t6_treat=F)
      for(p in 1:6){
        diag(K.r1[2:20,ind,p]) <- 1
        K.r1[1,p,p] <- 1
        if(p == 1) next
        K.r1[11,24+p,p] <- 1
      }
      allcoef <- t(apply(K.r1, 3, function(x){x %*% coef}))
      allvar <- apply(K.r1, 3, function(x){x %*% vcovvlm(mod) %*% t(x)}, simplify=F)
      
      r1_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
      coef = r1_meta$coefficients
      
      coef_mat[,i] <- coef
      se_mat[,i] <- diag(r1_meta$vcov)
      
      #recalibration of event rates adapted to train
      X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
      lp <- X %*% c(coef)
      coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
      
      #ite estimation
      ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                     coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                     coef[10]*test$height+coef[11]+
                     coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                     coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                     coef[20]*test$height) -
        expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                coef[10]*test$height)
      
      #ite prediction interval
      se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) - 
                          (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, r1_meta$vcov)
      
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
    se_mat.20 = mean(se_mat[20, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.sl.sim12 <- t(res.r1.sl12[1:26, ])
se.r1.sl.sim12 <- t(res.r1.sl12[27:47, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.sl12 = foreach(j = 1:m,
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
                        
                        mean_age <- c(52,56,64,70,77,78,82) #Age
                        sd_age <- c(4,2,1,3,4,6,2)
                        df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                        df$age <- df$age-mean(df$age)
                        
                        pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                        df$sex <- rbinom(n, 1, prob=pman[trial])
                        
                        mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                        sd_sbp <- c(9,11,5,12,9,10,16)
                        df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                        df$sbp <- df$sbp-mean(df$sbp)
                        
                        pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction
                        df$mi <- rbinom(n, 1, prob=pmi[trial])
                        
                        pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                        df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                        
                        psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                        df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                        
                        pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                        df$diab <- rbinom(n, 1, prob=pdiab[trial])
                        
                        plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                        df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                        
                        mean_height <- c(176,162,167,169,168,170,167) #Height
                        sd_height <- c(6,9,10,10,10,9,9)
                        df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                        df$height <- df$height-mean(df$height)
                        
                        df$treat <- rep(c(0,1), times = n/(2*k))
                        
                        b0 <- -1.4
                        b1 <- 0.03
                        b2 <- 0.7
                        b3 <- 0.02
                        b4 <- 0.82
                        b5 <- 0.8
                        b6 <- 0.7
                        b7 <- 0.1
                        b8 <- 0.33
                        b9 <- -0.02
                        b10 <- -0.3
                        b11 <- 0.015
                        b12 <- 0.04
                        b13 <- 0.1
                        b14 <- -0.008
                        
                        trial_eff <- (mean_age-60)/20
                        trial_eff <- trial_eff - mean(trial_eff)
                        
                        trt_het <- (mean_age-60)/100
                        trt_het <- trt_het - mean(trt_het)
                        
                        sbp_het <- 0.2 * pman
                        sbp_het <- sbp_het - mean(sbp_het)
                        
                        smoke_het <- 0.2 * pman
                        smoke_het <- smoke_het - mean(smoke_het)
                        
                        diab_het <- (mean_age-60)/60
                        diab_het <- diab_het - mean(diab_het)
                        
                        df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                          (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                          (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                        df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                          (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                          b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                          expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                  (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
                        ite_true <- split(c(df$ite), ceiling(seq_along(c(df$ite))/nt))
                        pout <- plogis(df$logodds)
                        df$death <- rbinom(n, 1, pout)
                        
                        c.ben <- c()
                        c.ben.se <- c()
                        a <- c()
                        b <- c()
                        mse.fs <- c()
                        in_int <- c()
                        coef_mat <- matrix(NA, nrow = 20, ncol = k)
                        se_mat <- matrix(NA, nrow = 20, ncol = k)
                        
                        testerror = try({
                          for(i in 1:k){
                            test <- df[df$trial==i,]
                            train <- df[df$trial!=i,]
                            
                            test0 <- test
                            test0$treat <- 0
                            test1 <- test
                            test1$treat <- 1
                            
                            #applying model to train
                            mod <- glm(death~-1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+
                                                                 height)*treat*factor(trial),
                                       data = train,family = binomial,control = glm.control(maxit = 100000))
                            coef <- mod$coefficients
                            coef <- coef[!is.na(coef)]
                            
                            K.fs <- array(0, dim=c(20, length(coef), 6))
                            ind <- c(rep(F,6),rep(T,19),rep(F,length(coef)-25))
                            if(i==1){vec <- 3:6} else {vec <- 2:7}
                            for(p in 1:6){
                              diag(K.fs[2:20,ind,p]) <- 1
                              K.fs[1,p,p] <- 1
                              if(p == 1) next
                              K.fs[2,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
                              K.fs[3,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
                              K.fs[4,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
                              K.fs[5,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
                              K.fs[6,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
                              K.fs[7,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
                              K.fs[8,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
                              K.fs[9,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
                              K.fs[10,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
                              K.fs[11,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":treat"))-2,p] <- 1
                              K.fs[12,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"age:treat"))-2,p] <- 1
                              K.fs[13,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(sex)1:treat"))-2,p] <- 1
                              K.fs[14,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":sbp:treat"))-2,p] <- 1
                              K.fs[15,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1:treat"))-2,p] <- 1
                              K.fs[16,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1:treat"))-2,p] <- 1
                              K.fs[17,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(smoke)1:treat"))-2,p] <- 1
                              K.fs[18,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"factor(diab)1:treat"))-2,p] <- 1
                              K.fs[19,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"lvhn1:treat"))-2,p] <- 1
                              K.fs[20,which(names(coef) %in% paste0("factor(trial)",vec[!vec==i],"height:treat"))-2,p] <- 1
                            }
                            allcoef <- t(apply(K.fs, 3, function(x){x %*% coef}))
                            
                            X <- model.matrix(death~-1+factor(trial)+(age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+
                                                                        factor(diab)+lvhn1+height)*treat*factor(trial),data = train)
                            idx <- which(colnames(X) %in% names(coef))
                            
                            X2 <- X[,idx]
                            allvar <- apply(K.fs, 3, function(x){x %*% cov(X2) %*% t(x)}, simplify=F)
                            
                            fs_meta = mvmeta(allcoef~1,S=allvar, method = "reml",control = mvmeta.control(maxit = 100000))
                            coef = fs_meta$coefficients
                            
                            coef_mat[,i] <- coef
                            se_mat[,i] <- diag(fs_meta$vcov)
                            
                            #recalibration of event rates adapted to train
                            X <- model.matrix(as.formula("death~(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*treat"), data=train)
                            lp <- X %*% c(coef)
                            coef[1] <- coef[1]+glm(death ~ offset(lp), data=train, family="binomial")$coef[1]
                            
                            #ite estimation
                            ite <- expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                           coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                           coef[10]*test$height+coef[11]+
                                           coef[12]*test$age+coef[13]*test$sex+coef[14]*test$sbp+coef[15]*test$mi+
                                           coef[16]*test$stroke+coef[17]*test$smoke+coef[18]*test$diab+coef[19]*test$lvhn1+
                                           coef[20]*test$height) -
                              expit(coef[1]+coef[2]*test$age+coef[3]*test$sex+coef[4]*test$sbp+coef[5]*test$mi+
                                      coef[6]*test$stroke+coef[7]*test$smoke+coef[8]*test$diab+coef[9]*test$lvhn1+
                                      coef[10]*test$height)
                            
                            #ite prediction interval
                            se <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20))) -
                                                (exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef, fs_meta$vcov)
                            
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
                          se_mat.20 = mean(se_mat[20, ])),
                          seed = j,
                          success = TRUE)
                      }
stopCluster(cl)
res.fs.sl.sim12 <- t(res.fs.sl12[1:26, ])
se.fs.sl.sim12 <- t(res.fs.sl12[27:47, ])

#### t-learner ####

#### naive model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.na.tl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train0,family = binomial,
                                       control = glm.control(maxit=100000))
                           mod1 <- glm(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
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
                           
                           #ite estimation
                           ite <- predict(mod1, newdata=test, type="response") -
                             predict(mod0, newdata=test, type="response")
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, vcov(mod1))
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.na.tl.sim12 <- t(res.na.tl12[1:26, ])
se.na.tl.sim12 <- t(res.na.tl12[27:47, ])

#### random effects model ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.re.tl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+(1|trial),
                                         data = train0,family = binomial,
                                         control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           mod1 <- glmer(death~age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height+(1|trial),
                                         data = train1,family = binomial,
                                         control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=100000)))
                           coef0 <- mod0@beta
                           coef1 <- mod1@beta
                           coef_mat0[,i] <- summary(mod0)$coef[ ,"Estimate"]
                           coef_mat1[,i] <- summary(mod1)$coef[ ,"Estimate"]
                           se_mat0[,i] <- summary(mod0)$coef[ ,"Std. Error"]
                           se_mat1[,i] <- summary(mod1)$coef[ ,"Std. Error"]
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
                           lp0 <- X0 %*% coef0
                           lp1 <- X1 %*% coef1
                           
                           coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
                           coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
                           
                           #ite estimation
                           ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                                          coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                                          coef1[10]*test$height) -
                             expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                                     coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                                     coef0[10]*test$height)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, vcov(mod0))
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, vcov(mod1))
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.re.tl.sim12 <- t(res.re.tl12[1:26, ])
se.re.tl.sim12 <- t(res.re.tl12[27:47, ])

#### stratified intercept ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.si.tl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","msm","mvmeta")) %dopar% {
                       set.seed(j)
                       n <- 700
                       k <- 7
                       nt <- n/k
                       
                       trial <- rep(1:k, each = nt)
                       df <- as.data.frame(trial)
                       
                       mean_age <- c(52,56,64,70,77,78,82) #Age
                       sd_age <- c(4,2,1,3,4,6,2)
                       df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
                       df$age <- df$age-mean(df$age)
                       
                       pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
                       df$sex <- rbinom(n, 1, prob=pman[trial])
                       
                       mean_sbp <- c(186,182,170,185,190,188,197) #SBP
                       sd_sbp <- c(9,11,5,12,9,10,16)
                       df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
                       df$sbp <- df$sbp-mean(df$sbp)
                       
                       pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
                       df$mi <- rbinom(n, 1, prob=pmi[trial])
                       
                       pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
                       df$stroke <- rbinom(n, 1, prob=pstroke[trial])
                       
                       psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
                       df$smoke <- rbinom(n, 1, prob=psmoke[trial])
                       
                       pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
                       df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
                       
                       plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
                       df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
                       
                       mean_height <- c(176,162,167,169,168,170,167) #Height
                       sd_height <- c(6,9,10,10,10,9,9)
                       df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
                       df$height <- df$height-mean(df$height)
                       
                       df$treat <- rep(c(0,1), times = n/(2*k))
                       
                       b0 <- -1.4
                       b1 <- 0.03
                       b2 <- 0.7
                       b3 <- 0.02
                       b4 <- 0.82
                       b5 <- 0.8
                       b6 <- 0.7
                       b7 <- 0.1
                       b8 <- 0.33
                       b9 <- -0.02
                       b10 <- -0.3
                       b11 <- 0.015
                       b12 <- 0.04
                       b13 <- 0.1
                       b14 <- -0.008
                       
                       trial_eff <- (mean_age-60)/20
                       trial_eff <- trial_eff - mean(trial_eff)
                       
                       trt_het <- (mean_age-60)/100
                       trt_het <- trt_het - mean(trt_het)
                       
                       sbp_het <- 0.2 * pman
                       sbp_het <- sbp_het - mean(sbp_het)
                       
                       smoke_het <- 0.2 * pman
                       smoke_het <- smoke_het - mean(smoke_het)
                       
                       diab_het <- (mean_age-60)/60
                       diab_het <- diab_het - mean(diab_het)
                       
                       df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
                         (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
                       df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                         (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                                         b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
                         expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                                 (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
                       coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat0 <- matrix(NA, nrow = 10, ncol = k)
                       se_mat1 <- matrix(NA, nrow = 10, ncol = k)
                       
                       testerror = try({
                         for(i in 1:k){
                           test <- df[df$trial==i,]
                           train0 <- df0[df0$trial!=i,]
                           train1 <- df1[df1$trial!=i,]
                           
                           #applying the model to train
                           mod0 <- glm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train0,family = binomial,
                                       control = glm.control(maxit=100000))
                           mod1 <- glm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+factor(diab)+lvhn1+height,
                                       data = train1,family = binomial,
                                       control = glm.control(maxit=100000))
                           coef0 <- mod0$coefficients
                           coef1 <- mod1$coefficients
                           coef0[is.na(coef0)] <- 0
                           coef1[is.na(coef1)] <- 0
                           
                           K.si0 <- array(0, dim=c(10, 15, 6))
                           ind <- c(rep(F,6),rep(T,9))
                           for(p in 1:6){
                             diag(K.si0[2:10,ind,p]) <- 1
                             K.si0[1,p,p] <- 1
                           }
                           allcoef0 <- t(apply(K.si0, 3, function(x){x %*% coef0}))
                           allvar0 <- apply(K.si0, 3, function(x){x %*% vcov(mod0) %*% t(x)}, simplify=F)
                           
                           si_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
                           coef0 = si_meta0$coefficients
                           
                           coef_mat0[,i] <- coef0
                           se_mat0[,i] <- diag(si_meta0$vcov)
                           
                           K.si1 <- array(0, dim=c(10, 15, 6))
                           ind <- c(rep(F,6),rep(T,9))
                           for(p in 1:6){
                             diag(K.si1[2:10,ind,p]) <- 1
                             K.si1[1,p,p] <- 1
                           }
                           allcoef1 <- t(apply(K.si1, 3, function(x){x %*% coef1}))
                           allvar1 <- apply(K.si1, 3, function(x){x %*% vcov(mod1) %*% t(x)}, simplify=F)
                           
                           si_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
                           coef1 = si_meta1$coefficients
                           
                           coef_mat1[,i] <- coef1
                           se_mat1[,i] <- diag(si_meta1$vcov)
                           
                           #recalibration of event rates adapted to train
                           X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
                           X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
                           lp0 <- X0 %*% c(coef0)
                           lp1 <- X1 %*% c(coef1)
                           
                           #ite estimation
                           ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                                          coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                                          coef1[10]*test$height) -
                             expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                                     coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                                     coef0[10]*test$height)
                           
                           #ite prediction interval
                           se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, si_meta0$vcov)
                           se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, si_meta1$vcov)
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
                           coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
                           coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
                           coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
                           coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
                           coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
                           coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
                           coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
                           coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
                           coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
                           coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
                           se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
                           se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
                           se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
                           se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
                           se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
                           se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
                           se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
                           se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
                           se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
                           se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
                         seed = j,
                         success = TRUE)
                     }
stopCluster(cl)
res.si.tl.sim12 <- t(res.si.tl12[1:26, ])
se.si.tl.sim12 <- t(res.si.tl12[27:47, ])

#### rank-1 ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.r1.tl12 = foreach(j = 1:m,
                     .final = function(l) {do.call("cbind", lapply(
                       l[sapply(l,`[[`, "success")],
                       function(subl) c(subl$res, subl$seed)))},
                     .packages = c("tidyverse","Hmisc","lme4","magrittr","gtools","VGAM","msm","mvmeta")
) %dopar% {
  set.seed(j)
  n <- 700
  k <- 7
  nt <- n/k
  
  trial <- rep(1:k, each = nt)
  df <- as.data.frame(trial)
  
  mean_age <- c(52,56,64,70,77,78,82) #Age
  sd_age <- c(4,2,1,3,4,6,2)
  df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
  df$age <- df$age-mean(df$age)
  
  pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
  df$sex <- rbinom(n, 1, prob=pman[trial])
  
  mean_sbp <- c(186,182,170,185,190,188,197) #SBP
  sd_sbp <- c(9,11,5,12,9,10,16)
  df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
  df$sbp <- df$sbp-mean(df$sbp)
  
  pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction 
  df$mi <- rbinom(n, 1, prob=pmi[trial])
  
  pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
  df$stroke <- rbinom(n, 1, prob=pstroke[trial])
  
  psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
  df$smoke <- rbinom(n, 1, prob=psmoke[trial])
  
  pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
  df$diab <- rbinom(n, 1, prob=pdiab[trial]) 
  
  plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
  df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
  
  mean_height <- c(176,162,167,169,168,170,167) #Height
  sd_height <- c(6,9,10,10,10,9,9)
  df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
  df$height <- df$height-mean(df$height)
  
  df$treat <- rep(c(0,1), times = n/(2*k))
  
  b0 <- -1.4
  b1 <- 0.03
  b2 <- 0.7
  b3 <- 0.02
  b4 <- 0.82
  b5 <- 0.8
  b6 <- 0.7
  b7 <- 0.1
  b8 <- 0.33
  b9 <- -0.02
  b10 <- -0.3
  b11 <- 0.015
  b12 <- 0.04
  b13 <- 0.1
  b14 <- -0.008
  
  trial_eff <- (mean_age-60)/20
  trial_eff <- trial_eff - mean(trial_eff)
  
  trt_het <- (mean_age-60)/100
  trt_het <- trt_het - mean(trt_het)
  
  sbp_het <- 0.2 * pman
  sbp_het <- sbp_het - mean(sbp_het)
  
  smoke_het <- 0.2 * pman
  smoke_het <- smoke_het - mean(smoke_het)
  
  diab_het <- (mean_age-60)/60
  diab_het <- diab_het - mean(diab_het)
  
  df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
    (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
  df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                    (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                    b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
    expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
            (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
  coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
  coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
  se_mat0 <- matrix(NA, nrow = 10, ncol = k)
  se_mat1 <- matrix(NA, nrow = 10, ncol = k)
  
  testerror = try({
    for(i in 1:k){
      test <- df[df$trial==i,]
      train0 <- df0[df0$trial!=i,]
      train1 <- df1[df1$trial!=i,]
      
      #applying the model to train
      mod0 <- rrvglm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height, family = binomialff, data = train0,Rank = 1)
      mod1 <- rrvglm(death~-1+factor(trial)+age+factor(sex)+sbp+factor(mi)+factor(stroke)+factor(smoke)+diab+lvhn1+height, family = binomialff, data = train1,Rank = 1)
      
      coef0 <- mod0@coefficients
      coef1 <- mod1@coefficients
      coef0[is.na(coef0)] <- 0
      coef1[is.na(coef1)] <- 0
      
      K.r10 <- array(0, dim=c(10, 15, 6))
      ind <- c(rep(F,6),rep(T,9))
      for(p in 1:6){
        diag(K.r10[2:10,ind,p]) <- 1
        K.r10[1,p,p] <- 1
      }
      allcoef0 <- t(apply(K.r10, 3, function(x){x %*% coef0}))
      allvar0 <- apply(K.r10, 3, function(x){x %*% vcovvlm(mod0) %*% t(x)}, simplify=F)
      
      r1_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
      coef0 = r1_meta0$coefficients
      
      coef_mat0[,i] <- coef0
      se_mat0[,i] <- diag(r1_meta0$vcov)
      
      K.r11 <- array(0, dim=c(10, 15, 6))
      ind <- c(rep(F,6),rep(T,9))
      for(p in 1:6){
        diag(K.r11[2:10,ind,p]) <- 1
        K.r11[1,p,p] <- 1
      }
      allcoef1 <- t(apply(K.r11, 3, function(x){x %*% coef1}))
      allvar1 <- apply(K.r11, 3, function(x){x %*% vcovvlm(mod1) %*% t(x)}, simplify=F)
      
      r1_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
      coef1 = r1_meta1$coefficients
      
      coef_mat1[,i] <- coef1
      se_mat1[,i] <- diag(r1_meta1$vcov)
      
      #recalibration of event rates adapted to train
      X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
      X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
      lp0 <- X0 %*% c(coef0)
      lp1 <- X1 %*% c(coef1)
      
      #ite estimation
      ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                     coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                     coef1[10]*test$height) -
        expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                coef0[10]*test$height)
      
      #ite prediction interval
      se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, r1_meta0$vcov)
      se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, r1_meta1$vcov)
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
  })
  if(any(class(testerror) == "try-error")) return(list(success = FALSE))
  list(
    res = c(c.ben = mean(c.ben),
            c.ben.se = mean(c.ben.se),
            a = mean(a),
            b = mean(b),
            mse = mean(mse.r1),
            in_int = mean(in_int),
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
            coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
            coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
            coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
            coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
            coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
            coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
            coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
            coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
            coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
            coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
            se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
            se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
            se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
            se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
            se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
            se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
            se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
            se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
            se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
            se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
    seed = j,
    success = TRUE)
}
stopCluster(cl)
res.r1.tl.sim12 <- t(res.r1.tl12[1:26, ])
se.r1.tl.sim12 <- t(res.r1.tl12[27:47, ])

#### fully stratified ####
m <- 1000
cl = makeCluster(detectCores())
registerDoParallel(cl)
res.fs.tl12 = foreach(
  j = 1:m,
  .final = function(l) do.call("cbind",
                               lapply(
                                 l[sapply(l, `[[`, "success")],
                                 function(subl) c(subl$res, subl$seed))),
  .packages = c("tidyverse", "Hmisc", "lme4", "magrittr", "gtools", "VGAM","msm","mvmeta")) %dopar%
  {
    set.seed(j)
    n <- 700
    k <- 7
    nt <- n/k
    
    trial <- rep(1:k, each = nt)
    df <- as.data.frame(trial)
    
    mean_age <- c(52,56,64,70,77,78,82) #Age
    sd_age <- c(4,2,1,3,4,6,2)
    df$age <- round(rnorm(n,mean=mean_age[trial], sd=sd_age[trial]),0)
    df$age <- df$age-mean(df$age)
    
    pman <- c(0.8,0.4,0.5,0.6,0.5,0.7,0.5) #Sex
    df$sex <- rbinom(n, 1, prob=pman[trial])
    
    mean_sbp <- c(186,182,170,185,190,188,197) #SBP
    sd_sbp <- c(9,11,5,12,9,10,16)
    df$sbp <- round(rnorm(n,mean=mean_sbp[trial], sd=sd_sbp[trial]),0)
    df$sbp <- df$sbp-mean(df$sbp)
    
    pmi <- c(0.1,0.005,0.01,0.02,0.05,0.01,0.04) #Previous myocardial infarction
    df$mi <- rbinom(n, 1, prob=pmi[trial])
    
    pstroke <- c(0.002,0.06,0.02,0.02,0.001,0.008,0.04) #Previous stroke
    df$stroke <- rbinom(n, 1, prob=pstroke[trial])
    
    psmoke <- c(0.5,0.2,0.3,0.4,0.3,0.25,0.3) #Current smoker
    df$smoke <- rbinom(n, 1, prob=psmoke[trial])
    
    pdiab <- c(0.03,0.001,0.002,0.07,0.003,0.01,0.002) #Diabetes
    df$diab <- rbinom(n, 1, prob=pdiab[trial])
    
    plvhn1 <- c(0.13,0.11,0.05,0.25,0.05,0.06,0.04) #LVHN1 (left ventricular hypertrophy detected by electrocardiography)
    df$lvhn1 <- rbinom(n, 1, prob=plvhn1[trial])
    
    mean_height <- c(176,162,167,169,168,170,167) #Height
    sd_height <- c(6,9,10,10,10,9,9)
    df$height <- round(rnorm(n,mean=mean_height[trial], sd=sd_height[trial]),0)
    df$height <- df$height-mean(df$height)
    
    df$treat <- rep(c(0,1), times = n/(2*k))
    
    b0 <- -1.4
    b1 <- 0.03
    b2 <- 0.7
    b3 <- 0.02
    b4 <- 0.82
    b5 <- 0.8
    b6 <- 0.7
    b7 <- 0.1
    b8 <- 0.33
    b9 <- -0.02
    b10 <- -0.3
    b11 <- 0.015
    b12 <- 0.04
    b13 <- 0.1
    b14 <- -0.008
    
    trial_eff <- (mean_age-60)/20
    trial_eff <- trial_eff - mean(trial_eff)
    
    trt_het <- (mean_age-60)/100
    trt_het <- trt_het - mean(trt_het)
    
    sbp_het <- 0.2 * pman
    sbp_het <- sbp_het - mean(sbp_het)
    
    smoke_het <- 0.2 * pman
    smoke_het <- smoke_het - mean(smoke_het)
    
    diab_het <- (mean_age-60)/60
    diab_het <- diab_het - mean(diab_het)
    
    df$logodds <- b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
      (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10*df$treat+trt_het[df$trial]*df$treat+
      (b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke)*df$treat
    df$ite <- expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
                      (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height+b10+trt_het[df$trial]+
                      b11*df$age+b12*df$sex+b13*df$sbp+b14*df$smoke) -
      expit(b0+trial_eff[df$trial]+b1*df$age+b2*df$sex+(sbp_het[df$trial]+b3)*df$sbp+b4*df$mi+b5*df$stroke+(smoke_het[df$trial]+b6)*df$smoke+
              (diab_het[df$trial]+b7)*df$diab+b8*df$lvhn1+b9*df$height)
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
    coef_mat0 <- matrix(NA, nrow = 10, ncol = k)
    coef_mat1 <- matrix(NA, nrow = 10, ncol = k)
    se_mat0 <- matrix(NA, nrow = 10, ncol = k)
    se_mat1 <- matrix(NA, nrow = 10, ncol = k)
    
    testerror = try({
      for(i in 1:k){
        test <- df[df$trial==i,]
        train0 <- df0[df0$trial!=i,]
        train1 <- df1[df1$trial!=i,]
        
        #applying the model to train
        mod0 <- glm(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),
                    data = train0,family = binomial,control = glm.control(maxit = 100000))
        mod1 <- glm(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),
                    data = train1,family = binomial,control = glm.control(maxit = 100000))
        coef0 <- mod0$coefficients
        coef1 <- mod1$coefficients
        coef0 <- coef0[!is.na(coef0)]
        coef1 <- coef1[!is.na(coef1)]
        
        K.fs0 <- array(0, dim=c(10, length(coef0), 6))
        ind <- c(rep(F,6),rep(T,9),rep(F,length(coef0)-15))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs0[2:10,ind,p]) <- 1
          K.fs0[1,p,p] <- 1
          if(p == 1) next
          K.fs0[2,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs0[3,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs0[4,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs0[5,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs0[6,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs0[7,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs0[8,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs0[9,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs0[10,which(names(coef0) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
        }
        allcoef0 <- t(apply(K.fs0, 3, function(x){x %*% coef0}))
        
        X0 <- model.matrix(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),data = train0)
        idx0 <- which(colnames(X0) %in% names(coef0))
        X0.2 <- X0[,idx0]
        allvar0 <- apply(K.fs0, 3, function(x){x %*% cov(X0.2) %*% t(x)}, simplify=F)
        
        fs_meta0 = mvmeta(allcoef0~1,S=allvar0, method = "reml")
        coef0 = fs_meta0$coefficients
        
        coef_mat0[,i] <- coef0
        se_mat0[,i] <- diag(fs_meta0$vcov)
        
        K.fs1 <- array(0, dim=c(10, length(coef1), 6))
        ind <- c(rep(F,6),rep(T,9),rep(F,length(coef1)-15))
        if(i==1){vec <- 3:6} else {vec <- 2:7}
        for(p in 1:6){
          diag(K.fs1[2:10,ind,p]) <- 1
          K.fs1[1,p,p] <- 1
          if(p == 1) next
          K.fs1[2,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":age"))-2,p] <- 1
          K.fs1[3,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(sex)1"))-2,p] <- 1
          K.fs1[4,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":sbp"))-2,p] <- 1
          K.fs1[5,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(mi)1"))-2,p] <- 1
          K.fs1[6,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(stroke)1"))-2,p] <- 1
          K.fs1[7,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(smoke)1"))-2,p] <- 1
          K.fs1[8,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":factor(diab)1"))-2,p] <- 1
          K.fs1[9,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":lvhn1"))-2,p] <- 1
          K.fs1[10,which(names(coef1) %in% paste0("factor(trial)",vec[!vec==i],":height"))-2,p] <- 1
        }
        allcoef1 <- t(apply(K.fs1, 3, function(x){x %*% coef1}))
        
        X1 <- model.matrix(death ~ -1+factor(trial)+(age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height)*factor(trial),data = train1)
        idx1 <- which(colnames(X1) %in% names(coef1))
        X1.2 <- X1[,idx1]
        allvar1 <- apply(K.fs1, 3, function(x){x %*% cov(X1.2) %*% t(x)}, simplify=F)
        
        fs_meta1 = mvmeta(allcoef1~1,S=allvar1, method = "reml")
        coef1 = fs_meta1$coefficients
        
        coef_mat1[,i] <- coef1
        se_mat1[,i] <- diag(fs_meta1$vcov)
        
        #recalibration of event rates adapted to train
        X0 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train0)
        X1 <- model.matrix(as.formula("death~age+sex+sbp+mi+stroke+smoke+diab+lvhn1+height"), data=train1)
        lp0 <- X0 %*% c(coef0)
        lp1 <- X1 %*% c(coef1)
        
        coef0[1] <- coef0[1]+glm(death ~ offset(lp0), data=train0, family="binomial")$coef[1]
        coef1[1] <- coef1[1]+glm(death ~ offset(lp1), data=train1, family="binomial")$coef[1]
        
        #ite estimation
        ite <- expit(coef1[1]+coef1[2]*test$age+coef1[3]*test$sex+coef1[4]*test$sbp+coef1[5]*test$mi+
                       coef1[6]*test$stroke+coef1[7]*test$smoke+coef1[8]*test$diab+coef1[9]*test$lvhn1+
                       coef1[10]*test$height) -
          expit(coef0[1]+coef0[2]*test$age+coef0[3]*test$sex+coef0[4]*test$sbp+coef0[5]*test$mi+
                  coef0[6]*test$stroke+coef0[7]*test$smoke+coef0[8]*test$diab+coef0[9]*test$lvhn1+
                  coef0[10]*test$height)
        
        #ite prediction interval
        se0 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef0, fs_meta0$vcov)
        se1 <- deltamethod(~(exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)/(1+exp(x1+x2+x3+x4+x5+x6+x7+x8+x9+x10))), coef1, fs_meta1$vcov)
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
    if (any(class(testerror) == "try-error")) return(list(success = FALSE))
    
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
        coef_mat.4 = mean(coef_mat0[4, ]),
        coef_mat.5 = mean(coef_mat0[5, ]),
        coef_mat.6 = mean(coef_mat0[6, ]),
        coef_mat.7 = mean(coef_mat0[7, ]),
        coef_mat.8 = mean(coef_mat0[8, ]),
        coef_mat.9 = mean(coef_mat0[9, ]),
        coef_mat.10 = mean(coef_mat0[10, ]),
        coef_mat.11 = mean(coef_mat1[1, ]) - mean(coef_mat0[1, ]),
        coef_mat.12 = mean(coef_mat1[2, ]) - mean(coef_mat0[2, ]),
        coef_mat.13 = mean(coef_mat1[3, ]) - mean(coef_mat0[3, ]),
        coef_mat.14 = mean(coef_mat1[4, ]) - mean(coef_mat0[4, ]),
        coef_mat.15 = mean(coef_mat1[5, ]) - mean(coef_mat0[5, ]),
        coef_mat.16 = mean(coef_mat1[6, ]) - mean(coef_mat0[6, ]),
        coef_mat.17 = mean(coef_mat1[7, ]) - mean(coef_mat0[7, ]),
        coef_mat.18 = mean(coef_mat1[8, ]) - mean(coef_mat0[8, ]),
        coef_mat.19 = mean(coef_mat1[9, ]) - mean(coef_mat0[9, ]),
        coef_mat.20 = mean(coef_mat1[10, ]) - mean(coef_mat0[10, ]),
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
        se_mat.11 = mean(se_mat1[1, ]) - mean(se_mat0[1, ]),
        se_mat.12 = mean(se_mat1[2, ]) - mean(se_mat0[2, ]),
        se_mat.13 = mean(se_mat1[3, ]) - mean(se_mat0[3, ]),
        se_mat.14 = mean(se_mat1[4, ]) - mean(se_mat0[4, ]),
        se_mat.15 = mean(se_mat1[5, ]) - mean(se_mat0[5, ]),
        se_mat.16 = mean(se_mat1[6, ]) - mean(se_mat0[6, ]),
        se_mat.17 = mean(se_mat1[7, ]) - mean(se_mat0[7, ]),
        se_mat.18 = mean(se_mat1[8, ]) - mean(se_mat0[8, ]),
        se_mat.19 = mean(se_mat1[9, ]) - mean(se_mat0[9, ]),
        se_mat.20 = mean(se_mat1[10, ]) - mean(se_mat0[10, ])),
      seed = j,
      success = TRUE
    )
  }
stopCluster(cl)
res.fs.tl.sim12 <- t(res.fs.tl12[1:26, ])
se.fs.tl.sim12 <- t(res.fs.tl12[27:47, ])

save(res.na.sl.sim12,se.na.sl.sim12,
     res.re.sl.sim12,se.re.sl.sim12,
     res.si.sl.sim12,se.si.sl.sim12,
     res.r1.sl.sim12,se.r1.sl.sim12,
     res.fs.sl.sim12,se.fs.sl.sim12,
     res.na.tl.sim12,se.na.tl.sim12,
     res.re.tl.sim12,se.re.tl.sim12,
     res.si.tl.sim12,se.si.tl.sim12,
     res.r1.tl.sim12,se.r1.tl.sim12,
     res.fs.tl.sim12,se.fs.tl.sim12,
     file = "res_scenario12.Rdata")