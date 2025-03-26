#### Code from J.E.Pustejovsky 
#https://jepusto.com/posts/Bug-in-nlme-with-fixed-sigma/

library(metafor)
library(nlme)
library(glmmTMB)


##############################################################################################
# Basic random effects model
##############################################################################################

bcg_example <- function(method = "REML", constant_var = FALSE) {
  
  data(dat.bcg)
  dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  dat$g <- 1
  dat$trial <- as.factor(dat$trial)
  
  v_bar <- mean(dat$vi)
  if (constant_var) dat$vi <- v_bar
  
  # random-effects model using rma.uni()
  LOR_uni_fit <- rma(yi, vi, data=dat, method = method)
  LOR_uni <- with(LOR_uni_fit, 
                  data.frame(f = "rma.uni", 
                             logLik = logLik(LOR_uni_fit),
                             df = attr(logLik(LOR_uni_fit), "df"),
                             est = as.numeric(b), 
                             se = se, 
                             tau = sqrt(tau2)))
  
  # random-effects model using rma.mv()
  LOR_mv_fit <- rma.mv(yi, vi, random = ~ 1 | trial, data=dat, method = method)
  LOR_mv <- with(LOR_mv_fit, 
                 data.frame(f = "rma.mv", 
                            logLik = logLik(LOR_mv_fit),
                            df = attr(logLik(LOR_mv_fit), "df"),#model-based degrees of freedom, i.e. the number of (estimated) parameters
                            est = as.numeric(b), 
                            se = se, 
                            tau = sqrt(sigma2)))
  
  # random-effects model using lme()
  if (constant_var) {
    LOR_lme_fit <- lme(yi ~ 1, data = dat, method = method, 
                       random = ~ 1 | trial,
                       control = lmeControl(sigma = sqrt(v_bar)))
    tau <- sqrt(as.numeric(coef(LOR_lme_fit$modelStruct$reStruct, unconstrained = FALSE)) * v_bar) 
  } else {
    LOR_lme_fit <- lme(yi ~ 1, data = dat, method = method, 
                       random = ~ 1 | trial,
                       weights = varFixed(~ vi),
                       control = lmeControl(sigma = 1))
    tau <- sqrt(as.numeric(coef(LOR_lme_fit$modelStruct$reStruct, unconstrained = FALSE)))
  }
  LOR_lme <- data.frame(f = "lme", 
                        logLik = logLik(LOR_lme_fit),
                        df = attr(logLik(LOR_lme_fit), "df"),  
                        est = as.numeric(fixef(LOR_lme_fit)), 
                        se = as.numeric(sqrt(vcov(LOR_lme_fit))), 
                        tau = tau)
  
  
  # random-effects model using glmmTMB() and equalto
  if (method=="REML") REML <- TRUE else REML <- FALSE
  
  if (constant_var) {
    VCV <- diag(v_bar, nrow = nrow(dat))
    LOR_glmm_fit <- glmmTMB(yi ~ 1 + equalto(0 + trial|g,VCV),
                            REML = REML,
                            data=dat)
    tau <- sigma(LOR_glmm_fit) 
  } else {
    VCV <- diag(dat$vi)
    LOR_glmm_fit <- glmmTMB(yi ~ 1 + equalto(0 + trial|g,VCV),
                            REML = REML,
                            data=dat)
    tau <- sigma(LOR_glmm_fit) 
  }
  LOR_glmm <- data.frame(f = "glmmTMB", 
                        logLik = logLik(LOR_glmm_fit)[1],
                        df = attr(logLik(LOR_glmm_fit), "df"),                         
                        est = as.numeric(unlist(fixef(LOR_glmm_fit))[[1]]), 
                        se = as.numeric(sqrt(vcov(LOR_glmm_fit)[[1]])), 
                        tau = tau)
  
  
  rbind(LOR_uni, LOR_mv, LOR_glmm, LOR_lme)
  
}

bcg_example("REML", constant_var = FALSE)
bcg_example("REML", constant_var = TRUE)

bcg_example("ML", constant_var = FALSE)
bcg_example("ML", constant_var = TRUE)




##############################################################################################
# Bi-variate random effects model
##############################################################################################

bcg_bivariate <- function(method = "REML", constant_var = FALSE) {
  
  data(dat.bcg)
  dat_long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
  levels(dat_long$group) <- c("exp", "con")
  dat_long$group <- relevel(dat_long$group, ref="con")
  dat_long$g <- 1
  dat_long$ID <- as.factor(1:nrow(dat_long))
  dat_long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat_long)
  
  v_bar <- mean(dat_long$vi)
  if (constant_var) dat_long$vi <- v_bar
  
  
  # bivariate random-effects model using rma.mv()
  bv_rma_fit <- rma.mv(yi, vi, mods = ~ group, 
                       random = ~ group | study, 
                       struct = "UN", method = method,
                       data=dat_long)
  bv_rma <- with(bv_rma_fit, data.frame(f = "rma.mv",
                                        logLik = logLik(bv_rma_fit),
                                        df = attr(logLik(bv_rma_fit), "df"), 
                                        est1 = as.numeric(b)[1], 
                                        se1 = se[1], 
                                        est2 = as.numeric(b)[2],
                                        se2 = se[2],
                                        tau1 = sqrt(tau2[1]),
                                        tau2 = sqrt(tau2[2])))
  
  # bivariate random-effects model using lme()
  if (constant_var) {
    bv_lme_fit <- lme(yi ~ group, data = dat_long, method = method, 
                      random = ~ group | study,
                      control = lmeControl(sigma = sqrt(v_bar)))
    tau_sq <- colSums(coef(bv_lme_fit$modelStruct$reStruct, unconstrained = FALSE) * matrix(c(1,0,0, 1,2,1), 3, 2)) * v_bar
    
  } else {
    bv_lme_fit <- lme(yi ~ group, data = dat_long, method = method, 
                      random = ~ group | study,
                      weights = varFixed(~ vi),
                      control = lmeControl(sigma = 1))
    tau_sq <- colSums(coef(bv_lme_fit$modelStruct$reStruct, unconstrained = FALSE) * matrix(c(1,0,0, 1,2,1), 3, 2))
  }
  bv_lme <- data.frame(f = "lme",
                       logLik = logLik(bv_lme_fit),
                       df = attr(logLik(bv_lme_fit), "df"), 
                       est1 = as.numeric(fixef(bv_lme_fit))[1],
                       se1 = as.numeric(sqrt(diag(vcov(bv_lme_fit))))[1],
                       est2 = as.numeric(fixef(bv_lme_fit))[2],
                       se2 = as.numeric(sqrt(diag(vcov(bv_lme_fit))))[2],
                       tau1 = sqrt(tau_sq[1]),
                       tau2 = sqrt(tau_sq[2]))
  
  
  
  # bivariate random-effects model using glmmTMB() and equalto
  if (method=="REML") REML <- TRUE else REML <- FALSE
  
  if (constant_var) {
    VCV <- diag(v_bar, nrow = nrow(dat_long))
    bv_glmm_fit <- glmmTMB(yi ~ group + (group|study) + equalto(0 + ID|g,VCV), 
                            REML = REML,
                            data=dat_long)
    tau <- c(exp(bv_glmm_fit$fit$par)[2:3]) 
    ## 1st theta is random intercept study(intercept), 2nd is random slope study(groupexp), and 3rd is the covariance term 

  } else {
    VCV <- diag(dat_long$vi)
    bv_glmm_fit <- glmmTMB(yi ~ group + (group|study) + equalto(0 + ID|g,VCV),
                            REML = REML,
                            data=dat_long)
    tau <- c(exp(bv_glmm_fit$fit$par)[2:3]) 
  }
  bv_glmm <- data.frame(f = "glmmTMB", 
                        logLik = logLik(bv_glmm_fit)[1],
                        df = attr(logLik(bv_glmm_fit), "df"), 
                        est1 = as.numeric(unlist(fixef(bv_glmm_fit))[[1]]), 
                        se1 = as.numeric(sqrt(diag(vcov(bv_glmm_fit)[[1]])))[1], 
                        est2 = as.numeric(unlist(fixef(bv_glmm_fit))[[2]]), 
                        se2 = as.numeric(sqrt(diag(vcov(bv_glmm_fit)[[1]])))[2], 
                        tau1 = tau[1],
                        tau2 = tau[2])
  
  
  
  rbind(bv_rma, bv_glmm, bv_lme)
  
}

bcg_bivariate("REML", constant_var = FALSE)
bcg_bivariate("REML", constant_var = TRUE)

bcg_bivariate("ML", constant_var = FALSE) ##non-convergence for glmmTMB
bcg_bivariate("ML", constant_var = TRUE) ##very different tau1 and ta2 for glmmTMB




##############################################################################################
# Three-level random-effects model
##############################################################################################


Konstantopoulos <- function(method = "REML", constant_var = FALSE) {
  
  dat <- get(data(dat.konstantopoulos2011))
  v_bar <- mean(dat$vi)
  dat$g <- 1
  dat$obs <- as.factor(1:nrow(dat))
  if (constant_var) dat$vi <- v_bar
  
  # multilevel random-effects model using rma.mv()
  ml_rma_fit <- rma.mv(yi, vi, random = ~ 1 | district/school, data=dat, method = method)
  
  ml_rma <- with(ml_rma_fit, 
                 data.frame(f = "rma.mv", 
                            logLik = logLik(ml_rma_fit),
                            df = attr(logLik(ml_rma_fit), "df"), 
                            est = as.numeric(b), 
                            se = se, 
                            tau1 = sqrt(sigma2[1]), 
                            tau2 = sqrt(sigma2[2])))
  
  # multilevel random-effects model using lme()
  if (constant_var) {
    ml_lme_fit <- lme(yi ~ 1, data = dat, method = method, 
                      random = ~ 1 | district / school,
                      control = lmeControl(sigma = sqrt(v_bar)))
    tau <- sqrt(as.numeric(coef(ml_lme_fit$modelStruct$reStruct, unconstrained = FALSE)) * v_bar)
    
  } else {
    ml_lme_fit <- lme(yi ~ 1, data = dat, method = method, 
                      random = ~ 1 | district / school,
                      weights = varFixed(~ vi),
                      control = lmeControl(sigma = 1))
    tau <- sqrt(as.numeric(coef(ml_lme_fit$modelStruct$reStruct, unconstrained = FALSE)))
    
  }  
  ml_lme <- data.frame(f = "lme",
                       logLik = logLik(ml_lme_fit),
                       df = attr(logLik(ml_lme_fit), "df"), 
                       est = as.numeric(fixef(ml_lme_fit)),
                       se = as.numeric(sqrt(diag(vcov(ml_lme_fit)))),
                       tau1 = tau[2],
                       tau2 = tau[1])
  
  
  # multilevel random-effects model using glmmTMB()
  if (method=="REML") REML <- TRUE else REML <- FALSE
  
  if (constant_var) {
    VCV <- diag(v_bar, nrow = nrow(dat))
    ml_glmm_fit <- glmmTMB(yi ~ 1 + equalto(0 + obs|g,VCV) + (1|district),
                            REML = REML,
                            data=dat)
    tau <- c(exp(ml_glmm_fit$fit$par)["theta"], sigma(ml_glmm_fit))
  } else {
    VCV <- diag(dat$vi)
    ml_glmm_fit <- glmmTMB(yi ~ 1 + equalto(0 + obs|g,VCV) + (1|district),
                           REML = REML,
                           data=dat)
    tau <- c(exp(ml_glmm_fit$fit$par)["theta"], sigma(ml_glmm_fit))
  }
  ml_glmm <- data.frame(f = "glmmTMB", 
                         logLik = logLik(ml_glmm_fit)[1],
                         df = attr(logLik(ml_glmm_fit), "df"), 
                         est = as.numeric(unlist(fixef(ml_glmm_fit))[[1]]), 
                         se = as.numeric(sqrt(vcov(ml_glmm_fit)[[1]])), 
                         tau1 = tau[1],
                         tau2 = tau[2])
  
  
  
  rbind(ml_rma, ml_glmm, ml_lme)
  
}

Konstantopoulos("REML", constant_var = FALSE)
Konstantopoulos("REML", constant_var = TRUE)

Konstantopoulos("ML", constant_var = FALSE)
Konstantopoulos("ML", constant_var = TRUE)
