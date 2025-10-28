########## Function to simulate Standardized Mean Difference (SMD) ############
# simulates Hedge's g (Hedges, 1981) and 
# Cohen's d sampling variances (this is escalc in metafor does)

simulate_smd <- function(k = 20, mu = 0.3, tau2 = 0.16,
                         n_T = 50, n_C = 50, mu_C = 100,
                         sigma = 10, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  theta <- rnorm(k, mean = mu, sd = sqrt(tau2))
  
  # set up vectors
  yi <- vi <- numeric(k)
  mTi <- mCi <- numeric(k)
  sdTi <- sdCi <- numeric(k)
  mu_Ti <- numeric(k) #for true treated means
  
  for (i in seq_len(k)) {
    muT_true <- mu_C + theta[i] * sigma   #true mean of treated group in study i 
    xT <- rnorm(n_T, mean = muT_true, sd = sigma)
    xC <- rnorm(n_C, mean = mu_C,  sd = sigma)
    
    mT <- mean(xT); mC <- mean(xC)
    s2T <- var(xT); s2C <- var(xC)
    
    # pooled SD
    sp <- sqrt(((n_T - 1) * s2T + (n_C - 1) * s2C) / (n_T + n_C - 2))
    
    # get Cohen's d
    d <- (mT - mC) / sp
    
    # compute Hedges' g (corrected SMD) and its sampling variance
    J  <- 1 - 3 / (4 * (n_T + n_C - 2) - 1)
    g  <- J * d
    vd <- ( (n_T + n_C) / (n_T * n_C) + (g^2) / (2 * (n_T + n_C - 2)) )##compute Cohen's d not Hedge's d variance
    
    yi[i] <- g
    vi[i] <- vd
    mTi[i] <- mT; mCi[i] <- mC
    sdTi[i] <- sqrt(s2T); sdCi[i] <- sqrt(s2C)
    mu_Ti[i] <- muT_true
  }
  
  data.frame(
    measure = "SMD",
    study= seq_len(k),
    theta= theta,
    yi= yi,   #Hedge's g
    vi= vi,   #var(d)
    n_T= n_T, n_C=n_C,
    mu_T= mu_Ti, mu_C= mu_C,  #true means   
    mT= mTi, mC= mCi, #sample means 
    sdT= sdTi, sdC=sdCi,   #sample SDs 
    sigma= sigma,
    stringsAsFactors = FALSE
  )
}

# # example
# out_smd <- simulate_smd(seed = 1)
# head(out_smd)
# 
# ## escalc cross-check (Cohen's d, uncorrected)
# library(metafor)
# esc <- escalc(
#   measure = "SMD",
#   m1i = out_smd$mT, sd1i = out_smd$sdT, n1i = out_smd$n_T,
#   m2i = out_smd$mC, sd2i = out_smd$sdC, n2i = out_smd$n_C
# )
# 
# head(esc)
# 
# # check
# c(max_abs_diff_yi = max(abs(out_smd$yi - esc$yi)),
#   max_abs_diff_vi = max(abs(out_smd$vi - esc$vi)))
# 
# head(esc)
# 
# #save dataset
# save(out_smd, file=here("data", "dat_smd.rdata"))



########## Function to simulate log response ratios (lnRR) ############


# function
simulate_lnrr <- function(
    k = 20, mu = log(1.3), tau2 = 0.16,
    n_T = 50, n_C = 50,
    mu_C = 100, sigma = 10,
    seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  theta <- rnorm(k, mu, sqrt(tau2))
  
  yi <- vi <- numeric(k)
  mTi <- mCi <- sdTi <- sdCi <- numeric(k)
  mu_Ti <- numeric(k)
  
  for (i in seq_len(k)) {
    muT_true <- mu_C * exp(theta[i]) #true means of treated group for study i
    
    # resample if a sample mean is noot positive (ln undefined)
    repeat {
      xT <- rnorm(n_T, mean = mu_Ti, sd = sigma)
      xC <- rnorm(n_C, mean = mu_C,  sd = sigma)
      mT <- mean(xT); mC <- mean(xC)
      if (mT > 0 && mC > 0) break
    }
    
    sdT <- sd(xT); sdC <- sd(xC)
    
    lnRR    <- log(mT / mC)
    v_lnRR  <- (sdT^2) / (n_T * mT^2) + (sdC^2) / (n_C * mC^2)
    
    yi[i] <- lnRR; vi[i] <- v_lnRR
    mTi[i] <- mT; mCi[i] <- mC
    sdTi[i] <- sdT; sdCi[i] <- sdC
    mu_Ti[i] <- muT_true
  }
  
  data.frame(
    measure = "lnRR",
    study = seq_len(k),
    theta = theta,
    yi = yi, vi = vi,
    n_T = n_T, n_C = n_C,
    mT = mTi, mC = mCi, #sample means
    sdT = sdTi, sdC = sdCi,  #sample SDs
    mu_C = mu_C, mu_T=mu_Ti, #true means
    sigma = sigma,
    stringsAsFactors = FALSE
  )
}

# # Example
# out_lnrr <- simulate_lnrr(seed=1)
# head(out_lnrr)
# 
# ## escalc cross-check
# library(metafor)
# esc <- escalc(measure = "ROM",
#               m1i = out_lnrr$mT, sd1i = out_lnrr$sdT, n1i = out_lnrr$n_T,
#               m2i = out_lnrr$mC, sd2i = out_lnrr$sdC, n2i = out_lnrr$n_C)
# head(esc)
# 
# # check
# c(max_abs_diff_yi = max(abs(out_lnrr$yi - esc$yi)),
#   max_abs_diff_vi = max(abs(out_lnrr$vi - esc$vi)))
# 
# 
# #save dataset
# save(out_lnrr, file=here("data", "dat_lnrr.rdata"))


########## Function to simulate log odds ratios (OR) ############


simulate_or <- function(k = 20, mu = log(1.5), tau2 = 0.16,
                        n_T = 50, n_C = 50, p_C = 0.20, cc = 0.5, 
                        seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  logit <- function(p) log(p/(1-p))
  invlogit <- function(x) 1/(1+exp(-x))
  
  theta <- rnorm(k, mu, sqrt(tau2))
  alpha <- logit(p_C)
  pT <- invlogit(alpha + theta)
  
  y <- v <- numeric(k)
  a <- b <- c <- d <- numeric(k)
  
  for (i in seq_len(k)) {
    a[i] <- rbinom(1, n_T, pT[i])
    b[i] <- n_T - a[i]
    c[i] <- rbinom(1, n_C, p_C)
    d[i] <- n_C - c[i]
    
    if (any(c(a[i],b[i],c[i],d[i]) == 0)) {
      a[i] <- a[i]+cc; b[i] <- b[i]+cc; c[i] <- c[i]+cc; d[i] <- d[i]+cc
    }
    
    logOR <- log((a[i]*d[i])/(b[i]*c[i]))
    v_logOR <- 1/a[i] + 1/b[i] + 1/c[i] + 1/d[i]
    
    y[i] <- logOR; v[i] <- v_logOR
  }
  
  out <- data.frame(measure="OR", study=seq_len(k),
                    theta=theta,
                    yi=y, vi=v,
                    n_T=n_T, n_C=n_C,
                    p_C=p_C,
                    a=a, b=b, 
                    c=c, d=d,
                    stringsAsFactors=FALSE)
  return(out)
}

# # Example
# out_or <- simulate_or(seed=1)
# head(out_or)
# 
# ## escalc cross-check
# library(metafor)
# esc <- escalc(measure="OR", ai=out_or$a, bi=out_or$b, ci=out_or$c, di=out_or$d)
# head(esc)
# 
# #save dataset
# save(out_or, file=here("data", "dat_OR.rdata"))




########## Function to simulate log incidence rate ratios (IRR) ############


simulate_irr <- function(k = 20, mu = log(1.5), tau2 = 0.16,
                         t_T = 50, t_C = 50, lambda_C = 0.10,
                         cc = 0.5, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  theta     <- rnorm(k, mean = mu, sd = sqrt(tau2))
  lambda_T  <- lambda_C * exp(theta)
  
  # set up vectors
  yi <- vi <- numeric(k)
  yT_raw <- yC_raw <- numeric(k)
  
  for (i in seq_len(k)) {
    yT_raw[i] <- rpois(1, lambda_T[i] * t_T)
    yC_raw[i] <- rpois(1, lambda_C * t_C)
    
    # if count is zero add a continuity-correction (=0.5; same as metafor default) ##### NOTE THIS IN SIM STUDY
    yTcc <- if (yT_raw[i] == 0) yT_raw[i] + cc else yT_raw[i]
    yCcc <- if (yC_raw[i] == 0) yC_raw[i] + cc else yC_raw[i]
    
    yi[i] <- log( (yTcc / t_T) / (yCcc / t_C) )
    vi[i] <- 1 / yTcc + 1 / yCcc
  }
  
  data.frame(
    measure = "IRR",
    study   = seq_len(k),
    theta   = theta,
    yi = yi, vi = vi,
    t_T = rep(t_T, k), t_C = rep(t_C, k),
    lambda_C = lambda_C, lambda_T = lambda_T,
    yT = yT_raw, yC = yC_raw,
    cc = cc,
    stringsAsFactors = FALSE
  )
}

# # Example
# out_irr <- simulate_irr(seed = 1)
# head(out_irr)
# 
# ## escalc cross-check 
# library(metafor)
# esc <- escalc(measure = "IRR",
#               x1i= out_irr$yT, t1i=out_irr$t_T,
#               x2i= out_irr$yC, t2i=out_irr$t_C,
#               add= out_irr$cc[1], to= "only0")
# head(esc)
# 
# 
# # check
# c(max_abs_diff_yi = max(abs(out_irr$yi - esc$yi)),
#   max_abs_diff_vi = max(abs(out_irr$vi - esc$vi)))
# 
# 
# #save dataset
# save(out_irr, file=here("data", "dat_IRR.rdata"))
