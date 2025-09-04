########## Simulation study for Standardized Mean Difference (SMD) ############
# simulates Hedge's g (Hedges, 1981) and corresponding sampling variances


# function
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

# example
out_smd <- simulate_smd(seed = 1)
head(out_smd)

## escalc cross-check (Cohen's d, uncorrected)
library(metafor)
esc <- escalc(
  measure = "SMD",
  m1i = out_smd$mT, sd1i = out_smd$sdT, n1i = out_smd$n_T,
  m2i = out_smd$mC, sd2i = out_smd$sdC, n2i = out_smd$n_C
)

head(esc)

# check
c(max_abs_diff_yi = max(abs(out_smd$yi - esc$yi)),
  max_abs_diff_vi = max(abs(out_smd$vi - esc$vi)))

head(esc)

#save dataset
save(out_smd, file=here("data", "dat_smd.rdata"))

