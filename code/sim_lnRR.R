########## Simulation study for log response ratios (lnRR) ############



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

# Example
out_lnrr <- simulate_lnrr(seed=1)
head(out_lnrr)

## escalc cross-check
library(metafor)
esc <- escalc(measure = "ROM",
              m1i = out_lnrr$mT, sd1i = out_lnrr$sdT, n1i = out_lnrr$n_T,
              m2i = out_lnrr$mC, sd2i = out_lnrr$sdC, n2i = out_lnrr$n_C)
head(esc)

# check
c(max_abs_diff_yi = max(abs(out_lnrr$yi - esc$yi)),
  max_abs_diff_vi = max(abs(out_lnrr$vi - esc$vi)))


#save dataset
save(out_lnrr, file=here("data", "dat_lnrr.rdata"))
