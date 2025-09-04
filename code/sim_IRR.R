########## Simulation study for log incidence rate ratios (IRR) ############



# function
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

# Example
out_irr <- simulate_irr(seed = 1)
head(out_irr)

## escalc cross-check 
library(metafor)
esc <- escalc(measure = "IRR",
              x1i= out_irr$yT, t1i=out_irr$t_T,
              x2i= out_irr$yC, t2i=out_irr$t_C,
              add= out_irr$cc[1], to= "only0")
head(esc)


# check
c(max_abs_diff_yi = max(abs(out_irr$yi - esc$yi)),
  max_abs_diff_vi = max(abs(out_irr$vi - esc$vi)))


#save dataset
save(out_irr, file=here("data", "dat_IRR.rdata"))
