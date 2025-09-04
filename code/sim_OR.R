########## Simulation study for log odds ratios (OR) ############


# function 
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

# Example
out_or <- simulate_or(seed=1)
head(out_or)

## escalc cross-check
library(metafor)
esc <- escalc(measure="OR", ai=out_or$a, bi=out_or$b, ci=out_or$c, di=out_or$d)
head(esc)

#save dataset
save(out_or, file=here("data", "dat_OR.rdata"))
