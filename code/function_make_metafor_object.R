
convert_equalto_to_metafor <- function(fit.equalto) {
  fit.metafor <- list()
  
  # Extract coefficients
  fit.metafor$b <- matrix(fit.equalto$fit$parfull["beta"], nrow = 1, dimnames = list("intrcpt", NULL))
  fit.metafor$beta <- fit.metafor$b
  fit.metafor$se <- sqrt(fit.equalto$sdr$cov.fixed["betadisp", "betadisp"])
  fit.metafor$zval <- fit.metafor$b / fit.metafor$se
  fit.metafor$pval <- 2 * (1 - pnorm(abs(fit.metafor$zval)))
  fit.metafor$ci.lb <- fit.metafor$b - 1.96 * fit.metafor$se
  fit.metafor$ci.ub <- fit.metafor$b + 1.96 * fit.metafor$se
  
  # Variance components
  fit.metafor$vb <- fit.equalto$sdr$cov.fixed
  fit.metafor$sigma2 <- c(fit.equalto$fit$par["betadisp"], fit.equalto$fit$par["theta"])
  fit.metafor$tau2 <- fit.metafor$sigma2[1]
  fit.metafor$rho <- 0  # Assume zero unless provided elsewhere
  fit.metafor$gamma2 <- 0  # Assume zero unless explicitly modelled
  fit.metafor$phi <- 0  # Assume zero unless specified
  
  # Heterogeneity statistics
  fit.metafor$QE <- NA  # Not available directly from glmmTMB output
  fit.metafor$QEdf <- NA
  fit.metafor$QEp <- NA
  
  # Model structure
  fit.metafor$k <- nrow(fit.equalto$frame)
  fit.metafor$p <- 1
  fit.metafor$int.only <- TRUE
  fit.metafor$intercept <- TRUE
  
  # Data components
  fit.metafor$yi <- fit.equalto$frame$y
  fit.metafor$vi <- rep(NA, length(fit.metafor$yi))  # Variances not directly extracted
  fit.metafor$X <- matrix(1, nrow = length(fit.metafor$yi), ncol = 1, dimnames = list(NULL, "intrcpt"))
  fit.metafor$V <- diag(fit.metafor$vi)
  
  # Metadata
  fit.metafor$method <- "REML"
  fit.metafor$weighted <- TRUE
  fit.metafor$test <- "z"
  fit.metafor$digits <- rep(4, 9)
  
  class(fit.metafor) <- c("rma.mv", "rma")
  return(fit.metafor)
}
