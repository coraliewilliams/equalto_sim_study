# -------------------
# Run simultations 
# author: Coralie Williams 

##### Performance measures to record
# - run time in seconds
# - warning/error messages
# - errors or warnings
# - fixed effect estimate (mu_est)
# - fixed effect standarde error (mu_se)
# - fixed effect p-value (pvalue)
# - tau2 estimate (tau2_est)


suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(furrr)
  library(metafor)
  library(glmmTMB)
  library(data.table)
  library(progressr)
})


# ---- load functions ---------------------------------------------------------
source("R/02_functions.R")


# ---- parallel setup ----
Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
plan(multisession, workers = max(1, future::availableCores() - 1)) ##use all cores but 1
handlers("rstudio", global = TRUE)
options(progressr.enable = TRUE)
options(future.globals.maxSize = 2 * 1024^3) 


# ---- functions to record time + erors ----------------------------------------------------------------

time_fit <- function(expr) {
  
  t0   <- proc.time()[["elapsed"]]
  warn <- character(); err <- NULL
  
  fit <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") }
    ),
    error = function(e) { err <<- conditionMessage(e); structure(list(), class = "try-error") }
  )
  
  dt <- proc.time()[["elapsed"]] - t0
  list(fit = fit,
       time = dt,
       warn = if (length(warn)) paste(warn, collapse = " | ") else NA_character_,
       error = err)
}



# Fit model and extract info  --------------------------------------------------------------


# --- Random effect models (metafor + glmmTMB) 
fit_rma_uni <- function(dat) {
  tf <- time_fit(rma(yi, vi, data = dat, control = list(REMLf = FALSE)))
  if (inherits(tf$fit, "try-error")) return(NULL)
  m <- tf$fit
  data.table(model="rma.uni",
             runtime=tf$time,
             error=ifelse(is.null(tf$error), NA, tf$error),
             warn=tf$warn, 
             logLik=as.numeric(logLik(m)),
             est=as.numeric(m$b),
             se=as.numeric(m$se),
             pvalue=as.numeric(m$pval),
             tau2=as.numeric(m$tau2))
}

fit_glmmTMB_RE <- function(dat) {
  dat <- dat %>% mutate(g = 1L, id = factor(seq_len(n())))
  V <- diag(dat$vi)
  tf <- time_fit(glmmTMB(yi ~ 1 + equalto(0 + id | g, V), REML = TRUE, data = dat))
  if (inherits(tf$fit, "try-error")) return(NULL)
  m <- tf$fit
  co <- as.data.frame(summary(m)$coefficients$cond)
  data.table(model="glmmTMB_RE",
             runtime=tf$time,
             error=ifelse(is.null(tf$error), NA, tf$error), 
             warn=tf$warn,
             logLik=as.numeric(logLik(m)[1]),
             est=co$Estimate[1],
             se=co$`Std. Error`[1],
             pvalue=co$`Pr(>|z|)`[1],
             tau2=sigma(m)^2)
}



# --- OR: binomial–normal models (metafor + two glmmTMB) 
fit_OR_models <- function(dat_or) {
  out <- list()
  # metafor binomial-normal (UM.RS, Laplace)
  tf1 <- time_fit(rma.glmm(measure="OR", ai=a, bi=b, ci=c, di=d,
                           data=dat_or, model="UM.RS", method="ML",
                           nAGQ=1, drop00=FALSE))
  if (!inherits(tf1$fit, "try-error")) {
    m <- tf1$fit
    out[[length(out)+1]] <- data.table(model="rma.glmm_OR",
                                       runtime=tf1$time,
                                       error=ifelse(is.null(tf1$error), NA, tf1$error),
                                       warn=tf1$warn,
                                       logLik=as.numeric(logLik(m)),
                                       est=as.numeric(m$b),
                                       se=as.numeric(m$se),
                                       pvalue=as.numeric(m$pval),
                                       tau2=as.numeric(m$tau2))
  }
  
  # long format for glmmTMB
  n       <- c(dat_or$a+dat_or$b, dat_or$c+dat_or$d)
  event   <- c(dat_or$a, dat_or$c)
  study   <- rep(dat_or$study, 2)
  treat   <- c(rep(1, nrow(dat_or)), rep(0, nrow(dat_or)))
  control <- 1 - treat
  treat12 <- treat - 0.5
  or_long <- data.table(study, n, event, treat, control, treat12)

  # glmmTMB variant 1
  tf2 <- time_fit(glmmTMB(cbind(event, n-event) ~ factor(treat) + (treat12-1|study) + (1|study),
                          family=binomial("logit"),
                          dispformula = ~ 0, 
                          REML=FALSE,
                          data=or_long))
  if (!inherits(tf2$fit, "try-error")) {
    m <- tf2$fit
    est <- summary(m)$coefficients$cond[2,]
    out[[length(out)+1]] <- data.table(model="glmmTMB.binomial.1", runtime=tf2$time,
                                   error=ifelse(is.null(tf2$error), NA, tf2$error), warn=tf2$warn,
                                   logLik=as.numeric(logLik(m)[1]),
                                   est=as.numeric(unlist(fixef(m))[2]),
                                   se=as.numeric(sqrt(vcov(m)$cond[2,2])),
                                   pvalue=summary(m)$coefficients$cond[2, "Pr(>|z|)"],
                                   tau2 = as.numeric(VarCorr(m)$cond$study[1]))
  }
  
  
  # glmmTMB variant 2
  tf3 <- time_fit(glmmTMB(cbind(event, n-event) ~ factor(treat) + (control + treat - 1 | study),
                          family=binomial("logit"),
                          dispformula = ~ 0, 
                          REML=FALSE,
                          data=or_long))
  if (!inherits(tf3$fit, "try-error")) {
    m <- tf3$fit
    est <- summary(m)$coefficients$cond[2,]
    VC  <- VarCorr(m)$cond$study
    tau2_combo <- as.numeric(VC[1,1] + VC[2,2] - 2*VC[1,2])
    out[[length(out)+1]] <- data.table(model="glmmTMB.binomial.2", runtime=tf3$time,
                                   error=ifelse(is.null(tf3$error), NA, tf3$error), warn=tf3$warn,
                                   logLik=as.numeric(logLik(m)[1]),
                                   est=as.numeric(unlist(fixef(m))[2]),
                                   se=as.numeric(sqrt(vcov(m)$cond[2,2])),
                                   pvalue=summary(m)$coefficients$cond[2, "Pr(>|z|)"],
                                   tau2 = tau2_combo)
  }
  
  bind_rows(out)
}



# --- IRR: poisson–normal models (metafor + glmmTMB) 
fit_IRR_models <- function(dat_irr) {
  out <- list()
  
  # metafor poisson-normal (UM.RS, Laplace)
  tf1 <- time_fit(rma.glmm(measure="IRR",
                           x1i=yT, t1i=t_T, x2i=yC, t2i=t_C,
                           data=dat_irr, model="UM.RS", method="ML",
                           nAGQ=1, drop00=FALSE))
  if (!inherits(tf1$fit, "try-error")) {
    m <- tf1$fit
    out[[length(out)+1]] <- data.table(model="rma.glmm_IRR",
                                       runtime=tf1$time,error=ifelse(is.null(tf1$error), NA, tf1$error),
                                       warn=tf1$warn,
                                       logLik=as.numeric(logLik(m)),
                                       est=as.numeric(m$b),
                                       se=as.numeric(m$se),
                                       pvalue=as.numeric(m$pval),
                                       tau2=as.numeric(m$tau2))
  }
  
  # long format for glmmTMB
  irr_long <- data.table(
    study = c(rep(dat_irr$study,2)),
    trt   = c(rep(-0.5, nrow(dat_irr)), rep(0.5, nrow(dat_irr))), #-0.5 = C, 0.5 = T => to match metafor contrast coding
    t     = c(dat_irr$t_C, dat_irr$t_T),
    y     = c(dat_irr$yC,   dat_irr$yT))
  
  tf2 <- time_fit(glmmTMB(y ~ factor(trt) + offset(log(t)) + (0+trt|study) + (1|study),
                          family=poisson("log"), REML=FALSE, data=irr_long))
  if (!inherits(tf2$fit, "try-error")) {
    m <- tf2$fit
    co <- as.data.frame(summary(m)$coefficients$cond)
    out[[length(out)+1]] <- data.table(model="glmmTMB.poisson",
                                       runtime=tf2$time,
                                       error=ifelse(is.null(tf2$error), NA, tf2$error),
                                       warn=tf2$warn,
                                       logLik=as.numeric(logLik(m)[1]),
                                       est=co$Estimate[2],
                                       se=co$`Std. Error`[2],
                                       pvalue=co$`Pr(>|z|)`[2],
                                       tau2=as.numeric(VarCorr(m)$cond$study[1]))
  }
  
  rbindlist(out, use.names = TRUE)
}




# ---- Sim function --------------------------------------------------------------

sim_func <- function(row_df, dat_dir) {
  
  # parse parameters
  par <- row_df %>%
    mutate(
      tau2 = as.numeric(tau2), mu = as.numeric(mu),
      baseline_value = as.numeric(baseline_value),
      k = as.integer(k), n_t = as.integer(n_t), n_c = as.integer(n_c),
      t_t = as.integer(t_t), t_c = as.integer(t_c),
      sd_t = as.numeric(sd_t), sd_c = as.numeric(sd_c),
      mu_c = as.numeric(mu_c),
      seed = as.integer(seed), replicate = as.integer(replicate)
    ) %>% as.list()
  
  if (!is.na(par$seed)) set.seed(par$seed)
  
  # simulate dataset
  dat <- switch(par$measure,
                "SMD"  = simulate_smd(k = par$k, mu = par$mu, tau2 = par$tau2,
                                      n_T = par$n_t, n_C = par$n_c, mu_C = par$mu_c, sigma = par$sd_t),
                "lnRR" = simulate_lnrr(k = par$k, mu = par$mu, tau2 = par$tau2,
                                       n_T = par$n_t, n_C = par$n_c, mu_C = par$mu_c, sigma = par$sd_t),
                "OR"   = simulate_or(k = par$k, mu = par$mu, tau2 = par$tau2,
                                     n_T = par$n_t, n_C = par$n_c, p_C = par$baseline_value),
                "IRR"  = simulate_irr(k = par$k, mu = par$mu, tau2 = par$tau2,
                                      t_T = par$t_t, t_C = par$t_c, lambda_C = par$baseline_value),
                stop("Unknown measure: ", par$measure)
  )
  
  # save sim dataset 
  dat_path <- file.path(dat_dir, sprintf("simdat_%05d.rds", par$param_id))
  saveRDS(dat, dat_path, compress = FALSE)
  
  # fit models
  fits <- list(
    fit_rma_uni(dat),
    fit_glmmTMB_RE(dat),
    if (par$measure == "OR")  fit_OR_models(dat)  else NULL,
    if (par$measure == "IRR") fit_IRR_models(dat) else NULL
  )
  fits <- fits[!vapply(fits, is.null, logical(1))]
  
  # combine model results + annotate
  out <- data.table::rbindlist(fits, use.names = TRUE, fill = TRUE)
  out$scenario  <- par$scenario
  out$replicate <- par$replicate
  out$param_id  <- par$param_id
  out$measure   <- par$measure
  out$dat_file  <- dat_path
  
  out
}


# ## manual testing
# row_df <- example_row
#  
# ### test one row
# example_row <- PARAM_GRID[200, , drop = FALSE]
# example_res <- sim_func(example_row)
# print(example_res)




# ---- Fucntion sims in parallel ---------------------

run_one_grid <- function(PARAM_GRID, label, k_per_worker = 8L,
                         base_dat_dir = "results/raw/dat",
                         base_results_dir = "results/raw/results") {
  
  # check folders exist
  dat_dir     <- file.path(base_dat_dir, label)
  results_dir <- file.path(base_results_dir)
  if (!dir.exists(dat_dir))     stop("Missing data dir: ", dat_dir)
  if (!dir.exists(results_dir)) stop("Missing results dir: ", results_dir)
  
  # chunking
  n <- nrow(PARAM_GRID)
  W <- future::nbrOfWorkers()
  chunk_list <- chunk_indices_by_count(n, k_chunks = max(1L, W * k_per_worker))
  
  # parallel map
  pkg_needed <- c("glmmTMB", "metafor", "dplyr", "tibble", "data.table")
  progressr::with_progress({
    res_chunks <- furrr::future_map(
      chunk_list,
      ~ run_chunk(PARAM_GRID, .x, dat_dir = dat_dir),
      .progress = TRUE,
      .options  = furrr::furrr_options(seed = TRUE, packages = pkg_needed, globals = TRUE)
    )
    
    # bind + save all results in one csv
    results <- data.table::rbindlist(res_chunks, use.names = TRUE, fill = TRUE)
    out_csv <- file.path(results_dir, sprintf("sim_results_%s.csv", label))
    data.table::fwrite(results, out_csv)
  })
}



# ---- Read parameter grids ---------------------------------------------------------------
PARAM_GRID_SMD <- readr::read_csv("data/param_grid_smd.csv")
PARAM_GRID_lnRR <- readr::read_csv("data/param_grid_lnrr.csv")
PARAM_GRID_OR <- readr::read_csv("data/param_grid_OR.csv")
PARAM_GRID_IRR <- readr::read_csv("data/param_grid_IRR.csv")


# ---- Run sims and save output -----------------------------------------------
run_one_grid(PARAM_GRID_SMD, "SMD")
run_one_grid(PARAM_GRID_lnRR, "lnRR")
run_one_grid(PARAM_GRID_OR, "OR")
run_one_grid(PARAM_GRID_IRR, "IRR")


