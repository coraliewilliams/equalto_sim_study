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
  library(R.utils)
})


# ---- load functions ---------------------------------------------------------
source("R/02_functions.R")


# ---- parallel setup ----
Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
plan(multisession, workers = max(1, future::availableCores() - 2)) ##use all cores but 1
handlers("rstudio")
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

sim_func <- function(row_df, data_dir, res_dir) {
  
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
  data_file <- file.path(data_dir, sprintf("simdat_%05d.rds", par$param_id))
  tmpd <- paste0(data_file, ".tmp")
  saveRDS(dat, tmpd, compress = FALSE)
  file.rename(tmpd, data_file)
  
  # fit models
  fits <- list(
    fit_rma_uni(dat),
    fit_glmmTMB_RE(dat),
    if (par$measure == "OR")  fit_OR_models(dat)  else NULL,
    if (par$measure == "IRR") fit_IRR_models(dat) else NULL
  )
  fits <- fits[!vapply(fits, is.null, logical(1))]
  
  # combine model results + annotate
  res <- data.table::rbindlist(fits, use.names = TRUE, fill = TRUE)
  res$scenario  <- par$scenario
  res$replicate <- par$replicate
  res$param_id  <- par$param_id
  res$measure   <- par$measure
  res$dat_file  <- data_file
  
  res_file <- file.path(res_dir, sprintf("res_%05d.rds", par$param_id))
  tmpr <- paste0(res_file, ".tmp")
  saveRDS(res, tmpr, compress = FALSE)
  file.rename(tmpr, res_file)
  
  invisible(NULL)
}


# ## manual testing
# row_df <- example_row
# data_dir = "results/raw/SMD/data"
# res_dir = "results/raw/SMD/res"
# 
# ### test one row
# example_row <- PARAM_GRID_SMD[1716, , drop = FALSE]
# example_res <- sim_func(example_row, data_dir, res_dir)
# print(example_res)




# ---- Fucntions to run sims in parallel ---------------------

chunk_indices_by_count <- function(n, k_chunks) {
  split(seq_len(n), as.integer(cut(seq_len(n), breaks = k_chunks, labels = FALSE)))
}

run_chunk_df <- function(chunk_df, data_dir, res_dir, timeout_sec = 300L) {
  n <- nrow(chunk_df)
  for (i in seq_len(n)) {
    row_df <- chunk_df[i, , drop = FALSE]
    try({
      R.utils::withTimeout(
        sim_func(row_df, data_dir = data_dir, res_dir = res_dir),
        timeout = timeout_sec, onTimeout = "error"
      )
    }, silent = TRUE)
    if ((i %% 200L) == 0L) gc(FALSE)
  }
  invisible(NULL)
}


## single runner (parallel inside, sequential across grids)
run_one_grid <- function(PARAM_GRID, label, k_per_worker = 8L,
                         base_dir = "results/raw",
                         timeout_sec = 300L) {
  
  data_dir <- file.path(base_dir, label, "data")
  res_dir  <- file.path(base_dir, label, "res")
  if (!dir.exists(data_dir)) stop("Missing data dir: ", data_dir)
  if (!dir.exists(res_dir))  stop("Missing res dir: ", res_dir)
  
  n <- nrow(PARAM_GRID)
  W <- future::nbrOfWorkers()
  k_chunks  <- max(1L, min(n, W * k_per_worker))
  idx_split <- chunk_indices_by_count(n, k_chunks)
  chunk_dfs <- lapply(idx_split, function(ix) PARAM_GRID[ix, , drop = FALSE])
  
  pkg_needed <- c("glmmTMB", "metafor", "dplyr", "tibble", "data.table", "R.utils")
  
  invisible(
    furrr::future_map(
      chunk_dfs,
      ~ run_chunk_df(.x, data_dir = data_dir, res_dir = res_dir, timeout_sec = timeout_sec),
      .options = furrr::furrr_options(
        seed       = TRUE,
        packages   = pkg_needed,
        globals    = TRUE
      )
    )
  )
}



# ---- Read parameter grids ---------------------------------------------------------------
PARAM_GRID_SMD <- readr::read_csv("data/param_grid_smd.csv")
PARAM_GRID_lnRR <- readr::read_csv("data/param_grid_lnrr.csv")
PARAM_GRID_OR <- readr::read_csv("data/param_grid_OR.csv")
PARAM_GRID_IRR <- readr::read_csv("data/param_grid_IRR.csv")



# ---- Run sims and save output -----------------------------------------------
run_one_grid(PARAM_GRID_SMD,  "SMD", k_per_worker = 10L, base_dir = "results/raw")
run_one_grid(PARAM_GRID_lnRR, "lnRR", k_per_worker = 1L, base_dir = "results/raw")



run_one_grid(PARAM_GRID_OR[1:4000,], "OR", k_per_worker = 3L, base_dir = "results/raw")
run_one_grid(PARAM_GRID_IRR, "IRR", k_per_worker = 10L, base_dir = "results/raw")




#-----Sim data and result files concatenation -------------------------------


### function to concatenate all raw res and data files into one list
raw_conc <- function(label) {
  
  base_dir <- file.path("results", "raw", label)
  res_dir  <- file.path(base_dir, "res")
  dat_dir  <- file.path(base_dir, "data")
  
  # read all .rds files 
  res_files <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
  dat_files <- list.files(dat_dir, pattern = "\\.rds$", full.names = TRUE)
  
  # combine into lists (each file corresponds to one simulation replicate)
  res_list <- lapply(res_files, readRDS)
  dat_list <- lapply(dat_files, readRDS)
  
  # store in global environment as res_all_<label>
  assign(paste0("res_all_", label), res_list, envir = .GlobalEnv)
  
  # output file names
  out_res <- file.path("results", "raw", paste0("res_all_", label, ".rds"))
  out_dat <- file.path("results", "raw", paste0("simdat_all_", label, ".rds"))
  
  # save
  saveRDS(res_list, out_res)
  saveRDS(dat_list, out_dat)
  
  message("✓ combined and saved raw output: ", label)
}

### --- concatenate raw sim data and raw sim result files 
raw_conc("SMD")
raw_conc("lnRR")
raw_conc("OR")
raw_conc("SMD")


### ----- bind all results into one dataframe for analyses/plots 
library(data.table)

res_SMD <- rbindlist(res_all_SMD, use.names = TRUE, fill = TRUE)
write.csv(res_SMD, "results/res_smd.csv")

res_SMD <- rbindlist(res_all_lnRR, use.names = TRUE, fill = TRUE)
write.csv(res_SMD, "results/res_lnrr.csv")