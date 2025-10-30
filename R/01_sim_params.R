# -------------------
# Parameter grid for meta-analysis simulation study
# author: Coralie Williams 



suppressPackageStartupMessages({
  library(tidyverse)
})

# Sim settings and fixed parameters -----------------------------------------

REPS <- 10L   # replicates per unique condition
K    <- 20L   # number of studies per meta-analysis

# Fixed sample sizes (treatment and control groups)
N_T <- 50L; N_C <- 50L   # for SMD, lnRR, OR
T_T <- 50L; T_C <- 50L   # exposure times for IRR

# Fixed within-group SDs and baseline mean for continuous measures
SD_T <- 10; SD_C <- 10
MU_C <- 100


# Varying parameters per effect size measure -------------------------------
spec_es <- list(
  SMD = list(
    tau2_levels = c(0.04, 0.09, 0.25),
    mu_levels   = c(null = 0, moderate = 0.3),
    vary_baseline = FALSE
  ),
  
  lnRR = list(
    tau2_levels = c(0.04, 0.16, 0.36),
    mu_levels   = c(null = 0, moderate = log(1.3)),
    vary_baseline = FALSE
  ),
  
  OR = list(
    tau2_levels     = c(0.04, 0.16, 0.36),
    mu_levels       = c(null = 0, moderate = log(1.5)),
    vary_baseline   = TRUE,
    baseline_name   = "p_c",
    baseline_levels = c(moderate = 0.20, rare = 0.05)
  ),
  
  IRR = list(
    tau2_levels     = c(0.04, 0.16, 0.36),
    mu_levels       = c(null = 0, moderate = log(1.5)),
    vary_baseline   = TRUE,
    baseline_name   = "lambda_c",
    baseline_levels = c(moderate = 0.10, rare = 0.01)
  )
)





# -- Param grid function -------------------------------------------------------


make_param_grid <- function(reps = REPS, k = K, measures = c("SMD", "lnRR", "OR", "IRR")) {
  stopifnot(all(measures %in% names(spec_es)))
  
  grid <- purrr::map_dfr(measures, function(measure) {
    sp <- spec_es[[measure]]
    
    # pair mu with its label 
    mu_tbl <- tibble::tibble(
      mu       = unname(sp$mu_levels),
      mu_label = names(sp$mu_levels)
    )
    
    # cross μ pairs with measure and tau²
    base <- tidyr::crossing(
      mu_tbl,
      measure = measure,
      tau2    = sp$tau2_levels
    )
    
    # properly paired baseline 
    if (isTRUE(sp$vary_baseline)) {
      baseline_tbl <- tibble::tibble(
        baseline_label = names(sp$baseline_levels),
        baseline_value = unname(sp$baseline_levels)
      )
      base <- tidyr::crossing(base, baseline_tbl)
    } else {
      base <- dplyr::mutate(base,
                            baseline_label = NA_character_,
                            baseline_value = NA_real_
      )
    }
    
    # Append fixed design elements + replicate expansion
    base |>
      dplyr::mutate(
        k   = k,
        n_t = dplyr::if_else(measure %in% c("SMD", "lnRR", "OR"), N_T, NA_integer_),
        n_c = dplyr::if_else(measure %in% c("SMD", "lnRR", "OR"), N_C, NA_integer_),
        t_t = dplyr::if_else(measure == "IRR", T_T, NA_integer_),
        t_c = dplyr::if_else(measure == "IRR", T_C, NA_integer_),
        sd_t = dplyr::if_else(measure %in% c("SMD", "lnRR"), SD_T, NA_real_),
        sd_c = dplyr::if_else(measure %in% c("SMD", "lnRR"), SD_C, NA_real_),
        mu_c = dplyr::if_else(measure %in% c("SMD", "lnRR"), MU_C, NA_real_),
        replicate = list(seq_len(reps))
      ) |>
      tidyr::unnest(replicate)
  })
  
  grid |>
    dplyr::arrange(measure, tau2, mu_label, baseline_label, replicate) |>
    dplyr::mutate(
      param_id = dplyr::row_number(),
      seed = as.integer((param_id * 10L + 4728L) %% .Machine$integer.max),
      scenario = dplyr::case_when(
        is.na(baseline_label) ~ sprintf("%s_tau2=%.2f_%s", measure, tau2, mu_label),
        TRUE ~ sprintf("%s_tau2=%.2f_%s_%s", measure, tau2, mu_label, baseline_label)
      )
    ) |>
    dplyr::relocate(param_id, seed, scenario)
}



# param grid for SMD 
PARAM_GRID_SMD <- make_param_grid(reps=1000, measures="SMD")
write.csv(PARAM_GRID_SMD, "data/param_grid_smd.csv")

# param grid for lnRR
PARAM_GRID_LNRR<- make_param_grid(reps=1000, measures="lnRR")
write.csv(PARAM_GRID_LNRR, "data/param_grid_lnRR.csv")

# param grid for OR
PARAM_GRID_OR <- make_param_grid(reps=1000, measures="OR")
write.csv(PARAM_GRID_OR, "data/param_grid_OR.csv")

# param grid for IRR
PARAM_GRID_IRR <- make_param_grid(reps=1000, measures="IRR")
write.csv(PARAM_GRID_IRR, "data/param_grid_IRR.csv")






#### expect 12 scenarios/conditions --- for each SMD, lnRR
# 3 tau^2 x 2 mu = 6 conditions

#### expect 24 scenarios/conditions ----- for OR, IRR
# 3 tau^2 x 2 mu x 2 baseline group mean = 12 conditions

## total sims
# If REP=1000 ===> 12*1000 + 24*1000 = 36,000


# #### check output
# scenario_params <- PARAM_GRID |>
#   dplyr::group_by(scenario, measure, tau2, mu, mu_label,
#                   baseline_label, baseline_value,
#                   k, n_t, n_c, t_t, t_c, sd_t, sd_c, mu_c) |>
#   dplyr::summarise(n_rows = dplyr::n()) |>
#   dplyr::arrange(measure, tau2, mu_label, baseline_label)
# 
# scenario_params
