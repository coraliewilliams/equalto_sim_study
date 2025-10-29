library(targets)
library(tarchetypes)
library(readr)
library(dplyr)
library(purrr)

# Load functions automatically from /R
tar_option_set(packages = c(
  "dplyr", "purrr", "glmmTMB", "metafor", "tibble"
))

source("R/03_run_sims.R")
source("R/02_functions.R")

list(
  # Load parameter grid only once
  tar_target(
    param_grid,
    readr::read_csv("data/param_grid_smd.csv"),
    format = "file" # ensures tracking changes
  ),
  
  # Run each row of param grid → returns one tibble per scenario/rep
  tar_target(
    sim_results_row,
    sim_func(param_grid[row, , drop = FALSE]),
    pattern = map(row = seq_len(nrow(param_grid)))
  ),
  
  # Combine all results into a single data frame
  tar_target(
    sim_results_all,
    bind_rows(sim_results_row)
  ),
  
  # Write to file once complete
  tar_target(
    save_results,
    {
      dir.create("results", showWarnings = FALSE)
      write_csv(sim_results_all, "results/raw/sim_smd_output.csv")
      "results/raw/sim_smd_output.csv"
    },
    format = "file"
  )
)
