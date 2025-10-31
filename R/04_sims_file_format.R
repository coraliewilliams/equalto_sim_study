
### ----- Function to concatenate all raw res and data files into one list

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

### concatenate raw sim data and raw sim result files
raw_conc("SMD")
raw_conc("lnRR")
raw_conc("OR")
raw_conc("SMD")


### bind all results into one dataframe for analyses/plots 
library(data.table)

res_SMD <- rbindlist(res_all_SMD, use.names = TRUE, fill = TRUE)
write.csv(res_SMD, "results/res_smd.csv")

res_SMD <- rbindlist(res_all_lnRR, use.names = TRUE, fill = TRUE)
write.csv(res_SMD, "results/res_lnrr.csv")
