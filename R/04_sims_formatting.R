
### function to concatenate all raw res and data files into one list

raw_conc <- function(label) {
  
  # paths
  base_dir <- file.path("results", "raw", label)
  res_dir  <- file.path(base_dir, "res")
  dat_dir  <- file.path(base_dir, "data")
  
  # read all .rds files in each subfolder
  res_files <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
  dat_files <- list.files(dat_dir, pattern = "\\.rds$", full.names = TRUE)
  
  # combine into lists (each element corresponds to one simulation replicate)
  res_list <- lapply(res_files, readRDS)
  dat_list <- lapply(dat_files, readRDS)
  
  # output file names
  out_res <- file.path("results", "raw", paste0("res_all_", label, ".rds"))
  out_dat <- file.path("results", "raw", paste0("simdat_all_", label, ".rds"))
  
  # save
  saveRDS(res_list, out_res)
  saveRDS(dat_list, out_dat)
  
  message("✓ Combined and saved: ", label)
}




### concatenate raw files
raw_conc("SMD")
raw_conc("lnRR")
raw_conc("OR")
raw_conc("SMD")




### function to read all single files and row bind to 