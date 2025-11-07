library(ggplot2);library(devtools);library(cowplot);library(ggdark);theme_set(dark_theme_bw())
library(ggdist);library(tidyverse);library(dplyr);library(patchwork); library(scales);
library(latex2exp); library(xtable); library(grid);library(readr)


# load result files
res_smd <- read_csv("results/res_smd.csv")
res_lnrr <- read_csv("results/res_lnrr.csv")
res_OR <- read_csv("results/res_OR.csv")
res_IRR <- read_csv("results/res_IRR.csv")


################# 1. Convergence #########################

res_smd$conv <- ifelse(is.na(res_smd$warn)|is.na(res_smd$error),0,1)
res_lnrr$conv <- ifelse(is.na(res_lnrr$warn)|is.na(res_lnrr$error),0,1)
res_OR$conv <- ifelse(is.na(res_OR$warn)|is.na(res_OR$error),0,1)
res_IRR$conv <- ifelse(is.na(res_IRR$warn)|is.na(res_IRR$error),0,1)



table(res_smd$model, res_smd$conv)
table(res_lnrr$model, res_lnrr$conv)
table(res_OR$model, res_OR$conv)
table(res_smd$model, res_smd$conv)


################## 2. Formatting ##########################





################# 3. Derive sim results ####################




################ 4. Plots/tables - overall mean ###############

plot(res_smd$est[which(res_smd$model=="rma.uni")],
     res_smd$est[which(res_smd$model=="glmmTMB_RE")])


plot(res_lnrr$est[which(res_lnrr$model=="rma.uni")],
     res_lnrr$est[which(res_lnrr$model=="glmmTMB_RE")])
