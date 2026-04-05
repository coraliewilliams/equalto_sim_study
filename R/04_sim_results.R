library(ggplot2);library(devtools);library(cowplot);
library(ggdist);library(tidyverse);library(dplyr);library(patchwork); library(scales);
library(latex2exp); library(xtable); library(grid);library(readr); library(data.table)

#### For plots
theme_update(plot.title = element_text(colour = "black"))

# color-blind friendly palette (Okabe–Ito)
col.m <- c(
  SMD  = "#0072B2",  # blue
  lnRR = "#009E73",  # green
  OR   = "#D55E00",  # vermilion
  IRR  = "#CC79A7"   # purple
)
col.m.pastel <- alpha(col.m, 0.4)

# load simulation result files (one for each effect size measure)
res_smd <- read_csv("results/res_smd.csv")
res_lnrr <- read_csv("results/res_lnrr.csv")
res_OR <- read_csv("results/res_OR.csv")
res_IRR <- read_csv("results/res_IRR.csv")


################# 1. Convergence #########################

# #### There are zero "error" messages for all simulated models
# res_smd[!is.na(res_smd$error),]
# res_lnrr[!is.na(res_lnrr$error),]
# res_OR[!is.na(res_OR$error),]
# res_IRR[!is.na(res_IRR$error),]
# #### Only IRR and OR have warning messages
# res_smd[!is.na(res_smd$warn),]
# res_lnrr[!is.na(res_lnrr$warn),]
# res_OR[!is.na(res_OR$warn),]
# res_IRR[!is.na(res_IRR$warn),]

# # overview of type of warning messages / which loglikelihoods were NA
# table(res_OR$warn)
# table(res_IRR$warn)
# res_OR[is.na(res_OR$logLik),]
# res_OR[is.na(res_IRR$logLik),]


res_smd$conv <- is.na(res_smd$warn)&is.na(res_smd$error)
res_lnrr$conv <- is.na(res_lnrr$warn)&is.na(res_lnrr$error)
res_OR$conv <- is.na(res_OR$warn)&is.na(res_OR$error)
res_IRR$conv <- is.na(res_IRR$warn)&is.na(res_IRR$error)


table(res_smd$model, res_smd$conv)/6000*100
table(res_lnrr$model, res_lnrr$conv)/6000*100
table(res_OR$model, res_OR$conv)/12000*100
table(res_IRR$model, res_IRR$conv)/12000*100

# not captured in the above output:
## -- for OR - 1 rma.uni_OR model did not return anything (Fisher score didn't converge)
## -- for IRR - 171 rma.glmm_IRR models did not return anything ("step size truncated due to divergence")
# ---> see 03_run_sims.R script CHECKS at the end (line 429) for more details



################## 2. Formatting ##########################

# true mu value
res_smd$mu <- ifelse(grepl("_null", res_smd$scenario), 0, 
                     ifelse(grepl("_moderate", res_smd$scenario), 0.3, NA_real_))
res_lnrr$mu <- ifelse(grepl("_null", res_lnrr$scenario), 0, 
                     ifelse(grepl("_moderate", res_lnrr$scenario), log(1.3), NA_real_))
res_OR$mu <- ifelse(grepl("_null", res_OR$scenario), 0, 
                     ifelse(grepl("_moderate", res_OR$scenario), log(1.5), NA_real_))
res_IRR$mu <- ifelse(grepl("_null", res_IRR$scenario), 0, 
                     ifelse(grepl("_moderate", res_IRR$scenario), log(1.5), NA_real_))

# true tau2 value
m <- regexpr("tau2=([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?\\d+)?)", res_smd$scenario, perl = TRUE)
res_smd$true_tau2 <- ifelse(m > 0, as.numeric(sub("tau2=([^ ]+).*", "\\1", regmatches(res_smd$scenario, m))), NA_real_)

m <- regexpr("tau2=([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?\\d+)?)", res_lnrr$scenario, perl = TRUE)
res_lnrr$true_tau2 <- ifelse(m > 0, as.numeric(sub("tau2=([^ ]+).*", "\\1", regmatches(res_lnrr$scenario, m))), NA_real_)

m <- regexpr("tau2=([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?\\d+)?)", res_OR$scenario, perl = TRUE)
res_OR$true_tau2 <- ifelse(m > 0, as.numeric(sub("tau2=([^ ]+).*", "\\1", regmatches(res_OR$scenario, m))), NA_real_)

m <- regexpr("tau2=([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?\\d+)?)", res_IRR$scenario, perl = TRUE)
res_IRR$true_tau2 <- ifelse(m > 0, as.numeric(sub("tau2=([^ ]+).*", "\\1", regmatches(res_IRR$scenario, m))), NA_real_)


### Filter out all sim repetitions if at least one model had a convergence warning or didn't compute

# remove all sim repetitions where rma.glmm didn't converge (no result) for proper comparison
# this is for the plot of agreement Poisson-Normal
res_IRR2 <- res_IRR |>
  group_by(param_id) |>
  filter(n() >= 4) |>
  ungroup()
# remove sim repetitions where rma.uni didn't converge (no result) for proper comparison
res_OR2 <- res_OR |> 
  filter(param_id != "11046") 


## Bind all results from effect size measures
res_all <- bind_rows(
  res_smd,
  res_lnrr,
  res_OR2,
  res_IRR2
)

## Keep only models with no warnings
# i.e. filter out non-converge/warnings for all models but grouping them by name
group1 <- c("rma.uni", "glmmTMB_RE")
group2 <- c("rma.glmm_OR", "glmmTMB.binomial.1", "glmmTMB.binomial.2", 
            "rma.glmm_IRR", "glmmTMB.poisson")

res <- res_all |>
  mutate(model_group = dplyr::case_when(
      model %in% group1 ~ "group1_NN",
      model %in% group2 ~ "group2_GLMM",
      TRUE ~ NA_character_ )) |>
  filter(!is.na(model_group)) |>
  group_by(param_id, scenario, measure, model_group) |>
  filter(!any(conv == FALSE, na.rm = TRUE)) |>
  ungroup()


#### Summaries for latex tables of convergence for Supplementary Materials 
# All sim results (Table S2)
tab.all <- as.data.frame(table(res_all$conv, res_all$model, res_all$measure))
res_all_summary <- tab.all |> 
  filter((Var2 == "rma.uni" | Var2 == "glmmTMB_RE" |
            Var2 == "rma.glmm_OR" | Var2 == "glmmTMB.binomial.1" | Var2 == "glmmTMB.binomial.2" | 
            Var2 == "rma.glmm_IRR" |  Var2 == "glmmTMB.poisson" ) & Var1 == TRUE) |> 
  dplyr::select(Var2, Var3, Freq) |> 
  pivot_wider(names_from = Var3, values_from = Freq) |> 
  mutate(model = factor(Var2, levels=c("rma.uni", "glmmTMB_RE", "rma.glmm_OR",
                                       "glmmTMB.binomial.1", "glmmTMB.binomial.2",
                                       "rma.glmm_IRR", "glmmTMB.poisson"))) |> 
  arrange(model)


# Converged sim results (Table S3)
tab.all <- as.data.frame(table(res$conv, res$model, res$measure))
res_all_summary <- tab.all |> 
  filter((Var2 == "rma.uni" | Var2 == "glmmTMB_RE" |
            Var2 == "rma.glmm_OR" | Var2 == "glmmTMB.binomial.1" | Var2 == "glmmTMB.binomial.2" | 
            Var2 == "rma.glmm_IRR" |  Var2 == "glmmTMB.poisson" ) & Var1 == TRUE) |> 
  dplyr::select(Var2, Var3, Freq) |> 
  pivot_wider(names_from = Var3, values_from = Freq) |> 
  mutate(model = factor(Var2, levels=c("rma.uni", "glmmTMB_RE", "rma.glmm_OR",
                                       "glmmTMB.binomial.1", "glmmTMB.binomial.2",
                                       "rma.glmm_IRR", "glmmTMB.poisson"))) |> 
  arrange(model)





################# 3. Derive sim results on converged cases ####################


# derive Wald z-test CI - forgot to store CI - and save pvalue for plots
zcrit  <- qnorm(1 - 0.05/2)
res$z <- res$est / res$se #H0 assumes true mean is zero
res$p_wald <- 2 * pnorm(-abs(res$z))
res$mu_ci_lb <- res$est - zcrit*res$se #ci lower bound
res$mu_ci_ub <- res$est + zcrit*res$se #ci upper bound
res$mu_ci_width <- res$mu_ci_ub - res$mu_ci_lb #ci width
# coverage mu estimate for Wald test
res$cov_mu <- res$mu_ci_lb <= res$mu & res$mu <= res$mu_ci_ub

# derive t-test p-value and CI ---> not used but in case I get a question about it
df <- 19 - 1 #all sim datasets were nrow=nstudy=20
t_stat <- (res$est-0) / res$se
res$p_t <- 2 * pt(-abs(t_stat), df = df)
res$mu_ci_lb_t <- res$est - qt(1 - 0.05/2, df = df) * res$se
res$mu_ci_ub_t <- res$est + qt(1 - 0.05/2, df = df) * res$se
res$mu_ci_width_t <- res$mu_ci_ub_t - res$mu_ci_lb_t
# coverage mu estimate for t-test
res$cov_mu_t <- res$mu >= res$mu_ci_lb_t & res$mu <= res$mu_ci_ub_t


# sampling variance
res_sample_var <- res |> 
  summarise(mu_S2 = sum((mu - mean(mu))^2) / (n() - 1),
            tau2_S2 = sum((tau2 - mean(tau2))^2) / (n() - 1),
            .by = c(measure, model, scenario))

# RMSE dataset for overall mean mu 
rmse.dat.mu <- res |> 
  summarise(rmse = sqrt(mean((mu - est)^2)),
            .by = c(measure, model, scenario))

# RMSE dataset for overall mean mu 
rmse.dat.tau2 <- res |> 
  summarise(rmse = sqrt(mean((true_tau2 - tau2)^2)),
            .by = c(measure, model, scenario))

# coverage
cov.dat <- res |> 
  summarise(cov_prop = mean(cov_mu, na.rm = TRUE),
            n = n(),
            .by = c(measure, model, scenario))

# CI width dataset for overall mean mu 
ci.dat <- res |> 
  summarise(ci_width = mean(mu_ci_width, na.rm = TRUE),
            n = n(),
            .by = c(measure, model, scenario))

# Type I error
typeI.dat <- res |>
  filter(mu == 0) |>
  summarise(typeI = mean(pvalue < 0.05, na.rm = TRUE),
            n = sum(!is.na(p_wald)),
            .by = c(measure, model, scenario))

# Power 
power.dat <- res |>
  filter(mu != 0) |> 
  summarise(power = 1 - mean(p_wald > 0.05, na.rm = TRUE),
            n = n(),
            .by = c(measure, model, scenario))






############### 4. Plots/Tables - runtime ###############

# table mean runtime
res_runtime <- res |>
  group_by(measure, model) |>
  summarise(mean_runtime = mean(runtime, na.rm = TRUE),
            sd_runtime= sd(runtime, na.rm = TRUE),
            n = n(),
            .groups = "drop"
  )
print(xtable(res_runtime, digits = c(0, 0, 0, 2, 2, 0)),
      include.rownames = FALSE)



# plot for SMD and lnRR
p_rt1 <- res |>
  filter(measure %in% c("SMD", "lnRR"),
         model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = runtime, colour = measure, fill = measure)) +
  geom_boxplot() +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = "time(sec)", title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for oR
p_rt2 <- res |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = runtime, colour = measure, fill = measure)) +
  geom_boxplot() +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "time(sec)", title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_rt3 <- res |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = runtime, colour = measure, fill = measure)) +
  geom_boxplot() +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "time(sec)", title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for runtime figure
plot_rt <- p_rt1 + p_rt2 + p_rt3 +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_runtime.pdf", width = 5, height = 12)






# plot for log10 SMD and lnRR
p_logrt1 <- res |>
  filter(measure %in% c("SMD", "lnRR"),
         model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = runtime+0.001, colour = measure, fill = measure)) +
  geom_boxplot() +
  scale_y_log10() + 
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = "time(sec)", title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for oR
p_logrt2 <- res |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = runtime+0.001, colour = measure, fill = measure)) +
  geom_boxplot() +
  scale_y_log10() + 
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "time(sec)", title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_logrt3 <- res |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = runtime+0.001, colour = measure, fill = measure)) +
  geom_boxplot() +
  scale_y_log10() + 
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "time(sec)", title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for runtime figure
plot_logrt <- p_logrt1 + p_logrt2 + p_logrt3 +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_logruntime.pdf", width = 5, height = 12)
ggsave(filename = "results/figures/Figure_logruntime.png", width = 5, height = 12)






################ 5. Plots - converged results agreements ###############

#### NORMAL-NORMAL
# (1) mu estimate dataset - package agreement
mu_hat <- res |> 
  filter(model %in% c("rma.uni", "glmmTMB_RE")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(),
         measure = factor(measure, levels = c("SMD", "lnRR", "OR", "IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, est) |> 
  pivot_wider(names_from = model, values_from = est) |> 
  drop_na(`rma.uni`, glmmTMB_RE) 

# plot
plot_mu <- ggplot(mu_hat, aes(x = `rma.uni`, y = glmmTMB_RE, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ measure, scales="free", ncol=4) +
  labs(
    x = TeX("$\\hat{\\mu}$ rma.uni"),
    y = TeX("$\\hat{\\mu}$ glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# (2) mu SE dataset - package difference
mu_se_hat <- res |> 
  filter(model %in% c("rma.uni", "glmmTMB_RE")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(),
         measure = factor(measure, levels = c("SMD", "lnRR", "OR", "IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, se) |> 
  pivot_wider(names_from = model, values_from = se) |> 
  drop_na(`rma.uni`, glmmTMB_RE) 

# plot
plot_mu_se <- ggplot(mu_se_hat, aes(x = `rma.uni`, y = glmmTMB_RE, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ measure, scales="free", ncol=4) +
  labs(
    x = TeX("$\\hat{\\mu}$ standard error rma.uni"),
    y = TeX("$\\hat{\\mu}$ standard error glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# (3) tau2 estimate dataset - package agreement
tau2_hat <- res_all |> 
  filter(model %in% c("rma.uni", "glmmTMB_RE")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(),
         measure = factor(measure, levels = c("SMD", "lnRR", "OR", "IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, tau2) |> 
  pivot_wider(names_from = model, values_from = tau2) |> 
  drop_na(`rma.uni`, glmmTMB_RE) 

plot_tau2 <- ggplot(tau2_hat, aes(x = `rma.uni`, y = glmmTMB_RE, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ measure, scales="free", ncol=4) +
  labs(
    x = TeX("$\\hat{\\tau^2}$ rma.uni"),
    y = TeX("$\\hat{\\tau^2}$ glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

### patchwork for figure of estimates
plot_est <- plot_mu + plot_mu_se + plot_tau2 + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_NormalNormal_converged.pdf", width = 11, height = 9)
ggsave(filename = "results/figures/Figure_NormalNormal_converged.png", width = 11, height = 9)



#### BINOMIAL-NORMAL
# mu estimate plot - package agreement
mu_hat_binom_plot <- res |> 
  filter(model %in% c("glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("OR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, est) |>
  pivot_wider(names_from = model, values_from = est) |>
  drop_na(rma.glmm_OR, `glmmTMB.binomial.1`, `glmmTMB.binomial.2`) |> 
  #stack the two glmmTMB models for faceting
  pivot_longer(cols = c(`glmmTMB.binomial.1`, `glmmTMB.binomial.2`),
               names_to = "glmm_model",
               values_to = "est_glmm") |> 
  ggplot(aes(x = rma.glmm_OR, y = est_glmm, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ glmm_model, ncol = 2, scales = "free") +
  labs(
    x = TeX("$\\hat{\\mu}$ rma.glmm[OR]"),
    y = TeX("$\\hat{\\mu}$ glmmTMB"),
    title = "Binomial-Normal models agreement"
  ) +
  theme_bw() +
  theme(legend.position = "none")


# mu SE plot - package agreement
mu_se_binom_plot <- res |> 
  filter(model %in% c("glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("OR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, se) |>
  pivot_wider(names_from = model, values_from = se) |>
  drop_na(rma.glmm_OR, `glmmTMB.binomial.1`, `glmmTMB.binomial.2`) |> 
  #stack the two glmmTMB models for faceting
  pivot_longer(cols = c(`glmmTMB.binomial.1`, `glmmTMB.binomial.2`),
               names_to = "glmm_model",
               values_to = "est_glmm") |> 
  ggplot(aes(x = rma.glmm_OR, y = est_glmm, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ glmm_model, ncol = 2, scales = "free") +
  labs(
    x = TeX("$\\hat{\\mu}$ standard error rma.glmm[OR]"),
    y = TeX("$\\hat{\\mu}$ standard error glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# tau2 estimate plot - package aggreement
tau2_hat_binom_plot <- res |> 
  filter(model %in% c("glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("OR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, tau2) |>
  pivot_wider(names_from = model, values_from = tau2) |>
  drop_na(rma.glmm_OR, `glmmTMB.binomial.1`, `glmmTMB.binomial.2`) |> 
  #stack the two glmmTMB models for faceting
  pivot_longer(cols = c(`glmmTMB.binomial.1`, `glmmTMB.binomial.2`),
               names_to = "glmm_model",
               values_to = "est_glmm") |> 
  ggplot(aes(x = rma.glmm_OR, y = est_glmm, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ glmm_model, ncol = 2, scales = "free") +
  labs(
    x = TeX("$\\hat{\\tau^2}$ rma.glmm[OR]"),
    y = TeX("$\\hat{\\tau^2}$ glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

### patchwork for figure of estimates
plot_est_binom <- mu_hat_binom_plot + mu_se_binom_plot + tau2_hat_binom_plot + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_binom_converged.pdf", width = 7, height = 10)


#### POISSON-NORMAL
# mu estimate plot - package agreement
mu_hat_pois_plot <- res |> 
  filter(model %in% c("glmmTMB.poisson", "rma.glmm_IRR")) |>
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, est) |>
  pivot_wider(names_from = model, values_from = est) |>
  drop_na(`rma.glmm_IRR`, `glmmTMB.poisson`) |>
  ggplot(aes(x = `rma.glmm_IRR`, y = glmmTMB.poisson, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  labs(
    x = TeX("$\\hat{\\mu}$ rma.glmm[IRR]"),
    y = TeX("$\\hat{\\mu}$ glmmTMB"),
    title = "Poisson–Normal models agreement"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# mu SE plot - package agreement
mu_se_pois_plot <- res |> 
  filter(model %in% c("glmmTMB.poisson", "rma.glmm_IRR")) |>
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, se) |>
  pivot_wider(names_from = model, values_from = se) |>
  drop_na(`rma.glmm_IRR`, `glmmTMB.poisson`) |>
  ggplot(aes(x = `rma.glmm_IRR`, y = glmmTMB.poisson, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  labs(
    x = TeX("$\\hat{\\mu}$ standard error rma.glmm[IRR]"),
    y = TeX("$\\hat{\\mu}$ standard error glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# tau2 estimate plot - package agreement
tau2_hat_pois_plot <- res |> 
  filter(model %in% c("glmmTMB.poisson", "rma.glmm_IRR")) |>
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, tau2) |>
  pivot_wider(names_from = model, values_from = tau2) |>
  drop_na(`rma.glmm_IRR`, `glmmTMB.poisson`) |>
  ggplot(aes(x = `rma.glmm_IRR`, y = glmmTMB.poisson, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  labs(
    x = TeX("$\\hat{\\tau^2}$ rma.glmm[IRR]"),
    y = TeX("$\\hat{\\tau^2}$ glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


### patchwork for figure of estimates
plot_est_poisson <- mu_hat_pois_plot + mu_se_pois_plot + tau2_hat_pois_plot + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))


ggsave(filename = "results/figures/Figure_poisson_converged.pdf", width = 5, height = 10)




################ 6. Plots - converged results differences ###############

#### NORMAL-NORMAL
# (1) mu estimate dataset - package differences
plot_mu_diff <- ggplot(
  mu_hat |> mutate(diff = `rma.uni` - glmmTMB_RE),
  aes(x = `rma.uni`, y = diff, colour = measure)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  #scale_y_continuous(limits=c(-0.0001, 0.0001))+
  scale_colour_manual(values = col.m) +
  facet_wrap(~ measure, scales = "free_x", ncol = 4) +
  labs(
    x = TeX("$\\hat{\\mu}$ rma.uni"),
    y = TeX("$\\hat{\\mu}$ difference (rma $-$ glmmTMB)"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# (2) mu SE dataset - package difference
plot_mu_se_diff <- ggplot(
  mu_se_hat |> mutate(diff = `rma.uni` - glmmTMB_RE),
  aes(x = `rma.uni`, y = diff, colour = measure)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  #scale_y_continuous(limits=c(-0.0001, 0.0001))+
  scale_colour_manual(values = col.m) +
  facet_wrap(~ measure, scales = "free_x", ncol = 4) +
  labs(
    x = TeX("$\\hat{\\mu}$ SE rma.uni"),
    y = TeX("$\\hat{\\mu}$ SE difference (rma $-$ glmmTMB)"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# (3) tau2 estimate dataset - package difference
plot_tau2_diff <- ggplot(
  tau2_hat |> mutate(diff = `rma.uni` - glmmTMB_RE),
  aes(x = `rma.uni`, y = diff, colour = measure)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  #scale_y_continuous(limits=c(-0.0001, 0.0001))+
  scale_colour_manual(values = col.m) +
  facet_wrap(~ measure, scales = "free_x", ncol = 4) +
  labs(
    x = TeX("$\\hat{\\tau^2}$ rma.uni"),
    y = TeX("$\\hat{\\tau^2}$ difference (rma $-$ glmmTMB)"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

### patchwork for figure of estimates difference format
plot_est_diff <- plot_mu_diff + plot_mu_se_diff + plot_tau2_diff + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_NormalNormal_converged_diff.png", width = 11, height = 9)





#### BINOMIAL-NORMAL
# mu estimate plot - package difference
mu_hat_binom_plot <- res |> 
  filter(model %in% c("glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("OR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, est) |>
  pivot_wider(names_from = model, values_from = est) |>
  drop_na(rma.glmm_OR, `glmmTMB.binomial.1`, `glmmTMB.binomial.2`) |> 
  #stack the two glmmTMB models for faceting
  pivot_longer(cols = c(`glmmTMB.binomial.1`, `glmmTMB.binomial.2`),
               names_to = "glmm_model",
               values_to = "est_glmm") |> 
  ggplot(aes(x = rma.glmm_OR, y = est_glmm, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ glmm_model, ncol = 2, scales = "free") +
  labs(
    x = TeX("$\\hat{\\mu}$ rma.glmm[OR]"),
    y = TeX("$\\hat{\\mu}$ glmmTMB"),
    title = "Binomial-Normal models agreement"
  ) +
  theme_bw() +
  theme(legend.position = "none")


# mu SE plot - package difference
mu_se_binom_plot <- res |> 
  filter(model %in% c("glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("OR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, se) |>
  pivot_wider(names_from = model, values_from = se) |>
  drop_na(rma.glmm_OR, `glmmTMB.binomial.1`, `glmmTMB.binomial.2`) |> 
  #stack the two glmmTMB models for faceting
  pivot_longer(cols = c(`glmmTMB.binomial.1`, `glmmTMB.binomial.2`),
               names_to = "glmm_model",
               values_to = "est_glmm") |> 
  ggplot(aes(x = rma.glmm_OR, y = est_glmm, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ glmm_model, ncol = 2, scales = "free") +
  labs(
    x = TeX("$\\hat{\\mu}$ standard error rma.glmm[OR]"),
    y = TeX("$\\hat{\\mu}$ standard error glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# tau2 estimate plot - package aggreement
tau2_hat_binom_plot <- res |> 
  filter(model %in% c("glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR")) |> 
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("OR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, tau2) |>
  pivot_wider(names_from = model, values_from = tau2) |>
  drop_na(rma.glmm_OR, `glmmTMB.binomial.1`, `glmmTMB.binomial.2`) |> 
  #stack the two glmmTMB models for faceting
  pivot_longer(cols = c(`glmmTMB.binomial.1`, `glmmTMB.binomial.2`),
               names_to = "glmm_model",
               values_to = "est_glmm") |> 
  ggplot(aes(x = rma.glmm_OR, y = est_glmm, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  facet_wrap(~ glmm_model, ncol = 2, scales = "free") +
  labs(
    x = TeX("$\\hat{\\tau^2}$ rma.glmm[OR]"),
    y = TeX("$\\hat{\\tau^2}$ glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

### patchwork for figure of estimates
plot_est_binom <- mu_hat_binom_plot + mu_se_binom_plot + tau2_hat_binom_plot + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_binom_converged.pdf", width = 7, height = 10)


#### POISSON-NORMAL
# mu estimate plot - package agreement
mu_hat_pois_plot <- res |> 
  filter(model %in% c("glmmTMB.poisson", "rma.glmm_IRR")) |>
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, est) |>
  pivot_wider(names_from = model, values_from = est) |>
  drop_na(`rma.glmm_IRR`, `glmmTMB.poisson`) |>
  ggplot(aes(x = `rma.glmm_IRR`, y = glmmTMB.poisson, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  labs(
    x = TeX("$\\hat{\\mu}$ rma.glmm[IRR]"),
    y = TeX("$\\hat{\\mu}$ glmmTMB"),
    title = "Poisson–Normal models agreement"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# mu SE plot - package agreement
mu_se_pois_plot <- res |> 
  filter(model %in% c("glmmTMB.poisson", "rma.glmm_IRR")) |>
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, se) |>
  pivot_wider(names_from = model, values_from = se) |>
  drop_na(`rma.glmm_IRR`, `glmmTMB.poisson`) |>
  ggplot(aes(x = `rma.glmm_IRR`, y = glmmTMB.poisson, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  labs(
    x = TeX("$\\hat{\\mu}$ standard error rma.glmm[IRR]"),
    y = TeX("$\\hat{\\mu}$ standard error glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


# tau2 estimate plot - package agreement
tau2_hat_pois_plot <- res |> 
  filter(model %in% c("glmmTMB.poisson", "rma.glmm_IRR")) |>
  group_by(measure, model) |> 
  mutate(.idx = row_number(), measure = factor(measure, levels = c("IRR"))) |> 
  ungroup() |> 
  select(measure, model, .idx, tau2) |>
  pivot_wider(names_from = model, values_from = tau2) |>
  drop_na(`rma.glmm_IRR`, `glmmTMB.poisson`) |>
  ggplot(aes(x = `rma.glmm_IRR`, y = glmmTMB.poisson, colour = measure)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_colour_manual(values = col.m) +
  labs(
    x = TeX("$\\hat{\\tau^2}$ rma.glmm[IRR]"),
    y = TeX("$\\hat{\\tau^2}$ glmmTMB"),
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")


### patchwork for figure of estimates
plot_est_poisson <- mu_hat_pois_plot + mu_se_pois_plot + tau2_hat_pois_plot + 
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))


ggsave(filename = "results/figures/Figure_poisson_converged.pdf", width = 5, height = 10)





################ 7. Plots - RMSE for mu and tau2 ###############

######### mu RMSE ---
# plot for SMD and lnRR
p_rmsea <- rmse.dat.mu |>
  filter(measure %in% c("SMD", "lnRR"), model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = rmse, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = TeX("$\\hat{\\mu}$ RMSE"), title="lnRR and SMD") +
  theme(legend.position = "none")


# plot for OR
p_rmseb <- rmse.dat.mu |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = rmse, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = TeX("$\\hat{\\mu}$ RMSE"), title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_rmsec <- rmse.dat.mu |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = rmse, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = TeX("$\\hat{\\mu}$ RMSE"), title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for runtime figure
plot_rmse <- p_rmsea + p_rmseb + p_rmsec +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_RMSE_mu.pdf", width = 5, height = 11)
ggsave(filename = "results/figures/Figure_RMSE_mu.png", width = 5, height = 11)



######### tau2 RMSE ---
# plot for SMD and lnRR
p_rmsea <- rmse.dat.tau2 |>
  filter(measure %in% c("SMD", "lnRR"), model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = rmse, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = TeX("$\\hat{\\tau^2}$ RMSE"), title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for OR
p_rmseb <- rmse.dat.tau2 |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = rmse, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = TeX("$\\hat{\\tau^2}$ RMSE"), title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_rmsec <- rmse.dat.tau2 |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = rmse, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = TeX("$\\hat{\\tau^2}$ RMSE"), title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for runtime figure
plot_rmse <- p_rmsea + p_rmseb + p_rmsec +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_RMSE_tau2.pdf", width = 5, height = 11)
ggsave(filename = "results/figures/Figure_RMSE_tau2.png", width = 5, height = 11)














################ 8. Plots - inference overall mean ###############

#### 95% Coverate rate --------
# plot for SMD and lnRR
p_cova <- cov.dat |>
  filter(measure %in% c("SMD", "lnRR"), model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = cov_prop, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = "Coverage", title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for OR
p_covb <- cov.dat |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = cov_prop, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "Coverage", title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_covc <- cov.dat |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = cov_prop, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "Coverage", title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for runtime figure
plot_cov <- p_cova + p_covb + p_covc +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_coverage_mu.pdf", width = 5, height = 11)
ggsave(filename = "results/figures/Figure_coverage_mu.png", width = 5, height = 11)


#### 95% CI widths --------
# plot for SMD and lnRR
p_cia <- ci.dat |>
  filter(measure %in% c("SMD", "lnRR"), model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = ci_width, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = "CI width (95%)", title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for OR
p_cib <- ci.dat |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = ci_width, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "CI width (95%)", title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_cic <- ci.dat |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = ci_width, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "CI width (95%)", title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for runtime figure
plot_ci <- p_cia + p_cib + p_cic +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_ci_width_mu.pdf", width = 5, height = 11)
ggsave(filename = "results/figures/Figure_ci_width_mu.png", width = 5, height = 11)




#### Type I error rate -----
# plot for SMD and lnRR
p_typeIa <- typeI.dat |>
  filter(measure %in% c("SMD", "lnRR"), model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = typeI, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = "Type I error", title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for OR
p_typeIb <- typeI.dat |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = typeI, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "Type I error", title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_typeIc <- typeI.dat |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = typeI, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "Type I error", title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for figure
plot_typeI <- p_typeIa + p_typeIb + p_typeIc +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_typeIerror.pdf", width = 5, height = 11)
ggsave(filename = "results/figures/Figure_typeIerror.png", width = 5, height = 11)




#### Power rate -------
# plot for SMD and lnRR
p_powera <- power.dat |>
  filter(measure %in% c("SMD", "lnRR"), model %in% c("rma.uni", "glmmTMB_RE")) |> 
  ggplot(aes(x = model, y = power, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  facet_wrap(~ measure, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(y = "Power", title="lnRR and SMD") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for OR
p_powerb <- power.dat |>
  filter(measure %in% c("OR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.binomial.1",  "glmmTMB.binomial.2", "rma.glmm_OR"))) |> 
  ggplot(aes(x = model, y = power, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "Power", title="OR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# plot for IRR
p_powerc <- power.dat |>
  filter(measure %in% c("IRR")) |>
  mutate(model = factor(model, levels = c("glmmTMB_RE", "rma.uni",
                                          "glmmTMB.poisson", "rma.glmm_IRR"))) |> 
  ggplot(aes(x = model, y = power, colour = measure, fill = measure)) +
  #coord_flip() +
  stat_halfeye() +
  geom_boxplot(width = 0.4) +
  scale_colour_manual(values = col.m) +
  scale_fill_manual(values = col.m.pastel) +
  theme_bw() +
  labs(y = "Power", title="IRR") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

### pathwork for figure
plot_power <- p_powera + p_powerb + p_powerc +
  plot_layout(ncol=1, guides = "collect") +
  plot_annotation(title = "", theme = theme(plot.background = element_rect(fill = "white", colour = NA)))

ggsave(filename = "results/figures/Figure_power.pdf", width = 5, height = 11)
ggsave(filename = "results/figures/Figure_power.png", width = 5, height = 11)

