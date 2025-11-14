library(glmmTMB)
library(metafor)
dat <- dat.assink2016
dat$id <- as.factor(dat$id)
dat$g <- rep(1, nrow(dat))
head(dat, 9)
# study esid id      yi     vi pubstatus year deltype g
#     1    1  1  0.9066 0.0740         1  4.5 general 1
#     1    2  2  0.4295 0.0398         1  4.5 general 1
#     1    3  3  0.2679 0.0481         1  4.5 general 1
#     1    4  4  0.2078 0.0239         1  4.5 general 1
#     1    5  5  0.0526 0.0331         1  4.5 general 1
#     1    6  6 -0.0507 0.0886         1  4.5 general 1
#     2    1  7  0.5117 0.0115         1  1.5 general 1
#     2    2  8  0.4738 0.0076         1  1.5 general 1
#     2    3  9  0.3544 0.0065         1  1.5 general 1


V <- diag(dat$vi)
round(V[1:5,1:5], 3)
#       [,1] [,2]  [,3]  [,4]  [,5]
# [1,] 0.074 0.00 0.000 0.000 0.000
# [2,] 0.000 0.04 0.000 0.000 0.000
# [3,] 0.000 0.00 0.048 0.000 0.000
# [4,] 0.000 0.00 0.000 0.024 0.000
# [5,] 0.000 0.00 0.000 0.000 0.033


fit.rma.tmb <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g,V),
                       dispformula=~0,
                       data=dat,
                       REML=TRUE)
### construct variance-covariance matrix 
### assuming rho = 0.6 for effect sizes within studies
VCV <- vcalc(vi=vi, cluster=study, obs=id, data=dat, rho=0.6)
round(VCV[1:5,1:5], 3)
#       [,1]  [,2]  [,3]  [,4]  [,5]
# [1,] 0.074 0.033 0.036 0.025 0.030
# [2,] 0.033 0.040 0.026 0.019 0.022
# [3,] 0.036 0.026 0.048 0.020 0.024
# [4,] 0.025 0.019 0.020 0.024 0.017
# [5,] 0.030 0.022 0.024 0.017 0.033

fit.lm.tmb <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g,VCV),
                       data=dat,
                       REML=TRUE)


# equivalent ml in metafor
fit.rma.mv <- rma.mv(yi, VCV, 
                     random = list(~1|study, ~1|id),
                     data=dat)


summary(fit.lm.tmb)
summary(fit.rma.mv)




# Family: gaussian  ( identity )
# Formula:          yi ~ 1 + (1 | study) + equalto(0 + id | g, VCV)
# Dispersion:          ~0
# Data: dat
# 
#       AIC      BIC   logLik deviance df.resid 
#     237.7    242.9   -116.9    233.7       98 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev. Corr                                                                                                                                 
#   study  (Intercept) 0.1925   0.43875                                            
#   g      id1         0.0740   0.27203                                 
#          id2         0.0398   0.19950  0.00                            
#          id3         0.0481   0.21932  0.00 0.00                       
#          id4         0.0239   0.15460  0.00 0.00 0.00                   
#          id5         0.0331   0.18193  0.00 0.00 0.00 0.00          
#          ..            ..        ..     ..   ..   ..   ..
# [ reached getOption("max.print") -- omitted 95 rows ]
# Number of obs: 100, groups:  study, 17; g, 1
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.4032     0.1111   3.629 0.000285 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1