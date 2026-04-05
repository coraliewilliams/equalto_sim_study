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



# -------------------- Random effect models
fit.rm <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g,V),
                       dispformula=~0,
                       data=dat,
                       REML=TRUE)


fit.rm2 <- glmmTMB(yi ~ 1 + (1|study),
                       dispformula= ~0 + id,
                       map = list(betadisp = factor(rep(NA, nrow(dat)))),
                       start = list(betadisp = log(sqrt(dat$vi))),
                       data=dat,
                       REML=TRUE)





### construct variance-covariance matrix 
### assuming rho = 0.6 for effect sizes within studies
VCV <- vcalc(vi=vi, cluster=study, obs=id, data=dat, rho=0.6)
round(VCV[1:5,1:5], 3)
#       [,1] [,2]  [,3]  [,4]  [,5]
# [1,] 0.074 0.00 0.000 0.000 0.000
# [2,] 0.000 0.04 0.000 0.000 0.000
# [3,] 0.000 0.00 0.048 0.000 0.000
# [4,] 0.000 0.00 0.000 0.024 0.000
# [5,] 0.000 0.00 0.000 0.000 0.033 




# -------------------- Multilevel models
fit.ml.tmb <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g,VCV) + (1|id),
                      dispformula=~0,
                      data=dat,
                      REML=TRUE)

fit.ml.tmb <- glmmTMB(yi ~ 1 + (1|study) + equalto(0 + id|g,VCV),
                      data=dat,
                      REML=TRUE)


# equivalent ml in metafor
fit.rma.mv <- rma.mv(yi, VCV, 
                     random = list(~1|study, ~1|id),
                     data=dat)


summary(fit.ml.tmb)

# Family: gaussian  ( identity )
# Formula:          yi ~ 1 + (1 | study) + equalto(0 + id | g, VCV)
# Data: dat
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 156.1     164.0     -75.1     150.1        97 
# 
# Random effects:
#   
# Conditional model:
#   Groups   Name        Variance Std.Dev. Corr                                            
#   study    (Intercept) 0.08073  0.28414                                                  
#   g        id1         0.07400  0.27203                                                  
#            id2         0.03980  0.19950  0.60                                            
#            id3         0.04810  0.21932  0.60 0.60                                       
#            id4         0.02390  0.15460  0.60 0.60 0.60                                  
#            id5         0.03310  0.18193  0.60 0.60 0.60 0.60                             
#            id6         0.08860  0.29766  0.60 0.60 0.60 0.60 0.60                        
#            id7         0.01150  0.10724  0.00 0.00 0.00 0.00 0.00 0.00                   
#            id8         0.00760  0.08718  0.00 0.00 0.00 0.00 0.00 0.00 0.60              
#            id9         0.00650  0.08062  0.00 0.00 0.00 0.00 0.00 0.00 0.60 0.60         
#            id10        0.33250  0.57663  0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00    
#            ...         ...      ...      ...  ...  ...  ...  ...  ...  ...  ...  ...  ...
#  Residual             0.15454  0.39312                                                  
# Number of obs: 100, groups:  study, 17; g, 1
# 
# Dispersion estimate for gaussian family (sigma^2): 0.155 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.36775    0.09726   3.781 0.000156 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   


summary(fit.rma.mv)

# Multivariate Meta-Analysis Model (k = 100; method: REML)
# 
#   logLik  Deviance       AIC       BIC      AICc   
# -72.7667  145.5334  151.5334  159.3188  151.7861   
# 
# Variance Components:
#   
#             estim    sqrt  nlvls  fixed  factor 
# sigma^2.1  0.0807  0.2841     17     no   study 
# sigma^2.2  0.1545  0.3931    100     no      id 
# 
# Test for Heterogeneity:
#   Q(df = 99) = 745.4385, p-val < .0001
# 
# Model Results:
#   
#   estimate      se    zval    pval   ci.lb   ci.ub      
# 0.3678  0.0965  3.8097  0.0001  0.1786  0.5570  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1