library(MASS)
library(metafor)
library(insight) ##to get dfs 

#remove.packages("glmmTMB") #if needed (reload R afterwards)
#build package from branch: https://github.com/coraliewilliams/glmmTMB/tree/equalto_covstruc
#remotes::install_github("coraliewilliams/glmmTMB", ref="equalto_covstruc", subdir="glmmTMB")
library(glmmTMB)
# check available covstruc (should see equalto in number 13)
(glmmTMB:::.valid_covstruct)


########## Fit a multilevel meta-analysis ##############
# y ~ b0 + u + m + e
# 
# u ~ (0, s2_u.I) - study level
# m ~ (0, s2_m.I) - observation level (effect sizes)
# e ~ (0, VCV)    - residuals/errors (VCV is a block diag with effect sizes correlated from the same study)


##### Simulate some data

# set up indices
set.seed(123)
k.studies <- 100 # number of studies
k.per.study <- 10 # number of effect sizes per study
study <- rep(seq_len(k.studies), times=k.per.study) # study id 
k <- length(study) # total number of effect sizes
id <- seq_len(k) # id of between study effect sizes 
es.id <- unlist(lapply(k.per.study, seq_len)) # id of within study effect sizes


# fixed effect coeff true value
b0 <- 0

# simulate sampling errors variances
vi <- rbeta(k, 2, 20)

# get VCV matrix assuming within-study correlation of rho=0.5 
VCV <- matrix(0, nrow = k, ncol = k) 
for (i in 2:k) {
  for (j in 1:i) {
    if (study[i] == study[j]) {
      VCV[i,j] <- 0.5 * sqrt(vi[i] * vi[j]) 
    }
  }
}
VCV[upper.tri(VCV)] <- t(VCV)[upper.tri(VCV)] # fill in upper diagonal
diag(VCV) <- vi # fill in diagonal


# random effects
sigma2.u <- 0.2
sigma2.m <- 0.3
u <- rnorm(k.studies, 0, sqrt(sigma2.u))[study]
m <- rnorm(k, 0, sqrt(sigma2.m))
#e <- rnorm(k, 0, sqrt(vi))[study] #without within-study correlation (diag VCV)
e <- mvrnorm(n = 1, mu = rep(0, length(vi)), Sigma = VCV) #excluding within-study correlation


# compute y
y <- b0 + u + m + e


# combine into dataframe
dat <- data.frame(y = y,
                  vi = vi,
                  study = study,
                  id = id,
                  es.id = es.id,
                  obs = as.factor(id), # unique ID for propto and equalto (really important it is a factor variable)
                  g = 1 # group ID for propto and equalto
)

##############################################

#Fit metafor model -------------------------------------------------------
ptm <- proc.time()
fit.metafor <- rma.mv(y, VCV,
                      random = list(~1 | study, ~1 | id),
                      test = "t",
                      dfs = "contain",
                      data=dat) #used REML by default
metafor_time <- (proc.time() - ptm)[3]

fit.metafor$sigma2[1] # u (study) 
fit.metafor$sigma2[2] # m (effect size- within study)


#Fit glmmTMB model with equalto --------------------------------------------
ptm <- proc.time()
fit.equalto <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + obs|g,VCV),
                       data=dat,
                       control = glmmTMBControl(parallel = 8),
                       REML=T)
equalto_time <- (proc.time() - ptm)[3]

(exp(fit.equalto$fit$par[[2]]))^2 # u (study) (late)
sigma(fit.equalto)^2 # m (effect size- within study)
# e (sampling error) **this should be zero 



##### Compare output 
metafor <- data.frame(model = "metafor", 
                      time = as.numeric(metafor_time),
                      est = fit.metafor$b[[1]], 
                      se = fit.metafor$se[[1]], 
                      zval = fit.metafor$zval,
                      sigma.u = sqrt(fit.metafor$sigma2[1]),
                      sigma.s = sqrt(fit.metafor$sigma2[2]),
                      sigma.total = sum(sqrt((fit.metafor$sigma2))),
                      logLik = logLik(fit.metafor)[1],
                      df = attr(logLik(fit.metafor), "df"), ## model-based 
                      df.res = df.residual(fit.metafor)) ## residual df

equalto <- data.frame(model = "equalto",
                      time = as.numeric(equalto_time),
                      est = unlist(fixef(fit.equalto))[[1]],
                      se = as.numeric(sqrt(vcov(fit.equalto)[[1]])),
                      zval = summary(fit.equalto)$coefficients$cond[3],
                      sigma.u = exp(fit.equalto$fit$par[[2]]),
                      sigma.s = sigma(fit.equalto),
                      sigma.total = sum(exp(fit.equalto$fit$par)),
                      logLik = logLik(fit.equalto),
                      df = attr(logLik(fit.equalto), "df"), ## model-based df
                      df.res = df.residual(fit.equalto)) ## residual df


output <- rbind(metafor, equalto)
output




#--------------------------------------------------------------------------------------
# Use microbenchmark to compare 
library(microbenchmark)

# run 100 times
bench_fit <- microbenchmark(
  metafor = rma.mv(y, VCV, random = list(~1 | study, ~1 | id),
                   test = "t",
                   dfs = "contain",
                   data=dat),
  glmmTMB = glmmTMB(y ~ 1 + (1|study) + equalto(0 + obs|g,VCV),
                    data=dat,
                    REML=T)
)

bench_fit
plot(bench_fit)



#-------------------------------------------------------------------------------------
#### test speed (should be under 0.1 if blas is installed)
system.time(A<-solve(diag(1000+.1)))


#-------------------------------------------------------------------------------------
###### Profiling equalto (to see which function is taking the most time)
Rprof(tmp <- tempfile())
fit.equalto <- glmmTMB(y ~ 1 + (1|study) + equalto(0 + obs|g,VCV),
                       data=dat,
                       REML=T)
Rprof()
summaryRprof(tmp)
unlink(tmp)




#Fit glmmTMB model with propto ---------------------------------------------
colnames(VCV) <- 1:nrow(dat)
row.names(VCV) <- 1:nrow(dat)
n.thetas <- nrow(dat) + (nrow(dat) * (nrow(dat)-1)) / 2 + 1 + 1 
# number of thetas = (log-sds) + (corr parameters) + 1 (log-lambda) +  1 (random effect)

fit.propto <- glmmTMB(y ~ 1 + (1|study) + propto(0 + obs|g,VCV),
                      #start = list(theta=c(rep(0, n.thetas-1), #default of start is 0 for non-fixed values
                      #                     rep(log(1), 1))),   #set lambda (prop) parameter to 1 => log(1/2)=1
                      #map = list(theta=factor(c(1:(n.thetas-1), NA))), #fix the last theta parameter
                      data=dat,
                      REML=T)

(exp(fit.propto$fit$par[[2]]))^2 # u (study) (late)
sigma(fit.propto)^2 # m (effect size- within study)
(exp(fit.propto$fit$par[[3]]))^2 # e (sampling error) **this should be zero **where propto inputs vi information; 

# get total variance
sum((exp(fit.propto$fit$par))^2)

# check glmmTMB convergence
fit.propto$sdr$pdHess
diagnose(fit.propto)

