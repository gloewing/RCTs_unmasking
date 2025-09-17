# mediation analysis vignette
library(lmtp)
# for more information see: 
#     1) https://muse.jhu.edu/pub/56/article/883479/pdf
#     2) https://beyondtheate.com/


expit <- function(x){ exp(x) / (1 + exp(x))}
CDE <- 3
n <- 100 # sample size

set.seed(1)

# treatment
trt.probability <- 0.5 # 50/50 trt/control split
A <- rbinom(n = n, # sample size
            size = 1, # binary treatment
            prob = trt.probability) # psilocybin vs. control

# baseline covariates
male_female <- rbinom(n = n, size = 1, prob = 0.5) # biological sex
age <- sample(20:50, size = n, replace = TRUE) / 50 # standardized age
psych_history <- runif(n) # psychedelic history PRIOR to study in [0,1]
id <- 1:n

#--------------------------
# post-treatment variables
#--------------------------
set.setting <- A * abs(rnorm(n)) + rbinom(n, size = 1, prob = 0.5) # set/setting variable randomly drawn

# belief
beta.B <- rep(1, 3) #rnorm(3) # coefficients for belief model
belief.prob <- expit( beta.B[1] * A * set.setting + 
                        beta.B[2] * A + 
                        beta.B[3] * set.setting) # arbitrary belief model
B <- sapply(belief.prob, function(x) rbinom(1, 1, prob = x)) # draw belief as binary

# combine variables as all causes of E
E.causes <- as.matrix(data.frame(A = A,
                                 B = B,
                                 SS = set.setting,
                                 mf = male_female,
                                 age = age,
                                 psych_h = psych_history)) # combine into matrix

# expectancy
# intercept.E <- -3 # arbitary model intercept
beta.E <- rep(1, ncol(E.causes)) #rnorm(ncol(E.causes)) # draw expectancy model coefficients
exp_prob <- expit(E.causes %*% beta.E ) # probability of expectancy for each individual
E <- sapply(exp_prob, function(x) rbinom(1, 1, prob = x)) # draw expectancy as binary

# outcome model
Y.mat <- as.matrix(cbind(E.causes, E))
beta.Y <- rnorm(ncol(Y.mat))
CDE.idx <- which(colnames(Y.mat) == "A")
beta.Y[CDE.idx] <- CDE # set CDE
Y <- Y.mat %*% beta.Y + rnorm(n, mean = 0, sd = 0.25)

df <- data.frame(A = A,
                 E = E, 
                 B = B,
                 SS = set.setting,
                 mf = male_female,
                 age = age,
                 psych_h = psych_history, 
                 Y = Y,
                 dummy = 1, # necessary trick for first timepoint confounders with randomization of A
                 id = id)

# package
# no baseline confounders because the first timepoint's "treatment" is treatment/control

# confounders
# Z <- c("B", "SS", "mf", "age", "psych_h") # confounders of E -> Y path
# L <- list(NULL, Z) # no time-varying confounders for first timepoint so first element of L is empty
# A <- c("A", "E")
# lmtp::create_node_list(trt = A, baseline = W, time_vary = L, tau = 2)

# static treatment policy for multivariate exposure
# info: https://beyondtheate.com/10_R_multivariate.html#multivariate-shift-functions

###############################
# set four regimes of interest
###############################
# A == 1, E == 0
d_A1_E0 <- function(data, trt) {
  rep( ifelse(trt == "A", 1, 0),
       nrow(data)
       )
}

# A == 0, E == 1
d_A0_E1 <- function(data, trt) {
  rep( ifelse(trt == "E", 1, 0),
       nrow(data)
  )
  }

# A == 1, E == 1
d_A1_E1 <- function(data, trt) {
  return( rep(1, nrow(data)) ) 
}

# A == 0, E == 0
d_A0_E0 <- function(data, trt) {
  return( rep(0, nrow(data)) ) 
}


# --------------------------
# Fit dataset
# --------------------------
# for available learners
# https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html

# fit dataset
# fit_sdr <- lmtp_sdr(
#   data = df, 
#   trt = c("A"), # second timepoint's "treatment"/exposure is expectancy, E 
#   outcome = "Y", 
#   outcome_type = "continuous",
#   baseline = NULL, # no time-varying confounders: Baseline confounders are always included in Ht according to documentation
#   time_vary = list(c("B", "SS", "mf", "age", "psych_h")), # confounders of E -> Y path
#   k = Inf, # whether to use all previous history variables pooled/concatenated
#   shift = d1,
#   id = "id",
#   mtp = FALSE,
#   folds = 2,
#   learners_trt = c("SL.glm", "SL.xgboost"), 
#   learners_outcome = c("SL.xgboost", "SL.glm")
# )

# should work but doesnt
# fit_sdr <- lmtp_sdr(
#   data = df,
#   trt = c("A", "E"), # second timepoint's "treatment"/exposure is expectancy, E
#   outcome = "Y",
#   outcome_type = "continuous",
#   baseline = NULL, # no time-varying confounders: Baseline confounders are always included in Ht according to documentation
#   time_vary = list(NULL, c("B", "SS", "psych_h", "age", "mf")), # confounders of E -> Y path
#   k = 0, # do not accumulate confounders from prior timepoints in history 
#   shift = d1,
#   id = "id",
#   mtp = FALSE,
#   folds = 2,
#   learners_trt = c("SL.glm", "SL.xgboost"),
#   learners_outcome = c("SL.xgboost", "SL.glm")
# )



X <- c("psych_h", "age", "mf") # baseline variables
Z <- c("SS", "B") # post-treatment confounders

# d_A1_E0
fit_sdr_A1E0 <- lmtp_sdr(
  data = df,
  trt = c("A", "E"), # second timepoint's "treatment"/exposure is expectancy, E
  outcome = "Y",
  outcome_type = "continuous",
  baseline = NULL, # no time-varying confounders: Baseline confounders are always included in Ht according to documentation
  time_vary = list(X, Z), # confounders of E -> Y path
  k = Inf, # whether to use all previous history variables pooled/concatenated
  shift = d_A1_E0,
  id = "id",
  mtp = FALSE,
  folds = 2,
  learners_trt = c("SL.glm", "SL.xgboost"),
  learners_outcome = c("SL.xgboost", "SL.glm")
)

# d_A0_E0
fit_sdr_A0E0 <- lmtp_sdr(
  data = df,
  trt = c("A", "E"), # second timepoint's "treatment"/exposure is expectancy, E
  outcome = "Y",
  outcome_type = "continuous",
  baseline = NULL, # no time-varying confounders: Baseline confounders are always included in Ht according to documentation
  time_vary = list(X, Z), # confounders of E -> Y path
  k = Inf, # whether to use all previous history variables pooled/concatenated
  shift = d_A0_E0,
  id = "id",
  mtp = FALSE,
  folds = 2,
  learners_trt = c("SL.glm", "SL.xgboost"),
  learners_outcome = c("SL.xgboost", "SL.glm")
)

CDE.idx <- which(colnames(Y.mat) == "A")
beta.Y[CDE.idx]
lmtp_contrast(fit_sdr_A1E0, ref = fit_sdr_A0E0, type = "additive")




# E == 1
# d_A1_E0
fit_sdr_A1E1 <- lmtp_sdr(
  data = df,
  trt = c("A", "E"), # second timepoint's "treatment"/exposure is expectancy, E
  outcome = "Y",
  outcome_type = "continuous",
  baseline = NULL, # no time-varying confounders: Baseline confounders are always included in Ht according to documentation
  time_vary = list(X, Z), # confounders of E -> Y path
  k = Inf, # whether to use all previous history variables pooled/concatenated
  shift = d_A1_E1,
  id = "id",
  mtp = FALSE,
  folds = 2,
  learners_trt = c("SL.glm", "SL.xgboost"),
  learners_outcome = c("SL.xgboost", "SL.glm")
)

# d_A0_E0
fit_sdr_A0E1 <- lmtp_sdr(
  data = df,
  trt = c("A", "E"), # second timepoint's "treatment"/exposure is expectancy, E
  outcome = "Y",
  outcome_type = "continuous",
  baseline = NULL, # no time-varying confounders: Baseline confounders are always included in Ht according to documentation
  time_vary = list(X, Z), # confounders of E -> Y path
  k = Inf, # whether to use all previous history variables pooled/concatenated
  shift = d_A0_E1,
  id = "id",
  mtp = FALSE,
  folds = 2,
  learners_trt = c("SL.glm", "SL.xgboost"),
  learners_outcome = c("SL.xgboost", "SL.glm")
)

lmtp_contrast(fit_sdr_A1E1, ref = fit_sdr_A0E1, type = "additive")

