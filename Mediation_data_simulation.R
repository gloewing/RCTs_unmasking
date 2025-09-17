# mediation analysis vignette data simulation

CDE <- 3 # causal effect size
n <- 250 # sample size # 250
set.seed(123) # random seed # 123
save.wd <- "~/Desktop/NIMH Research/Causal/Causal_Psychedelics/Paper/Mediation_Vignette"
  
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
expit <- function(x){ exp(x) / (1 + exp(x))}
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
                 id = id)

setwd(save.wd)
write.csv(df, "CDE_data.csv", row.names = FALSE)
