TT = 50
N = 50
grid.star <- c(-2/5, 2/4, 3/4) 
time_cutoff <- 5 # cut off timepoints before this

div.factor <- 2
p0 <- 0.5

###################################################
# actual Monte-Carlo Simulation
sim.subject <- function(TT, g1, g0, p0, risk0, risk_slope, seed, arm){
  set.seed(seed)
  Y <- X <- Y.culm <- X.culm <- vector(length = TT)
  #rand.eff <- runif(1, -0.1, 0.1)
  
  # timepoint 1
  X[1] <- X.culm[1] <- rbinom(1, 1, prob = p0) # risk variable
  Y[1] <- Y.culm[1] <- rbinom(1, size = 1, 
                              prob = (risk0 - arm*risk_slope) * X[1]) # outcome
  
  for(t in 2:TT){
    X[t] <- rbinom(1, size = 1, 
                   prob = g0 + g1 * max(Y[ 1:(t-1) ]) )# risk
    Y[t] <- rbinom(1, size = 1, 
                   prob = (risk0 - arm*risk_slope) * X[t] ) # outcome
    Y.culm[t] <- mean(Y[1:t])
    X.culm[t] <- mean(X[1:t])
  }
  
  return( data.frame(Y = Y, 
                     X = X, 
                     Y.culm = Y.culm,
                     X.culm = X.culm,
                     participant = seed, 
                     arm = arm,
                     time = 1:TT) 
  )
  
}

# control group
control_arm <- 0
control.group <- lapply(1:N, 
                        function(i) sim.subject(TT = TT, 
                                                g1 = grid.star[1], 
                                                g0 = grid.star[2], 
                                                p0 = grid.star[2],
                                                risk0=grid.star[3], 
                                                risk_slope=grid.star[3]/div.factor, 
                                                seed = i, 
                                                arm = control_arm) )

control.group <- do.call(rbind, control.group)

# treatment group
treatment_arm <- 1
treatment.group <- lapply(1:N, 
                          function(i) sim.subject(TT = TT, 
                                                  g1 = grid.star[1], 
                                                  g0 = grid.star[2], 
                                                  p0 = grid.star[2],
                                                  risk0=grid.star[3], 
                                                  risk_slope=grid.star[3]/div.factor, 
                                                  seed = i, 
                                                  arm = treatment_arm) )

treatment.group <- do.call(rbind, treatment.group)

# join datasets
dat <- rbind(treatment.group, control.group)
# ensure particpants have unique IDs
dat$participant[dat$arm == 0] <- dat$participant[dat$arm == 0] + N


# figures
myColors2 <- c("red", "blue", "#E69F00", "darkgreen", "darkgrey")

plot1 <- dat %>% as_tibble() %>% 
  dplyr::mutate(Arm = ifelse(arm == 1, "Treat", "Control"),
                id = as.factor(participant),
                trial = time) %>%
  group_by(Arm, trial) %>%
  summarise(Y_mean = mean(Y.culm),
            se = Y_mean*(1-Y_mean) / sqrt(n()),
            upper = Y_mean + 1.96*se,
            lower = Y_mean - 1.96*se) %>%
  dplyr::filter(trial > time_cutoff) %>%
  ggplot(aes( y = Y_mean, x = trial, color = Arm)) +
  geom_ribbon(aes(x = trial, ymax = upper, ymin = lower, color = Arm, fill = Arm),
              alpha = 0.4) +
  geom_line(aes(color = Arm),
            alpha = 1, linewidth = 1) +
  ylab(TeX('Culmulative Outcome Value:  $\\frac{1}{t} \\sum_{j=1}^t Y_j$')) + # paste0('Culumlative Outcome Value: ',  
  xlab(paste0('Timepoint: ', TeX('$t$'))) + 
  scale_fill_manual(values = myColors2) +
  scale_color_manual(values = myColors2) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(0.75), face="bold"),
         axis.text=element_text(color="black", size=rel(1)),
         axis.title = element_text(color="black", size=rel(1.25)),
         legend.key.size = unit(1.5, "line"), # added in to increase size
         legend.text = element_text(color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(color="black", size = rel(1))) +
  labs(title = "Efficacy-Expectancy Feedback") +
  coord_cartesian(ylim = c(0, 0.5))



############################################
# no feedback
grid.star[1] <- 0 # no feedback
# control group
control_arm <- 0
control.group <- lapply(1:N, 
                        function(i) sim.subject(TT = TT, 
                                                g1 = grid.star[1], 
                                                g0 = grid.star[2], 
                                                p0 = grid.star[2],
                                                risk0=grid.star[3], 
                                                risk_slope=grid.star[3]/div.factor, 
                                                seed = i, 
                                                arm = control_arm) )

control.group <- do.call(rbind, control.group)

# treatment group
treatment_arm <- 1
treatment.group <- lapply(1:N, 
                          function(i) sim.subject(TT = TT, 
                                                  g1 = grid.star[1], 
                                                  g0 = grid.star[2], 
                                                  p0 = grid.star[2],
                                                  risk0=grid.star[3], 
                                                  risk_slope=grid.star[3]/div.factor, 
                                                  seed = i, 
                                                  arm = treatment_arm) )

treatment.group <- do.call(rbind, treatment.group)

# join datasets
dat <- rbind(treatment.group, control.group)
# ensure particpants have unique IDs
dat$participant[dat$arm == 0] <- dat$participant[dat$arm == 0] + N


# figures
myColors2 <- c("red", "blue", "#E69F00", "darkgreen", "darkgrey")

plot0 <- dat %>% as_tibble() %>% 
  dplyr::mutate(Arm = ifelse(arm == 1, "Treat", "Control"),
                id = as.factor(participant),
                trial = time) %>%
  group_by(Arm, trial) %>%
  summarise(Y_mean = mean(Y.culm),
            se = Y_mean*(1-Y_mean) / sqrt(n()),
            upper = Y_mean + 1.96*se,
            lower = Y_mean - 1.96*se) %>%
  dplyr::filter(trial > time_cutoff) %>%
  ggplot(aes( y = Y_mean, x = trial, color = Arm)) +
  geom_ribbon(aes(x = trial, ymax = upper, ymin = lower, color = Arm, fill = Arm),
              alpha = 0.4) +
  geom_line(aes(color = Arm),
            alpha = 1, linewidth = 1) +
  ylab(TeX('Culmulative Outcome Value:  $\\frac{1}{t} \\sum_{j=1}^t Y_j$')) + # paste0('Culumlative Outcome Value: ',  
  xlab(paste0('Timepoint: ', TeX('$t$'))) + 
  scale_fill_manual(values = myColors2) +
  scale_color_manual(values = myColors2) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(0.75), face="bold"),
         axis.text=element_text(color="black", size=rel(1)),
         axis.title = element_text(color="black", size=rel(1.25)),
         legend.key.size = unit(1.5, "line"), # added in to increase size
         legend.text = element_text(color="black", size = rel(1)), 
         legend.title = element_text(face="bold", color="black", size = rel(1)),
         strip.text = element_text(color="black", size = rel(1))) +
  labs(title = "No Feedback") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  coord_cartesian(ylim = c(0, 0.5))
  

library(patchwork)

plots <- plot1 + 
            plot0 + 
            plot_layout(guides = "collect") & theme(legend.position = "bottom")
setwd("/Users/loewingergc/Desktop/NIMH Research/Causal/Causal_Psychedelics/Paper/Figures/Drafts")
ggsave( plot = plots,
        paste0("effi_exp_feedback.pdf"),
        width = 8,
        height = 6)
