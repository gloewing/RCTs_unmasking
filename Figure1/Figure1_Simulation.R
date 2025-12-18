library(patchwork)
library(ggplot2)
library(dplyr)
wd <- "/Users/loewingergc/Desktop/NIMH Research/Causal/Causal_Psychedelics/Paper/Figures/Drafts/Belief-Conditional"

#-----------------------------------------------------------
## Unmasking + unmeasured confounding (BELIEF ONLY) scenario
#-----------------------------------------------------------
N <- 5000 
set.seed(1)

gamma_0 <- -3.5
gamma_1 <- 3.5 # -3.5
gamma_2 <- 5.5 # 5.5
beta_0 <- 3#30
beta_1 <- -3 ### controlled direct effect of treatment 3
beta_2 <- -3 # -3
beta_3 <- -25 # -25

A <- rbinom(N, 1, 0.5)
U <- runif(N, -0.1,1.5)
B <- rbinom(N, 1, plogis(gamma_0 + gamma_1 * A + gamma_2 * U))
Y <- rnorm(N, beta_0 + beta_1 * A + beta_2 * B + beta_3 * U, sd = 1)

# numerically check
w_b1a1 <- (mean(A == 1 & B == 1) / mean(A == B)) / 
  ((mean(A == 1 & B == 1) / mean(A == B)) + (mean(A == 1 & B == 0) / mean(A != B)))
w_b0a1 <- (mean(A == 1 & B == 0) / mean(A != B)) / 
  ((mean(A == 1 & B == 1) / mean(A == B)) + (mean(A == 1 & B == 0) / mean(A != B)))
w_b1a0 <- (mean(A == 0 & B == 1) / mean(A != B)) /
  ((mean(A == 0 & B == 0) / mean(A == B)) + (mean(A == 0 & B == 1) / mean(A != B)))
w_b0a0 <- (mean(A == 0 & B == 0) / mean(A == B)) /
  ((mean(A == 0 & B == 0) / mean(A == B)) + (mean(A == 0 & B == 1) / mean(A != B)))

mean(Y[B == 1 & A == 1]) - mean(Y[B == 1 & A == 0]) # belief = 1-conditional effect
mean(Y[B == 0 & A == 1]) - mean(Y[B == 0 & A == 0]) # belief = 0-conditional effect


mean(Y[A == 1]) - mean(Y[A == 0]) # belief = 1-conditional effect

## Correct Guessing Rate (CGR) - adjusted ATE [in this case, a mixture of belief-conditional effects]
mean(Y[B == 1 & A == 1]) * w_b1a1 + 
  mean(Y[B == 0 & A == 1]) * w_b0a1 -
  mean(Y[B == 1 & A == 0]) * w_b1a0 -
  mean(Y[B == 0 & A == 0]) * w_b0a0
mean(Y[A == 1]) - mean(Y[A == 0]) # marginal ATE

myColors2 <- c("#4285F4","#808080", "#E69F00", "darkgreen", "darkgrey")
df <- data.frame(Y = Y,
                 B = as.factor(B),
                 A = as.factor(A)) %>%
  mutate(A = ifelse(A == 1, "Treatment", "Control"))

################
# Marginal ATE
################
# Create the plot with 95% CIs
summary_p1 <- df %>%
  mutate(A=forcats::fct_relevel(A,c("Treatment", "Control"))) %>%
           group_by( A) %>%
  summarise(
    mean_Y = mean(Y),
    se_Y = 1.96 * sd(Y) / sqrt(n()),
    .groups = "drop"
  )

# Create the barplot with standard error bars
p1 <- ggplot(summary_p1, aes(x = factor(A), y = mean_Y, fill = A)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7,
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_Y - se_Y, ymax = mean_Y + se_Y), 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.8) +
  labs(
    # title = "Treatment Effect",
    x = "Treatment Arm (A)",
    y = "LS mean change from baseline in \nCAPS-5 total severity score",
    fill = "Treatment Arm (A)"
  ) +
  scale_fill_manual(values = myColors2) +
  scale_color_manual(values = myColors2) +
  #scale_x_discrete(labels = c("Control", "Treatment")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

################################
# Belief Conditional ATE
################################
# Create the plot with 95% CIs

summary_p2 <- df %>%
  mutate(A=forcats::fct_relevel(A,c("Treatment", "Control"))) %>%
  group_by(B, A) %>%
  summarise(
    mean_Y = mean(Y),
    se_Y = 1.96 * sd(Y) / sqrt(n()),
    .groups = "drop"
  )

# Create the barplot with standard error bars
p2 <- ggplot(summary_p2, aes(x = factor(B), y = mean_Y, fill = A)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), 
           width = 0.7,
           color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = mean_Y - se_Y, ymax = mean_Y + se_Y), 
                position = position_dodge(width = 0.8), 
                width = 0.25, linewidth = 0.8) +
  labs(
    # title = "Outcome Stratified by Belief",
    x = "Belief (B)",
    y = "LS mean change from baseline in \nCAPS-5 total severity score",
    fill = "Treatment Arm (A)"
  ) +
  scale_fill_manual(values = myColors2) +
  scale_color_manual(values = myColors2) +
  scale_x_discrete(labels = c("Belief", "No Belief")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())



################################
# Belief Conditional Effect
################################
# belief conditional effects in sample
lm.mod <- lm(Y ~ A*B, data = df)
model.sum <- summary(lm.mod)

# conditional on belief = 0 
belief0 <- coef(lm.mod)[2]
belief0.se <- model.sum$coefficients[2,2]

# conditional on belief = 1
L <- c(0, 1, 0, 1)  # Corresponds todifference 
cov_matrix <- vcov(lm.mod) # Get the covariance matrix of the coefficients
belief1 <- sum(L * coef(lm.mod)) # Compute the contrast estimate
belief1.se <- sqrt(t(L) %*% cov_matrix %*% L) # Compute the standard error of the contrast

# data.frame for plot
df.mod <- data.frame(b1 = belief0,
                     b1.se = belief0.se,
                     b0 = belief1,
                     b0.se = belief1.se)

# Reformat the data for ggplot:
plot_df <- data.frame( group = c("B", "A"),
                       value = as.numeric(df.mod[1, c(1,3)]),
                       se = as.numeric(df.mod[1, c(2,4)]) )

p3=ggplot(plot_df, aes(x = group, y = value)) +
  geom_point() + # aes(color = group)
  geom_errorbar(aes(ymin = value - 1.96* se, 
                    ymax = value + 1.96* se), #,color = group), 
                width = 0.2) +
  # labs(x = "Group", y = "Value") +
  labs(
    title = "Estimates Stratified by Belief",
    x = "Belief (B)",
    y = "Treatment Effect: Treatment vs. Control Arms",
    fill = "Treatment Arm (A)" # Customize legend title.
  ) +
  scale_fill_manual(values = myColors2) +
  scale_color_manual(values = myColors2) +
  scale_x_discrete(labels = c("Belief", "No Belief")) + # Customize x-axis labels for binary variable.
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  theme( axis.title.x = element_blank(),
         plot.title = element_text(hjust = 0.5) ) +
  geom_hline(yintercept = beta_1, linetype = "dashed", color = "gray35") +
  coord_cartesian( ylim = c(-1, 2.5) ) +
  annotate(
    "text",
    x = 1.5,             # Horizontally centered between "A" (1) and "B" (2)
    y = beta_1 + 0.1,     # A bit above the dashed line
    label = "True Causal Effect",
    color = "gray35",
    vjust = 0)


plots2 <- p1 + 
  p2 + 
  plot_layout(axis_titles = "collect", guides = "collect") & 
  theme(legend.position = "bottom") 


setwd(wd)
ggsave( plot = plots,
        paste0("belief_conditional2.pdf"),
        width = 9,
        height = 5)
