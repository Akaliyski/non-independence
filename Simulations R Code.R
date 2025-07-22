############################################################
# Simulation Code for "Non-Independence of Nations: Revisiting a Centuries-Old Methodological Challenge"
# Code Author: Plamen Akaliyski 
# 
# Purpose:
# This script generates simulated data to illustrate the effects 
# of confounding bias in cross-national research.
# Specifically, it compares regression estimates of X on Y with 
# and without controlling for a confounding variable Z under different
# causal models. 
#
# Last updated: 21 July 2025
############################################################


# Set your working directory if needed:
# setwd("your/path/here")

# -------------------------------
# Load Required Libraries
# -------------------------------

library(patchwork)
library(ggplot2)
library(dplyr)
library(cowplot)


# -------------------------------
# Model 1: Confounding
# -------------------------------

### 1. Simulation Functions for Model 1 ###
simulate_model1_correct <- function(n, a, b, c) {
  # n: sample size
  # a: effect of Z -> X (confounding influence on X)
  # b: true causal effect of X -> Y
  # c: effect of Z -> Y (confounding influence on Y)
  # (Here we set a = c)
  
  Z <- rnorm(n)
  X <- a * Z + rnorm(n, sd = sqrt(1 - a^2))
  
  # Compute variance of (bX + cZ) = b^2 + c^2 + 2*b*c^2.
  var_eY <- 1 - (b^2 + c^2 + 2 * b * c^2)
  if (var_eY < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  Y <- b * X + c * Z + rnorm(n, sd = sqrt(var_eY))
  
  # Do not re-standardize X, Y, Z to preserve the intended covariance structure.
  model_noZ   <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model1_replicated <- function(n, a, b, c, reps = 1000) {
  out_noZ   <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  
  for (i in seq_len(reps)) {
    sim <- simulate_model1_correct(n, a, b, c)
    out_noZ[i, ]   <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
sample_sizes <- c(30, 60, 120)         # Three sample sizes.
conf_levels  <- c(0.3, 0.6)              # Confounding levels (for a and c).
b_levels     <- c(0, 0.3, 0.5)            # True effect of X -> Y.

results_list <- list()
for (n in sample_sizes) {
  for (conf in conf_levels) {
    for (b_val in b_levels) {
      key <- paste0("n", n, "_conf", conf * 100, "_b", b_val * 100)
      results_list[[key]] <- simulate_model1_replicated(n, a = conf, b = b_val, c = conf, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for (nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val    <- as.numeric(gsub("n", "", parts[1]))
    conf_val <- as.numeric(gsub("conf", "", parts[2])) / 100
    b_val    <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz   <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Without controlling for Z", Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Controlling for Z", Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Convert Model to a factor with the requested order:
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect and Confound.
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound)))
)

# Convert SampleSize to a factor for shape mapping.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor (for setting x-axis positions).
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound to a factor using sprintf to force one decimal.
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.3", "0.6"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5),
                width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(y = "Estimated Coefficient for X",
       title = "Model 1",
       size = "Confounding") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.6")))

# Add one dashed segment per true effect group.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
m1 <- p + geom_text(data = hline_df,
                   aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                   vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(m1)

### 7. Save the Plot ###
ggsave("Model1.png", p, width = 12, height = 6, units = "in")










# -------------------------------
# Model 2
# -------------------------------

### 1. Simulation Functions for Model 2 ###
simulate_model2_correct <- function(n, b, c) {
  # n: sample size
  # b: true effect of X -> Y
  # c: effect of Z -> Y (confounding influence on Y); in Model 2, X and Z are independent.
  # X ~ N(0,1) and Z ~ N(0,1). Therefore, Var(bX + cZ) = b^2 + c^2.
  
  X <- rnorm(n)
  Z <- rnorm(n)
  
  var_eY <- 1 - (b^2 + c^2)
  if (var_eY < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  Y <- b * X + c * Z + rnorm(n, sd = sqrt(var_eY))
  
  # Do not re-standardize in order to preserve the intended variance.
  model_noZ   <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model2_replicated <- function(n, b, c, reps = 1000) {
  out_noZ   <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  for (i in seq_len(reps)) {
    sim <- simulate_model2_correct(n, b, c)
    out_noZ[i, ]   <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
# For Model 2, we use:
#  - sample sizes: 30, 60, 120
#  - confounding levels (c): 0.3, 0.5
#  - true effect (b): 0, 0.3, 0.5
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.3, 0.6)   # These are values for c.
b_levels     <- c(0, 0.3, 0.5) # Values for b.

results_list <- list()
for (n in sample_sizes) {
  for (c_val in conf_levels) {
    for (b_val in b_levels) {
      key <- paste0("n", n, "_c", c_val*100, "_b", b_val*100)
      results_list[[key]] <- simulate_model2_replicated(n, b = b_val, c = c_val, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for (nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val    <- as.numeric(gsub("n", "", parts[1]))
    c_val    <- as.numeric(gsub("c", "", parts[2])) / 100  # c is the confounding effect in Model 2.
    b_val    <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz   <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = c_val, TrueEffect = b_val,
                           Model = "Without controlling for Z", Beta = noz[1],
                           Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = c_val, TrueEffect = b_val,
                           Model = "Controlling for Z", Beta = withz[1],
                           Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Switch the order: "Without controlling for Z" should come first and "Controlling for Z" second.
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect (b) and Confound (c).
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", c=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", c=", results_df$Confound)))
)

# Convert SampleSize to factor.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor (for x-axis placement).
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound (c) to a factor using sprintf to force one decimal.
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.3", "0.5"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
# This produces three dashed segments: one for b = 0, one for b = 0.3, one for b = 0.5.
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5), width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(y = "Estimated Coefficient for X",
       title = "Model 2",
       size = "Confounding (Z→Y)") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.6")))

# Add one dashed segment per unique b group.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
m2 <- p + geom_text(data = hline_df,
                   aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                   vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(m2)

### 7. Save the Plot ###
ggsave("Model2.png", p, width = 12, height = 6, units = "in")









# -------------------------------
# Model 4
# -------------------------------

### 1. Simulation Functions for Model 4 ###
simulate_model4_correct <- function(n, a, b) {
  # n: sample size
  # a: instrument strength; effect of Z -> X.
  # b: true effect of X -> Y.
  #
  # Generate Z ~ N(0,1)
  Z <- rnorm(n)
  # Generate X influenced by Z: X = a * Z + error, with error sd such that Var(X)=1.
  X <- a * Z + rnorm(n, sd = sqrt(1 - a^2))
  
  # Generate Y solely as a function of X: Y = b * X + error,
  # with error sd = sqrt(1 - b^2) so that Var(Y)=1.
  var_eY <- 1 - b^2
  if(var_eY < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  Y <- b * X + rnorm(n, sd = sqrt(var_eY))
  
  # Fit two models:
  # Model without control: Y ~ X.
  # Model with control: Y ~ X + Z.
  model_noZ <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(m) {
    ci <- confint(m)["X", ]
    beta <- coef(m)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model4_replicated <- function(n, a, b, reps = 1000) {
  mat_noZ <- matrix(NA, nrow = reps, ncol = 3)
  mat_withZ <- matrix(NA, nrow = reps, ncol = 3)
  
  for(i in seq_len(reps)) {
    sim <- simulate_model4_correct(n, a, b)
    mat_noZ[i, ] <- sim$noZ
    mat_withZ[i, ] <- sim$withZ
  }
  
  list(noZ = colMeans(mat_noZ),
       withZ = colMeans(mat_withZ))
}

### 2. Run Simulation Over All Conditions ###
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.3, 0.6)  # a values.
b_levels     <- c(0, 0.3, 0.5)

results_list <- list()
for(n in sample_sizes) {
  for(a_val in conf_levels) {
    for(b_val in b_levels) {
      key <- paste0("n", n, "_a", a_val*100, "_b", b_val*100)
      results_list[[key]] <- simulate_model4_replicated(n, a = a_val, b = b_val, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for(nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val <- as.numeric(gsub("n", "", parts[1]))
    a_val <- as.numeric(gsub("a", "", parts[2])) / 100  # a is instrument strength, used as "Confound"
    b_val <- as.numeric(gsub("b", "", parts[3])) / 100
    
    # Reverse the order: "Controlling for Z" comes second now becomes:
    # We want model order: "Without controlling for Z" appears first, then "Controlling for Z".
    # (That is, in each pair, without control comes first.)
    withz <- res_list[[nm]]$withZ
    noz <- res_list[[nm]]$noZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = a_val, TrueEffect = b_val,
                           Model = "Without controlling for Z",
                           Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = a_val, TrueEffect = b_val,
                           Model = "Controlling for Z",
                           Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Ensure the Model order is "Without controlling for Z" then "Controlling for Z".
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect (b) and Confound (a).
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", a=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", a=", results_df$Confound)))
)

# Convert SampleSize to a factor.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor (for x-axis placement).
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound to a factor (using sprintf) for dot size mapping.
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.3", "0.5"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5), width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(x = "Scenario",
       y = "Estimated Coefficient for X",
       title = "Model 4",
       size = "Confounding") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.5")))

# Add one dashed segment per unique b group (TrueEffect), with legend for True Effect Size.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
m4 <- p + geom_text(data = hline_df,
                   aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                   vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(m4)

### 7. Save the Plot ###
ggsave("Model4.png", p, width = 12, height = 6, units = "in")









# -------------------------------
# Model 5
# -------------------------------

### 1. Simulation Functions for Model 5 ###
simulate_model5_correct <- function(n, conf, b) {
  # n: sample size
  # conf: confounding level (affects both U's effect and Z's effect on X)
  # b: true effect of X -> Y
  #
  # U ~ N(0,1)
  # Z ~ N(0,1)
  # X = conf * U + conf * Z + error_X, with:
  #      Var(error_X) = 1 - 2*(conf^2)   (to make Var(X)=1)
  U <- rnorm(n)
  Z <- rnorm(n)
  error_X <- rnorm(n, sd = sqrt(1 - 2 * conf^2))
  X <- conf * U + conf * Z + error_X
  
  # Y = b * X + conf * U + error_Y, with:
  #      Var(Y) = b^2 + conf^2 + 2*b*(conf^2) + Var(error_Y) = 1,
  # so set Var(error_Y) = 1 - (b^2 + conf^2 + 2 * b * conf^2)
  error_Y_var <- 1 - (b^2 + conf^2 + 2 * b * conf^2)
  if (error_Y_var < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  error_Y <- rnorm(n, sd = sqrt(error_Y_var))
  Y <- b * X + conf * U + error_Y
  
  # Fit two regression models:
  # Model 1 (without Z): Y ~ X
  # Model 2 (with Z): Y ~ X + Z
  model_noZ <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model5_replicated <- function(n, conf, b, reps = 1000) {
  out_noZ <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  for (i in seq_len(reps)) {
    sim <- simulate_model5_correct(n, conf, b)
    out_noZ[i, ] <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
# We vary:
#   - Sample sizes: 30, 60, 120.
#   - Confounding levels (conf): 0.3 and 0.5.
#   - True effect (b): 0, 0.3, 0.5.
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.3, 0.6)
b_levels     <- c(0, 0.3, 0.5)

results_list <- list()
for(n in sample_sizes) {
  for(conf in conf_levels) {
    for(b_val in b_levels) {
      key <- paste0("n", n, "_conf", conf*100, "_b", b_val*100)
      results_list[[key]] <- simulate_model5_replicated(n, conf = conf, b = b_val, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for(nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val <- as.numeric(gsub("n", "", parts[1]))
    conf_val <- as.numeric(gsub("conf", "", parts[2])) / 100  # confounding level.
    b_val <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    # Build the data frame so that "Without controlling for Z" appears first.
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Without controlling for Z",
                           Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Controlling for Z",
                           Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Keep the order of the models so that "Without controlling for Z" comes first.
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect (b) and Confound (conf).
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound)))
)

# Convert SampleSize to factor.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor for x-axis placement.
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound to a factor using sprintf so that levels are exactly "0.3" and "0.5".
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.3", "0.5"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5), width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(x = "Scenario", 
       y = "Estimated Coefficient for X",
       title = "Model 5",
       size = "Confounding") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.6")))

# Add one dashed segment per unique true effect group.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
m5 <- p + geom_text(data = hline_df,
                   aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                   vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(m5)

### 7. Save the Plot ###
ggsave("Model5.png", p, width = 12, height = 6, units = "in")








# Combining the 4 plots together

# Remove axis titles from the individual plots so they do not repeat.
m1_clean <- m1 & theme(axis.title = element_blank())
m2_clean <- m2 & theme(axis.title = element_blank())
m4_clean <- m4 & theme(axis.title = element_blank())
m5_clean <- m5 & theme(axis.title = element_blank())

# Combine the plots: top row = m1 and m2, bottom row = m4 and m5;
# collect the legends and position them at the bottom.
combined <- (m1_clean + m2_clean) / (m4_clean + m5_clean) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Now add overall x and y axis labels.
# Create a combined plot using cowplot's draw_plot and add labels.
final_plot <- ggdraw() +
  draw_plot(combined, x = 0, y = 0, width = 1, height = 1) +
  draw_label(NULL, x = 0.5, y = 0.02, hjust = 0.5, size = 12) +
  draw_label("Estimated Coefficient for X", x = 0.01, y = 0.5, hjust = 0.5, angle = 90, size = 12)

# Display the final plot
print(final_plot)

# Save the final combined plot.
ggsave("combined_models.png", final_plot, width = 14, height = 10, units = "in")














### The same models but with higher confounding (and total variance increased to 2)
# For Supplementary Material Figure S1

# -------------------------------
# Model 1: Confounding
# -------------------------------

### 1. Simulation Functions for Model 1 ###
simulate_model1_correct <- function(n, a, b, c, totvar = 2) {
  # n: sample size
  # a: effect of Z -> X (confounding influence on X)
  # b: true causal effect of X -> Y
  # c: effect of Z -> Y (confounding influence on Y)
  # totvar: total variance for Y (instead of fixing to 1)
  
  Z <- rnorm(n)
  X <- a * Z + rnorm(n, sd = sqrt(1 - a^2))
  
  # Replace 1 with totvar in the error variance calculation:
  var_eY <- totvar - (b^2 + c^2 + 2 * b * c^2)
  if (var_eY < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  Y <- b * X + c * Z + rnorm(n, sd = sqrt(var_eY))
  
  # Do NOT re-standardize X, Y, Z to preserve the intended covariance structure.
  model_noZ   <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model1_replicated <- function(n, a, b, c, reps = 1000, totvar = 2) {
  out_noZ   <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  
  for (i in seq_len(reps)) {
    sim <- simulate_model1_correct(n, a, b, c, totvar)
    out_noZ[i, ]   <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
sample_sizes <- c(30, 60, 120)         # Three sample sizes.
conf_levels  <- c(0.75, 0.85)              # Note: now you can use a confounding level 
b_levels     <- c(0, 0.3, 0.5)            # True effect of X -> Y.

results_list <- list()
for (n in sample_sizes) {
  for (conf in conf_levels) {
    for (b_val in b_levels) {
      key <- paste0("n", n, "_conf", conf * 100, "_b", b_val * 100)
      results_list[[key]] <- simulate_model1_replicated(n, a = conf, b = b_val, c = conf,
                                                        reps = 1000, totvar = 2)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for (nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val    <- as.numeric(gsub("n", "", parts[1]))
    conf_val <- as.numeric(gsub("conf", "", parts[2])) / 100
    b_val    <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz   <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Without controlling for Z", Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Controlling for Z", Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Convert Model to a factor with the requested order:
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect and Confound.
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound)))
)

# Convert SampleSize to a factor for shape mapping.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor (for setting x-axis positions).
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound to a factor using sprintf to force one decimal.
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.75", "0.85"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5),
                width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(y = "Estimated Coefficient for X",
       title = "Model 1",
       size = "Confounding") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.75", "0.85")))

# Add one dashed segment per true effect group.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
hm1 <- p + geom_text(data = hline_df,
                    aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                    vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(hm1)

### 7. Save the Plot ###
ggsave("H_Model1.png", p, width = 12, height = 6, units = "in")













# -------------------------------
# Model 2
# -------------------------------

library(ggplot2)
library(dplyr)

### 1. Simulation Functions for Model 2 ###
simulate_model2_correct <- function(n, b, c) {
  # n: sample size
  # b: true effect of X -> Y
  # c: effect of Z -> Y (confounding influence on Y); in Model 2, X and Z are independent.
  # X ~ N(0,1) and Z ~ N(0,1). Therefore, Var(bX + cZ) = b^2 + c^2.
  
  X <- rnorm(n)
  Z <- rnorm(n)
  
  var_eY <- 2 - (b^2 + c^2)
  if (var_eY < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  Y <- b * X + c * Z + rnorm(n, sd = sqrt(var_eY))
  
  # Do NOT re-standardize in order to preserve the intended variance.
  model_noZ   <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model2_replicated <- function(n, b, c, reps = 1000) {
  out_noZ   <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  for (i in seq_len(reps)) {
    sim <- simulate_model2_correct(n, b, c)
    out_noZ[i, ]   <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
# For Model 2, we use:
#  - sample sizes: 30, 60, 120
#  - confounding levels (c): 0.3, 0.5
#  - true effect (b): 0, 0.3, 0.5
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.75, 0.85)   # These are values for c.
b_levels     <- c(0, 0.3, 0.5) # Values for b.

results_list <- list()
for (n in sample_sizes) {
  for (c_val in conf_levels) {
    for (b_val in b_levels) {
      key <- paste0("n", n, "_c", c_val*100, "_b", b_val*100)
      results_list[[key]] <- simulate_model2_replicated(n, b = b_val, c = c_val, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for (nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val    <- as.numeric(gsub("n", "", parts[1]))
    c_val    <- as.numeric(gsub("c", "", parts[2])) / 100  # c is the confounding effect in Model 2.
    b_val    <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz   <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = c_val, TrueEffect = b_val,
                           Model = "Without controlling for Z", Beta = noz[1],
                           Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = c_val, TrueEffect = b_val,
                           Model = "Controlling for Z", Beta = withz[1],
                           Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Switch the order: "Without controlling for Z" should come first and "Controlling for Z" second.
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect (b) and Confound (c).
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", c=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", c=", results_df$Confound)))
)

# Convert SampleSize to factor.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor (for x-axis placement).
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound (c) to a factor using sprintf to force one decimal.
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.75", "0.85"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
# This produces three dashed segments: one for b = 0, one for b = 0.3, one for b = 0.5.
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5), width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(y = "Estimated Coefficient for X",
       title = "Model 2",
       size = "Confounding (Z→Y)") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.75", "0.85")))

# Add one dashed segment per unique b group.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
hm2 <- p + geom_text(data = hline_df,
                    aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                    vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(hm2)

### 7. Save the Plot ###
ggsave("H_Model2.png", p, width = 12, height = 6, units = "in")













# -------------------------------
# Model 4
# -------------------------------

library(ggplot2)
library(dplyr)

### 1. Simulation Functions for Model 4 ###
simulate_model4_correct <- function(n, a, b) {
  # n: sample size
  # a: instrument strength; effect of Z -> X.
  # b: true effect of X -> Y.
  #
  # Generate Z ~ N(0,1)
  Z <- rnorm(n)
  # Generate X influenced by Z: X = a * Z + error, with error sd such that Var(X)=1.
  X <- a * Z + rnorm(n, sd = sqrt(1 - a^2))
  
  # Generate Y solely as a function of X: Y = b * X + error,
  # with error sd = sqrt(1 - b^2) so that Var(Y)=1.
  var_eY <- 2 - b^2
  if(var_eY < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  Y <- b * X + rnorm(n, sd = sqrt(var_eY))
  
  # Fit two models:
  # Model without control: Y ~ X.
  # Model with control: Y ~ X + Z.
  model_noZ <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(m) {
    ci <- confint(m)["X", ]
    beta <- coef(m)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model4_replicated <- function(n, a, b, reps = 1000) {
  mat_noZ <- matrix(NA, nrow = reps, ncol = 3)
  mat_withZ <- matrix(NA, nrow = reps, ncol = 3)
  
  for(i in seq_len(reps)) {
    sim <- simulate_model4_correct(n, a, b)
    mat_noZ[i, ] <- sim$noZ
    mat_withZ[i, ] <- sim$withZ
  }
  
  list(noZ = colMeans(mat_noZ),
       withZ = colMeans(mat_withZ))
}

### 2. Run Simulation Over All Conditions ###
# We vary:
#   - Sample sizes: 30, 60, 120.
#   - Instrument strength (a): 0.3 and 0.5.
#   - True effect (b): 0, 0.3, and 0.5.
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.75, 0.85)  # a values.
b_levels     <- c(0, 0.3, 0.5)

results_list <- list()
for(n in sample_sizes) {
  for(a_val in conf_levels) {
    for(b_val in b_levels) {
      key <- paste0("n", n, "_a", a_val*100, "_b", b_val*100)
      results_list[[key]] <- simulate_model4_replicated(n, a = a_val, b = b_val, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for(nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val <- as.numeric(gsub("n", "", parts[1]))
    a_val <- as.numeric(gsub("a", "", parts[2])) / 100  # a is instrument strength, used as "Confound"
    b_val <- as.numeric(gsub("b", "", parts[3])) / 100
    
    # Reverse the order: "Controlling for Z" comes second now becomes:
    # We want model order: "Without controlling for Z" appears first, then "Controlling for Z".
    # (That is, in each pair, without control comes first.)
    withz <- res_list[[nm]]$withZ
    noz <- res_list[[nm]]$noZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = a_val, TrueEffect = b_val,
                           Model = "Without controlling for Z",
                           Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = a_val, TrueEffect = b_val,
                           Model = "Controlling for Z",
                           Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Ensure the Model order is "Without controlling for Z" then "Controlling for Z".
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect (b) and Confound (a).
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", a=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", a=", results_df$Confound)))
)

# Convert SampleSize to a factor.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor (for x-axis placement).
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound to a factor (using sprintf) for dot size mapping.
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.75", "0.85"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5), width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(x = "Scenario",
       y = "Estimated Coefficient for X",
       title = "Model 4",
       size = "Confounding") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.5")))

# Add one dashed segment per unique b group (TrueEffect), with legend for True Effect Size.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
hm4 <- p + geom_text(data = hline_df,
                    aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                    vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(hm4)

### 7. Save the Plot ###
ggsave("H_Model4.png", p, width = 12, height = 6, units = "in")









# -------------------------------
# Model 5
# -------------------------------

### 1. Simulation Functions for Model 5 ###
simulate_model5_correct <- function(n, conf, b) {
  # n: sample size
  # conf: confounding level (affects both U's effect and Z's effect on X)
  # b: true effect of X -> Y
  #
  # U ~ N(0,1)
  # Z ~ N(0,1)
  # X = conf * U + conf * Z + error_X, with:
  #      Var(error_X) = 1 - 2*(conf^2)   (to make Var(X)=1)
  U <- rnorm(n)
  Z <- rnorm(n)
  error_X <- rnorm(n, sd = sqrt(2 - 2 * conf^2))
  X <- conf * U + conf * Z + error_X
  
  # Y = b * X + conf * U + error_Y, with:
  #      Var(Y) = b^2 + conf^2 + 2*b*(conf^2) + Var(error_Y) = 1,
  # so set Var(error_Y) = 1 - (b^2 + conf^2 + 2 * b * conf^2)
  error_Y_var <- 2 - (2 * b^2 + conf^2 + 2 * b * conf^2)
  if (error_Y_var < 0) stop("Negative error variance in Y generation. Adjust parameters.")
  error_Y <- rnorm(n, sd = sqrt(error_Y_var))
  Y <- b * X + conf * U + error_Y
  
  # Fit two regression models:
  # Model 1 (without Z): Y ~ X
  # Model 2 (with Z): Y ~ X + Z
  model_noZ <- lm(Y ~ X)
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model5_replicated <- function(n, conf, b, reps = 1000) {
  out_noZ <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  for (i in seq_len(reps)) {
    sim <- simulate_model5_correct(n, conf, b)
    out_noZ[i, ] <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
# We vary:
#   - Sample sizes: 30, 60, 120.
#   - Confounding levels (conf): 0.3 and 0.5.
#   - True effect (b): 0, 0.3, 0.5.
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.75, 0.85)
b_levels     <- c(0, 0.3, 0.5)

results_list <- list()
for(n in sample_sizes) {
  for(conf in conf_levels) {
    for(b_val in b_levels) {
      key <- paste0("n", n, "_conf", conf*100, "_b", b_val*100)
      results_list[[key]] <- simulate_model5_replicated(n, conf = conf, b = b_val, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for(nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val <- as.numeric(gsub("n", "", parts[1]))
    conf_val <- as.numeric(gsub("conf", "", parts[2])) / 100  # confounding level.
    b_val <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    # Build the data frame so that "Without controlling for Z" appears first.
    df <- rbind(df,
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Without controlling for Z",
                           Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Confound = conf_val, TrueEffect = b_val,
                           Model = "Controlling for Z",
                           Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
# Keep the order of the models so that "Without controlling for Z" comes first.
results_df$Model <- factor(results_df$Model, levels = c("Without controlling for Z", "Controlling for Z"))

# Create a "Scenario" factor combining TrueEffect (b) and Confound (conf).
results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", conf=", results_df$Confound)))
)

# Convert SampleSize to factor.
results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))

# Create a numeric version of the Scenario factor for x-axis placement.
results_df$ScenarioNum <- as.numeric(results_df$Scenario)

# Convert Confound to a factor using sprintf so that levels are exactly "0.3" and "0.5".
results_df$Confound2 <- factor(sprintf("%.1f", results_df$Confound), levels = c("0.75", "0.85"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Confound2)) +
  geom_hline(yintercept = c(0, 0.25, 0.50, 0.75, 1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5), width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(x = "Scenario", 
       y = "Estimated Coefficient for X",
       title = "Model 5",
       size = "Confounding") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.6")))

# Add one dashed segment per unique true effect group.
p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y, linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect Size",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

# Add text labels at the midpoint of each dashed segment.
hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
hm5 <- p + geom_text(data = hline_df,
                    aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                    vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(hm5)

### 7. Save the Plot ###
ggsave("H_Model5.png", p, width = 12, height = 6, units = "in")





# Combined plot 
# Remove axis titles from the individual plots so they do not repeat.
hm1_clean <- hm1 & theme(axis.title = element_blank())
hm2_clean <- hm2 & theme(axis.title = element_blank())
hm4_clean <- hm4 & theme(axis.title = element_blank())
hm5_clean <- hm5 & theme(axis.title = element_blank())

# Combine the plots: top row = m1 and m2, bottom row = m4 and m5;
# collect the legends and position them at the bottom.
combined <- (hm1_clean + hm2_clean) / (hm4_clean + hm5_clean) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Now add overall x and y axis labels.
# Create a combined plot using cowplot's draw_plot and add labels.
final_plot <- ggdraw() +
  draw_plot(combined, x = 0, y = 0, width = 1, height = 1) +
  draw_label(NULL, x = 0.5, y = 0.02, hjust = 0.5, size = 12) +
  draw_label("Estimated Coefficient for X", x = 0.01, y = 0.5, hjust = 0.5, angle = 90, size = 12)

# Display the final plot
print(final_plot)

# Save the final combined plot.
ggsave("H_combined_models.png", final_plot, width = 14, height = 10, units = "in")










# -------------------------------
# Model 6 (Collider Variation with Z = coll * (Y - X))
# -------------------------------

# In this model the process is identical to the original Model 1 except that 
# Z is generated as a descendant of X and Y:
#
#   X ~ N(0,1)
#   Y = b * X + u,     where u ~ N(0, sqrt(1 - b^2)) so that Var(Y)=1.
#   Z = coll * (Y - X) + error_Z.
#
# The systematic part of Z is: coll * (Y - X) = coll * [(b - 1) * X + u].
# Its theoretical variance is:  V_sys = coll^2 * [ (b - 1)^2 + (1 - b^2) ].
#
# Instead of forcing Var(Z)=1, we adjust the variance:
#   targetZ = 1   if V_sys <= 1,
#   targetZ = V_sys + 0.01   if V_sys > 1.
# Then error_Z ~ N(0, sqrt(targetZ - V_sys)).
#
# With this design, Y ~ X (without Z) recovers the true coefficient b,
# but controlling for Z (Y ~ X + Z) introduces collider bias.
# Under reasonable parameter values the induced bias is positive.


simulate_model1_collider_diff <- function(n, b, coll) {
  # n: sample size
  # b: true causal effect of X on Y.
  # coll: weight used in generating Z.
  
  # 1. Generate X ~ N(0,1)
  X <- rnorm(n)
  
  # 2. Generate Y = b * X + u, with u ~ N(0, sqrt(1 - b^2)) so that Var(Y)=1.
  u <- rnorm(n, sd = sqrt(1 - b^2))
  Y <- b * X + u
  
  # 3. Generate Z as a descendant:
  #    Z = coll * (Y - X) + error_Z.
  # Note that Y - X = (b - 1) * X + u.
  # The theoretical variance of the systematic part is:
  #    V_sys = coll^2 * [ (b - 1)^2 * Var(X) + Var(u) ]
  #          = coll^2 * [ (b - 1)^2 + (1 - b^2) ].
  V_sys <- coll^2 * ( (b - 1)^2 + (1 - b^2) )
  
  # Set target variance for Z:
  targetZ <- ifelse(V_sys <= 1, 1, V_sys + 0.01)
  noise_var_Z <- targetZ - V_sys
  if(noise_var_Z < 0) stop("Negative noise variance for Z; adjust parameters.")
  error_Z <- rnorm(n, sd = sqrt(noise_var_Z))
  
  Z <- coll * (Y - X) + error_Z
  
  # 4. Fit two regression models:
  #    Model without controlling for Z: Y ~ X.
  model_noZ <- lm(Y ~ X)
  #    Model with controlling for Z: Y ~ X + Z.
  model_withZ <- lm(Y ~ X + Z)
  
  extract_info <- function(mod) {
    ci <- confint(mod)["X", ]  # corrected from "colliderint"
    beta <- coef(mod)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(noZ = extract_info(model_noZ),
       withZ = extract_info(model_withZ))
}

simulate_model1_collider_diff_replicated <- function(n, b, coll, reps = 1000) {
  out_noZ <- matrix(NA, nrow = reps, ncol = 3)
  out_withZ <- matrix(NA, nrow = reps, ncol = 3)
  
  for(i in seq_len(reps)) {
    sim <- simulate_model1_collider_diff(n, b, coll)
    out_noZ[i, ] <- sim$noZ
    out_withZ[i, ] <- sim$withZ
  }
  
  list(noZ = colMeans(out_noZ),
       withZ = colMeans(out_withZ))
}

### 2. Run Simulation Over All Conditions ###
# We vary:
#   - Sample sizes: 30, 60, 120.
#   - True effect (b): 0, 0.3, 0.5.
#   - Collider weight (coll): 0.3 and 0.6.
sample_sizes <- c(30, 60, 120)
b_levels <- c(0, 0.3, 0.5)
collider_levels <- c(0.3, 0.6)

results_list <- list()
for(n in sample_sizes) {
  for(coll in collider_levels) {
    for(b_val in b_levels) {
      key <- paste0("n", n, "_coll", coll * 100, "_b", b_val * 100)
      results_list[[key]] <- simulate_model1_collider_diff_replicated(n, b = b_val, coll = coll, reps = 1000)
    }
  }
}

### 3. Extract the Results into a Data Frame ###
extract_results <- function(res_list) {
  df <- data.frame()
  for(nm in names(res_list)) {
    parts <- strsplit(nm, "_")[[1]]
    n_val <- as.numeric(gsub("n", "", parts[1]))
    coll_val <- as.numeric(gsub("coll", "", parts[2])) / 100
    b_val <- as.numeric(gsub("b", "", parts[3])) / 100
    
    noz <- res_list[[nm]]$noZ
    withz <- res_list[[nm]]$withZ
    
    df <- rbind(df,
                data.frame(SampleSize = n_val, Collider = coll_val, TrueEffect = b_val,
                           Model = "Without controlling for Z", 
                           Beta = noz[1], Lower = noz[2], Upper = noz[3]),
                data.frame(SampleSize = n_val, Collider = coll_val, TrueEffect = b_val,
                           Model = "Controlling for Z", 
                           Beta = withz[1], Lower = withz[2], Upper = withz[3])
    )
  }
  df
}

results_df <- extract_results(results_list)
print(results_df)

### 4. Prepare Data for Plotting ###
results_df$Model <- factor(results_df$Model,
                           levels = c("Without controlling for Z", "Controlling for Z"))

results_df$Scenario <- factor(
  paste0("b=", results_df$TrueEffect, ", coll=", results_df$Collider),
  levels = sort(unique(paste0("b=", results_df$TrueEffect, ", coll=", results_df$Collider)))
)

results_df$SampleSize <- factor(results_df$SampleSize, levels = c(30, 60, 120))
results_df$ScenarioNum <- as.numeric(results_df$Scenario)
results_df$Collider2 <- factor(sprintf("%.1f", results_df$Collider), levels = c("0.3", "0.6"))

### 5. Compute Horizontal Dashed Reference Segments (Group by TrueEffect) ###
hline_df <- results_df %>%
  group_by(TrueEffect) %>%
  summarize(
    xmin = min(ScenarioNum) - 0.4,
    xmax = max(ScenarioNum) + 0.4,
    y = unique(TrueEffect)
  ) %>%
  ungroup()

### 6. Build the ggplot2 Plot ###
p <- ggplot(results_df, aes(x = Scenario, y = Beta,
                            color = Model, shape = SampleSize,
                            size = Collider2)) +
  geom_hline(yintercept = c(0,0.25,0.50,0.75,1), color = "grey90") +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.5),
                width = 0, size = 0.5) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5, size = 14)) +
  labs(y = "Estimated Coefficient for X",
       title = "Model 6",
       size = "Collider Weight") +
  scale_color_manual(values = c("Without controlling for Z" = "#337ab7",
                                "Controlling for Z" = "#5cb85c")) +
  scale_size_manual(values = setNames(c(2.8, 3.0), c("0.3", "0.6")))

p <- p + geom_segment(data = hline_df,
                      aes(x = xmin, xend = xmax, y = y, yend = y,
                          linetype = factor(y)),
                      color = "gray40", size = 0.7, inherit.aes = FALSE) +
  scale_linetype_manual(name = "True Effect (b)",
                        values = setNames(rep("dashed", 3), c("0", "0.3", "0.5")),
                        labels = paste("b=", c(0, 0.3, 0.5), sep = ""))

hline_df <- hline_df %>% mutate(xmid = (xmin + xmax) / 2)
p <- p + geom_text(data = hline_df,
                   aes(x = xmid, y = y, label = paste("b=", y, sep = "")),
                   vjust = -1, color = "gray40", size = 3.5, inherit.aes = FALSE)

print(p)

### 7. Save the Plot ###
ggsave("Collider_Model.png", p, width = 12, height = 6, units = "in")








# -------------------------------
# Model 7: Mediated Confounding
# -------------------------------

library(ggplot2)
library(dplyr)

### 1. One‐run simulation ###
simulate_model7_three <- function(n, conf, b) {
  # 1a. Z ~ N(0,1)
  Z  <- rnorm(n)
  # 1b. M = conf*Z + e_M  (Var(M)=1)
  e_M <- rnorm(n, sd = sqrt(1 - conf^2));  M <- conf*Z + e_M
  # 1c. X = conf*M + e_X  (Var(X)=1)
  e_X <- rnorm(n, sd = sqrt(1 - conf^2));  X <- conf*M + e_X
  # 1d. Y = b*X + conf*M + e_Y  (Var(Y)=1)
  var_eY <- 1 - (b^2 + conf^2 + 2*b*conf^2)
  if(var_eY < 0) stop("Negative error variance for Y.")
  e_Y <- rnorm(n, sd = sqrt(var_eY));  Y <- b*X + conf*M + e_Y
  
  # 2. Fit the three models:
  m0 <- lm(Y ~ X)       # unadjusted
  mZ <- lm(Y ~ X + Z)   # proxy Z
  mM <- lm(Y ~ X + M)   # true mediator M
  
  extract <- function(m) {
    ci   <- confint(m)["X", ]
    beta <- coef(m)["X"]
    c(beta = beta, lower = ci[1], upper = ci[2])
  }
  
  list(unadj = extract(m0),
       viaZ  = extract(mZ),
       viaM  = extract(mM))
}

### 2. Replication ###
simulate_model7_three_rep <- function(n, conf, b, reps = 1000) {
  out0 <- outZ <- outM <- matrix(NA, nrow = reps, ncol = 3)
  for(i in seq_len(reps)) {
    sim <- simulate_model7_three(n, conf, b)
    out0[i,] <- sim$unadj
    outZ[i,] <- sim$viaZ
    outM[i,] <- sim$viaM
  }
  list(unadj = colMeans(out0),
       viaZ  = colMeans(outZ),
       viaM  = colMeans(outM))
}

### 3. Sweep over scenarios ###
sample_sizes <- c(30, 60, 120)
conf_levels  <- c(0.3, 0.6)
b_levels     <- c(0, 0.3, 0.5)

results <- list()
for(n in sample_sizes) {
  for(conf in conf_levels) {
    for(b in b_levels) {
      key <- paste0("n",n,"_conf",conf*100,"_b",b*100)
      results[[key]] <- simulate_model7_three_rep(n, conf, b, reps = 1000)
    }
  }
}

### 4. Gather into data.frame ###
df <- do.call(rbind, lapply(names(results), function(key) {
  parts    <- strsplit(key, "_")[[1]]
  n_val    <- as.numeric(sub("n","",parts[1]))
  conf_val <- as.numeric(sub("conf","",parts[2]))/100
  b_val    <- as.numeric(sub("b","",parts[3]))/100
  r       <- results[[key]]
  rbind(
    data.frame(SampleSize=n_val, Conf=conf_val, TrueEffect=b_val,
               Model="Unadjusted",        Beta=r$unadj[1], Lower=r$unadj[2], Upper=r$unadj[3]),
    data.frame(SampleSize=n_val, Conf=conf_val, TrueEffect=b_val,
               Model="Controlling for Z", Beta=r$viaZ[1],  Lower=r$viaZ[2],  Upper=r$viaZ[3]),
    data.frame(SampleSize=n_val, Conf=conf_val, TrueEffect=b_val,
               Model="Controlling for M", Beta=r$viaM[1],  Lower=r$viaM[2],  Upper=r$viaM[3])
  )
}))
df$Model      <- factor(df$Model,
                        levels=c("Unadjusted","Controlling for Z","Controlling for M"))
df$Scenario   <- factor(paste0("b=",df$TrueEffect,", conf=",df$Conf),
                        levels=sort(unique(paste0("b=",df$TrueEffect,", conf=",df$Conf))))
df$SampleSize <- factor(df$SampleSize, levels=c(30,60,120))
df$ScenarioNum<- as.numeric(df$Scenario)
df$Conf2      <- factor(sprintf("%.1f",df$Conf), levels=c("0.3","0.6"))

hline_df <- df %>% group_by(TrueEffect) %>%
  summarize(xmin=min(ScenarioNum)-.4, xmax=max(ScenarioNum)+.4, y=unique(TrueEffect)) %>%
  ungroup()

### 5. Plot ###
p <- ggplot(df, aes(x=Scenario,y=Beta,
                    color=Model, shape=SampleSize, size=Conf2)) +
  geom_hline(yintercept=c(0,0.25,0.5,0.75,1), color="grey90") +
  geom_point(position=position_dodge(.5), size=2) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper),
                position=position_dodge(.5), width=0, size=.5) +
  geom_segment(data=hline_df,
               aes(x=xmin,xend=xmax,y=y,yend=y,linetype=factor(y)),
               color="gray40",size=.7,inherit.aes=FALSE) +
  geom_text(data=hline_df,
            aes(x=(xmin+xmax)/2,y=y,label=paste("b=",y,sep="")),
            vjust=-1,color="gray40",size=3.5,inherit.aes=FALSE) +
  scale_color_manual(values=c("Unadjusted"="#999999",
                              "Controlling for Z"="#337ab7",
                              "Controlling for M"="#5cb85c")) +
  scale_size_manual(values=c("0.3"=2.8,"0.6"=3.0)) +
  scale_linetype_manual(name="True Effect (b)",
                        values=setNames(rep("dashed",3),c("0","0.3","0.5")),
                        labels=paste("b=",c(0,0.3,0.5),sep="")) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill="white",color=NA),
    plot.background  = element_rect(fill="white",color=NA),
    panel.grid.major = element_line(color="grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle=45,hjust=1,size=10),
    plot.title       = element_text(hjust=0.5,size=14)
  ) +
  labs(
    title="Model 7",
    y="Estimated Coefficient for X",
    size="Conf Strength"
  )

print(p)
ggsave("Model7_Mediated_Confounding.png", p, width=12, height=6, units="in")



