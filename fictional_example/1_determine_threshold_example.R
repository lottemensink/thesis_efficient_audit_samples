library(hitandrun)
library(quantreg)
library(tidyverse)
source("../functions/function_generate_yz_probs.R")
source("../functions/function_calc_deviance.R")
source("../functions/function_generate_w_probs.R")
source("../functions/function_prop_w_wy.R")
set.seed(123)

# load data and put it into suitable structure
data <- read.csv("example_data.csv")
colnames(data) <- c("X", "Y", "Z", "freq")
data$X <- as.numeric(factor(data$X))
data$Y <- as.numeric(factor(data$Y))
data$Z <- ifelse(data$Z == "Yes", 1, 0)

### STEP 1 ###
# Extract relevant data characteristics

# population size
pop_size <- sum(data$freq)
# sampling fraction of initial audit
sampling_fraction <- sum(data[data$Z == 1, ]$freq)/pop_size
# maximum expected measurement error
max_me <- 0.25
# relation between Y and X (including calculation of cramer's V)
YX <- aggregate(data$freq, by = list(data$X, data$Y), FUN = sum, drop = TRUE)
YX$x <- YX$x/pop_size
colnames(YX) <- c("X", "Y", "prob")
V <- sqrt(chisq.test(x = YX[, "prob"] * 10000)$statistic/(10000*(3-1)))

# specify number of XYZ distributions to generate
n_xyz <- 1000
# specify number of XYWZ distributions to generate per XYZ distribution
n_xywz <- 50
# specify maximum amount of tolerable bias
tolerable_bias <- 0.05
  
### STEP 2 ###
# Generate a large number of different sets of probabilities P(Y = y, Z = 1) 
# that add up to the right marginal probability

YZ_base <- rep((1 - sampling_fraction)/3, 3)
# generate P(Y = y, Z = 1) 
YZ_probs <- generate_yz_probs(sampling_fraction, n_xyz)

### STEP 3 ###
# Generate a large number of YZ distributions and corresponding data sets

XYZ <- data[, 1:3]
XYZ_distributions <- list()

for (i in 1:n_xyz){
  
  # add column for the probabilities P(X = x, Y = y, Z = z)
  XYZ[, "prob"] <- rep(NA, nrow(XYZ))
  
  for (j in 1:nrow(YX)){
    
    # create index for Y and X
    x_select <- YX[j, 1]
    y_select <- YX[j, 2]
    
    # add the right probabilities to the XYZ distribution
    XYZ[XYZ$Y == y_select & XYZ$X == x_select & XYZ$Z == 0, "prob"] <- 
      YZ_base[y_select] * YX[YX$X == x_select & YX$Y == y_select, "prob"] * 3
    XYZ[XYZ$Y == y_select & XYZ$X == x_select & XYZ$Z == 1, "prob"] <- 
      YZ_probs[i, y_select] * YX[YX$X == x_select & YX$Y == y_select, "prob"] * 3
  }
  
  freq <- rmultinom(1, size = pop_size, prob = XYZ$prob)
  XYZ_data <- cbind(XYZ, freq)
  
  XYZ_distributions[[i]] <- XYZ_data

}

### STEP 4 ###
# Calculate the deviance threshold for every XYZ distribution

deviance <- c()

for (i in 1:n_xyz){
  
  tab <- XYZ_distributions[[i]][, c(1:3,5)]
  
  tot  <- aggregate(tab$freq, by = tab[ , "X", drop = FALSE], FUN = sum)
  rtot <- aggregate(tab$freq, by = tab[ , c("X","Y")], FUN = sum)
  cst <- 2 * sum(tot$x * log(tot$x), na.rm=TRUE) - 2 * sum(rtot$x * log(rtot$x), na.rm=TRUE)
  
  deviance[i] <-calc_deviance(tab, cst)
}

### STEP 5 ###
# Determine the bias in the target estimand for every generated data set
# Here, we choose P(W = w) and P(Y = y | W = w as target estimands)

# dataframes to store bias in p(W = w) for every data set
error_pw1 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_pw2 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_pw3 <- matrix(nrow = n_xyz, ncol = n_xywz)
# dataframes to store bias in p(Y = y | W = w) for every data set
error_py1w1 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py2w1 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py3w1 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py1w2 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py2w2 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py3w2 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py1w3 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py2w3 <- matrix(nrow = n_xyz, ncol = n_xywz)
error_py3w3 <- matrix(nrow = n_xyz, ncol = n_xywz)

# expand xyz distributions to xywz distributions and determine bias
for (i in 1:n_xyz){
  
  XYZ_data <- XYZ_distributions[[i]]
  
  # generate 50 distributions of XYWZ for every XYZ distribution
  for (j in 1:n_xywz){
    
    # data frame to store XYWZ
    XYWZ <- rbind(XYZ_data[, 1:3], XYZ_data[, 1:3], XYZ_data[, 1:3])
    XYWZ$W <- rep(c(1,2,3), each = nrow(XYZ_data))
    XYWZ$freq <- rep(NA, nrow(XYWZ))
    
    # generate different P(W = w) for every possible combination of Y and X
    for (k in 1:nrow(YX)){
      
      # create index for Y and X
      x_select <- YX[k, 1]
      y_select <- YX[k, 2]
      
      # generate P(W = w) randomly taking maximum measurement error into account
      w_probs <- generate_w_probs(max_me, y_select)
      
      # use probabilities of W and multinomial distribution to generate XYWZ
      XYWZ[XYWZ$X == x_select & XYWZ$Y == y_select & XYWZ$Z == 0, "freq"] <-
        rmultinom(1, size = as.numeric(XYZ_data[XYZ_data$X == x_select & 
          XYZ_data$Y == y_select & XYZ_data$Z == 0, "freq"]), prob = w_probs)
      
      XYWZ[XYWZ$X == x_select & XYWZ$Y == y_select & XYWZ$Z == 1, "freq"] <- 
        rmultinom(1, size = as.numeric(XYZ_data[XYZ_data$X == x_select & 
          XYZ_data$Y == y_select & XYZ_data$Z == 1, "freq"]), prob = w_probs)
      
    }
    
    # calculate true probs using the distribution the data was generated from
    pw_true <- aggregate(freq ~ W, data = XYWZ, sum)
    pw_true$prob <- pw_true$freq/sum(XYWZ$freq)
    pyw_true <- aggregate(freq ~ Y + W, data = XYWZ, sum)
    pyw_true$prob <- pyw_true$freq/sum(XYWZ$freq)
    pyw_true[1:3,"prob"] <- pyw_true[1:3,"prob"]/sum(pyw_true[1:3,"prob"])
    pyw_true[4:6,"prob"] <- pyw_true[4:6,"prob"]/sum(pyw_true[4:6,"prob"])
    pyw_true[7:9,"prob"] <- pyw_true[7:9,"prob"]/sum(pyw_true[7:9,"prob"])
    
    # subset the audit sample 
    audit <- XYWZ[XYWZ$Z == 1, ]
    
    # calculate P(W = w) and P(Y = y | W = w) in the audit sample 
    probs <- prop_w_wy(audit, XYWZ)
    pw_audit <- probs[[1]][, "prob"]
    pyw_audit <- probs[[2]][, "prob"]
    
    # calculate bias in P(W = w) for every audit sample 
    error_pw1[i, j] <- (pw_audit - pw_true[, "prob"])[1]
    error_pw2[i, j] <- (pw_audit - pw_true[, "prob"])[2]
    error_pw3[i, j] <- (pw_audit - pw_true[, "prob"])[3]
    
    # calculate bias in P(Y = y | W = w) for every audit sample 
    error_py1w1[i, j] <- (pyw_audit - pyw_true[, "prob"])[1]
    error_py2w1[i, j] <- (pyw_audit - pyw_true[, "prob"])[2]
    error_py3w1[i, j] <- (pyw_audit - pyw_true[, "prob"])[3]
    error_py1w2[i, j] <- (pyw_audit - pyw_true[, "prob"])[4]
    error_py2w2[i, j] <- (pyw_audit - pyw_true[, "prob"])[5]
    error_py3w2[i, j] <- (pyw_audit - pyw_true[, "prob"])[6]
    error_py1w3[i, j] <- (pyw_audit - pyw_true[, "prob"])[7]
    error_py2w3[i, j] <- (pyw_audit - pyw_true[, "prob"])[8]
    error_py3w3[i, j] <- (pyw_audit - pyw_true[, "prob"])[9]
  }
}

bias_pw1 <- rowMeans(error_pw1)
bias_pw2 <- rowMeans(error_pw2)
bias_pw3 <- rowMeans(error_pw3)

bias_py1w1 <- rowMeans(error_py1w1)
bias_py2w1 <- rowMeans(error_py2w1)
bias_py3w1 <- rowMeans(error_py3w1)
bias_py1w2 <- rowMeans(error_py1w2)
bias_py2w2 <- rowMeans(error_py2w2)
bias_py3w2 <- rowMeans(error_py3w2)
bias_py1w3 <- rowMeans(error_py1w3)
bias_py2w3 <- rowMeans(error_py2w3)
bias_py3w3 <- rowMeans(error_py3w3)

bias_dev_data <- as.data.frame(cbind(deviance, bias_pw1, bias_pw2, bias_pw3, 
    bias_py1w1, bias_py2w1, bias_py3w1, bias_py1w2, bias_py2w2, bias_py3w2, 
    bias_py1w3, bias_py2w3, bias_py3w3))

### STEP 6 ###
# Model the relation between maximum expected absolute bias and deviance

data_pw <- bias_dev_data %>% pivot_longer(cols = c("bias_pw1", "bias_pw2", 
              "bias_pw3"), values_to = "bias_pw") %>%  
  select(deviance, bias_pw)


data_pyw <- bias_dev_data %>% pivot_longer(cols = c("bias_py1w1", "bias_py2w1", 
              "bias_py3w1", "bias_py1w2", "bias_py2w2", "bias_py3w2", "bias_py1w3", 
              "bias_py2w3", "bias_py3w3",), values_to = "bias_pyw") %>% 
  select(deviance, bias_pyw)

# add a root term for the deviance to the bias vs. deviance data
data_pw$dev_root <- sqrt(data_pw$deviance)

# fit polynomial quantile regression model for absolute bias
# (assume intercept to be zero)
quanfit_pw <- rq(abs(bias_pw) ~ deviance + dev_root -1, data = data_pw, tau = 1)

# add a root term for the deviance to the bias vs. deviance data
data_pyw$dev_root <- sqrt(data_pyw$deviance)

# fit polynomial quantile regression model for absolute bias
# (assume intercept to be zero)
quanfit_pyw <- rq(abs(bias_pyw) ~ deviance + dev_root -1, data = data_pyw, tau = 1)

### STEP 7 ###
# Determine the deviance threshold 

# extract coefficient (assume intercept to be 0)
b1_pw <- as.numeric(quanfit_pw$coefficients[1])
b2_pw <- as.numeric(quanfit_pw$coefficients[2])

# function to solve for bias
f <- function(x) {
  b1_pw*x + b2_pw*sqrt(x) - tolerable_bias
}

# calculate deviance threshold by solving function
dev_threshold_pw <- uniroot(f, c(0, max(data_pw$deviance)))$root

# extract coefficient (assume intercept to be 0)
b1_pyw <- as.numeric(quanfit_pyw$coefficients[1])
b2_pyw <- as.numeric(quanfit_pyw$coefficients[2])

# function to solve for bias
f <- function(x) {
  b1_pyw*x + b2_pyw*sqrt(x) - tolerable_bias
}

# calculate deviance threshold by solving function
dev_threshold_pyw <- uniroot(f, c(0, max(data_pyw$deviance)))$root

min(dev_threshold_pw, dev_threshold_pyw)

# make illustrative plots of bias and deviance relation for P(W = w)

library(ggplot2)
library(ggpubr)

devbias <- ggplot(data_pw, aes(x = deviance, y = bias_pw)) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  xlim(0,600) +
  theme(legend.title = element_blank(), text = element_text(family = "serif", size = 16)) +
  ggtitle("Relation between Deviance and Bias \n ")

modelrel <- ggplot(data_pw, aes(x = deviance, y = abs(bias_pw))) +
  geom_point() +
  geom_line(aes(y = fitted(quanfit_pw), color = "Maximum predicted bias"), 
            linetype = "solid", size = 1) +
  labs(x = "Deviance", y = "Bias") +
  scale_color_manual(values = c("Maximum predicted bias" = "red")) +
  theme_minimal() +
  xlim(0,600) +
  theme(legend.title = element_blank(), text = element_text(family = "serif", size = 16)) +
  ggtitle("Polynomial Quantile Regression Model \n for Absolute Bias")

grid <- ggarrange(devbias, modelrel, nrow = 1, legend = "bottom", common.legend = TRUE)
ggsave("plot_relation.png", grid, width = 10, height = 5, bg = "white")


