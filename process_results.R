# load results, functions and packages
load("workspaces/results_simulation_part1.RData")
load("workspaces/results_simulation_part2.RData")
load("workspaces/conditions.RData")

source("functions/function_calc_dev_threshold.R")

library(tidyverse)
library(quantreg)
library(ggpubr)

# amounts of maximum tolerable bias that were used to determine deviance threshold
max_bias <- c(0.05, 0.02, 0.01)

# initialize the vectors for all the results to save 
dev_thresholds <- c()
mean_finalaudit <- c()
mean_dplus <- c()
mean_dmin <- c()
bias_pw_before <- c()
bias_pw <- c()
bias_pyw_before <- c()
bias_pyw <- c()
bias_exc <- c()

# extract relevant results per condition
for (i in 1:nrow(sim_cons)){
  
  bias_dev_data <- results[[i]]
  
  # transform the data to make it suitable for deviance threshold calculation
  data_pw <- bias_dev_data %>% pivot_longer(cols = c("bias_pw1", "bias_pw2", "bias_pw3"), 
                                            values_to = "bias_pw") %>%  select(deviance, bias_pw)
  data_pyw <- bias_dev_data %>% pivot_longer(cols = c("bias_py1w1", "bias_py2w2", "bias_py3w3"), 
                                             values_to = "bias_pyw") %>% select(deviance, bias_pyw)
  data_pynw <- bias_dev_data %>% pivot_longer(cols = c("bias_py2w1", "bias_py3w1", "bias_py1w2", 
                                                       "bias_py3w2", "bias_py1w3", "bias_py2w3"), values_to = "bias_pynw") %>%
    select(deviance, bias_pynw)
  
  # calculate bias in P(W = w) before applying the sample size approach
  bias_pw_before[i] <- ((mean(abs(bias_dev_data$bias_pw1)) +
    mean(abs(bias_dev_data$bias_pw2)) + mean(abs(bias_dev_data$bias_pw3)))/3)
  
  # calculate bias in P(Y = y | W = w) before applying the sample size approach
  bias_pyw_before[i] <- ((mean(abs(bias_dev_data$bias_py1w1)) + mean(abs(bias_dev_data$bias_py1w2)) 
    + mean(abs(bias_dev_data$bias_py1w3))+ mean(abs(bias_dev_data$bias_py2w1)) + mean(abs(bias_dev_data$bias_py2w2)) + 
    mean(abs(bias_dev_data$bias_py2w3))+ mean(abs(bias_dev_data$bias_py3w1)) + mean(abs(bias_dev_data$bias_py3w2)) + 
    mean(abs(bias_dev_data$bias_py3w3)))/9)
  
  # extract relevant results for every amount of maximum tolerable bias
  for (j in 1:length(max_bias)){
    
    sim_data <- as.data.frame(results_part2[[(j-1)*12 + i]])
    sim_data <- na.omit(sim_data)
    
    # calculate deviance thresholds for every estimand and select the most conservative deviance threshold
    dev_threshold_pw <- calc_dev_threshold(data_pw, max_bias = max_bias[j], tau = 1)
    dev_threshold_pyw <- calc_dev_threshold(data_pyw, max_bias = max_bias[j], tau = 1)
    dev_threshold_pynw <- calc_dev_threshold(data_pynw, max_bias = max_bias[j], tau = 1)
    
    dev_thresholds[(j-1)*12 + i] <- min(dev_threshold_pw, dev_threshold_pyw, dev_threshold_pynw)
    
    # extract final audit sample size characteristics
    mean_finalaudit[(j-1)*12 + i] <- mean(sim_data$auditsize)
    mean_dplus[(j-1)*12 + i] <- mean(sim_data$dplus)
    mean_dmin[(j-1)*12 + i] <- mean(sim_data$dmin)
    
    # calculate bias in P(W = w) in final audit sample
    bias_pw[(j-1)*12 + i] <- ((mean(abs(sim_data$bias_pw1)) +
      mean(abs(sim_data$bias_pw2)) + mean(abs(sim_data$bias_pw3)))/3)

    # calculate bias in P(Y = y | W = w) in final audit sample
    bias_pyw[(j-1)*12 + i] <- ((mean(abs(sim_data$bias_py1w1)) + mean(abs(sim_data$bias_py1w2)) 
      + mean(abs(sim_data$bias_py1w3))+ mean(abs(sim_data$bias_py2w1)) + mean(abs(sim_data$bias_py2w2)) + 
      mean(abs(sim_data$bias_py2w3))+ mean(abs(sim_data$bias_py3w1)) + mean(abs(sim_data$bias_py3w2)) + 
      mean(abs(sim_data$bias_py3w3)))/9)
    
    # calculate proportion of cases in which maximum tolerable bias was exceeded
    bias_exc[(j-1)*12 + i] <- sum(abs(sim_data[, 8:19] > max_bias[j]))/(nrow(sim_data)*12)
  }
  
}

### creating the tables for the results section ###

# table with deviance threshold for every condition
tab_dev_thresholds <- as.data.frame(cbind(sim_cons, 
                        dev_thresholds_5 = round(dev_thresholds[1:12], 2),
                        dev_thresholds_2 = round(dev_thresholds[13:24], 2), 
                        dev_thresholds_1 = round(dev_thresholds[25:36], 2)))

# table with average final audit size, added units and removed units per condition
tab_auditsize <- as.data.frame(cbind(sim_cons, 
                         final_5 = round(mean_finalaudit[1:12], 0),
                         final_2 = round(mean_finalaudit[13:24], 0), 
                         final_1 = round(mean_finalaudit[25:36], 0),
                         dplus_5 = round(mean_dplus[1:12], 0), 
                         dplus_2 = round(mean_dplus[13:24], 0), 
                         dplus_1 = round(mean_dplus[25:36], 0),
                         dmin_5 = round(mean_dmin[1:12], 0),
                         dmin_2 = round(mean_dmin[13:24], 0), 
                         dmin_1 = round(mean_dmin[25:36], 0)))
  
# table with average absolute bias in P(W = w) before and after for every condition
tab_biaspw <- as.data.frame(cbind(sim_cons, bias_pw5 = round(bias_pw[1:12], 3),
                         bias_pw2 = round(bias_pw[13:24], 3), 
                         bias_pw1 = round(bias_pw[25:36], 3)))

# table with average absolute bias in P(Y = y | W = w) before and after for 
# every condition
tab_biaspyw <- as.data.frame(cbind(sim_cons, bias_pyw_before,
                         bias_pyw5 = round(bias_pyw[1:12], 3),
                         bias_pyw2 = round(bias_pyw[13:24], 3), 
                         bias_pyw1 = round(bias_pyw[25:36], 3)))

# table with proportion of data set in which observed bias in the audit sample 
# exceeded maximum tolerable bias
tab_biasexc <- as.data.frame(cbind(sim_cons, exc_5 = round(bias_exc[1:12], 3),
                        exc_2 = round(bias_exc[13:24], 3), 
                        exc_1 = round(bias_exc[25:36], 3)))

### creating the plots for the results section ###

# extract the data for the most desirable and least desirable condition
most_des <- results[[3]]
least_des <- results[[9]]

# plot bias in P(W = w) for most desirable condition
most_pw <- most_des %>% 
  pivot_longer(cols = c("bias_pw1", "bias_pw2", "bias_pw3"), 
               values_to = "bias_pw") %>%
  select(deviance, bias_pw)
pw_plot1 <- ggplot(aes(x = deviance, y = bias_pw), data = most_pw) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank(), text = element_text(family = "Times New Roman", size = 19),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  ggtitle("Bias in P(W = w)")

# plot bias in P(Y = y | W = w) for most desirable condition when y = w
most_pyw <- most_des %>% 
  pivot_longer(cols = c("bias_py1w1", "bias_py2w2", "bias_py3w3"), 
               values_to = "bias_pyw") %>%
  select(deviance, bias_pyw)
pyw_plot1 <- ggplot(aes(x = deviance, y = bias_pyw), data = most_pyw) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank(), text = element_text(family = "Times New Roman", size = 19),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  ggtitle("Bias in P(Y = y | W = w) for y = w")

# plot bias in P(Y = y | W = w) for most desirable condition when y != w
most_pynw <- most_des %>% 
  pivot_longer(cols = c("bias_py2w1", "bias_py3w1", "bias_py1w2", "bias_py3w2", 
                        "bias_py1w3", "bias_py2w3"), values_to = "bias_pynw") %>%
  select(deviance, bias_pynw)
pynw_plot1 <- ggplot(aes(x = deviance, y = bias_pynw), data = most_pynw) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank(), text = element_text(family = "Times New Roman", size = 19),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  ggtitle("Bias in P(Y = y | W = w) for y ??? w")


# plot bias in P(W = w) for least desirable condition
least_pw <- least_des %>% 
  pivot_longer(cols = c("bias_pw1", "bias_pw2", "bias_pw3"), 
               values_to = "bias_pw") %>%
  select(deviance, bias_pw)
pw_plot2 <- ggplot(aes(x = deviance, y = bias_pw), data = least_pw) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank(), text = element_text(family = "Times New Roman", size = 19),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  ggtitle("Bias in P(W = w)")

# plot bias in P(Y = y | W = w) for least desirable condition when y = w
least_pyw <- least_des %>% 
  pivot_longer(cols = c("bias_py1w1", "bias_py2w2", "bias_py3w3"), 
               values_to = "bias_pyw") %>%
  select(deviance, bias_pyw)
pyw_plot2 <- ggplot(aes(x = deviance, y = bias_pyw), data = least_pyw) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank(), text = element_text(family = "Times New Roman", size = 19),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  ggtitle("Bias in P(Y = y | W = w) for y = w")

# plot bias in P(Y = y | W = w) for least desirable condition when y != w
least_pynw <- least_des %>% 
  pivot_longer(cols = c("bias_py2w1", "bias_py3w1", "bias_py1w2", "bias_py3w2", 
                        "bias_py1w3", "bias_py2w3"), values_to = "bias_pynw") %>%
  select(deviance, bias_pynw)
pynw_plot2 <- ggplot(aes(x = deviance, y = bias_pynw), data = least_pynw) +
  geom_point() +
  labs(x = "Deviance", y = "Bias") +
  theme_minimal() +
  theme(legend.title = element_blank(), text = element_text(family = "Times New Roman", size = 19),
        plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  ggtitle("Bias in P(Y = y | W = w) for y ??? w")

# create grid with all the plots
grid <- ggarrange(pw_plot1, pyw_plot1, pynw_plot1, pw_plot2, pyw_plot2, pynw_plot2,
                   nrow = 2, ncol = 3)
