load("workspaces/results_simulation_part1.RData")
load("workspaces/parameters_conditions.RData")

source("functions/function_calc_dev_threshold.R")
source("functions/function_calc_deviance.R")
source("functions/function_prop_w_wy.R")
source("functions/function_generate_starts.R")
source("functions/function_apply_constraints.R")
source("functions/function_apply_approach.R")
source("functions/function_generate_audit.R")

library(tidyverse)
library(quantreg)
library(nloptr)
library(extraDistr)
library(hitandrun)

set.seed(123)

# before we apply the sample size procedure, we want to generate new XYWZ distributions
# with corresponding data sets to prevent problems with overfitting. For this, 
# we re-use the code of the simulation study part 1 with a few minor modifications

# create list to store the generated distributions
generated_XYWZ <- vector(mode = "list", length = nrow(sim_cons))
for (i in 1:nrow(sim_cons)){
  generated_XYWZ[[i]] <- vector(mode = "list", length = n_xyz)
}

# create grid with distribution X, Y and Z
variables <- list(Y = Y, X = X, Z = Z)
XYZ <- expand.grid(variables)

# loop over the different conditions and construct the relation between 
# bias and deviance in every condition
for (i in 1:nrow(sim_cons)){
  
  # specify YX in the selected condition
  YX <- YX_cons[[sim_cons[i,3]]]
  
  # generate YZ probs for the selected condition 
  # that is, taking measurement error into account
  YZ_base <- rep((1 - sim_cons[i,2])/3, 3)
  YZ_probs <- generate_yz_probs(sim_cons[i, 2], n_xyz)
  
  for (j in 1:n_xyz){
    
    # print condition 
    cat(sprintf('Condition %d XYZ distribution %d \n', i, j))
    
    # add column for the probabilities P(X = x, Y = y, W = w)
    XYZ[, "prob_XYZ"] <- rep(NA, nrow(XYZ))
    for (k in 1:nrow(YX_combinations)){
      
      # create index for Y and X
      x_select <- YX_combinations[k, 1]
      y_select <- YX_combinations[k, 2]
      
      # add the right probabilities to the XYWZ distribution
      XYZ[XYZ$Y == y_select & XYZ$X == x_select & XYZ$Z == 0, "prob_XYZ"] <- 
        YZ_base[y_select] * YX[x_select, y_select] * 3
      XYZ[XYZ$Y == y_select & XYZ$X == x_select & XYZ$Z == 1, "prob_XYZ"] <- 
        YZ_probs[j, y_select] * YX[x_select, y_select] * 3
    }
    
    
    # now generating 1 data set corresponding to this distribution
    XYZ_data <- cbind(XYZ[, 1:4], rmultinom(1, size = pop_size, prob = XYZ[, "prob_XYZ"]))
    colnames(XYZ_data) <- c("Y", "X", "Z", "prob_XYZ", "freq_XYZ")

    # data frame to store XYWZ
    XYWZ <- rbind(XYZ[, 1:3], XYZ[, 1:3], XYZ[, 1:3])
    XYWZ$W <- rep(c(1,2,3), each = nrow(XYZ_data))
    XYWZ$freq <- rep(NA, nrow(XYWZ))
    
    # generate different P(W = w) for every possible combination of Y and X
    for (k in 1:nrow(YX_combinations)){
      
      # create index for Y and X
      x_select <- YX_combinations[k, 1]
      y_select <- YX_combinations[k, 2]
      
      # generate P(W = w)
      w_probs <- generate_w_probs(sim_cons[i, 1], y_select)
      
      # use probabilities of W and multinomial distribution to generate XYWZ
      XYWZ[XYWZ$X == x_select & XYWZ$Y == y_select & XYWZ$Z == 0, "freq"] <-
        rmultinom(1, size = as.numeric(XYZ_data[XYZ_data$X == x_select & 
                                                  XYZ_data$Y == y_select & XYZ_data$Z == 0, "freq_XYZ"]), prob = w_probs)
      
      XYWZ[XYWZ$X == x_select & XYWZ$Y == y_select & XYWZ$Z == 1, "freq"] <- 
        rmultinom(1, size = as.numeric(XYZ_data[XYZ_data$X == x_select & 
                                                  XYZ_data$Y == y_select & XYZ_data$Z == 1, "freq_XYZ"]), prob = w_probs)
    }
    
    generated_XYWZ[[i]][[j]] <- XYWZ
        
      
  }
  
}

# we want to test various amounts of maximum tolerable bias in every condition, 
# therefore, we expand the conditions with max_bias
sim_cons_bias <- rbind(cbind(sim_cons, max_bias = 0.05),
                       cbind(sim_cons, max_bias = 0.02),
                       cbind(sim_cons, max_bias = 0.01))

# create list to save variances
varlist <- vector(mode = "list", length = nrow(sim_cons_bias))

for (i in 1:nrow(sim_cons_bias)){
  varlist[[i]] <- vector(mode = "list", length = n_xyz)
}

# set the right options for the nloptr function 
opts <- nl.opts()
opts$algorithm <- "NLOPT_LD_AUGLAG"
opts$tol_constraints_ineq <- c(1e-2, 1e-2)
opts$maxeval <- 10000
opts$local_opts <- list("algorithm" = "NLOPT_LD_SLSQP", "eval_grad_f" = NULL, 
                        xtol_rel = 1e-2, "maxeval" = 10000, "tol_constraints_ineq" = 1e-2,
                        ftol_abs = 1e-6, ftol_rel = 1e-6, NLOPT_INT = c(1,2))

# create infrastructure to store all the results 
results_part2 <- list()
data_results_part2 <- matrix(NA, nrow = n_xyz, ncol = 19)
colnames(data_results_part2) <- 
  c("conv", "dev_before", "dev_after", "dev_after_rounded", "dplus", "dmin", 
    "auditsize", "bias_pw1", "bias_pw2", "bias_pw3", "bias_py1w1", "bias_py2w1", 
    "bias_py3w1", "bias_py1w2", "bias_py2w2", "bias_py3w2", "bias_py1w3", 
    "bias_py2w3", "bias_py3w3")

# initialize vectors for storing bias in probabilities for every audit sample
error_pw1 <- c()
error_pw2 <- c()
error_pw3 <- c()

error_py1w1 <- c()
error_py2w1 <- c()
error_py3w1 <- c()
error_py1w2 <- c()
error_py2w2 <- c()
error_py3w2 <- c()
error_py1w3 <- c()
error_py2w3 <- c()
error_py3w3 <- c()

for (i in 1:nrow(sim_cons_bias)){
  
  index <- as.numeric(rownames(sim_cons[sim_cons$wy_cons == sim_cons_bias[i,"wy_cons"] & 
                    sim_cons$sf_cons == sim_cons_bias[i,"sf_cons"] &
                    sim_cons$yx_cons == sim_cons_bias[i, "yx_cons"],]))
  
  # select the bias vs. deviance data for this condition for all the probs
  # calculate deviance thresholds for all the probs 
  bias_dev_data <- results[[index]]
  
  data_pw <- bias_dev_data %>% pivot_longer(cols = c("bias_pw1", "bias_pw2", "bias_pw3"), 
                 values_to = "bias_pw") %>%  select(deviance, bias_pw)
  
  dev_threshold_pw <- calc_dev_threshold(data_pw, max_bias = sim_cons_bias[i,4], tau = 1)
  
  data_pyw <- bias_dev_data %>% pivot_longer(cols = c("bias_py1w1", "bias_py2w2", "bias_py3w3"), 
                 values_to = "bias_pyw") %>% select(deviance, bias_pyw)
  
  dev_threshold_pyw <- calc_dev_threshold(data_pyw, max_bias = sim_cons_bias[i,4], tau = 1)
  
  data_pynw <- bias_dev_data %>% pivot_longer(cols = c("bias_py2w1", "bias_py3w1", "bias_py1w2", 
                 "bias_py3w2", "bias_py1w3", "bias_py2w3"), values_to = "bias_pynw") %>%
                 select(deviance, bias_pynw)

  dev_threshold_pynw <- calc_dev_threshold(data_pynw, max_bias = sim_cons_bias[i,4], tau = 1)

  # select the lowest deviance threshold to limit all types of max bias 
  dev_threshold <- min(dev_threshold_pw, dev_threshold_pyw, dev_threshold_pynw)

  # select the XYWZ distributions that correspond to the current condition
  XYWZ_distributions <- generated_XYWZ[[index]]
  
  # use sampling fraction audit sample as minimum sample size
  min_sample <- 0
  
  # loop over the 1000 XYWZ distributions generated in the condition
  for (j in 1:n_xyz){
    
    cat(sprintf('Condition %d XYZ distribution %d \n', i, j))
    flush.console()

    # select XYWZ distribution
    XYWZ <- XYWZ_distributions[[j]]

    # create tab
    tab <- aggregate(XYWZ[, "freq"], by = list(X = XYWZ$X, 
                                         Y = XYWZ$Y, 
                                         Z = XYWZ$Z), 
                     FUN  = sum, drop = TRUE)
    colnames(tab) = c("X","Y","Z","freq")
    
    # calculate true probs using the distribution the data was generated from
    pw_true <- aggregate(freq ~ W, data = XYWZ, sum)
    pw_true$prob <- pw_true$freq/sum(XYWZ$freq)
    pyw_true <- aggregate(freq ~ Y + W, data = XYWZ, sum)
    pyw_true$prob <- pyw_true$freq/sum(XYWZ$freq)
    pyw_true[1:3,"prob"] <- pyw_true[1:3,"prob"]/sum(pyw_true[1:3,"prob"])
    pyw_true[4:6,"prob"] <- pyw_true[4:6,"prob"]/sum(pyw_true[4:6,"prob"])
    pyw_true[7:9,"prob"] <- pyw_true[7:9,"prob"]/sum(pyw_true[7:9,"prob"])
    
    # the constant term of the deviance does not depend on the procedure 
    # and can already be computed
    tot  <- aggregate(tab$freq, by = tab[ , "X", drop = FALSE], FUN = sum)
    rtot <- aggregate(tab$freq, by = tab[ , c("X","Y")], FUN = sum)
    cst <- 2 * sum(tot$x * log(tot$x), na.rm=TRUE) - 2 * sum(rtot$x * log(rtot$x), na.rm=TRUE)
    
    # calculate deviance before
    data_results_part2[j, "dev_before"] <- calc_deviance(tab, cst)
    
    # generate starting values and apply sample size approach
    starts <- generate_starts(tab)
    sol <- apply_approach(starts, dev_threshold, min_sample, tab, opts, cst)
    
    # store and round solution 
    solution <- sol$solution
    deltaplus <- solution[1:(length(solution)/2)]
    deltaplusr <- ifelse(deltaplus > 0.01, ceiling(deltaplus), round(deltaplus, 0))
    deltamin <- solution[(1+(length(solution)/2)):length(solution)]
    deltaminr <- ifelse(deltamin > 0.01, ceiling(deltamin), round(deltamin, 0))
    
    
    # create extended frequency table, with original freq and the new freq from solution
    tab_ext <- cbind(tab, freq_new = tab$freq +
                       c(-deltaplus, -deltamin) + 
                       c(deltamin, deltaplus),
                     freq_newr = tab$freq + 
                       c(-deltaplusr, -deltaminr) + 
                       c(deltaminr, deltaplusr))
    
    # store deviance after (rounded)
    data_results_part2[j, "dev_after"] <- calc_deviance(tab_ext[, c(1:3, 5)], cst)
    data_results_part2[j, "dev_after_rounded"] <- calc_deviance(round(tab_ext[, c(1:3, 6)], 0), cst)
    
    # list for storing variances estimated from each audit sample 
    est_vars <- list()
    
    for (k in 1:n_xywz){
      # generate audit sample
      audit <- generate_audit(XYWZ, tab_ext)
      
      # calculate P(W = w) in the audit sample 
      probs <- prop_w_wy(audit, XYWZ)
      pw_audit <- probs[[1]][, "prob"]
      pyw_audit <- probs[[2]][, "prob"]
      
      # calculate bias in P(W = w) for every audit sample 
      error_pw1[k] <- (pw_audit - pw_true[, "prob"])[1]
      error_pw2[k] <- (pw_audit - pw_true[, "prob"])[2]
      error_pw3[k] <- (pw_audit - pw_true[, "prob"])[3]
      
      # calculate bias in P(Y = y|W = w) for every audit sample 
      error_py1w1[k] <- unname(pyw_audit[1] - pyw_true[1, "prob"])
      error_py2w1[k] <- unname(pyw_audit[2] - pyw_true[2, "prob"])
      error_py3w1[k] <- unname(pyw_audit[3] - pyw_true[3, "prob"])
      error_py1w2[k] <- unname(pyw_audit[4] - pyw_true[4, "prob"])
      error_py2w2[k] <- unname(pyw_audit[5] - pyw_true[5, "prob"])
      error_py3w2[k] <- unname(pyw_audit[6] - pyw_true[6, "prob"])
      error_py1w3[k] <- unname(pyw_audit[7] - pyw_true[7, "prob"])
      error_py2w3[k] <- unname(pyw_audit[8] - pyw_true[8, "prob"])
      error_py3w3[k] <- unname(pyw_audit[9] - pyw_true[9, "prob"])
      
      # save the estimated variances from the audit sample in a list
      est_vars[[k]] <- c(probs[[1]][, "var"], probs[[2]][, "var"])
    }
    
    # store convergence
    data_results_part2[j, "conv"] <- sol$status
    # store number of added units
    data_results_part2[j, "dplus"] <- sum(deltaplusr)
    # store number of removed units
    data_results_part2[j, "dmin"] <- sum(deltaminr)
    # store audit sample size 
    data_results_part2[j, "auditsize"] <- sum(tab_ext[tab_ext$Z == 1, ]$freq_newr)
    
    # store bias in P(W = w)
    data_results_part2[j, "bias_pw1"] <- mean(error_pw1)
    data_results_part2[j, "bias_pw2"] <- mean(error_pw2)
    data_results_part2[j, "bias_pw3"] <- mean(error_pw3)
    
    # store bias in P(Y = y | W = w)
    data_results_part2[j, "bias_py1w1"] <- mean(error_py1w1)
    data_results_part2[j, "bias_py2w1"] <- mean(error_py2w1)
    data_results_part2[j, "bias_py3w1"] <- mean(error_py3w1)
    data_results_part2[j, "bias_py1w2"] <- mean(error_py1w2)
    data_results_part2[j, "bias_py2w2"] <- mean(error_py2w2)
    data_results_part2[j, "bias_py3w2"] <- mean(error_py3w2)
    data_results_part2[j, "bias_py1w3"] <- mean(error_py1w3)
    data_results_part2[j, "bias_py2w3"] <- mean(error_py2w3)
    data_results_part2[j, "bias_py3w3"] <- mean(error_py3w3)
    
    # store the results with respect to variance 
    varlist[[i]][[j]] <- list(var(error_pw1), var(error_pw2), var(error_pw3), 
                              var(error_py1w1), var(error_py2w1), var(error_py3w1),
                              var(error_py1w2), var(error_py2w2), var(error_py3w2),
                              var(error_py1w3), var(error_py2w3), var(error_py3w3),
                              est_vars)

  }

  results_part2[[i]] <- data_results_part2
  save(results_part2, file = "workspaces/results_simulation_part2.RData")
  
}

save(results_part2, file = "workspaces/results_simulation_part2.RData")
