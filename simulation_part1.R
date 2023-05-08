# load the required functions
source("functions/function_calc_deviance.R")
source("functions/function_prop_w_wy.R")
source("functions/function_generate_w_probs.R")
source("functions/function_generate_yz_probs.R")

# load packages
library(hitandrun)
library(tidyverse)

# set seed for replicability
set.seed(123)

# specify all conditions, investigate main effects and full factorial
conditions_wy <- list(wy_cons = c(1/2, 1/3, 1/4), sf_cons = 0.1, yx_cons = 3)
conditions_sf <- list(wy_cons = 1/4, sf_cons = c(0.003, 0.01, 0.03, 0.1), yx_cons = 3)
conditions_yx <- list(wy_cons = 1/4, sf_cons = 0.1, yx_cons = c(1, 2, 3))
conditions_ff <- list(wy_cons = c(1/2, 1/4), sf_cons = c(0.003, 0.1), yx_cons = c(1, 3))

# create a grid with all the conditons
sim_cons <- rbind(expand.grid(conditions_wy), expand.grid(conditions_sf), 
                  expand.grid(conditions_yx), expand.grid(conditions_ff))
# remove duplicate conditions and save to RData to reuse in part 2
sim_cons <- distinct(sim_cons)

# specify the number of categories in X, Y, W and Z
X <- c(1, 2, 3)
Y <- c(1, 2, 3)
W <- c(1, 2, 3)
Z <- c(0, 1)

# specify the values for the different YX relations, generated based on cramer's V
YX_cons <- list()

# cramer's V = 0.1
YX_cons[[1]] <- matrix(c(.4/3, .3/3, .3/3,
                         .3/3, .4/3, .3/3,
                         .3/3, .3/3, .4/3), ncol = 3, nrow = 3, byrow = TRUE)

# cramer's V = 0.2
YX_cons[[2]] <- matrix(c(7/45, 4/45, 4/45,
               4/45, 7/45, 4/45,
               4/45, 4/45, 7/45), ncol = 3, nrow = 3, byrow = TRUE)

# cramer's V = 0.4
YX_cons[[3]] <- matrix(c(.6/3, .2/3, .2/3,
                         .2/3, .6/3, .2/3,
                         .2/3, .2/3, .6/3), ncol = 3, nrow = 3, byrow = TRUE)

# create a grid with all possible combinations between X and Y
YX_combinations <- expand.grid(list(X = X, Y = Y))

# set the number of XYZ & XYWZ distributions to generate per condition
n_xyz <- 1000
n_xywz <- 50

# set population size 
pop_size <- 100000

save.image(file = "workspaces/parameters_conditions.RData")

# create list to store the generated distributions
generated_XYWZ <- vector(mode = "list", length = nrow(sim_cons))
for (i in 1:nrow(sim_cons)){
  generated_XYWZ[[i]] <- vector(mode = "list", length = n_xywz*n_xyz)
}

# vector to store the deviance for every XYZ distribution
deviance <- matrix(nrow = n_xyz, ncol = nrow(sim_cons))
# vectors to store bias in p(W = w) for every data set
error_pw1 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_pw2 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_pw3 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
# vectors to store bias in p(Y = y | W = w) for every data set
error_py1w1 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py2w1 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py3w1 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py1w2 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py2w2 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py3w2 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py1w3 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py2w3 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))
error_py3w3 <- matrix(nrow = n_xyz*n_xywz, ncol = nrow(sim_cons))

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
    
    # the constant term of the deviance does not depend on the procedure 
    # and can already be computed
    tot  <- aggregate(XYZ_data$freq, by = XYZ_data[ , "X", drop = FALSE], FUN = sum)
    rtot <- aggregate(XYZ_data$freq, by = XYZ_data[ , c("X","Y")], FUN = sum)
    cst <- 2 * sum(tot$x * log(tot$x), na.rm=TRUE) - 2 * sum(rtot$x * log(rtot$x), na.rm=TRUE)
    
    # store deviance
    deviance[j, i] <- calc_deviance(XYZ_data[, c(1:3, 5)], cst)
    
    # generate 50 distributions of XYWZ for every XYZ distribution
    for (k in 1:n_xywz){
      # data frame to store XYWZ
      XYWZ <- rbind(XYZ[, 1:3], XYZ[, 1:3], XYZ[, 1:3])
      XYWZ$W <- rep(c(1,2,3), each = nrow(XYZ_data))
      XYWZ$freq <- rep(NA, nrow(XYWZ))
      
      # generate different P(W = w) for every possible combination of Y and X
      for (m in 1:nrow(YX_combinations)){
        
        # create index for Y and X
        x_select <- YX_combinations[m, 1]
        y_select <- YX_combinations[m, 2]
        
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
      
      # calculate P(W = w) in the audit sample 
      probs <- prop_w_wy(audit, XYWZ)
      pw_audit <- probs[[1]][, "prob"]
      pyw_audit <- probs[[2]][, "prob"]
      
      # calculate bias in P(W = w) for every audit sample 
      error_pw1[(j-1)*n_xywz + k, i] <- (pw_audit - pw_true[, "prob"])[1]
      error_pw2[(j-1)*n_xywz + k, i] <- (pw_audit - pw_true[, "prob"])[2]
      error_pw3[(j-1)*n_xywz + k, i] <- (pw_audit - pw_true[, "prob"])[3]
      
      # calculate bias in P(Y = y | W = w) for every audit sample 
      error_py1w1[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[1]
      error_py2w1[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[2]
      error_py3w1[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[3]
      error_py1w2[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[4]
      error_py2w2[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[5]
      error_py3w2[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[6]
      error_py1w3[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[7]
      error_py2w3[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[8]
      error_py3w3[(j-1)*n_xywz + k, i] <- (pyw_audit - pyw_true[, "prob"])[9]
      
      # save the generated XYWZ distribution
      generated_XYWZ[[i]][[(j-1)*n_xywz + k]] <- XYWZ
      
    }
  
  }
  
}

save(generated_XYWZ, file = "workspaces/generated_XYWZ.RData")

# get the results into appropriate dataframes

results <- list()
data <- matrix(ncol = 13, nrow = n_xyz)

for (i in 1:nrow(sim_cons)){
  data[, 1] <- deviance[, i]
  
  data[, 2] <- colMeans(matrix(error_pw1[, i], n_xywz))
  data[, 3] <- colMeans(matrix(error_pw2[, i], n_xywz))
  data[, 4] <- colMeans(matrix(error_pw3[, i], n_xywz))
  
  data[, 5] <- colMeans(matrix(error_py1w1[, i], n_xywz))
  data[, 6] <- colMeans(matrix(error_py2w1[, i], n_xywz))
  data[, 7] <- colMeans(matrix(error_py3w1[, i], n_xywz))
  data[, 8] <- colMeans(matrix(error_py1w2[, i], n_xywz))
  data[, 9] <- colMeans(matrix(error_py2w2[, i], n_xywz))
  data[, 10] <- colMeans(matrix(error_py3w2[, i], n_xywz))
  data[, 11] <- colMeans(matrix(error_py1w3[, i], n_xywz))
  data[, 12] <- colMeans(matrix(error_py2w3[, i], n_xywz))
  data[, 13] <- colMeans(matrix(error_py3w3[, i], n_xywz))
  
  colnames(data) <- c("deviance", "bias_pw1", "bias_pw2", "bias_pw3", "bias_py1w1",
                      "bias_py2w1", "bias_py3w1", "bias_py1w2", "bias_py2w2",
                      "bias_py3w2", "bias_py1w3", "bias_py2w3", "bias_py3w3")
  
  results[[i]] <- as.data.frame(data)
}

save(results, file = "workspaces/results_simulation_part1.RData")

# extra object with results that also contains individual bias in every audit

extra_results_ind <- list()
data <- matrix(ncol = 25, nrow = n_xywz*n_xyz)

for (i in 1:nrow(sim_cons)){
  data[, 1] <- rep(deviance[, i], each = n_xywz)
  
  data[, 2] <- rep(colMeans(matrix(error_pw1[, i], n_xywz)), each = n_xywz)
  data[, 3] <- rep(colMeans(matrix(error_pw2[, i], n_xywz)), each = n_xywz)
  data[, 4] <- rep(colMeans(matrix(error_pw3[, i], n_xywz)), each = n_xywz)
  
  data[, 5] <- rep(colMeans(matrix(error_py1w1[, i], n_xywz)), each = n_xywz)
  data[, 6] <- rep(colMeans(matrix(error_py2w1[, i], n_xywz)), each = n_xywz)
  data[, 7] <- rep(colMeans(matrix(error_py3w1[, i], n_xywz)), each = n_xywz)
  data[, 8] <- rep(colMeans(matrix(error_py1w2[, i], n_xywz)), each = n_xywz)
  data[, 9] <- rep(colMeans(matrix(error_py2w2[, i], n_xywz)), each = n_xywz)
  data[, 10] <- rep(colMeans(matrix(error_py3w2[, i], n_xywz)), each = n_xywz)
  data[, 11] <- rep(colMeans(matrix(error_py1w3[, i], n_xywz)), each = n_xywz)
  data[, 12] <- rep(colMeans(matrix(error_py2w3[, i], n_xywz)), each = n_xywz)
  data[, 13] <- rep(colMeans(matrix(error_py3w3[, i], n_xywz)), each = n_xywz)
  
  
  data[, 14] <- error_pw1[, i]
  data[, 15] <- error_pw2[, i]
  data[, 16] <- error_pw3[, i]
  
  data[, 17] <- error_py1w1[, i]
  data[, 18] <- error_py2w1[, i]
  data[, 19] <- error_py3w1[, i]
  data[, 20] <- error_py1w2[, i]
  data[, 21] <- error_py2w2[, i]
  data[, 22] <- error_py3w2[, i]
  data[, 23] <- error_py1w3[, i]
  data[, 24] <- error_py2w3[, i]
  data[, 25] <- error_py3w3[, i]
  
  colnames(data) <- c("deviance", "bias_pw1", "bias_pw2", "bias_pw3", "bias_py1w1",
                      "bias_py2w1", "bias_py3w1", "bias_py1w2", "bias_py2w2",
                      "bias_py3w2", "bias_py1w3", "bias_py2w3", "bias_py3w3",
                      "error_pw1", "error_pw2", "error_pw3", "error_py1w1", 
                      "error_py2w1", "error_py3w1", "error_py1w2", "error_py2w2", 
                      "error_py3w2", "error_py1w3", "error_py2w3", "error_py3w3")
  
  extra_results_ind[[i]] <- as.data.frame(data)
}

save(extra_results_ind, file = "workspaces/results_simulation_part1_extra.RData")
