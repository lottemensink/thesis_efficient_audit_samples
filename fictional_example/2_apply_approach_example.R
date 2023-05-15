# load nloptr library for optimization problem
library(nloptr)

# load functions
source("../functions/function_generate_starts.R")
source("../functions/function_apply_constraints.R")

# load data and get it into the right format
data <- read.csv("example_data.csv")
colnames(data) <- c("X", "Y", "Z", "freq")
data$X <- as.numeric(factor(data$X))
data$Y <- as.numeric(factor(data$Y))
data$Z <- ifelse(data$Z == "Yes", 1, 0)

# set the user-specified constraints
dev_threshold <- 7.6
min_sample <- 0

# generate starting values, 
starting_values <- generate_starts(data)

# get apply_constraints function in the right form for the nloptr function
hin <- function(x) apply_constraints(x, bound = dev_threshold, min_sample = min_sample, cst = cst)
.hin <- match.fun(hin)
hin <- function(x) (-1) * .hin(x)
hinjac <- function(x) nl.jacobian(x, hin)

# set the right options for the nloptr function 
opts <- nl.opts()
opts$algorithm <- "NLOPT_LD_AUGLAG"
opts$tol_constraints_ineq <- c(1e-2, 1e-2)
opts$maxeval <- 10000
opts$local_opts <- list("algorithm" = "NLOPT_LD_SLSQP", "eval_grad_f" = NULL, 
                        xtol_rel = 1e-2, "maxeval" = 10000, "tol_constraints_ineq" = 1e-2,
                        ftol_abs = 1e-6, ftol_rel = 1e-6, NLOPT_INT = c(1,2))

# apply the sample size approach trying different starting values until
# convergence is reached
for (i in 1:nrow(starting_values)){
  
  # apply the constrained minimization procedure using the nloptr function
  sol <- nloptr(
    # specify starting values, these differ over iterations until algorithm converges
    x0 = starting_values[i, ],
    # specify target function: the sum of delta+ and delta- values
    eval_f = function(x) { sum(x) },
    # specify the derivative of the target function
    eval_grad_f = function(x) { rep(1, nrow(data)) },
    # specify lower bounds for the solution: the delta values cannot be below 0
    lb = rep(0, nrow(data)),
    # specify upper bounds for the solution: there can't be more units moved 
    # from Z = 0 to Z = 1 then are initially in Z = 0, and vice versa
    ub = tab$freq[c(which(data$Z == 0), which(data$Z == 1))],
    # non-linear constraint on the deviance and minimum sample
    eval_g_ineq = hin,
    eval_jac_g_ineq = hinjac,
    # set the prespecified options
    opts = opts)
  
  print(sol$status)
  
  # break loop when convergence has been reached
  if(sol$status > 0 & sol$status != 5){break}
  
}

# inspect solution
solution <- sol$solution
print(solution)

# solution provided by the algorithm is not rounded
# unless the algorithm suggests a value very close to zero, we round upwards
# this way, we can obtain the delta plus (to include in audit) and delta min
# (to exclude from audit) values
deltaplus <- solution[1:(length(solution)/2)]
deltaplus <- ifelse(deltaplus > 0.01, ceiling(deltaplus), round(deltaplus, 0))
print(deltaplus)
deltamin <- solution[(1+(length(solution)/2)):length(solution)]
deltamin <- ifelse(deltamin > 0.01, ceiling(deltamin), round(deltamin, 0))
print(deltamin)

