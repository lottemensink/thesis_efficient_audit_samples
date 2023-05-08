apply_approach <- function(starts, dev_threshold, min_sample, tab, opts, cst){
  
  # nloptr() function can only be used with numeric, tables not integers
  # convert table to numeric:
  tab[] <- lapply(tab, as.numeric)
  
  # get apply_constraints function in the right form for the nloptr function
  hin <- function(x) apply_constraints(x, bound = dev_threshold, min_sample = min_sample, cst = cst)
  .hin <- match.fun(hin)
  hin <- function(x) (-1) * .hin(x)
  hinjac <- function(x) nl.jacobian(x, hin)
  
  for (i in 1:nrow(starts)){
    
    # apply the constrained minimization procedure using the nloptr function
    sol <- nloptr(
      # specify starting values, these differ over iterations until algorithm converges
      x0 = starts[i, ],
      # specify target function: the sum of delta+ and delta- values
      eval_f = function(x) { sum(x) },
      # specify the derivative of the target function
      eval_grad_f = function(x) { rep(1, nrow(tab)) },
      # specify lower bounds for the solution: the delta values cannot be below 0
      lb = rep(0, nrow(tab)),
      # specify upper bounds for the solution: there can't be more units moved 
      # from Z = 0 to Z = 1 then are initially in Z = 0, and vice versa
      ub = tab$freq[c(which(tab$Z == 0), which(tab$Z == 1))],
      # non-linear constraint on the deviance and minimum sample
      eval_g_ineq = hin,
      eval_jac_g_ineq = hinjac,
      # set the prespecified options
      opts = opts)
    
    print(sol$status)
    
    if(sol$status > 0 & sol$status != 5){break}
  }
  
  return(sol)
  
}
