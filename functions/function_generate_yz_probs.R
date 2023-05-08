generate_yz_probs <- function(sim_cons_sf, n_xyz){
  
  # set up the constraints for 3 random numbers that add up to the sampling fraction
  n <- 3
  lower <- 0
  upper <- sim_cons_sf
  marg_prob <- sim_cons_sf
  constr <- list(constr = rbind(-diag(n), diag(n), rep(1, n), rep(-1, n)),
                 dir = rep("<=", 2*n+2),
                 rhs = c(rep(lower, n), rep(upper, n), marg_prob, -marg_prob))
  
  # generate P(Y = y, Z = 1) 
  YZ_probs <- as.data.frame(hitandrun(constr, n.samples = n_xyz))
  
  return(YZ_probs)
  
}