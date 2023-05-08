generate_w_probs <- function(sim_cons_wy, y_select){
  # set constraints to create 3 probabilities P(X = x, Y = y, W = w) 
  # that add up to the right marginal probability P(X = x, Y = y)
  n <- 3
  lower <- rep(0, 3)
  lower[y_select] <- sim_cons_wy - 1
  upper <- 1
  marg_prob <- 1
  constr <- list(constr = rbind(-diag(n), diag(n), rep(1, n), rep(-1, n)),
                 dir = rep("<=", 2*n+2),
                 rhs = c(lower, rep(upper, n), marg_prob, -marg_prob))
  
  # generate a probability that adds up to right marginal
  w_probs <- as.data.frame(hitandrun(constr, n.samples = 1))
  colnames(w_probs) <- c("W1", "W2", "W3")
  
  return(w_probs)
  
}
