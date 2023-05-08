# specify function for the constraints. This is a function of the 
# solution x, the user specified deviance threshold and the minimum sample 
# size. 
apply_constraints <- function(x, bound, min_sample, cst) {
  deltaplus <- x[1:(length(x)/2)]
  deltamin <- x[(1+(length(x)/2)):length(x)]
  # compute the new frequency table under solution x
  tab$freq <- tab$freq +
    c(-deltaplus, -deltamin) + 
    c(deltamin, deltaplus)
  # calculate deviance under solution x using new frequency table
  ktot <- aggregate(tab$freq, by = tab[ , c('X','Z')], FUN = sum)
  ind <- which(tab$freq > 0)
  ktot_ind <- which(ktot$x > 0)
  dev <- cst + 2 * sum(tab$freq[ind] * log(tab$freq[ind])) -
    2 * sum(ktot$x[ktot_ind] * log(ktot$x[ktot_ind]))
  ineq1 <- bound - dev
  # the constraint should be in the form hin(x) >= 0, therefore we return
  # the threshold - calculated deviance under the solution
  sol_sample <- sum(tab$freq[which(tab$Z == 1)])
  ineq2 <- sol_sample - min_sample
  return(c(ineq1, ineq2))
}
