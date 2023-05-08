prop_w_wy <- function(audit, data){
  
  ### P(W = w) ###
  
  # calculate P(X = x) in the population
  x_pop <- aggregate(freq ~ X, data = data, sum)
  x_pop$p_x <- x_pop$freq / sum(x_pop$freq)
  
  # calculate P(W = w, X = x) in the audit sample
  wx_audit <- aggregate(freq ~ X + W, data = audit, sum)
  # calculate P(X = x) in the audit sample
  x_audit <- aggregate(freq ~ X, data = audit, sum)
  # merge both probabilities
  audit_merged <- merge(wx_audit, x_audit, by = "X", suffixes = c('.wx','.x'))
  # calculate P(W = w | X = x) in the audit sample using the calculated probabilities
  audit_merged$p_wcx <- audit_merged$freq.wx / audit_merged$freq.x
  
  # add P(W = w | X = x) from audit sample to P(X = x) from population
  wx_merged <- merge(audit_merged, x_pop, by = "X")
  # calculate P(W = w) using both probabilities
  wx_merged$p_w <- wx_merged$p_wcx * wx_merged$p_x
  
  wx_merged$prob_w <- wx_merged$p_wcx * wx_merged$p_x
  wx_merged$var_w  <- (wx_merged$p_x^2 / wx_merged$freq.x) * 
    wx_merged$p_wcx * (1 - wx_merged$p_wcx)
  
  w_probs           <- matrix(NA, length(unique(audit$W)), 2)
  colnames(w_probs) <- c("prob", "var")
  rownames(w_probs) <- 1:nrow(w_probs)
  
  for (w in 1:nrow(w_probs)){
    
    # proportions
    w_probs[w, 1] <- sum(wx_merged$prob_w[wx_merged$W == w])
    # variances
    w_probs[w, 2] <-  sum(wx_merged$var_w[wx_merged$W == w])
  }
  
  ### P(Y = y | W = w) ###
  
  # calculate P(X = x, Y = y, W = w) in audit
  wxy_audit <- aggregate(freq ~ X + Y + W, data = audit, sum)
  # merge this probability with already calculated probabilities
  wxy_merged <- merge(wxy_audit, wx_merged[ , !(names(wx_merged) %in% c("p_w"))],
                      by = c('W','X'), all = TRUE, suffixes = (c(".wxy", ".xpop")))
  # calculate P(Y = y, W = w | X = x) in audit 
  wxy_merged$p_wycx <- wxy_merged$freq.wxy / wxy_merged$freq.x 
  # multiply P(Y = y, W = w | X = x) in audit by P(X = x) in population
  wxy_merged$nom <- wxy_merged$p_x * wxy_merged$p_wycx
  
  # calculate the three terms for the variance 
  wxy_merged$term1_var <- (wxy_merged$p_x^2 / wxy_merged$freq.x) * wxy_merged$p_wycx * (1 - wxy_merged$p_wycx)
  wxy_merged$term2_var <- (wxy_merged$p_x^2 / wxy_merged$freq.x) * wxy_merged$p_wcx * (1 - wxy_merged$p_wcx)
  wxy_merged$term3_var <- (wxy_merged$p_x^2 / wxy_merged$freq.x) * wxy_merged$p_wycx * (1 - wxy_merged$p_wcx)
  
  
  # store the proportions
  yw_probs <- matrix(NA, 9, 2)
  colnames(yw_probs) <- c("prob","var")
  rownames(yw_probs) <- c("Y1W1","Y2W1","Y3W1","Y1W2","Y2W2","Y3W2","Y1W3","Y2W3","Y3W3")
  
  for (k in 1:nrow(yw_probs)){
    # obtain second character in string, representing x category
    y <- as.integer(substr(rownames(yw_probs)[k],2,2))
    # obtain fourth character in string, representing w category
    w <- as.integer(substr(rownames(yw_probs)[k],4,4))
    
    # calculate proportions
    yw_probs[k,1] <- sum(wxy_merged$nom[wxy_merged$Y == y & wxy_merged$W == w]) / w_probs[w]
    
    # variances
    yw_probs[k,2] <- ( sum(wxy_merged$term1_var[wxy_merged$Y == y & wxy_merged$W == w]) +
                         (yw_probs[k,1])^2 * sum(wxy_merged$term2_var[wxy_merged$Y == y & wxy_merged$W == w]) +
                         (-2 * yw_probs[k,1]) * sum(wxy_merged$term3_var[wxy_merged$Y == y & wxy_merged$W == w]) ) / (w_probs[w, 1])^2
    
  }
  
  return(list(w_probs, yw_probs))
}
