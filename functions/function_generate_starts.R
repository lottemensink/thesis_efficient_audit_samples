generate_starts <- function(tab){
  # calculate n_i++, totals from covariate X 
  tot <- aggregate(tab$freq, by = tab[, "X", drop = FALSE], FUN = sum)
  names(tot) <- c("X", "tot")
  # calculate n_ij+, totals for combinations X and Y
  rtot <- aggregate(tab$freq, by = tab[, c("X", "Y"), drop = FALSE], FUN = sum)
  names(rtot) <- c("X", "Y", "rtot")
  # calculate n_i+k, totals for combinations X and Z
  ktot <- aggregate(tab$freq, by = tab[,c("X", "Z")], FUN = sum)
  names(ktot) <- c("X", "Z", "ktot")
  
  # merge tot, rtot and ktot with tab
  tab0 <- merge(tab, tot, by = "X")
  tab0 <- merge(tab0, rtot, by = c("X", "Y"))
  tab0 <- merge(tab0, ktot, by = c("X", "Z"))
  
  # compute individual contribution of every XY combination to the deviance
  tab0$contr <- log(tab0$freq) + log(tab0$tot) - log(tab0$rtot) - log(tab0$ktot)

  # order tab0 in the same way that tab is ordered
  tab0 <- tab0[order(tab0$Z, tab0$Y, tab0$X),]

  # if the contribution is negative, the current freq is too small. In this
  # case, we want additional units in the sample for this XY combination. Hence, 
  # we set deltaplus to 1 and deltamin to zero
  # if the contribution is positive, the current freq is too large. In this case, 
  # we want to remove units from the sample for this XY combination. Hence, 
  # we set deltaplus to 0 and deltamin to 1
  dplus <- (tab0$contr < 0)[tab0$Z == 1]
  dmin <- (tab0$contr > 0)[tab0$Z == 1]
  x0 <- rep(0, nrow(tab))
  x0[c(dplus,dmin)] <- 1
  
  # create matrix to fill with starting values
  x0_matrix <- matrix(ncol = 18, nrow = 15)
  
  # fill matrix with some generic starting values possibilities
  x0_matrix[1,] <- x0
  x0_matrix[2,] <- c(rep(1, nrow(tab)/2), rep(0, nrow(tab)/2))
  x0_matrix[3,] <- c(rep(0, nrow(tab)/2), rep(0, nrow(tab)/2))
  x0_matrix[4,] <- c(rep(1, nrow(tab)/2), as.integer(tab$freq[which(tab$Z == 1)] != 0))
  x0_matrix[5,] <- c(rep(0, nrow(tab)/2), as.integer(tab$freq[which(tab$Z == 1)] != 0))
  
  # also add some random starting values that fall between the bounds created by tab
  up_lim <- tab$freq[c(which(tab$Z == 0), which(tab$Z == 1))]
  starts <- rep(NA, 18)
  
  for (i in 1:10){
    for (j in 1:length(starts)){
      starts[j] <- runif(n = 1, min = 0, max = up_lim[j])
    }
    x0_matrix[i+5,] <- starts
  }
  
  
  
  return(x0_matrix)
}
