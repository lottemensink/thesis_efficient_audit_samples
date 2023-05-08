# function for calculating the deviance
calc_deviance <- function(tab, cst){
  
  ktot <- aggregate(tab[,4], by = tab[ , c('X','Z')], FUN = sum)
  ind <- which(tab[,4] != 0)
  ktot_ind <- which(ktot$x != 0)
  dev <- cst + 2 * sum(tab[,4][ind] * log(tab[,4][ind])) -
    2 * sum(ktot$x[ktot_ind] * log(ktot$x[ktot_ind]))
  
  return(dev)
}