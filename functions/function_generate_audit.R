generate_audit <- function(XYWZ, tab_ext){
  # obtain the change in freq per combination of X and Y, rounded to whole units
  XY_change <- cbind(tab_ext[tab_ext$Z == 1, c("X","Y")], 
                     round(tab_ext[tab_ext$Z == 1, "freq_newr"] - tab_ext[tab_ext$Z == 1, "freq"], 0))
  colnames(XY_change) <- c("X", "Y", "change")
  
  # create matrix to store the generated distribution of W 
  WYX_new <- matrix(NA, nrow(XY_change)*3, 4)
  colnames(WYX_new) <- c("X","Y","W","freq")
  
  # Loop over rows XY_change to generate W for every combination of X and Y
  for(k in 1:nrow(XY_change)){
    
    # obtain X category, Y category, and the number of units that need to be 
    # either removed or added for this particular combination of X and Y
    X_select <- paste0(XY_change[k, "X"])
    Y_select <- paste0(XY_change[k, "Y"])
    n_select <- XY_change[k, "change"]
    
    # fill in the values for X, Y and W in the new matrix that stores the 
    # generated distribution of W
    WYX_new[(1+((k-1)*3)):(3+((k-1)*3)), "X"] <- rep(X_select, 3)
    WYX_new[(1+((k-1)*3)):(3+((k-1)*3)), "Y"] <- rep(Y_select, 3)
    WYX_new[(1+((k-1)*3)):(3+((k-1)*3)), "W"] <- c(1,2,3)
    
    # if n_select > 0, units need to be added to the audit sample. We need the 
    # distribution of W for all units in the final audit sample. To do this, we 
    # generate the distribution of W for the units that need to be added to the 
    # sample, and add this to the distribution of W in the initial sample. 
    if (n_select > 0) {
      WYX_new[(1+((k-1)*3)):(3+((k-1)*3)), "freq"] <- 
        XYWZ[XYWZ$X == X_select & XYWZ$Y == Y_select & XYWZ$Z == 1, "freq"] +
        rmvhyper(nn = 1, n = XYWZ[XYWZ$X == X_select & XYWZ$Y == Y_select & XYWZ$Z == 0, "freq"], k = abs(n_select))
    }
    # if n_select < 0, units need to be removed from the audit sample. We need  
    # the distribution of W for all units in the final audit sample. To do this,  
    # we generate the distribution of W for the units that need to be removed 
    # from the sample, and substract this from the distribution of W in the 
    # initial sample. 
    else if (n_select < 0) {
      WYX_new[(1+((k-1)*3)):(3+((k-1)*3)), "freq"] <- 
        XYWZ[XYWZ$X == X_select & XYWZ$Y == Y_select & XYWZ$Z == 1, "freq"] -
        rmvhyper(nn = 1, n = XYWZ[XYWZ$X == X_select & XYWZ$Y == Y_select & XYWZ$Z == 1, "freq"], k = abs(n_select))
     
    }
    # if n_select == 0, no units need to be added or removed from the audit 
    # sample. The distribution of W in the final audit sample is equal to the 
    # ditribution of W in the initial audit sample. 
    else {
      WYX_new[(1+((k-1)*3)):(3+((k-1)*3)), "freq"] <- 
        XYWZ[XYWZ$X == X_select & XYWZ$Y == Y_select & XYWZ$Z == 1, "freq"]
    }
  }
  
  # convert the distribution of X, Y and W to a XYWZ frame and return
  WYX     <- as.data.frame(WYX_new)
  WYX$freq <- as.numeric(as.character(WYX$freq))
  
  return(WYX)
}

