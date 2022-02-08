
## JointNC:

## Arguments/Inputs:
# Y1: outcome of interest, assumed to be binary.
# Y2: secondary outcome.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T and to relate Y2 and T.

## Values/Outputs:
# beta_1.hat: final estimate of the treatment effect on Y1.
# sd.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.
# diagnostic.error: diagnostic if estimation fails.

## Details:
# JointNC estimates the primary outcome log-relative risk, and uses data on the secondary outcome to reduce bias 
# due to unmeasured confounders.

JointNC = function(Y1 = Y1, Y2 = Y2, T = T){
  n                     = length(Y1) # sample size
  
  # Parameter estimation
  # We assume log-linear functions relate the treatment and outcomes

  ret   <- list(beta_1.hat=NA, sd.beta_1.hat=NA, diagnostic.error=NA)
  sumT  <- sum(T)
  if (sumT == 0) {
    ret$diagnostic.error <- "No treated individuals"
    return(ret)
  }
  sumOneMinusT <- sum(1 - T)
  if (sumOneMinusT == 0) {
    diagnostic.error    <- "No untreated individuals"
    return(ret)
  }

  e.mu_1.tilde.hat      = sum((1 - T) * Y1) / sumOneMinusT # intercept for primary outcome Y1
  e.beta_1.tilde.hat    = sum(T * Y1) / sumT / e.mu_1.tilde.hat
  if (is.na(e.beta_1.tilde.hat)) {
    ret$diagnostic.error <- "No Y1 cases"
    return(ret)
  }
  if (e.beta_1.tilde.hat == 0){
    ret$diagnostic.error  <- "No treated Y1 cases"
    return(ret)
  }
  beta_1.tilde.hat      = log(e.beta_1.tilde.hat) # estimation of the treatment effect on Y1. 
  e.mu_2.tilde.hat      = sum((1 - T) * Y2) / sumOneMinusT # intercept for secondary outcome Y2
  e.beta_2.tilde.hat    = sum(T * Y2) / sumT / e.mu_2.tilde.hat
  if (is.na(e.beta_2.tilde.hat)){
    ret$diagnostic.error <- "No Y2 cases"
    return(ret)
  }
  if (e.beta_2.tilde.hat == 0){
    ret$diagnostic.error <- "No treated Y2 cases"
    return(ret)
  }
  beta_2.tilde.hat      = log(e.beta_2.tilde.hat) # estimation of the treatment effect on Y2.
  
  
      if(e.mu_1.tilde.hat == 1){
        diagnostic.error  = "All untreated individuals are Y1 cases"
        beta_1.hat        = NA
        Var_theta.hat     = matrix(NA, nrow = 4, ncol = 4)
      }else{
        if(e.beta_1.tilde.hat == (1 / e.mu_1.tilde.hat)){
          diagnostic.error  = "All treated individuals are Y1 cases"
          beta_1.hat        = NA
          Var_theta.hat     = matrix(NA, nrow = 4, ncol = 4)
        }else{
            if(e.beta_1.tilde.hat == 0){
              diagnostic.error  = "No treated Y1 cases"
              beta_1.hat        = NA
              Var_theta.hat     = matrix(NA, nrow = 4, ncol = 4)
            }else{
              if (isPlusInf(e.beta_1.tilde.hat)){
                diagnostic.error  = "No untreated Y1 cases"
                beta_1.hat        = NA
                Var_theta.hat     = matrix(NA, nrow = 4, ncol = 4)
              }else{
                  if(e.beta_2.tilde.hat == 0){
                    diagnostic.error  = "No treated Y2 cases"
                    beta_1.hat        = NA
                    Var_theta.hat     = matrix(NA, nrow = 4, ncol = 4)
                  }else{
                    if (isPlusInf(e.beta_2.tilde.hat)){
                      diagnostic.error  = "No untreated Y2 cases"
                      beta_1.hat        = NA
                      Var_theta.hat     = matrix(NA, nrow = 4, ncol = 4)
                    }else{
                      diagnostic.error      = NA
                      # using the estimated non-zero treatment effect on the secondary outcome to reduce confounding bias in the estimated treatment effect on the primary outcome
                      beta_1.hat            = beta_1.tilde.hat - beta_2.tilde.hat # de-biased treatment effect on Y1.
                      # Sandwich estimation of the variance
                      p1.hat                = e.mu_1.tilde.hat * (e.beta_1.tilde.hat * T + 1 - T)
                      p2.hat                = e.mu_2.tilde.hat * (e.beta_2.tilde.hat * T + 1 - T)
                      U.hat                 = rbind((Y1 - p1.hat) / (1 - p1.hat), T * (Y1 - p1.hat) / (1 - p1.hat), (Y2 - p2.hat), T * (Y2 - p2.hat))
                      drond_U.hat           = matrix(0, nrow = 4, ncol = 4)
                      drond_U.hat[1,1:2]    = c(sum(p1.hat * (Y1 - 1) / (1 - p1.hat)^2), sum(T * p1.hat * (Y1 - 1) / (1 - p1.hat)^2)) / n
                      drond_U.hat[2,1:2]    = c(sum(T * p1.hat * (Y1 - 1) / (1 - p1.hat)^2), sum(T^2 * p1.hat * (Y1 - 1) / (1 - p1.hat)^2)) / n
                      drond_U.hat[3,3:4]    = c(- sum(p2.hat), - sum(T * p2.hat)) / n
                      drond_U.hat[4,3:4]    = c(- sum(T * p2.hat), - sum(T^2 * p2.hat)) / n
                      inv                   = try(solve(drond_U.hat), silent=TRUE)
                      if (checkForTryError(inv)) {
                        inv              <- matrix(NA, nrow = 4, ncol = 4)
                        diagnostic.error <- "Sandwich covariance matrix is not invertible"
                      }
                      Var_theta.hat         = 1 / n * inv %*% (U.hat %*% t(U.hat) / n) %*% t(inv)
                  }}}}}}
  
  return(list(beta_1.hat = beta_1.hat, sd.beta_1.hat = sqrt(Var_theta.hat[2,2] + Var_theta.hat[4,4] - 2 * Var_theta.hat[2,4]), diagnostic.error = diagnostic.error))
}

## JointMH:

## Arguments/Inputs:
# Y1: outcome of interest, assumed to be binary.
# Y2: secondary outcome.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T and Y2 and T.
# W: categorical variable (or stratum indicator) used for adjustment.

## Values/Outputs:
# beta_1.hat: "de-biased" estimate of the treatment effect on Y1.
# sd.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.
# diagnostic.error: error diagnostic, if estimation fails.

## Details:
# JointMH estimates the primary outcome log-relative risk from stratification on the observed 
# covariate with Mantel-Haenszel-type weights, and uses data on the secondary outcome to reduce bias
# due to unmeasured confounders.

JointMH = function(Y1 = Y1, Y2 = Y2, T = T, W = NULL){
  if(is.null(W)){
    # in the absence of measured covariate, perform the joint approach with no covariates (JointNC)
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    n                 = length(Y1) # sample size
    
    if(length(W) != n){ # if W contains stratum indicator vectors, use the corresponding categorical variable instead
      if((nrow(W) != n)|(sum((rowSums(W == 0) + rowSums(W == 1)) == ncol(W)) != n)|(sum(rowSums(W) <= 1) != n)){
        return(list(beta_1.hat = NA, sd.beta_1.hat = NA, diagnostic.error = "W should be either a stratum indicator or a categorical variable", n.strata=NA))
      }else{
        if(sum(rowSums(W) == 1) != n){
          W = cbind(W, 1 - rowSums(W))
        }
        WW = as.factor(apply(W, 1, function(x) which(x == 1)))
        W = WW
      }
    }
    
    W                 = as.factor(W)
    levelsW           = levels(W)
    K                 = length(levelsW)
    T1                = T == 1
    T0                = T == 0
    Num_Y1 = Denom_Y1 = Num_Y2 = Denom_Y2 = pl = pl1 = pl0 = rep(NA, K)
    
    # Parameter estimation
    for(i in 1:K){
      l               = levelsW[i]
      Wl              = W == l
      WlT1            = Wl & T1
      WlT0            = Wl & T0
      nl              = sum(Wl) # number of individuals in the stratum
      nl1             = sum(WlT1) # number of individuals who received the treatment in the stratum
      nl0             = sum(WlT0) # number of individuals who did not received the treatment in the stratum
      pl[i]           = nl0 * nl1 / nl
      pl1[i]          = nl0 / nl 
      pl0[i]          = nl1 / nl
      Num_Y1[i]       = sum(Y1[WlT1]) # number of treated cases in the stratum, for the primary outcome Y1
      Denom_Y1[i]     = sum(Y1[WlT0]) # number of untreated cases in the stratum, for the primary outcome
      Num_Y2[i]       = sum(Y2[WlT1]) # number of treated cases in the stratum, for the secondary outcome Y2
      Denom_Y2[i]     = sum(Y2[WlT0]) # number of untreated cases in the stratum, for the secondary outcome
    }

    tmp <- is.na(pl)
    if (any(tmp)) {
      tmp      <- !tmp
      pl       <- pl[tmp]
      pl1      <- pl1[tmp]
      pl0      <- pl0[tmp]
      Num_Y1   <- Num_Y1[tmp]
      Denom_Y1 <- Denom_Y1[tmp] 
      Num_Y2   <- Num_Y2[tmp]
      Denom_Y2 <- Denom_Y2[tmp] 
    }
    
    tmp <- pl == 0
    m   <- sum(tmp)  
    if (m) {
      if (m == K){
        return(list(beta_1.hat = NA, sd.beta_1.hat = NA, diagnostic.error = "Check repartition of treated and untreated individuals in each stratum", n.strata=0))
      }
      # if in a given stratum all individuals are either treated or untreated, this stratum does not contribute
      tmp      <- !tmp
      pl       <- pl[tmp]
      pl1      <- pl1[tmp]
      pl0      <- pl0[tmp]
      Num_Y1   <- Num_Y1[tmp]
      Denom_Y1 <- Denom_Y1[tmp] 
      Num_Y2   <- Num_Y2[tmp]
      Denom_Y2 <- Denom_Y2[tmp] 
    }

    beta_1.tilde.hat  = log(sum(Num_Y1 * pl1) / sum(Denom_Y1 * pl0)) # estimated treatment effect on the primary outcome via stratification on W with MH weights
    beta_2.tilde.hat  = log(sum(Num_Y2 * pl1) / sum(Denom_Y2 * pl0)) # estimated treatment effect on the secondary outcome via stratification on W with MH weights
    # using the estimated non-zero treatment effect on the secondary outcome to reduce confounding bias in the estimated treatment effect on the primary outcome
    beta_1.hat        = beta_1.tilde.hat - beta_2.tilde.hat # de-biased treatment effect on Y1.
    
        if (isMinusInf(beta_1.tilde.hat)){
        diagnostic.error = "No treated Y1 cases in the contributing strata"
        beta_1.hat        = NA
        Var_theta.hat    = matrix(NA, nrow = 2, ncol = 2)
        }else{
          if (isPlusInf(beta_1.tilde.hat)){
            diagnostic.error = "No untreated Y1 cases in the contributing strata"
            beta_1.hat      = NA
            Var_theta.hat  = matrix(NA, nrow = 2, ncol = 2)
          }else{
            if(is.na(beta_1.tilde.hat)){
              diagnostic.error = "No Y1 cases in the contributing strata"
              beta_1.hat      = NA
              Var_theta.hat  = matrix(NA, nrow = 2, ncol = 2)
            }else{
              if (isMinusInf(beta_2.tilde.hat)){
                diagnostic.error = "No treated Y2 cases in the contributing strata"
                beta_1.hat        = NA
                Var_theta.hat    = matrix(NA, nrow = 2, ncol = 2)
              }else{
                if (isPlusInf(beta_2.tilde.hat)){
                  diagnostic.error = "No untreated Y2 cases in the contributing strata"
                  beta_1.hat      = NA
                  Var_theta.hat  = matrix(NA, nrow = 2, ncol = 2)
                }else{
                  if(is.na(beta_2.tilde.hat)){
                    diagnostic.error = "No Y2 cases in the contributing strata"
                    beta_1.hat      = NA
                    Var_theta.hat  = matrix(NA, nrow = 2, ncol = 2)
                  }else{
                    diagnostic.error = NA
                    # Sandwich variance estimation
                    U.hat             = rbind(pl1 * Num_Y1 - exp(beta_1.tilde.hat) * pl0 * Denom_Y1, pl1 * Num_Y2 - exp(beta_2.tilde.hat) * pl0 * Denom_Y2)
                    bar_U.hat         = rowMeans(U.hat)
                    drond_U.hat       = matrix(NA, nrow = 2, ncol = 2)
                    drond_U.hat[1,]   = rbind(- exp(beta_1.tilde.hat) * sum(Denom_Y1 * pl0), 0) / K
                    drond_U.hat[2,]   = rbind(0, - exp(beta_2.tilde.hat) * sum(Denom_Y2 * pl0)) / K
                    inv               = try(solve(drond_U.hat), silent=TRUE)
                    if (checkForTryError(inv)) {
                      inv              <- matrix(NA, nrow = 2, ncol = 2)
                      diagnostic.error <- "Sandwich covariance matrix is not invertible" 
                    }
                    Var_theta.hat     = 1 / (K) * inv %*% ((U.hat - bar_U.hat) %*% t(U.hat - bar_U.hat) / (K - 1)) %*% t(inv)
              }}}}}}
    
    return(list(beta_1.hat = beta_1.hat, sd.beta_1.hat = sqrt(Var_theta.hat[1,1] + Var_theta.hat[2,2] - 2*Var_theta.hat[1,2]), 
                diagnostic.error = diagnostic.error, n.strata=length(pl)))
  }
}


## SSJoint:

## Arguments/Inputs:
# Y1: outcome of interest, assumed to be binary.
# T: treatment of interest, assumed to be binary. A log-link is used to relate Y1 and T.
# Y2: secondary outcome.
# W: categorical covariates (= stratification variable with only a few strata) used for adjustment.

## Values/Outputs:
# beta_1.hat: "de-biased" estimate of the treatment effect on Y1.
# sd.beta_1.hat: sandwich estimate for the standard deviation of beta_1.hat.
# diagnostic.error: error diagnostic, if estimation fails.

## Details:
# SSJoint estimates the primary outcome log-relative risk within each stratum of the observed 
# covariate, and uses data on the secondary outcome to reduce bias due to unmeasured confounders. 
# These stratum-specific estimates are then efficiently combined through a weighted combination,
# where weights are inversely proportional to their respective variance.

SSJoint = function(Y1 = Y1, Y2 = Y2, T = T, W = NULL, minObsPerStrata=0){

  if(is.null(W)){
    # in the absence of measured covariates, perform the joint approach with no covariates (JointNC)
    JointNC(Y1 = Y1, Y2 = Y2, T = T)
  }else{
    n               = length(Y1) # sample size
    
    if(length(W) != n){ # if W contains stratum indicator vectors, use the corresponding categorical variable instead
      if((nrow(W) != n)|(sum((rowSums(W == 0) + rowSums(W == 1)) == ncol(W)) != n)|(sum(rowSums(W) <= 1) != n)){
        return(list(beta_1.hat = NA, sd.beta_1.hat = NA, diagnostic.error = "W should be either a stratum indicator or a categorical variable", n.strata=NA))
      }else{
        if(sum(rowSums(W) == 1) != n){
          W = cbind(W, 1 - rowSums(W))
        }
        WW = as.factor(apply(W, 1, function(x) which(x == 1)))
        W = WW
      }
    }
    
    W               = as.factor(W)
    levelsW         = levels(W)
    nlevelsW        = length(levelsW)
    beta_1.hat      = rep(NA, nlevelsW)
    sd.beta_1.hat   = beta_1.hat
    n.lowObs        = 0
    for(i in 1:nlevelsW){
      # apply function JointNC within each stratum of W
      l                = levelsW[i]
      tmp              = W == l
      if (sum(tmp) >= minObsPerStrata) {
        est              = JointNC(Y1 = Y1[tmp], Y2 = Y2[tmp], T = T[tmp])
        beta_1.hat[i]    = est$beta_1.hat
        sd.beta_1.hat[i] = est$sd.beta_1.hat
      } else {
        n.lowObs         = n.lowObs + 1
      }
    }
    n.strata0 <- length(beta_1.hat)
    tmp       <- is.na(beta_1.hat) | is.na(sd.beta_1.hat)
    m         <- sum(tmp)
    if (m){
      if (m == nlevelsW){
        return(list(beta_1.hat = NA, sd.beta_1.hat = NA, diagnostic.error = "Check the distribution of treated and untreated Y1 and Y2 cases in each stratum", n.strata=0))
      }
      beta_1.hat    <- beta_1.hat[!tmp]
      sd.beta_1.hat <- sd.beta_1.hat[!tmp]
    }
    n.strata        <- length(beta_1.hat)
    n.drop          <- n.strata0 - n.strata - n.lowObs
    if (n.drop) warning(paste0(n.drop, " strata have non-finite estimates and will be removed from the SS-Joint calculations."))

    pl              = 1 / sd.beta_1.hat^2 # weights inversely proportional to the variance of the estimate in each stratum
    sum_pl          = sum(pl)
    beta_1.hat      = sum(beta_1.hat * pl) / sum_pl # parameter estimation: weighted average of the stratum-specific estimates
    sd.beta_1.hat   = sqrt(sum((sd.beta_1.hat * pl)^2)) / sum_pl # variance estimation
    
    return(list(beta_1.hat = beta_1.hat, sd.beta_1.hat = sd.beta_1.hat, diagnostic.error = NA, n.strata=n.strata))
  }
}

