negativeControlOutcomeAdjustment <- function(Y1, Y2, Trt, W=NULL, method=c("Joint-MH", "Joint-NC", "SS-Joint"),
                                      minObsPerStrata=20) {

  tmp             <- checkData(Y1, Y2, Trt, W)
  Y1              <- tmp$Y1
  Y2              <- tmp$Y2
  Trt             <- tmp$Trt
  W               <- tmp[["W", exact=TRUE]]
  minObsPerStrata <- check_minObsPerStrata(minObsPerStrata)
  method          <- check_method(method, W, minObsPerStrata)  
  ret             <- mainFunc(Y1, Y2, Trt, W, method, minObsPerStrata)
  ret

} # END: negativecontroladjustment

mainFunc <- function(Y1, Y2, Trt, W, method, minObsPerStrata) {

  nmethod <- length(method)
  bv      <- "beta_1.hat"
  sev     <- "sd.beta_1.hat"
  errv    <- "error.message"
  nsv     <- "n.strata"
  vec     <- rep(NA, nmethod)
  ret     <- data.frame(method, vec, vec, vec, rep("", nmethod), stringsAsFactors=FALSE)
  colnames(ret) <- c("method", bv, sev, nsv, errv)

  ii <- 0
  if ("Joint-MH" %in% method) {
    ii  <- ii + 1
    tmp <- try(JointMH(Y1=Y1, Y2=Y2, T=Trt, W=W), silent=TRUE)
    ret <- updateRetObj(ret, ii, tmp, bv, sev, errv, nsv)
  }
  if ("Joint-NC" %in% method) {
    ii  <- ii + 1
    tmp <- try(JointNC(Y1=Y1, Y2=Y2, T=Trt), silent=TRUE)
    ret <- updateRetObj(ret, ii, tmp, bv, sev, errv, nsv)
  }
  if ("SS-Joint" %in% method) {
    ii  <- ii + 1
    tmp <- try(SSJoint(Y1=Y1, Y2=Y2, T=Trt, W=W, minObsPerStrata=minObsPerStrata), silent=TRUE)
    ret <- updateRetObj(ret, ii, tmp, bv, sev, errv, nsv)
  }

  ret

} # END: mainFunc

updateRetObj <- function(ret, row, obj, bv, sev, errv, nsv) {

  if (checkForTryError(obj)) {
    ret[row, errv] <- getErrorMsg(obj)
  } else {
    ret[row, bv]   <- obj$beta_1.hat
    ret[row, sev]  <- obj$sd.beta_1.hat
    ret[row, errv] <- obj$diagnostic.error

    n <- obj[["n.strata", exact=TRUE]]
    if (length(n) && is.finite(n)) ret[row, nsv] <- n 
  }

  ret

} # END: updateRetObj

getErrorMsg <- function(obj) {

  ret <- paste(as.character(obj), collapse=" ", sep="")
  ret <- gsub("\n", " ", ret, fixed=TRUE)
  ret <- gsub("\r", " ", ret, fixed=TRUE)

  ret

} # END: getErrorMsg

checkData <- function(Y1, Y2, Trt, W, min.nobs=3) {

  # Check for the correct type of objects
  Y1  <- check_vector(Y1, "Y1")
  n0  <- length(Y1)
  Y2  <- check_vector(Y2, "Y2", len=n0)
  Trt <- check_vector(Trt, "Trt", len=n0)

  # Remove missing/non-finite values
  ok <- is.finite(Y1) & is.finite(Y2) & is.finite(Trt)
  if (!all(ok)) {
    Y1  <- Y1[ok]
    Y2  <- Y2[ok]
    Trt <- Trt[ok]  
  }
  if (all(Y2 == 0)) stop("ERROR: all elements of Y2 are 0")

  n <- length(Y1)
  if (n < min.nobs) {
    msg <- paste0("ERROR: after removing missing/non-finite values, the data contains ",
                  n, " observations")
    stop(msg)
  }
  check_binary(Y1, "Y1")
  check_Y2(Y2)
  check_binary(Trt, "Trt")
  W <- check_W(W, ok)

  list(Y1=Y1, Y2=Y2, Trt=Trt, W=W)

} # END: checkData

check_vector <- function(x, name, len=0, numeric=1) {

  n <- length(x)
  if (!n) stop(paste0("ERROR: length(", name, ") = 0"))
  x <- try(as.vector(x), silent=TRUE)
  if (checkForTryError(x)) stop(paste0("ERROR: ", name, " cannot be coereced to a vector"))
  if (!is.vector(x)) stop(paste0("ERROR: ", name, " must be a vector"))
  #if (len && (len != length(x))) stop(paste0("ERROR: ", name, " must have length ", len))
  if (len && (len != length(x))) stop(paste0("ERROR: length(", name, ") != length(Y1)"))
  if (numeric && !is.numeric(x)) stop(paste0("ERROR: ", name, " must be numeric"))
  
  x

} # END: check_vector

checkForTryError <- function(obj) {
  ret <- "try-error" %in% class(obj)
  ret
}

check_binary <- function(x, name) {

  tmp <- x %in% 0:1
  if (!all(tmp)) stop(paste0("ERROR: ", name, " contains values that are not 0 or 1"))
  NULL

}

check_Y2 <- function(Y2) {

  if (any(Y2 < 0)) stop("ERROR: Y2 contains negative values")
  tmp <- (Y2 != floor(Y2)) | (Y2 != ceiling(Y2))
  if (any(tmp)) stop("ERROR: Y2 contains non-integer values")
  NULL

}

check_W <- function(W, subset) {

  # Must be called after checking Y1, Y2, T

  n <- length(W)
  if (!n) return(NULL)
  
  # n0 is the original correct length (or number of rows)
  n0 <- length(subset)

  if (is.factor(W) || is.vector(W)) {
    if (n != n0) stop(paste0("ERROR: W must have length ", n0))
    W <- W[subset] 
  } else if (is.matrix(W) || is.data.frame(W)) {
    if (nrow(W) != n0) stop(paste0("ERROR: W must have ", n0, " rows"))
    W <- W[subset, , drop=FALSE]
    if (is.matrix(W)) W <- split(W, rep(1:ncol(W), each=nrow(W)))
  } else {
    stop("ERROR: W must be a factor, vector, matrix or data frame of categorical covariates")
  }

  W <- interaction(W, drop=TRUE)

  if (length(levels(W)) < 2) {
    warning("W contains only one stratum and will be ignored")
    W <- NULL
  }

  W

}

check_method <- function(obj, W, minObsPerStrata) {

  valid0 <- c("Joint-MH", "Joint-NC", "SS-Joint")
  valid  <- tolower(valid0) 
  if (is.null(obj)) obj <- valid
  n <- length(obj)
  msg <- "ERROR: method must be NULL or contain values 'Joint-MH', 'Joint-NC', or 'SS-Joint'"
  if (!n) stop(msg) 
  if (!is.character(obj)) stop(msg) 
  obj <- unique(tolower(removeWhiteSpace(obj)))
  tmp <- !(obj %in% valid) 
  if (any(tmp)) stop(msg) 
 
  if (is.null(W)) {
    tmp <- obj %in% "joint-nc"
    obj <- obj[tmp]
    if (!length(obj)) stop("ERROR: only method='Joint-NC' can be used when W is NULL")
  }
  tmp <- valid %in% obj
  ret <- valid0[tmp]

  if (length(W) && (minObsPerStrata > 1) && ("SS-Joint" %in% ret)) {
    tmp <- table(W) < minObsPerStrata
    m   <- sum(tmp)  
    if (all(tmp)) {
      msg <- paste0("All strata contain less than ", minObsPerStrata, " observations. SS-Joint will not be called.") 
      warning(msg)
      ret <- ret[ret != "SS-Joint"]
      if (!length(ret)) stop("ERROR: no methods can be used by default. Try changing the method or minObsPerStrata options.")
    } else if (m) {
      msg <- paste0(m, " strata contain less than ", minObsPerStrata, " observations and will be removed from the SS-Joint calculations.") 
      warning(msg)
    }
  }

  ret

} # END: check_method

# Function to remove leading/trailing white space
removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

check_minObsPerStrata <- function(x) {

  nm <- "minObsPerStrata"
  if (!length(x)) x <- 20
  if (length(x) > 1) stop(paste0("ERROR: ", nm, " must be a single non-negative number"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a single non-negative number"))
  if (x < 0) stop(paste0("ERROR: ", nm, " must be non-negative"))
  x

}

isPlusInf <- function(x) {
  is.infinite(x) && (x > 0)
}

isMinusInf <- function(x) {
  is.infinite(x) && (x < 0)
}

