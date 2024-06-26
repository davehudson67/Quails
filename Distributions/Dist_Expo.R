# Truncated Exponential distribution:

## zL = lower truncation point
## zU = upper truncation point

#### Probability density function
dExpo <- function(x, r = 1.0, zL = NA, zU = NA, log = FALSE) {
  
    ## all parameters length = 1 or length(x)
    ntot <- length(x)
    if ((length(r) != 1 & length(r) != ntot)) {
      stop("Length of parameters must be = 1 or length(x)")
    }
    if ((length(zL) != 1 & length(zL) != ntot) | 
        (length(zU) != 1 & length(zU) != ntot)) {
      stop("Length of zL/zU must be = 1 or length(x)")
    }
    
    ## expand entries
    x <- cbind(x, zL, zU, r)
    colnames(x) <- NULL
    
    ## run checks (first column of checkPass corresponds to invalid
    ## inputs [so should return NA]; second column is x outside range
    ## [so should return 0 p.d.f. outside the range])
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 1)
    checkPass <- cbind(checkPass, (x[, 1] >= 0) & (is.na(x[, 2]) | x[, 1] >= x[, 2]) & (is.na(x[, 3]) | x[, 1] <= x[, 3]))
    
    ## extract vectors as needed
    zL <- x[, 2]
    zU <- x[, 3]
    r <- x[, 4]
    x <- x[, 1]
    
    ## log PDF
    logS <- rep(NA, length(x))
    logH <- rep(NA, length(x))
    zN_ind <- which(checkPass[, 1] & checkPass[, 2])
    logS[zN_ind] <- -r[zN_ind] * x[zN_ind]
    logH[zN_ind] <- log(r[zN_ind])
    logProb <- logH + logS
    
    #interval truncation
    zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2])
    if(length(zI_ind) > 0) {
        logS_zL <- -r[zI_ind] * zL[zI_ind]
        logS_zU <- -r[zI_ind] * zU[zI_ind]
        S_zL <-  exp(logS_zL)
        S_zU <- exp(logS_zU)
        logProb[zI_ind] <- logProb[zI_ind] - log(S_zL - S_zU)
    }
    
    #left truncation
    zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass[, 1] & checkPass[, 2])
    if(length(zL_ind) > 0) {
        logS_zL <- -r[zL_ind] * zL[zL_ind]
        logProb[zL_ind] <- logProb[zL_ind] - logS_zL
    }

    #right truncation
    zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2])
    if(length(zU_ind) > 0) {
        logS_zU <- -r[zU_ind] * zU[zU_ind]
        logProb[zU_ind] <- logProb[zU_ind] - log(1 - exp(logS_zU))
    }

    ## return correctly
    logProb[!checkPass[, 2]] <- -Inf
    logProb[!checkPass[, 1]] <- NA
    
    if(log) {
      return(logProb)
    } else {
      return(exp(logProb))
    }
}

## Cumulative distribution function (and survivor function)
pExpo <- function(q, r = 1.0, zL = NA, zU = NA, lower.tail = TRUE, log.p = FALSE) {
  
    ## all parameters length = 1 or length(x)
    ntot <- length(q)
    if ((length(r) != 1 & length(r) != ntot)) {
      stop("Length of parameters must be = 1 or length(q)")
    }
    if ((length(zL) != 1 & length(zL) != ntot) | 
        (length(zU) != 1 & length(zU) != ntot)) {
      stop("Length of zL/zU must be = 1 or length(q)")
    }
    
    ## expand entries
    x <- cbind(q, zL, zU, r)
    colnames(x) <- NULL
    
    ## run checks (first column of checkPass corresponds to invalid
    ## inputs [so should return NA]; second column is x below lower
    ## bound and third column is x above upper bound
    ## [so should return 0 c.d.f. before the lower bound and 1 above 
    ## the upper bound])
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 1)
    checkPass <- cbind(checkPass, (x[, 1] >= 0) & (is.na(x[, 2]) | x[, 1] >= x[, 2]))
    checkPass <- cbind(checkPass, is.na(x[, 3]) | x[, 1] <= x[, 3])
    
    ## extract vectors as needed
    zL <- x[, 2]
    zU <- x[, 3]
    r <- x[, 4]
    q <- x[, 1]
  
    ## log survivor function
    logS <-  -r * q
    S <- exp(logS)

    #interval truncation (same as right truncation for survivor function)
    zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
    if(length(zI_ind) > 0) {
        logS_zU <- -r[zI_ind] * zU[zI_ind]
        logS_zL <- -r[zI_ind] * zL[zI_ind]
        S_zU <- exp(logS_zU)
        S_zL <- exp(logS_zL)
        logS[zI_ind] <- log(S[zI_ind] - S_zU) - log(S_zL - S_zU)
    }

    #left truncation (same as no truncation for survivor function)
    zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
    if(length(zL_ind) > 0) {
        logS_zL <- -r[zL_ind] * zL[zL_ind]
        logS[zL_ind] <- logS[zL_ind] - logS_zL
    }

    #right truncation
    zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass[, 1] & checkPass[, 2] & checkPass[, 3])
    if(length(zU_ind) > 0) {
        logS_zU <- -r[zU_ind] * zU[zU_ind]
        S_zU <- exp(logS_zU)
        logS[zU_ind] <- log(S[zU_ind] - S_zU) - log(1 - S_zU)
    }

    ## return correctly
    logS[!checkPass[, 2]] <- 0
    logS[!checkPass[, 3]] <- -Inf
    logS[!checkPass[, 1]] <- NA
    
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
}

## Quantile function
qExpo <- function(p, r = 1.0, zL = NA, zU = NA, lower.tail = TRUE, log.p = FALSE) {
    
    ## all parameters length = 1 or length(x)
    ntot <- length(p)
    if ((length(r) != 1 & length(r) != ntot)) {
      stop("Length of parameters must be = 1 or length(p)")
    }
    if ((length(zL) != 1 & length(zL) != ntot) | 
        (length(zU) != 1 & length(zU) != ntot)) {
      stop("Length of zL/zU must be = 1 or length(p)")
    }
    
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    
    ## expand entries
    x <- cbind(p, zL, zU, r)
    colnames(x) <- NULL
    
    ## run checks (first column of checkPass corresponds to invalid
    ## inputs [so should return NA]; second column is p outside range
    ## [so should return NA if p not in (0, 1)]
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 1)
    checkPass <- cbind(checkPass, x[, 1] >= 0 & x[, 1] <= 1)
    
    whichPass <- which(checkPass[, 1] & checkPass[, 2])
    if(length(whichPass) > 0) {
        ## extract valid entries
        x <- x[whichPass, , drop = FALSE]
    
        ## extract vectors as needed
        zL <- x[, 2]
        zU <- x[, 3]
        r <- x[, 4]
        p <- x[, 1]
        
        out1 <- -log(1 - p) / r
        
        #interval truncation
        zI_ind <- which(!is.na(zL) & !is.na(zU))
        if(length(zI_ind) > 0) {
            S_zU <- exp(-r[zI_ind] * zU[zI_ind])
            CDF_zU <- 1 - S_zU
                  
            S_zL <-  exp(-r[zI_ind] * zL[zI_ind])
            CDF_zL <- 1 - S_zL
          
            p[zI_ind] <- (p[zI_ind] * (CDF_zU - CDF_zL)) + CDF_zL
            out1[zI_ind] <- -log(1 - p[zI_ind]) / r[zI_ind]
        }
      
        #left truncation
        zL_ind <- which(!is.na(zL) & is.na(zU))
        if(length(zL_ind) > 0) {
            S_zL <-  exp(-r[zL_ind] * zL[zL_ind])
            CDF_zL <- 1 - S_zL
                  
            p[zL_ind] <- p[zL_ind] * S_zL + CDF_zL
            out1[zL_ind] <- -log(1 - p[zL_ind]) / r[zL_ind]
        }
      
        #right truncation
        zU_ind <- which(is.na(zL) & !is.na(zU))
        if(length(zU_ind) > 0) {
            S_zU <-  exp(-r[zU_ind] * zU[zU_ind])
            CDF_zU <- 1 - S_zU
                  
            p[zU_ind] <- p[zU_ind] * CDF_zU
            out1[zU_ind] <- -log(1 - p[zU_ind]) / r[zU_ind]
        }
        out <- rep(NA, ntot)
        out[whichPass] <- out1
    } else {
        out <- rep(NA, ntot)
    }
    return(out)
}

## Random sampler
rExpo <-  function(n, r = 1.0, zL = NA, zU = NA, log = FALSE) {
    ## all parameters length = 1 or length(x)
    if(length(n) > 1) {
      n <- length(n)
      message("Since 'n' is a vector, 'length(n)' is taken to be the number of samples required.")
    }
    if ((length(r) != 1 & length(r) != n)) {
      stop("Length of parameters must be = 1 or n")
    }
    if ((length(zL) != 1 & length(zL) != n) | 
        (length(zU) != 1 & length(zU) != n)) {
      stop("Length of zL/zU must be = 1 or n")
    }
    
    ## expand entries
    x <- cbind(zL, zU, r)
    if(n == nrow(x)) {
      x <- cbind(1, x)
    } else {
      if(nrow(x) != 1) {
        stop("Error we're not catching")
      } else {
        x <- matrix(rep(c(1, x[1, ]), n), nrow = n, byrow = TRUE)
      }
    }
    colnames(x) <- NULL
    
    ## run checks (corresponds to invalid
    ## inputs [so should return NA])
    checkPass <- (is.na(x[, 2]) | is.na(x[, 3]) | x[, 2] < x[, 3])
    checkPass <- checkPass & (rowSums(x[, -c(1:3), drop = FALSE] > 0) == 1)
    
    ## extract vectors as needed
    zL <- x[, 2]
    zU <- x[, 3]
    r <- x[, 4]
    
    u <- runif(n, 0, 1)
    rs <- rep(NA, n)
    
    ## no truncation
    zN_ind <- which(checkPass)
    if(length(zN_ind) > 0) {
        rs[zN_ind] <- -log(1 - u[zN_ind]) / r[zN_ind]
    } else {
        return(rs)
    }
  
    #interval truncation  
    zI_ind <- which(!is.na(zL) & !is.na(zU) & checkPass)
    if(length(zI_ind) > 0) {
        S_zU <- exp(-r[zI_ind] * zU[zI_ind])
        CDF_zU <- 1 - S_zU
        S_zL <-  exp(-r[zI_ind] * zL[zI_ind])
        CDF_zL <- 1 - S_zL
        u[zI_ind] <- (u[zI_ind] * (CDF_zU - CDF_zL)) + CDF_zL
        rs[zI_ind] <- -log(1 - u[zI_ind]) / r[zI_ind]
    }
  
    #left truncation
    zL_ind <- which(!is.na(zL) & is.na(zU) & checkPass)
    if(length(zL_ind) > 0) {
        S_zL <-  exp(-r[zL_ind] * zL[zL_ind])
        CDF_zL <- 1 - S_zL
        u[zL_ind] <- u[zL_ind] * S_zL + CDF_zL
        rs[zL_ind] <- -log(1 - u[zL_ind]) / r[zL_ind]
    }
  
    #right truncation
    zU_ind <- which(is.na(zL) & !is.na(zU) & checkPass)
    if(length(zU_ind) > 0) {
        S_zU <-  exp(-r[zU_ind] * zU[zU_ind])
        CDF_zU <- 1 - S_zU
        u[zU_ind] <- u[zU_ind] * CDF_zU
        rs[zU_ind] <- -log(1 - u[zU_ind]) / r[zU_ind]
    }
  
    ## return correctly
    return(rs)
}

