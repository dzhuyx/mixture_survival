
library(numDeriv)
library(MASS)
library(rootSolve)
library(Rcpp)
sourceCpp("llh_fun.cpp")

expit <- function(x) {
      return(exp(x) / (1 + exp(x)))
}

logit <- function(x) {
      return(log(x / (1-x)))
}

# change input of hazard functions here 
h <- function(t, theta1) {
      out <- theta1[1] * t^(-theta1[2])
      return(out)
}

# # log likelihood
# # all parameters in theta are on original scale
# llh_fun <- function(Y, d, X0, X1, X2, alpha, beta1, beta2) {
#       n <- length(Y)
#       X0 <- matrix(unlist(X0), nrow = n)
#       X1 <- matrix(unlist(X1), nrow = n)
#       X2 <- matrix(unlist(X2), nrow = n)
#       
#       theta0 <- exp(X0 %*% alpha)
#       theta1 <- exp(X1 %*% beta1)
#       theta2 <- exp(X2 %*% beta2)
#       
#       temp1 <- sapply(1:n, function(i) {
#             if (d[i] == 1) {
#                   out <- d[i] * log(theta0[i] + h(Y[i], c(theta1[i], theta2[i])))
#             } else {
#                   out <- 0
#             }
#             return(out)
#       })
#       temp2 <- sapply(1:n, function(i) {
#             if (d[i] == 0) {
#                   temp <- 1 - (theta0[i] * Y[i] + theta1[i] * (Y[i] ^ (1- theta2[i]) / (1 - theta2[i]) - 1/(1-theta2[i])))
#                   out <- (1 - d[i]) * log(temp)
#             } else {
#                   out <- 0
#             }
#             # if (is.nan(log(temp))) print(i)
#             return(out)
#       })
#       
#       tempcheck <- sapply(1:n, function(i) {
#             temp <- 1 - (theta0[i] * Y[i] + theta1[i] * (Y[i] ^ (1- theta2[i]) / (1 - theta2[i]) - 1/(1-theta2[i])))
#             return(temp)
#       })
#       
#       if (any(temp2 > 0) | any(is.na(temp2)) | any(tempcheck < 0) | any(tempcheck > 1)) {
#             out <- max(-99999, -99*n)
#             return(out)
#       } else {
#             out <- sum(temp1, na.rm = T) + sum(temp2, na.rm = T)
#             out <- ifelse(is.nan(out) | is.infinite(out), max(-99999, -99*n), out)
#             return(out)
#       }
# }

llh_fun_nomis <- function(Y, d, X0, alpha) {
      n <- length(Y)
      X0 <- matrix(unlist(X0), nrow = n)
      theta0 <- exp(X0 %*% alpha)
      
      temp1 <- sapply(1:n, function(i) {
            if (d[i] == 1) {
                  out <- d[i] * log(theta0[i])
            } else {
                  out <- 0
            }
            return(out)
      })
      temp2 <- sapply(1:n, function(i) {
            if (d[i] == 0) {
                  temp <- 1 - theta0[i] * Y[i]
                  out <- (1 - d[i]) * log(temp)
            } else {
                  out <- 0
            }
            # if (is.nan(log(temp))) print(i)
            return(out)
      })
      
      tempcheck <- sapply(1:n, function(i) {
            temp <- 1 - theta0[i] * Y[i]
            return(temp)
      })
      if (any(temp2 > 0) | any(is.na(temp2)) | any(tempcheck < 0) | any(tempcheck > 1)) {
            out <- -Inf
            return(out)
      } else {
            out <- sum(temp1, na.rm = T) + sum(temp2, na.rm = T)
            out <- ifelse(is.nan(out) | is.infinite(out), -Inf, out)
            return(out)
      }
}


gof_mle_poly <- function(Y, d, X0, X1, X2, fit0, avec, t0 = 2) {
      p <- ncol(X0)
      p1 <- ncol(X1)
      p2 <- ncol(X2)
      n <- length(Y)
      
      X0 <- matrix(unlist(X0), nrow = n)
      X1 <- matrix(unlist(X1), nrow = n)
      X2 <- matrix(unlist(X2), nrow = n)
      
      alpha0 <- fit0$par[1:p]
      beta10 <- fit0$par[1:p1 + p]
      beta20 <- fit0$par[1:p2 + p + p1]
      
      theta00 <- c(exp(X0 %*% alpha0))
      theta10 <- c(exp(X1 %*% beta10))
      theta20 <- c(exp(X2 %*% beta20))

      # ttt <- theta00*Y + theta10 * (theta20^Y -1)/log(theta20)
      
      # compute invSigma for testing
      # avec <- c(0, quantile(Y[d==1], c(1:10)/10))
      Uvec <- rep(NA, length(avec) - 1)
      Utemp <- NULL
      Mmattemp <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j + 1]
            Uvec[j] <- sum(Y >= a1 & Y < a2 & d == 1)
            
            Mmattemp <- rbind(Mmattemp, do.call("c", lapply(1:n, function(i) {
                  tempfun <- function(t) {
                        out <- (theta00[i] + h(t+t0, c(theta10[i], theta20[i]))) /
                              (1 - theta00[i] * t - theta10[i] * ((t+t0) ^ (1- theta20[i]) / (1 - theta20[i]) - 
                                                                        t0^(1-theta20[i])/(1-theta20[i])))
                        return(out)
                  }
                  if (Y[i] < a1) {
                        out <- 0
                  } else {
                        out <- integrate(Vectorize(tempfun), a1, min(Y[i], a2))$value
                  }
                  return(out)
            })))
            
            Utemp <- rbind(Utemp, do.call("c", lapply(1:n, function(i) {
                  return(as.numeric(Y[i] >= a1 & Y[i] < a2 & d[i] == 1))
            })))
      }
      evec <- rowSums(Mmattemp)
      Mmat <- Utemp - Mmattemp
      
      theta <- fit0$par
      Amat <- diag(Uvec / n)
      Cmat <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j+1]
            Cmat <- cbind(Cmat, colMeans(do.call(rbind, lapply(1:n, function(i) {
                  alpha <- theta[1:p]
                  beta1 <- theta[1:p1+p]
                  beta2 <- theta[1:p2+p+p1]
                  
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  theta1 <- c(exp(X1[i, ] %*% beta1))
                  theta2 <- c(exp(X2[i, ] %*% beta2))
                  
                  ftemp <- theta0 + h(Y[i]+t0, c(theta1, theta2))
                  Stemp <- 1 - theta0*Y[i] - theta1 * ((Y[i]+t0) ^ (1- theta2) / (1 - theta2) - 
                                                             t0^(1-theta2)/(1-theta2)) 
                  
                  f_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                       
                        out <- theta0 + h(Y[i]+t0, c(theta1, theta2))
                        return(out)
                  }
                  
                  F_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                        
                        out <- theta0*Y[i] - theta1 * ((Y[i]+t0) ^ (1- theta2) / (1 - theta2) - 
                                                             t0^(1-theta2)/(1-theta2))
                        return(out)
                  }
                  
                  fder <- gradient(f = f_fun, x = theta)
                  Fder <- gradient(f = F_fun, x = theta)
                  
                  return((fder / ftemp + Fder / Stemp) * (Y[i] >= a1) * (Y[i] < a2) * d[i])
            }))))
            
      }
      Imat <- fit0$hessian / n
      
      lder <- sapply(1:n, function(i) {
            ftemp <- function(theta) {
                  alpha <- theta[1:p]
                  beta1 <- theta[1:p1+p]
                  beta2 <- theta[1:p2+p+p1]
                  
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  theta1 <- c(exp(X1[i, ] %*% beta1))
                  theta2 <- c(exp(X2[i, ] %*% beta2))
                  
                  ftemp <- theta0 + h(Y[i]+t0, c(theta1, theta2))
                  Stemp <- 1 - theta0*Y[i] - theta1 * ((Y[i]+t0) ^ (1- theta2) / (1 - theta2)-
                                                             t0^(1-theta2)/(1-theta2))
                  
                  f_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                        
                        out <- theta0 + h(Y[i]+t0, c(theta1, theta2))
                        return(out)
                  }
                  
                  F_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                       
                        out <- theta0*Y[i] - theta1 * ((Y[i]+t0) ^ (1- theta2) / (1 - theta2) - 
                                                             t0^(1-theta2)/(1-theta2))
                        return(out)
                  }
                  
                  fder <- gradient(f_fun, theta)
                  Fder <- gradient(F_fun, theta)
                  
                  if (d[i] == 1) {
                        out <- fder / ftemp
                  } else {
                        out <- -Fder / Stemp
                  }
                  return(out)
            }
            return(ftemp(fit0$par))
      })
      
      temp <- (Mmat - t(Cmat) %*% solve(Imat) %*% lder)
      Sigma <- cov(t(temp))
      invSigma <- solve(Sigma)
      
      chisqstat <- t(Uvec - evec) %*% invSigma %*% (Uvec - evec) / n
      pval <- pchisq(chisqstat, df = length(avec)-1, lower.tail = F)
      result <- list(chisqstat = chisqstat, pval = pval, Uvec = Uvec, evec = evec, Sigma = Sigma)
      return(result)
}




# extpoly
gof_mle_extpoly <- function(Y, d, X0, X1, X2, X3, fit0, avec) {
      p <- ncol(X0)
      p1 <- ncol(X1)
      p2 <- ncol(X2)
      p3 <- ncol(X3)
      n <- length(Y)
      
      X0 <- matrix(unlist(X0), nrow = n)
      X1 <- matrix(unlist(X1), nrow = n)
      X2 <- matrix(unlist(X2), nrow = n)
      X3 <- matrix(unlist(X3), nrow = n)
      
      alpha0 <- fit0$par[1:p]
      beta10 <- fit0$par[1:p1 + p]
      beta20 <- fit0$par[1:p2 + p + p1]
      beta30 <- fit0$par[1:p3 + p + p1 + p2]
      
      
      theta00 <- c(exp(X0 %*% alpha0))
      theta10 <- c(exp(X1 %*% beta10))
      theta20 <- c(exp(X2 %*% beta20))
      theta30 <- c(exp(X3 %*% beta30))
      
      
      
      # ttt <- theta00*Y + theta10 * (theta20^Y -1)/log(theta20)
      
      # compute invSigma for testing
      # avec <- c(0, quantile(Y[d==1], c(1:10)/10))
      Uvec <- rep(NA, length(avec) - 1)
      Utemp <- NULL
      Mmattemp <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j + 1]
            Uvec[j] <- sum(Y > a1 & Y <= a2 & d == 1)
            
            Mmattemp <- rbind(Mmattemp, do.call("c", lapply(1:n, function(i) {
                  tempfun <- function(t) {
                        out <- (theta00[i] + h(t+theta30[i], c(theta10[i], theta20[i]))) /
                              (1 - theta00[i] * t - theta10[i] * ((t+theta30[i]) ^ (1- theta20[i]) / (1 - theta20[i]) - 
                                                                        theta30[i]^(1-theta20[i])/(1-theta20[i])))
                        return(out)
                  }
                  if (Y[i] < a1) {
                        out <- 0
                  } else {
                        out <- integrate(Vectorize(tempfun), a1, min(Y[i], a2))$value
                  }
                  return(out)
            })))
            
            Utemp <- rbind(Utemp, do.call("c", lapply(1:n, function(i) {
                  return(Y[i] > a1 & Y[i] <= a2 & d[i] == 1)
            })))
      }
      evec <- rowSums(Mmattemp)
      Mmat <- Utemp - Mmattemp
      
      theta <- fit0$par
      Amat <- diag(Uvec / n)
      Cmat <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j+1]
            Cmat <- cbind(Cmat, colMeans(do.call(rbind, lapply(1:n, function(i) {
                  alpha <- theta[1:p]
                  beta1 <- theta[1:p1+p]
                  beta2 <- theta[1:p2+p+p1]
                  beta3 <- theta[1:p3+p+p1+p2]
                  
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  theta1 <- c(exp(X1[i, ] %*% beta1))
                  theta2 <- c(exp(X2[i, ] %*% beta2))
                  theta3 <- c(exp(X3[i, ] %*% beta3))
                  
                  
                  ftemp <- theta0 + h(Y[i]+theta3, c(theta1, theta2))
                  Stemp <- 1 - theta0*Y[i] - theta1 * ((Y[i]+theta3) ^ (1- theta2) / (1 - theta2) - 
                                                             theta3^(1-theta2)/(1-theta2)) 
                  
                  f_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        beta3 <- theta[1:p3+p+p1+p2]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                        theta3 <- c(exp(X3[i, ] %*% beta3))
                        
                        out <- theta0 + h(Y[i]+theta3, c(theta1, theta2))
                        return(out)
                  }
                  
                  F_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        beta3 <- theta[1:p3+p+p1+p2]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                        theta3 <- c(exp(X3[i, ] %*% beta3))
                        
                        out <- theta0*Y[i] - theta1 * ((Y[i]+theta3) ^ (1- theta2) / (1 - theta2) - 
                                                             theta3^(1-theta2)/(1-theta2))
                        return(out)
                  }
                  
                  fder <- gradient(f = f_fun, x = theta)
                  Fder <- gradient(f = F_fun, x = theta)
                  
                  return((fder / ftemp + Fder / Stemp) * (Y[i] > a1) * (Y[i] <= a2) * d[i])
            }))))
            
      }
      Imat <- fit0$hessian / n
      
      lder <- sapply(1:n, function(i) {
            ftemp <- function(theta) {
                  alpha <- theta[1:p]
                  beta1 <- theta[1:p1+p]
                  beta2 <- theta[1:p2+p+p1]
                  beta3 <- theta[1:p3+p+p1+p2]
                  
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  theta1 <- c(exp(X1[i, ] %*% beta1))
                  theta2 <- c(exp(X2[i, ] %*% beta2))
                  theta3 <- c(exp(X3[i, ] %*% beta3))
                  
                  ftemp <- theta0 + h(Y[i]+theta3, c(theta1, theta2))
                  Stemp <- 1 - theta0*Y[i] - theta1 * ((Y[i]+theta3) ^ (1- theta2) / (1 - theta2)-
                                                            theta3^(1-theta2)/(1-theta2))
                  
                  f_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        beta3 <- theta[1:p3+p+p1+p2]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                        theta3 <- c(exp(X3[i, ] %*% beta3))
                        
                        out <- theta0 + h(Y[i]+theta3, c(theta1, theta2))
                        return(out)
                  }
                  
                  F_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        beta3 <- theta[1:p3+p+p1+p2]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(exp(X2[i, ] %*% beta2))
                        theta3 <- c(exp(X3[i, ] %*% beta3))
                        
                        out <- theta0*Y[i] - theta1 * ((Y[i]+theta3) ^ (1- theta2) / (1 - theta2) - 
                                                             theta3^(1-theta2)/(1-theta2))
                        return(out)
                  }
                  
                  fder <- gradient(f_fun, theta)
                  Fder <- gradient(F_fun, theta)
                  
                  if (d[i] == 1) {
                        out <- fder / ftemp
                  } else {
                        out <- -Fder / Stemp
                  }
                  return(out)
            }
            return(ftemp(fit0$par))
      })
      
      temp <- (Mmat - t(Cmat) %*% solve(Imat) %*% lder)
      Sigma <- cov(t(temp))
      invSigma <- solve(Sigma)
      
      chisqstat <- t(Uvec - evec) %*% invSigma %*% (Uvec - evec) / n
      pval <- pchisq(chisqstat, df = length(avec)-1, lower.tail = F)
      result <- list(chisqstat = chisqstat, pval = pval, Uvec = Uvec, evec = evec, Sigma = Sigma)
      return(result)
}






gof_mle_exp <- function(Y, d, X0, X1, X2, fit0, avec) {
      p <- ncol(X0)
      p1 <- ncol(X1)
      p2 <- ncol(X2)
      n <- length(Y)
      
      X0 <- matrix(unlist(X0), nrow = n)
      X1 <- matrix(unlist(X1), nrow = n)
      X2 <- matrix(unlist(X2), nrow = n)
      
      alpha0 <- fit0$par[1:p]
      beta10 <- fit0$par[1:p1 + p]
      beta20 <- fit0$par[1:p2 + p + p1]
      
      theta00 <- c(exp(X0 %*% alpha0))
      theta10 <- c(exp(X1 %*% beta10))
      theta20 <- c(expit(X2 %*% beta20))
      # ttt <- theta00*Y + theta10 * (theta20^Y -1)/log(theta20)
      
      # compute invSigma for testing
      # avec <- c(0, quantile(Y[d==1], c(1:10)/10))
      Uvec <- rep(NA, length(avec) - 1)
      Utemp <- NULL
      Mmattemp <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j + 1]
            Uvec[j] <- sum(Y > a1 & Y < a2 & d == 1) + sum(Y == a2 & d == 1) 
            
            Mmattemp <- rbind(Mmattemp, do.call("c", lapply(1:n, function(i) {
                  tempfun <- function(t) {
                        out <- (theta00[i] + theta10[i] * theta20[i] ^ t) /
                              (1 - theta00[i] * t - theta10[i] / log(theta20[i]) * (theta20[i] ^ t - 1))
                        return(out)
                  }
                  if (Y[i] < a1) {
                        out <- 0
                  } else {
                        out <- integrate(Vectorize(tempfun), a1, min(Y[i], a2))$value
                  }
                  return(out)
            })))
            
            Utemp <- rbind(Utemp, do.call("c", lapply(1:n, function(i) {
                  return(Y[i] > a1 & Y[i] <= a2 & d[i] == 1)
            })))
      }
      evec <- rowSums(Mmattemp)
      Mmat <- Utemp - Mmattemp
      
      theta <- fit0$par
      Amat <- diag(Uvec / n)
      Cmat <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j+1]
            Cmat <- cbind(Cmat, colMeans(do.call(rbind, lapply(1:n, function(i) {
                  alpha <- theta[1:p]
                  beta1 <- theta[1:p1+p]
                  beta2 <- theta[1:p2+p+p1]
                  
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  theta1 <- c(exp(X1[i, ] %*% beta1))
                  theta2 <- c(expit(X2[i, ] %*% beta2))
                  
                  ftemp <- theta0 + theta1 * theta2 ^ Y[i]
                  Stemp <- 1 - theta0*Y[i] - theta1 / log(theta2) * (theta2^Y[i] - 1)
                  
                  f_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(expit(X2[i, ] %*% beta2))
                        
                        out <- theta0 + theta1 * theta2 ^ Y[i]
                        return(out)
                  }
                  
                  F_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(expit(X2[i, ] %*% beta2))
                        
                        out <- theta0*Y[i] - theta1 / log(theta2) * (theta2 ^ Y[i] - 1)
                        return(out)
                  }
                  
                  fder <- gradient(f = f_fun, x = theta)
                  Fder <- gradient(f = F_fun, x = theta)
                  
                  return((fder / ftemp + Fder / Stemp) * (Y[i] > a1) * (Y[i] <= a2) * d[i])
            }))))
            
      }
      Imat <- fit0$hessian / n
      
      lder <- sapply(1:n, function(i) {
            ftemp <- function(theta) {
                  alpha <- theta[1:p]
                  beta1 <- theta[1:p1+p]
                  beta2 <- theta[1:p2+p+p1]
                  
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  theta1 <- c(exp(X1[i, ] %*% beta1))
                  theta2 <- c(expit(X2[i, ] %*% beta2))
                  
                  ftemp <- theta0 + theta1 * theta2 ^ Y[i]
                  Stemp <- 1 - theta0*Y[i] - theta1 / log(theta2) * (theta2 ^ Y[i] - 1)
                  
                  f_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(expit(X2[i, ] %*% beta2))
                        
                        out <- theta0 + theta1 * theta2 ^ Y[i]
                        return(out)
                  }
                  
                  F_fun <- function(theta) {
                        alpha <- theta[1:p]
                        beta1 <- theta[1:p1+p]
                        beta2 <- theta[1:p2+p+p1]
                        
                        theta0 <- c(exp(X0[i, ] %*% alpha))
                        theta1 <- c(exp(X1[i, ] %*% beta1))
                        theta2 <- c(expit(X2[i, ] %*% beta2))
                        
                        out <- theta0*Y[i] - theta1 / log(theta2) * (theta2^Y[i] - 1)
                        return(out)
                  }
                  
                  fder <- gradient(f_fun, theta)
                  Fder <- gradient(F_fun, theta)
                  
                  if (d[i] == 1) {
                        out <- fder / ftemp
                  } else {
                        out <- -Fder / Stemp
                  }
                  return(out)
            }
            return(ftemp(fit0$par))
      })
      
      temp <- (Mmat - t(Cmat) %*% solve(Imat) %*% lder)
      Sigma <- cov(t(temp))
      invSigma <- solve(Sigma)
      
      chisqstat <- t(Uvec - evec) %*% invSigma %*% (Uvec - evec) / n
      pval <- pchisq(chisqstat, df = length(avec)-1, lower.tail = F)
      result <- list(chisqstat = chisqstat, pval = pval, Uvec = Uvec, evec = evec, Sigma = Sigma)
      return(result)
}



gof_mle_nomis <- function(Y, d, X0, fit0, avec) {
      alpha0 <- fit0$par
      n <- length(Y)
      theta00 <- c(exp(X0 %*% alpha0))
      X0 <- matrix(unlist(X0), nrow = n)
      
      # compute invSigma for testing
      # avec <- c(0, quantile(Y[d==1], c(1:10)/10))
      Uvec <- rep(NA, length(avec) - 1)
      Utemp <- NULL
      Mmattemp <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j + 1]
            Uvec[j] <- sum(Y >= a1 & Y < a2 & d == 1)
            
            Mmattemp <- rbind(Mmattemp, do.call("c", lapply(1:n, function(i) {
                  tempfun <- function(t) {
                        out <- theta00[i] /(1 - theta00[i] * t)
                        return(out)
                  }
                  if (Y[i] <= a1) {
                        out <- 0
                  } else {
                        out <- integrate(Vectorize(tempfun), a1, min(Y[i], a2))$value
                  }
                  return(out)
            })))
            
            Utemp <- rbind(Utemp, do.call("c", lapply(1:n, function(i) {
                  return(Y[i] >= a1 & Y[i] < a2 & d[i] == 1)
            })))
      }
      evec <- rowSums(Mmattemp)
      Mmat <- Utemp - Mmattemp
      
      theta <- fit0$par
      Amat <- diag(Uvec / n)
      Cmat <- NULL
      for (j in 1:(length(avec) - 1)) {
            a1 <- avec[j]
            a2 <- avec[j+1]
            Cmat <- cbind(Cmat, colMeans(do.call(rbind, lapply(1:n, function(i) {
                  alpha <- theta
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  
                  ftemp <- theta0
                  Stemp <- 1 - theta0*Y[i] 
                  fder <- X0[i, ] * theta0
                  Fder <-  Y[i] * X0[i, ] * theta0
                  return((fder / ftemp + Fder / Stemp) * (Y[i] >= a1) * (Y[i] < a2) * d[i])
            }))))
            
      }
      Imat <- fit0$hessian / n
      
      lder <- sapply(1:n, function(i) {
            ftemp <- function(theta) {
                  alpha <- theta
                  theta0 <- c(exp(X0[i, ] %*% alpha))
                  
                  ftemp <- theta0
                  Stemp <- 1 - theta0*Y[i]
                  
                  fder <- X0[i, ] * theta0
                  Fder <-  Y[i] * X0[i, ] * theta0
                  
                  if (d[i] == 1) {
                        out <- fder / ftemp
                  } else {
                        out <- -Fder / Stemp
                  }
                  return(out)
            }
            return(ftemp(fit0$par))
      })
      
      temp <- (Mmat - t(Cmat) %*% solve(Imat) %*% lder)
      Sigma <- cov(t(temp))
      invSigma <- solve(Sigma)
      
      chisqstat <- t(Uvec - evec) %*% invSigma %*% (Uvec - evec) / n
      pval <- pchisq(chisqstat, df = length(avec), lower.tail = F)
      result <- list(chisqstat = chisqstat, pval = pval, Uvec = Uvec, evec = evec, Sigma = Sigma)
      return(result)
}
