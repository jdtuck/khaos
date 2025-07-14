## sparse_khaos
sobol.sparse_khaos <- function(object){
  p <- ncol(object$vars)
  ST <- rep(NA, p)
  S <- rep(NA, p)

  # generate indices
  idx1 <- object$vars != 0
  idx <- list()
  for (i in 1:p){
    idx[[i]] <- which(idx1[,i])
  }

  dims <- 1:p
  for (i in 1:p){
    tmp <- which(idx1[,i])
    tmp2 <- setdiff(dims, i)
    idxt <- c()
    for (j in tmp2){
      idxt <- c(idxt, which(idx1[,j]))
    }
    idx[[i]] <- setdiff(idx[[i]], idxt)
  }

  # Sobol
  for (i in 1:p){
    S[i] <- sum(object$beta_hat[idx[[i]]]^2)/sum(object$beta_hat^2)
  }

  # Total Sobol
  for (i in 1:p){
    ST[i] <- sum(object$beta_hat[object$vars[,i] != 0]^2)/sum(object$beta_hat^2)
  }

  return(list(S=S,ST=ST))
}


## adaptive_khaos
soboladaptive_khaos <- function(object){
  mcmc_use = 1000
  nbasis = object$nbasis[mcmc_use]
  vars = object$vars[mcmc_use,1:nbasis,]
  coef = object$beta[mcmc_use,1:nbasis]
  p <- ncol(vars)
  ST <- rep(NA, p)
  S <- rep(NA, p)

  # generate indices
  idx1 <- !is.na(vars)
  idx <- list()
  for (i in 1:p){
    idx[[i]] <- which(idx1[,i])
  }

  dims <- 1:p
  for (i in 1:p){
    tmp <- which(idx1[,i])
    tmp2 <- setdiff(dims, i)
    idxt <- c()
    for (j in tmp2){
      idxt <- c(idxt, which(idx1[,j]))
    }
    idx[[i]] <- setdiff(idx[[i]], idxt)
  }

  # Sobol
  for (i in 1:p){
    S[i] <- sum(object$coef[idx[[i]]]^2)/sum(object$coeff^2)
  }

  # Total Sobol
  for (i in 1:p){
    ST[i] <- sum(coef[which(vars==i, arr.ind=T)[,1]]^2)/sum(coef^2)
  }

  return(list(S=S,ST=ST))
}

