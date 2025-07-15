#' @title Sparse Kahos Sensitivity Analysis
#'
#' @description Decomposes the variance of the Sparse Khaos model into variance due to main effects and total effects
#' @param object a fitted model output from the \code{sparse_khaos} function.
#' @param ... additional arguments passed
#' @details Performs analytical Sobol' decomposition for each MCMC iteration in mcmc.use (each corresponds to a MARS model), yeilding a posterior distribution of sensitivity indices.  Can obtain Sobol' indices as a function of one functional variable.
#' @return a list with two elements:
#'  \item{S}{a vector of sensitivity indices with number of entries the number of parameters.}
#'  \item{ST}{a vector of sensitivity indices with number of entries the number of parameters.}
#'
#' @keywords Sobol decomposition
#' @seealso \link{sparse_khaos} for model fitting and \link{predict.sparse_khaos} for prediction.
#' @export
#'
sobol.sparse_khaos <- function(object, ...){
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
    idx[[i]] <- which(vars[,i]==i)
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
    S[i] <- sum(coef[idx[[i]]]^2)/sum(coef^2)
  }

  # Total Sobol
  for (i in 1:p){
    ST[i] <- sum(coef[which(vars==i, arr.ind=T)[,1]]^2)/sum(coef^2)
  }

  return(list(S=S,ST=ST))
}

#' @title Sensitivity Analysis
#'
#' @description Decomposes the variance of the Sparse Khaos model into variance due to main effects and total effects
#' @param object a fitted model output
#' @export
sobol <- function(object, ...) {
  UseMethod("sobol")
}


#' @export
sobol.default <- function(object, ...) {
  cat("This is a generic function\n")
}
