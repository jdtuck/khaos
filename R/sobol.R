## sparse_khaos
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
  S[i] <- sum(object$coef[idx[[i]]]^2)
}

# Total Sobol
for (i in 1:p){
  ST[i] <- sum(object$coef[object$vars[,i] != 0]^2)/sum(object$coeff^2)
}
