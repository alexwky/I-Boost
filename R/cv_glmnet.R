get.lambda.series <- function(X, Y, alpha, offset=rep(0,nrow(X)), standardize=TRUE, family="cox", seq.length=100, max.min.ratio=NA, penalty.factor=rep(1,ncol(X))) {
  if (standardize) X <- apply(X, 2, FUN=function(x) if(sd(x) > 0) (x - mean(x)) / (sd(x)*sqrt(length(x)-1)) * sqrt(length(x)) else rep(0,length(x)))
  penalty.factor <- penalty.factor / sum(penalty.factor) * ncol(X)
  X <- t(t(X) / penalty.factor)
  n <- nrow(X)
  if (family == "cox") {
    X <- X[order(Y[,1],-Y[,2]),]
    offset <- offset[order(Y[,1],-Y[,2])]
    Y <- Y[order(Y[,1],-Y[,2])]
    totjump <- length(unique(Y[Y[,2] == 1,1]))
    A <- matrix(0, totjump, n)
    B <- matrix(0, totjump, n)
    all.eta <- exp(offset)
    cum.eta <- c(rev(cumsum(rev(all.eta))),0)
    eta.atrisk <- cum.eta[1]
    k <- 0
    w <- rep(0, n)
    z <- rep(0, n)
    Y.tmp <- Y[Y[,2]==1,1]
    nofail <- table(Y.tmp)[match(Y.tmp, names(table(Y.tmp)))]
    for (i in 1:n) {
      if (i == 1) temp <- 0 else temp <- Y[i-1,1]
      if (Y[i,2] == 1 & Y[i,1] > temp) {
        A[k+1,i:n] <- ((all.eta[i:n] * eta.atrisk - all.eta[i:n]^2) / eta.atrisk^2) * nofail[k+1]
        B[k+1,i:n] <- nofail[k+1] * all.eta[i:n] / eta.atrisk
        k <- k + 1
      }
      w[i] <- sum(A[1:k,i])
      z[i] <- (Y[i,2] - sum(B[1:k,i])) / (w[i] + as.integer(abs(w[i]) < 1e-8))
      eta.atrisk <- cum.eta[i+1]
    }
    tmp <- apply(X * w * z, 2, sum)
    lambda.max <- max(abs(tmp)) / (n * alpha)
  } else if (family == "binomial") {
    beta0 <- mean(4 * Y - 2)
    lambda.max <- max(apply(X, 2, FUN=function(x) abs(t(x)%*%(4*Y-2-beta0)))) * 0.25 / (n * alpha)
  } else if (family == "gaussian") {
    lambda.max <- max(apply(X, 2, FUN=function(x) abs(t(x)%*%Y))) / (n * alpha)
  }
  if (is.na(max.min.ratio)) lambda.min <- (0.0001 * (n >= ncol(X)) + 0.01 * (n < ncol(X))) * lambda.max else
    lambda.min <- max.min.ratio * lambda.max
  c(lambda.max * (lambda.min / lambda.max)^(seq(0,seq.length,length=seq.length)/seq.length))
}

eval.fit <- function(beta, reg.obj, offset=rep(0,nrow(reg.obj$X))) {
  offset <- offset[order(reg.obj$Y[,1],-reg.obj$Y[,2])]
  X <- reg.obj$X[order(reg.obj$Y[,1],-reg.obj$Y[,2]),]
  Y <- reg.obj$Y[order(reg.obj$Y[,1],-reg.obj$Y[,2])]
  Xb <- X %*% beta + offset
  at.risk <- matrix(0, nrow(X), nrow(X))
  for (i in 1:nrow(X)) at.risk[i,i:nrow(X)] <- 1
  for (i in 2:nrow(X))
    if (Y[i,1] == Y[i-1,1]) at.risk[i,] <- at.risk[i-1,]
  like <- apply(Xb - log(at.risk %*% exp(Xb)), 2, FUN=function(x,y) sum(x*y), y=Y[,2])
  c(like)
}

cum.order <- function(x,from) {
  a <- rep(0,length=length(x))
  for (i in 1:(length(x)-1))
    a[i] <- sum(x[from[i]:length(x)] < x[i])
  c(a)
}

cal.Cindex <- function(beta, reg.obj, offset=rep(0,nrow(reg.obj$X))) {
  offset <- offset[order(reg.obj$Y[,1],-reg.obj$Y[,2])]
  X <- reg.obj$X[order(reg.obj$Y[,1],-reg.obj$Y[,2]),]
  Y <- reg.obj$Y[order(reg.obj$Y[,1],-reg.obj$Y[,2])]
  Xb <- X %*% beta + offset
  km <- survfit(Surv(Y[,1],1-Y[,2])~1)
  estsurv <- stepfun(km$time, c(1,km$surv))
  GY <- estsurv(Y[,1])
  from <- cumsum(rep(1,length=length(Y[,1])))
  for (i in (length(Y[,1])-1):1)
    from[i] <- ifelse(Y[i,1] == Y[i+1,1], from[i+1], from[i])
  from <- from + 1
  denom <- sum(apply(cbind(Y[,2],GY,length(Y[,1])-from+1),1,FUN=function(x)
    c(ifelse(x[1] == 1, x[2]^(-2) * x[3], 0))))
  X.cum.order <- apply(Xb, 2, FUN=cum.order, from=from)
  num <- apply(X.cum.order, 2, FUN=function(x,y)
    sum(apply(cbind(y,x),1,FUN=function(x) c(ifelse(x[1] == 1, x[2]^(-2) * x[3], 0)))),
    y = cbind(Y[,2],GY))
  c(num / denom)
}

rescale.en <- function(select, naive=naive, alpha=alpha, lambda=lambda, all, error, offset, ...) {
  sub.fit <- glmnet(all$X[select,], all$Y[select], alpha=alpha, lambda=lambda, offset=offset[select], ...)
  df <- sub.fit$df
  if (naive) beta <- sub.fit$beta else {
    scale <- kronecker((1 + lambda*(1-alpha)/2),rep(1,ncol(all$X)))
    beta <- sub.fit$beta * as.numeric(scale)
  }
  if (error == "cox") {
    train.like <- eval.fit(beta, list(X=all$X[select,],Y=all$Y[select]), offset=offset[select])
    all.like <- eval.fit(beta, all, offset=offset)
    cvm <- all.like - train.like
    cvm[is.na(cvm)] <- -1e8
    c(list(cvm=cvm, beta=beta, df=df))
  } else if (error == "cindex") {
    c(list(cvm=cal.Cindex(beta, list(X=all[!select,],Y=all[!select]),offset=offset[!select]), beta=beta, df=df))
  }
}

#' Tuning parameter selection by cross-validation
#'
#' Performs k-fold cross-validation for glmnet in the same way as \code{cv.glmnet} in the package \code{glmnet} and outputs the cross-validation error for each fold.
#' @param X The matrix of all predictors. Each row represents a subject (with total of n rows) and each column represents a feature (with total of d columns). Each column of X should be standardized.
#' @param Y The survival time, represented by a Surv object.
#' @param foldid A vector of integers between 1 and \code{k} indicating  what fold each subject is in, where \code{k} is the number of folds.
#' @param alpha The second parameter in elastic net (\eqn{\alpha}).
#' @param error The choice of cross-validation error. \code{error="cox"} for deviance under the Cox model and \code{error="Cindex"} for the C-index.
#' @param offset A vector of offset terms.
#' @param setlambda If \code{setlambda=TRUE}, then a user-supplied vector of tuning parameters (\eqn{\lambda}) would be used; otherwise, the sequence would be chosen automatically.
#' @param lambda The user-supplied sequence of tuning parameters (\eqn{\lambda}).
#' @param ... Other arguments that can be passed to \code{glmnet}.
#' @export
#' @import glmnet survival stats
#' @seealso \code{glmnet}

cv_glmnet <- function(X, Y, foldid, alpha=1, error="cox", offset=rep(0,nrow(X)), setlambda=F, lambda=0, ...) {
  naive <- TRUE
  if (!setlambda) {
    init.fit <- glmnet(X, Y, alpha=alpha, offset=offset, ...)
    lambda <- init.fit$lambda
  }

  select.folds <- list()
  for (j in unique(foldid)) select.folds <- c(select.folds, list((foldid != j)))
  res <- lapply(select.folds, FUN=rescale.en, naive=naive, alpha=alpha, lambda=lambda, all=list(X=X,Y=Y), error=error, offset=offset, ...)

  cvm <- rep(0,length(res[[1]]$cvm))
  df <- rep(0,length(res[[1]]$df))
  all.beta <- list()
  all.cvm <- list()
  for (i in 1:length(res)) {
    short <- min(length(cvm),length(res[[i]]$cvm))
    cvm <- cvm[1:short] - res[[i]]$cvm[1:short] / length(res)
    df <- df[1:short] + res[[i]]$df[1:short] / length(res)
    all.beta <- c(all.beta, list(res[[i]]$beta))
    all.cvm <- c(all.cvm, list(-res[[i]]$cvm))
  }
  lambda <- lambda[1:length(cvm)]

  min.cvm <- min(cvm)
  lambda.min <- lambda[which.min(cvm)]
  min.pos <- as.numeric(which.min(cvm))
  c(list(cvm=cvm, lambda=lambda, df=df, min.cvm=min.cvm, lambda.min=lambda.min, min.pos=min.pos, all.beta=all.beta, all.cvm=all.cvm))
}
