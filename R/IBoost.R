get.permute.lambda <- function(X, Y, alpha=1, permN=100, standardize=T, family="cox", quan=0.5, offset=rep(0,nrow(X))) {
  n <- nrow(X)
  all.lambda <- rep(0, permN)
  for (k in 1:permN) {
    perm <- sample(n)
    all.lambda[k] <- get.lambda.series(X, Y[perm], alpha=alpha, standardize=standardize, family=family, offset=offset[perm])[1]
  }
  c(quantile(all.lambda,quan))
}


#' Regression Parameter Estimation Using I-Boost
#'
#' Takes in data matrices and perform I-Boost-CV or I-Boost-Permutation
#' @param X The matrix of all predictors. Each row represents a subject (with total of n rows) and each column represents a feature (with total of d columns). Each column of X should be standardized.
#' @param Y The survival time, represented by a \code{Surv} object.
#' @param data.type The list of indices representing the types of the predictors. Each element of the list is a vector of integers between 1 and d that corresponds to the column numbers of a type of predictors in X. The indices should be non-overlapping.
#' @param method The version of I-Boost to be used; “permute” for I-Boost-Permutation, and “CV” for I-Boost-CV.
#' @param iter.max The maximum number of iterations.
#' @param v The penalty factor for the current estimate at each iteration.
#' @param m.stop The stopping criterion; if the parameters are not updated consecutively for m.stop iterations, then the algorithm terminates.
#' @param alpha.series The set of values to be considered for the second tuning parameter in elastic net (\eqn{\alpha}) in cross-validation; only applicable when \code{method=”CV”}.
#' @param n.fold The number of folds in cross-validation for the selection of tuning parameters; only applicable when \code{method=”CV”}.
#' @param permN The number of permutation data sets; only applicable when \code{method="permute"}.
#' @param seed 	The initial random seed for partitioning the data set for cross-validation or generating permutation data sets.
#' @return A list of two elements, \code{beta} and \code{iter.no}. \code{beta} is a vector of estimated regression parameters. \code{iter.no} is the number of iterations used in the estimation.
#' @export
#' @import glmnet survival stats
#' @seealso \code{glmnet}, \code{survival}.
IBoost <- function(X,Y,data.type,method="permute",iter.max=2000,v=ifelse(method=="permute",0.01,0.01),m.stop=5,alpha.series=c(0.05,seq(0.1,1,by=0.1)),n.fold=5,permN=100,seed=12345678) {

  n <- nrow(X)
  p <- ncol(X)
  n.type <- length(data.type)
  f_m <- rep(0,n)
  beta.iter <- matrix(0,nrow=p,ncol=iter.max)
  rownames(beta.iter) <- colnames(X)
  k <- 1
  cum <- 0

  while (k <= iter.max & cum <= m.stop) {
    feature.beta <- matrix(0,p,n.type)
    min.cvm.type <- 1e8
    for (j in 1:n.type) {
      pos <- data.type[[j]][which(apply(X[,data.type[[j]]],2,sd)>0)]
      set.seed(seed+417*k)
      if (method == "permute") {
        best.alpha <- 1
        lambda <- get.lambda.series(X[,pos],Y,alpha=1,offset=as.numeric(f_m),standardize=FALSE)
        lambda.permute <- get.permute.lambda(X[,pos], Y, offset=as.numeric(f_m), quan=0.5, permN=permN, standardize=FALSE)
        best.lambda <- c(lambda[lambda > lambda.permute], lambda.permute)
        fit.en <- glmnet(X[,pos], Y, family="cox", offset=as.numeric(f_m), alpha=best.alpha, lambda=best.lambda)
        feature.beta[pos,j] <- fit.en$beta[,length(best.lambda)]
      } else if (method == "CV") {
        foldid <- rep(cumsum(rep(1,n.fold)), length=nrow(X))
        foldid <- sample(foldid, replace=F)
        min.cvm = 1e8
        for (a in alpha.series) {
          cvm <- rep(0, 51)
          lambda <- get.lambda.series(X[,pos],Y,alpha=a, offset=as.numeric(f_m),seq.length=50)
          lambda <- c(lambda[1] + (lambda[1]+lambda[2])/10, lambda)
          cv.fit.en <- cv_glmnet(X[,pos], Y, family="cox",
                                        offset=as.numeric(f_m), alpha=a, foldid=foldid, error="cox", setlambda=TRUE, lambda=lambda)
          short <- min(length(cvm),length(cv.fit.en$lambda))
          cvm <- cvm[1:short] + cv.fit.en$cvm[1:short]
          if (min(cvm) < min.cvm) {
            best.alpha = a
            best.lambda = cv.fit.en$lambda[1:which.min(cvm)]
            min.cvm <- min(cvm)
          }
        }
        if (min.cvm < min.cvm.type) {
          best.alpha.type <- best.alpha
          best.lambda.type <- best.lambda
          min.cvm.type <- min.cvm
          best.pos <- pos
        }
      }
    }
    if (method == "permute") {
      like <- eval.fit(feature.beta,list(X=X,Y=Y),offset=f_m)
      if (n.type == 1) this.beta <- feature.beta else this.beta <- feature.beta[,which.max(like)]
    } else if (method == "CV") {
      fit.en <- glmnet(X[,best.pos], Y, family="cox", offset=as.numeric(f_m), alpha=best.alpha.type, lambda=best.lambda.type)
      this.beta <- rep(0,p)
      this.beta[best.pos] <- fit.en$beta[,length(fit.en$lambda)]
    }
    beta.est <- this.beta * v
    Xb <- rep(0,nrow(X))
    for (jj in 1:nrow(X)) Xb[jj] <- sum(X[jj,] * beta.est)
    names(Xb) <- rownames(X)
    f_m <- f_m + Xb
    if (k > 1) {
      beta.iter[,k] <- beta.iter[,k-1] + beta.est
    } else {
      beta.iter[,1] <- beta.est
    }
    k <- k + 1
    if (sum(abs(beta.est)) == 0) cum <- cum + 1 else cum <- 0
  }
  beta.out <- beta.iter[,k-1]

  list(beta=beta.out,iter.no=k-1)

}
