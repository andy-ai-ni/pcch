#' SCAD-penalized Cox proportional hazards regression with a case-cohort design.
#'
#' @param time follow-up time.
#' @param status event indicator with 0=censored, 1=event.
#' @param covar covariate matrix (only numeric values are allowed).
#' @param sub.coh indicator of random subcohort membership with 0=not in subcohort, 1=in subcohort.
#' @param alpha sampling probability of the random subcohort. Required if \code{full.size} is not provided. Ignored if \code{full.size} is provided.
#' @param full.size full cohort size. Required if \code{alpha} is not provided.
#' @param no.pen.covar.position column positions of covariates that should not receive penalty. Default is NULL.
#' @param lambda.range range of lambda values used to create solution path. Default is [0.001, 0.5].
#' @param nlambda number of lambda values within \code{lambda.range}. Default is 100.
#' @param threshold convergence criterion. Default is 10^(-4).
#' @param maxit maximum number of iterations when fitting the model. Default is 200.
#' @description
#' \code{pcch.fit} performs SCAD-penalized Cox proportional hazards regression with a case-cohort design. The input parameters can be based on a dataset
#' that contains the full cohort (with the covariates of observations outside the case-cohort being missing) or a dataset that only contains the case-cohort
#' observations. The function produces a solution path for the given grid of tuning parameter values. It also selects the best tuning parameter and the
#' corresponding models based on generalized cross-validation (GCV) and BIC-type criterion. In current implementation, only SCAD penalty is available and
#' the number of covariates should not exceeds the number of events in the dataset.
#' @return An object of class \code{pcch.fit} containing solution path and fitted regression coefficients based on generalized cross-validation (GCV)
#' and BIC-type tuning parameter selection criteria.
#' @export
#' @examples
#' Busselton.cch=na.omit(Busselton.full)
#' 
#' ## use dataset containing the full cohort
#' pcch.fit1=pcch.fit(time=Busselton.full$time, status=Busselton.full$status, covar=as.matrix(Busselton.full[,4:35]), sub.coh=Busselton.full$sub.coh.value, 
#' full.size=nrow(Busselton.full))
#' 
#' ## use dataset containing only the case-cohort
#' pcch.fit2=pcch.fit(time=Busselton.cch$time, status=Busselton.cch$status, covar=as.matrix(Busselton.cch[,4:35]), sub.coh=Busselton.cch$sub.coh.value, 
#' full.size=nrow(Busselton.full))
#' @references Ni, A., Cai, J., and Zeng, D. (2016) "Variable selection for case-cohort studies with failure time outcome". \emph{Biometrika} \strong{103}, 3, pp. 547-562.

pcch.fit <-
function(time, status, covar, sub.coh, alpha=NULL,	full.size=NULL, no.pen.covar.position=NULL, lambda.range=c(0.001, 0.5), nlambda=100,
                  threshold=10^(-4), maxit=200){
  
  if(lambda.range[1]<0 | lambda.range[2]<0){
    stop("'lambda' must be positive real number")
  }
  if(lambda.range[1]>lambda.range[2]){
    stop("'lambda.range' must be a valid range")
  }
  if(is.null(alpha) & is.null(full.size)){
    stop("Either 'alpha' or 'full.size'is required")
  }
  if(ncol(covar)>=sum(status)){
    stop("In current implementation the number of covariates should be smaller than the number of cases in the sample")
  }
  
  if(!is.null(full.size) & is.null(alpha)){
    alpha=sum(sub.coh)/full.size
  }
  if(is.null(full.size) & !is.null(alpha)){
    full.size=round(sum(sub.coh)/alpha)
  }
  if(!is.null(full.size) & !is.null(alpha)){
    alpha=sum(sub.coh)/full.size
  }
  
  beta.names=colnames(covar)
  covar=as.matrix(covar)
  
  ## get initial values from unpenalized regression
  cat("generating initial values...")
  unp.fit=unpenalized.fit(time=time, status=status, covar=covar, sub.coh=sub.coh, alpha=alpha, full.size=full.size, threshold=threshold, maxit=maxit)
  if(class(unp.fit)=="numeric"){
    beta.ini=std.ini=rep(NA, ncol(covar))
    for(i in 1:ncol(covar)){
      unp.fiti=unpenalized.fit(time=time, status=status, covar=as.matrix(covar[,i]), sub.coh=sub.coh, alpha=alpha, full.size=full.size, threshold=threshold, maxit=maxit)
      beta.ini[i]=unp.fiti[[1]]
      std.ini[i]=unp.fiti[[2]]
    }
  }else{
    beta.ini=unp.fit[[1]]
    std.ini=unp.fit[[2]]
  }
  
  ## GCV and BIC based tuning parameter selection
  cat("done \n")
  cat("generating solution path and selecting final models...\n")
  cat("lambda value number: ")
  output.cv=tuning.selection(time=time, status=status, covar=covar, sub.coh=sub.coh, alpha=alpha, full.size=full.size, 
                             no.pen.covar.position=no.pen.covar.position, beta.ini=beta.ini, std.ini=std.ini, lambda.range=lambda.range, 
                             nlambda=nlambda, threshold=threshold, maxit=maxit)
  cat("\n done \n")
  
  out=NULL
  out$call=match.call()
  out$dim=ncol(covar)
  out$fit.gcv=cbind(output.cv[[7]], output.cv[[8]])
  rownames(out$fit.gcv)=beta.names
  colnames(out$fit.gcv)=c("coefficient","se")
  out$fit.bic=cbind(output.cv[[14]], output.cv[[15]])
  rownames(out$fit.bic)=beta.names
  colnames(out$fit.bic)=c("coefficient","se")	
  out$lambda.gcv=output.cv[[3]]
  out$lambda.bic=output.cv[[10]]
  out$lambda.seq=output.cv[[1]]
  out$solution.path=output.cv[[16]]

  if(output.cv[[4]]==1){
    cat("WARNING: the lambda selected by GCV method is at the lower limit of lambda range. Consider change the lambda range.\n")
  }	
  if(output.cv[[5]]==1){
    cat("WARNING: the lambda selected by GCV method is at the upper limit of lambda range. Consider change the lambda range.\n")
  }	
  if(output.cv[[11]]==1){
    cat("WARNING: the lambda selected by BIC method is at the lower limit of lambda range. Consider change the lambda range.\n")
  }
  if(output.cv[[12]]==1){
    cat("WARNING: the lambda selected by BIC method is at the upper limit of lambda range. Consider change the lambda range.\n")
  }	
  
  class(out)="pcch.fit"  
  return(out)
}
