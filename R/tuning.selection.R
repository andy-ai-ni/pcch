
tuning.selection=function(time, status, covar, sub.coh, alpha, full.size, no.pen.covar.position, beta.ini, std.ini, lambda.range, nlambda, threshold, maxit){

	n=length(time)
	l0=seq(lambda.range[1], lambda.range[2], length.out=nlambda)*log(log(n))
	aics=bics=rep(0, length(l0))
	betas=stds=matrix(NA, nrow=ncol(covar), ncol=nlambda)

	iter=0
	pfit=1
	while(class(pfit)=="numeric"){
  	for(k in 1:nlambda){
  		lambda0=l0[k]
  		pfit=penalized.fit(time=time, status=status, covar=covar, sub.coh=sub.coh, alpha=alpha, full.size=full.size,
  		                   no.pen.covar.position=no.pen.covar.position, lambda=lambda0, beta.ini=beta.ini, 
  		                   std.ini=std.ini, threshold=threshold, maxit=maxit)
  		if(class(pfit)=="numeric"){
  		  cat("The original initial values failed to converge. Trying different initial values \n")
  		  break
  		}
  		aics[k]=log(-pfit[[3]]/n)+2*pfit[[4]]/n   ### GCV
  		bics[k]=log(-pfit[[3]]/n)+log(n)*pfit[[4]]/n   ### BIC-type
  		betas[,k]=pfit[[1]]
  		stds[,k]=pfit[[2]]
  	  cat(paste(k, " "))		
  	}
	  beta.ini=0.9*beta.ini
	  iter=iter+1
	  if(iter==20){
	    stop("Penalized model failed to converge")
	  }
	}

	aic.min=min(aics)
	position.aic=match(aic.min, aics)
	lambda.aic=l0[position.aic]
	lambda.aic.at.min=1*(lambda.aic==min(l0))
	lambda.aic.at.max=1*(lambda.aic==max(l0))
	beta.aic=betas[,position.aic]
	std.aic=stds[,position.aic]
	
	bic.min=min(bics)	
	position.bic=match(bic.min, bics)
	lambda.bic=l0[position.bic]
	lambda.bic.at.min=1*(lambda.bic==min(l0))
	lambda.bic.at.max=1*(lambda.bic==max(l0))
	beta.bic=betas[,position.bic]
	std.bic=stds[,position.bic]

	return(list(l0, position.aic, lambda.aic, lambda.aic.at.min, lambda.aic.at.max, aic.min, beta.aic, std.aic, 
	            position.bic, lambda.bic, lambda.bic.at.min, lambda.bic.at.max, bic.min, beta.bic, std.bic, 
	            betas, stds, aics, bics))
}




		




















