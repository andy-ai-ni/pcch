
penalized.fit=function(time, status, covar, sub.coh, alpha, full.size, no.pen.covar.position, lambda, beta.ini, std.ini, threshold, maxit){

	d0=ncol(covar)
	n=length(time)
	nsub=sum(sub.coh)
	t0=time[status==1]
	N=length(t0)

	I0=ifelse(time%*%t(rep(1,N))>=rep(1,n)%*%t(t0), 1, 0)  ## risk sets (n by N)
	
	inv.alpha=t(rep(1/alpha,N))
	weight=status%*%t(rep(1,N))+sub.coh*(1-status)%*%inv.alpha  ## weight matrix, n by N
	wI0=weight*I0  ## weighted risk sets (n by N)
	
	covar.cch=as.matrix(covar[sub.coh==1 | status==1,])  ## keep only covariates of subjects in the case-cohort
	comp=1*(rowSums(is.na(covar))==0)  ## indicator for complete cases
	covar.cch.comp=as.matrix(covar[(sub.coh==1 | status==1) & comp==1,])
	time.cch.comp=time[(sub.coh==1 | status==1) & comp==1]
	status.cch.comp=status[(sub.coh==1 | status==1) & comp==1]
	sub.coh.cch.comp=sub.coh[(sub.coh==1 | status==1) & comp==1]
	I0.cch.comp=I0[(sub.coh==1 | status==1) & comp==1,]
	wI0.cch.comp=wI0[(sub.coh==1 | status==1) & comp==1,]
	n.cch.comp=length(time.cch.comp)

	epsilon=threshold  ## the accuracy 
	a=3.7
	th0=lambda*std.ini
	no.pen.list=rep(0, d0)
	no.pen.list[no.pen.covar.position]=1
	th0[no.pen.list==1]=0
	
	step=0
	delta0=1
	
	beta.ini2=beta.ini
	std.ini2=std.ini

	while(delta0>epsilon & step<=maxit){
  	step=step+1
  	dl=ddl=0
		beta0=beta.ini2[abs(beta.ini2)>th0]
		no.pen.list0=no.pen.list[abs(beta.ini2)>th0]
		
		if(length(beta0)>0){
			abeta=abs(beta0)
			std0=std.ini2[abs(beta.ini2)>th0]
			## compute first derivative of the penalty
			# if(penalty=="scad"){
			# 	pbeta0=(abeta<=lambda*std0)+((a*lambda*std0-abeta)/((a-1)*lambda*std0))*((abeta>lambda*std0)*(abeta<(a*lambda*std0)))
			# }
			# if(penalty=="lasso"){
			# 	pbeta0=rep(1, length(abeta))
			# }
			# if(penalty=="alasso"){
			# 	pbeta0=1/abs(beta.ini[abs(beta.ini2)>th0])
			# }
				
			pbeta0=(abeta<=lambda*std0)+((a*lambda*std0-abeta)/((a-1)*lambda*std0))*((abeta>lambda*std0)*(abeta<(a*lambda*std0)))
			pbeta0[no.pen.list0==1]=0  ## no penalty for specified variables
			
			covar0=as.matrix(covar.cch.comp[,abs(beta.ini2)>th0])
  		expx=exp(covar0%*%beta0)  ## n by 1 dimension
  		if(Inf %in% expx){
  		  return(1)
  		}else{
  		temp1=t(expx)%*%wI0.cch.comp+((t(expx)%*%wI0.cch.comp)==0)  ## n.cch.comp*S(0), 1 by N dimension
  		w=(expx%*%t(rep(1,N)))/(rep(1,n.cch.comp)%*%temp1)  ## n.cch.comp by N
  		temp20=t(w*wI0.cch.comp)  ## N by n
  		temp2=temp20%*%covar0  ## S(1)/S(0) N by d matrix, d=number of covariates
  		temp3=colSums(temp20)
			temp4=matrix(0, ncol=ncol(covar0), nrow=nrow(covar0))
  		temp4[status.cch.comp==1,]=covar0[status.cch.comp==1,]-temp2  ## "Z - Z-bar" 1st derivative, N by d
			dl=dl+colSums(temp4)
  		ddl=ddl-t(covar0)%*%diag(temp3)%*%covar0+t(temp2)%*%temp2  ## S(2)/S(0) - (S(1)/S(0))@2, d by d
			if(length(abeta)>1){
				temp5=dl-full.size*lambda*diag(std0*pbeta0/abeta)%*%beta0
				temp6=ddl-full.size*lambda*diag(std0*pbeta0/abeta)
			}
			if(length(abeta)==1){
				temp5=dl-full.size*lambda*std0*(pbeta0/abeta)*beta0		
				temp6=ddl-full.size*lambda*std0*pbeta0/abeta
			}

			singular.ddl=svd(temp6)$d
    	if((max(singular.ddl)/min(singular.ddl))<=10^4){
     		iddl=solve(temp6)
			}else{
			  # stop("Hessian matrix is near singular. Consider changing the design matrix or increasing sample size.")
			  iddl=ginv(ddl)
			}
  		beta0=beta0-iddl%*%temp5
  		delta0=max(abs(iddl%*%temp5))
			beta.est=rep(0,d0)
			beta.est[abs(beta.ini2)>th0]=beta0
			beta.ini2=beta.est
		  }  ## end else
		}  ## end if
		else{
			delta0=epsilon/2
		}
	} ## end while loop


	
	if(step>maxit){
		converge=0
	}
	if(step<=maxit & length(beta0)==0){
		converge=1
		beta.est=rep(0, d0)
		beta.std=rep(0, d0)
		xbeta=rep(0, n.cch.comp)
		expx=exp(xbeta)
	}
	if(step<=maxit & length(beta0)>0){
		converge=1

		### To compute SE of beta

		p=length(beta0)
		M=Mk=matrix(0, ncol=p, nrow=n.cch.comp)
		V1=matrix(0, ncol=n.cch.comp, nrow=p)

		xbeta=covar0%*%beta0
  	expx=exp(xbeta)  ## n by 1 dimension
		temp1=t(expx)%*%wI0.cch.comp+((t(expx)%*%wI0.cch.comp)==0)  ## n*S(0), 1 by N dimension
		w=(expx%*%t(rep(1,N)))/(rep(1,n.cch.comp)%*%temp1)  ## n by N
		temp20nowt=t(w*I0.cch.comp)  ## N by n
		temp20=t(w*wI0.cch.comp)  ## N by n		
  	temp2=temp20%*%covar0  ## S(1)/S(0) N by p matrix, d=number of covariates
		temp4=matrix(0, ncol=p, nrow=n.cch.comp)
  	temp4[status.cch.comp==1,]=covar0[status.cch.comp==1,]-temp2  ## "Z - Z-bar" 1st derivative, N by p
		for(i in 1:n.cch.comp){
  		if(sub.coh.cch.comp[i]!=0){
      		mtemp1=temp20nowt[,i]%*%t(rep(1,p))  ## N by p
      		mtemp2=rep(1,N)%*%t(covar0[i,])-temp2  ## N by p
      		mtemp3=colSums(mtemp1*mtemp2)  ## 1 by p
      		Mk[i,]=status.cch.comp[i]*temp4[i,]-mtemp3
			}
		}
  	M=M+Mk
  		
  	vtp7=matrix(0, nrow=n.cch.comp, ncol=p)
  	for(j in 1:N){
  	  vtp2=matrix(1,nrow=n.cch.comp,ncol=p)*temp1[j]
  	  Rj=matrix(0,nrow=n.cch.comp,ncol=p)
  	  for(i in 1:n.cch.comp){
  	    if(sub.coh.cch.comp[i]!=0){
  	      rjtp1=covar0[i,]-temp2[j,]
  	      rjtp2=I0.cch.comp[i,j]%*%expx[i]%*%t(rep(1,p))
  	      Rj[i,]=rjtp1*rjtp2
  	    }
  	  }
  	  vtp7=vtp7+Rj/vtp2
  	}
  	vtp8=(1-status.cch.comp)%*%t(rep(1,p))
  	vtp9=vtp8*vtp7
  	V1=V1+t(vtp9)		

  	Q=V=matrix(0,nrow=p,ncol=p)
  	for(i in 1:n.cch.comp){
  	  Q=Q+M[i,]%*%t(M[i,])/sum(sub.coh)
  	  V=V+V1[,i]%*%t(V1[,i])/sum(sub.coh)
  	}
  	A=-ddl/full.size
  	beta.std0=sqrt(diag(solve(A)%*%(Q+V*(full.size-sum(sub.coh))/sum(sub.coh))%*%solve(A))/full.size)	
		beta.std=rep(0,d0)
		beta.std[beta.est!=0]=beta.std0
	}

	## compute effective number of parameters

	if(!is.null(ncol(ddl))){
		if((max(singular.ddl)/min(singular.ddl))<=10^4){
     	effno=sum(diag(solve(temp6)%*%ddl))
		}else{
     	effno=sum(diag(ginv(temp6)%*%ddl))
		}
	}
	else{
		effno=0
	}

	lkhd=sum(xbeta[status.cch.comp==1])-sum(log(colSums((expx%*%t(rep(1,N)))*wI0.cch.comp)))
	
	return(list(beta.est, beta.std, lkhd, effno, converge))
}






