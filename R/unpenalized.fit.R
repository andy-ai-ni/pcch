
unpenalized.fit=function(time, status, covar, sub.coh, alpha, full.size, threshold, maxit){

	epsilon=threshold  ## the accuracy
	d0=ncol(covar)
	beta.est=rep(0, d0)
	n=length(time)
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
	
	step=0
	delta0=1

	while(delta0>epsilon & step<=maxit){
  		step=step+1
  		dl=ddl=0
  
  		xbeta=covar.cch.comp%*%beta.est
  		expx=exp(xbeta)  ## n.cch.comp by 1 dimension
  		if(Inf %in% expx){
  		  return(1)
  		}else{
  		temp1=t(expx)%*%wI0.cch.comp+((t(expx)%*%wI0.cch.comp)==0)  ## n.cch.comp*S(0), 1 by N dimension
  		w=(expx%*%t(rep(1,N)))/(rep(1,n.cch.comp)%*%temp1)  ## n.cch.comp by N
  		temp20=t(w*wI0.cch.comp)  ## N by n
  		temp2=temp20%*%covar.cch.comp  ## S(1)/S(0) N by d matrix, d=number of covariates
  		temp3=colSums(temp20)
		  temp4=matrix(0, ncol=ncol(covar.cch.comp), nrow=nrow(covar.cch.comp))
  		temp4[status.cch.comp==1,]=covar.cch.comp[status.cch.comp==1,]-temp2  ## "Z - Z-bar" 1st derivative, N by d
		  dl=dl+colSums(temp4) 
  		ddl=ddl-t(covar.cch.comp)%*%diag(temp3)%*%covar.cch.comp+t(temp2)%*%temp2  ## S(2)/S(0) - (S(1)/S(0))@2, d by d
		
		  singular.ddl=svd(ddl)$d
    	if((max(singular.ddl)/min(singular.ddl))<=10^4){
     		iddl=solve(ddl)
		  }else{
     		# stop("Hessian matrix is near singular. Consider changing the design matrix or increasing sample size.")
		    iddl=ginv(ddl)
		  }
  		beta.est=beta.est-iddl%*%dl
  		delta0=max(abs(iddl%*%dl))}
	}  ## end while loop


	### To compute SE of beta

	M=matrix(0, ncol=d0, nrow=n.cch.comp)
  V1=matrix(0, ncol=n.cch.comp, nrow=d0)
  Mk=matrix(0, ncol=d0, nrow=n.cch.comp)
  xbeta=covar.cch.comp%*%beta.est
  expx=exp(xbeta)  ## n.cch.comp by 1 dimension
  temp1=t(expx)%*%wI0.cch.comp+((t(expx)%*%wI0.cch.comp)==0)  ## n.cch.comp*S(0), 1 by N dimension
  w=(expx%*%t(rep(1,N)))/(rep(1,n.cch.comp)%*%temp1)  ## n.cch.comp by N
  temp20=t(w*wI0.cch.comp)  ## N by n
  temp20nowt=t(w*I0.cch.comp)  ## N by n.cch.comp
  temp2=temp20%*%covar.cch.comp  ## S(1)/S(0) N by d matrix, d=number of covariates
  temp4=matrix(0, ncol=ncol(covar.cch.comp), nrow=nrow(covar.cch.comp))
  temp4[status.cch.comp==1,]=covar.cch.comp[status.cch.comp==1,]-temp2  ## "Z - Z-bar" 1st derivative, N by d  
  
  for(i in 1:n.cch.comp){
    if(sub.coh.cch.comp[i]!=0){
      mtemp1=temp20nowt[,i]%*%t(rep(1,d0))
      mtemp2=rep(1,N)%*%t(as.matrix(covar.cch.comp[i,]))-temp2
      mtemp3=colSums(mtemp1*mtemp2)
      Mk[i,]=status.cch.comp[i]*temp4[i,]-mtemp3
    }
  }
  M=M+Mk

  # E1 = matrix(0, nrow=N, ncol=d0)
  # for(i in 1:n.cch.comp){
  #   if(sub.coh.cch.comp[i]!=0){
  #     ritemp1=rep(1,N)%*%t(as.matrix(covar.cch.comp[i,]))-temp2
  #     ritemp2=(I0.cch.comp[i,]*expx[i])%*%t(rep(1,d0))
  #     Ri=ritemp1*ritemp2
  #     E1=E1+((1-status.cch.comp[i])*sub.coh.cch.comp[i]/sum(sub.coh))*Ri
  #   }
  # }
  # E2=colSums(((1-status.cch.comp)%*%t(rep(1,N)))*I0.cch.comp)/full.size

  # vtp5=E2%*%t(rep(1,d0))
  vtp7=matrix(0, nrow=n.cch.comp, ncol=d0)
  for(j in 1:N){
    vtp2=matrix(1,nrow=n.cch.comp,ncol=d0)*temp1[j]
    # vtp3=rep(1,n.cch.comp)%*%t(E1[j,])
    # vtp4=I0.cch.comp[,j]%*%t(rep(1,d0))
    # vtp6=rep(1,n.cch.comp)%*%t(vtp5[j,])+(rep(1,n.cch.comp)%*%t(vtp5[j,])==0)
    Rj=matrix(0,nrow=n.cch.comp,ncol=d0)
    for(i in 1:n.cch.comp){
      if(sub.coh.cch.comp[i]!=0){
        rjtp1=covar.cch.comp[i,]-temp2[j,]
        rjtp2=I0.cch.comp[i,j]%*%expx[i]%*%t(rep(1,d0))
        Rj[i,]=rjtp1*rjtp2
      }
    }
    # vtp7=vtp7+(Rj-vtp3*vtp4/vtp6)/vtp2
    vtp7=vtp7+Rj/vtp2
  }
  vtp8=(1-status.cch.comp)%*%t(rep(1,d0))
  vtp9=vtp8*vtp7
  V1=V1+t(vtp9)
  
  Q=V=matrix(0,nrow=d0,ncol=d0)
  for(i in 1:n.cch.comp){
    Q=Q+M[i,]%*%t(M[i,])/sum(sub.coh)
    V=V+V1[,i]%*%t(V1[,i])/sum(sub.coh)
  }
  A=-ddl/full.size
  beta.std=sqrt(diag(solve(A)%*%(Q+V*(full.size-sum(sub.coh))/sum(sub.coh))%*%solve(A))/full.size)
	lkhd=sum(xbeta[status.cch.comp==1])-sum(log(colSums((expx%*%t(rep(1,N)))*wI0.cch.comp)))
	return(list(beta.est, beta.std, lkhd))
}



	




















