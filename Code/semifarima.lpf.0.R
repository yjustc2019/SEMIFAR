###################################semifarima.lpf - extension of semifar.lpf with an addition MA-part################################
###################################content is heavily simplified, e.g no more integer differencing and###############################
###################################other adjustments#################################################################################
semifarima.lpf = function(x, p, q, mse.RANGE,b0 = NULL, pg, kn, bb, IF, par.est = "fracdiff")
{
  alpha.SEMI<-0.01
  alpha.MLE<-0.05
  IFM=c("opt","naive","var") ###("opt", "naive", "var")
  IFM=IFM[IF]
  n.DATA<-length(x)  
  korder<-pg+1            
  ko21<-2*korder+1        
  ti<-(1:n.DATA)/n.DATA   
  
  #--- iteration settings
  
  nITER<-11
  icrit<-0
  convc=1   #### b0 converges 
  
  ################ In the following Bk=1/(beta^2) with beta=int(u^(p+1)K(u))du
  ################ All formulas for K(u) may be found in Table 5.7 of Mueller (988)
  ################ Rk is the kernel constant in the variance, see e.g. Table 1 of Feng and Heiler (2009)   
  if(korder==2){Kk<-1 #### Kk=(factorial(p+1))^2/(2*(p+1))
  if(kn==1){
    Bk<-9                      #### =3^2      
    Rk<-0.5}    
  if(kn==2){
    Bk<-25                    #### =5^2
    Rk<-0.6}    
  if(kn==3){
    Bk<-49                    #### =7^2
    Rk<-0.7143} 
  if(kn==4){
    Bk<-81                    #### =9^2
    Rk<-0.8159}	
  }
  
  if(korder==4){Kk<-72 #### Kk=(factorial(p+1))^2/(2*(p+1))
  if(kn==1){
    Bk<-(35/3)^2          #### =(5*7/3)^2         
    Rk<-1.125}
  if(kn==2){
    Bk<-21^2               #### =(7*9/3)^2 
    Rk<-1.250}
  if(kn==3){
    Bk<-33^2              #### =(9*11/3)^2 
    Rk<-1.4073}	
  if(kn==4){
    Bk<-(143/3)^2     #### =(11*13/3)^2
    Rk<-1.5549}	
  }
  
  
  ###########starting bandwidth if not selected by user############
  if(is.null(b0) & korder==2){b0=0.1}
  if(is.null(b0) & korder==4){b0=0.2}
  
  ########### smoothing iteration
  for(iITER in 1:nITER)
  {
    cat("iteration =",iITER,fill=TRUE) 
    
    n.DATA<-length(x)
    ti<-(1:n.DATA)/n.DATA                        
    
    if( icrit==0 )
    {
      g0<-smooth.lpf(x, 0, pg, kn, b0, bb)
      x.DETRENDED<-x-g0 
      ############ estimation of theta ##########
      p = p
      q = q
      if(iITER==1) 
      {
        result<-mle.farima.new(x.DETRENDED, p, q, 0.05, 0)  
        
        theta<-result$theta
        eta<-result$eta
        d<-result$d
        d<-result$d
        CI<-result$CI
        BIC.opt<-result$BIC.opt
        r<-result$r
        
      }
      
      if(iITER>1)
      {
        
        result<-mle.farima.new(x.DETRENDED, p, q, 0.05, 0)
        
        p<-result$p
        q<-result$q
        theta<-result$theta
        eta<-result$eta
        d<-result$d
        CI<-result$CI
        BIC<-result$BIC
        rB<-result$r
      }
      
      ############## Cf ##################
      ar = if(p > 0){ar = theta[3:(3 + (p-1))]}else{ar = 0} ####SL
      ma = if(q > 0){ma = theta[(-1 : -(p+2))]}else{ma = 0} ####SL
      var.pred=theta[1]
      Cf <- sum(c(1, -ma))^2 / (1 - sum(ar))^2 * var.pred / (2*pi) ###SL
      #browser()
      ############ inflation factors ###SL
      if(convc==1){                    
        if(korder==2 && IFM=="opt"){b2<-(b0)**((5-2*d)/(7-2*d))}
        if(korder==2 && IFM=="naive"){b2<-(b0)**((5-2*d)/(9-2*d))}    
        if(korder==2 && IFM=="var"){b2<-(b0)**((1/2))}
        if(korder==4 && IFM=="opt"){b2<-(b0)**((9-2*d)/(11-2*d))}
        if(korder==4 && IFM=="naive"){b2<-(b0)**((9-2*d)/(13-2*d))}  
        if(korder==4 && IFM=="var"){b2<-(b0)**((1/2))}
        
        b2<-min(0.49,b2) 
        
        if(iITER<=2)
          #{gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, bb)}
        {if(pg==1){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, 0)}
          if(pg==3){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, 0)}}
        else
          #{gk<-smooth.lpf(x, pg+1, pg+2, kn, b2, bb)}
        {if(pg==1){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, 0)}
          if(pg==3){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, 0)}}
        
        index<-max(1,trunc(mse.RANGE*n.DATA)):trunc((1-mse.RANGE)*n.DATA) 
        
        gkk<-gk[index]**2
        Intgk<-sum(gkk)/n.DATA
        
        b0old<-b0
        
        ##### The function kdf.t is proposed in Feng (2007) based on K(u) in Mueller (1988)
        if(d!=0)
        {
          Vc<-2*gamma(1-2*d)*sin(pi*d) # 
          
          if(korder==2 && kn==1){
            Rk<-(1/2)^2*Vc*kdf.t(0,0,d)
          }  
          if(korder==2 && kn==2){
            Rk<-(3/4)^2*Vc*(kdf.t(0,0,d)-
                              2*kdf.t(2,0,d)+
                              kdf.t(2,2,d))
          }
          if(korder==2 && kn==3){
            Rk<-(15/16)^2*Vc*(kdf.t(0,0,d)+
                                4*kdf.t(2,2,d)+
                                kdf.t(4,4,d)-
                                4*kdf.t(2,0,d)+
                                2*kdf.t(4,0,d)-
                                4*kdf.t(4,2,d))  
          }
          if(korder==2 && kn==4){
            Rk<-(35/32)^2*Vc*(kdf.t(0,0,d)+
                                9*kdf.t(2,2,d)+
                                9*kdf.t(4,4,d)+
                                kdf.t(6,6,d)-
                                6*kdf.t(2,0,d)+
                                6*kdf.t(4,0,d)-
                                2*kdf.t(6,0,d)-
                                18*kdf.t(4,2,d)+
                                6*kdf.t(6,2,d)-
                                6*kdf.t(6,4,d))  
          }
          
          if(korder==4 && kn==1){
            Rk<-(3/8)^2*Vc*(9*kdf.t(0,0,d)-
                              30*kdf.t(2,0,d)+
                              25*kdf.t(2,2,d))
          }
          if(korder==4 && kn==2){
            Rk<-(15/32)^2*Vc*(9*kdf.t(0,0,d)+
                                100*kdf.t(2,2,d)+
                                49*kdf.t(4,4,d)-
                                60*kdf.t(2,0,d)+
                                42*kdf.t(4,0,d)-
                                140*kdf.t(4,2,d))  
          }
          if(korder==4 && kn==3){
            Rk<-(105/64)^2*Vc*(kdf.t(0,0,d)+   
                                 25*kdf.t(2,2,d)+
                                 49*kdf.t(4,4,d)+
                                 9*kdf.t(6,6,d)-
                                 10*kdf.t(2,0,d)+
                                 14*kdf.t(4,0,d)-
                                 6*kdf.t(6,0,d)-
                                 70*kdf.t(4,2,d)+
                                 30*kdf.t(6,2,d)-
                                 42*kdf.t(6,4,d))
          }  
          if(korder==4 && kn==4){
            Rk<-(315/512)^2*Vc*(9*kdf.t(0,0,d)+   
                                  400*kdf.t(2,2,d)+
                                  42^2*kdf.t(4,4,d)+
                                  36^2*kdf.t(6,6,d)+
                                  11^2*kdf.t(8,8,d)-
                                  120*kdf.t(2,0,d)+
                                  252*kdf.t(4,0,d)-
                                  216*kdf.t(6,0,d)+
                                  66*kdf.t(8,0,d)-
                                  40*42*kdf.t(4,2,d)+
                                  40*36*kdf.t(6,2,d)-
                                  40*11*kdf.t(8,2,d)-
                                  84*36*kdf.t(6,4,d)+
                                  84*36*kdf.t(8,4,d)-
                                  72*11*kdf.t(8,6,d))
          }
          
          const1<-(Bk*Kk*Rk*(1-2*mse.RANGE)*(1-2*d))**(1/(ko21-2*d))
        }
        
        else
        {
          const1<-(Bk*Kk*Rk*(1-2*mse.RANGE)*2*pi)**(1/ko21)     #### Here, delta=0, but Kk depends on p.
        }
        
        const2<-(Cf/Intgk)**(1/(ko21-2*d))
        const3<-n.DATA**((2*d-1)/(ko21-2*d))
        
        b0<-const1*const2*const3      
        b0<-min(0.49,b0) 
        b0<-max(n.DATA**(-5/7),b0) # nicht kleiner als 0.005 oder 0.05?
      }  
      
      
      cat("Selected b0=", b0, fill=TRUE)
      
      deltaITER<-0.01*max(b0,b0old)
      if( (abs(b0-b0old)<deltaITER)&(iITER>3) ){
        icrit<-1
      }
      if( (abs(b0-b0old)>=deltaITER)&(iITER==(nITER-1)) ){
        convc=0
        b0=(b0+b0old)/2   #
      }
      if( iITER==nITER ){
        cat("Alg doesn't converge. b0 is the mean of the last two values!!!", fill=TRUE) 	       
      }
    } 
    if(icrit == 1)
    {
      cat("Algorithm converged after", iITER, "iterations with selected bandwidth b0 =", b0, fill = TRUE )
      break
    }
  }
  
  if(d==0)
  {
    nu<-pi
  }
  else
  {
    nu<-2**(2*d)*gamma(1-2*d)*sin(pi*d)
    nu<-nu/(d*(2*d+1))
  }
  
  eta <- theta[c(-1)]
  d <- theta[2]
  
  
  ###### confidence intervals for theta
  CI<-conf.farima(eta,p,q,n.DATA,alpha.MLE)$CI
  ####significance test trend
  VARg0<-(n*b0)**(2*d-1)*nu*Cf 
  CRITg0<-qnorm( (1-alpha.SEMI/2) )*sqrt(VARg0) 
  SIG<-1
  
  MEAN<-mean(x)
  g0.limits<-c(MEAN-CRITg0,MEAN+CRITg0)
  
  if( max(abs(g0-MEAN)) <= CRITg0 )
  {
    SIGtext<-"Trend g0 is not significant!"
    SIG<-0
  }  
  else
  {
    SIGtext<-"Trend g0 is significant!"
  } 
  result<-list(x = x,
               mse.RANGE=mse.RANGE, 
               n = n.DATA, 
               p = p, 
               q=q, 
               eta = eta, 
               n = n.DATA,
               theta=theta,
               eta=eta,
               d=d,
               BIC=BIC,
               g0=g0,
               b0=b0,
               Cf=Cf,
               nu=nu,
               r=r,
               CI=CI,
               SIG=SIG,
               SIGtext=SIGtext)
  drop(result)
  
}
################################mle.farima###############################################
mle.farima.new <- function(x, p, q, alpha, iCI)
{
  ndata<-length(x)
  
  p = p
  q = q
  result<-mle.farima(x, p, q)
  
  eta<-result$etaEST       
  loglik<-result$loglikEST 
  
  BICpq<--2*loglik+log(ndata)*(p+q)
  BIC<-BICpq   
  eta<-eta
  p<-p          
  q<-q
  
  
  # sigma2BIC and residuals rBIC: 
  # 1. fractional filter (m and delta); 2. ARMA filter
  
  d<-eta[1]
  
  r<-x
  n<-length(r)
  b<-ar.coef.farima(n,eta[1])$b
  r<-c(r*0,r)
  r<-filter(r,b,sides=1)[(n+1):(2*n)]
  
  if((p==0)&(q==0))
  {
    n<-length(r)
    sigma2<-sum(r**2)/n
  }
  else
  {
    if((p!=0)&(q!=0))
    {
      model=list(ar=c(eta[2:(p+1)]),ma=c(eta[(p+1):(p+q+1)]))
      result=arma.filt(r,model)
    }
    else
    {
      if(p==0)
      {
        model=list(ma=c(eta[(p+1):(p+q+1)]))
        result=arma.filt(r,model)
      }
      else
      {
        model=list(ar=c(eta[2:(p+1)]))
        result=arma.filt(r,model)
      }
    }
    sigma2<-arima(r,order=c(p,0,q))$sigma2            
    r<-r[c(-(1:p))]-result[-(1:p)]
  }
  
  theta<-c(sigma2, eta)
  # confidence intervals for eta
  
  r.conf<-conf.farima(theta[c(-1)],p,q,ndata,alpha)
  V<-r.conf$V
  CI<-r.conf$CI
  
  # result
  
  drop(list(x=x,
            alpha=0.05,
            iCI=iCI,
            p=p,
            q=q,
            theta=theta,
            eta=eta,
            d=d,
            CI=CI, 
            V=V,
            BIC=BIC,
            r=r))
}
#######################################mle.farima################################### 
mle.farima <- function(x, p, q)
{ 
  eta <- (1:(p+q+1))*0
  loglik <- 0
  r<-x
  n<-length(x)
  ###############fracdiff######################
  if(par.est == "fracdiff")
  {
    result=fracdiff(r,nar=p, nma=q,M=100,drange = c(0,0.5))
    eta<-as.numeric(coef(result))
    loglik<-result$log.likelihood
  }
  ###############arfima######################### 
  if(par.est == "arfima")
  {
    result=result=arfima(r, order = c(p,0,q))
    eta = c(result$modes[[1]]$dfrac,(as.numeric(coef(result)[1:(p+q)])))
    loglik<-result$modes[[1]]$loglik
  }
  etaEST<-eta      
  loglikEST<-loglik  
  
  drop(list(x = x, p=p, q=q, etaEST=etaEST, loglikEST=loglikEST))
  
}
##########################################smooth.lpf################################
smooth.lpf<-function(y,v,p,kn,b,bb)
{ 
  n<-length(y)          
  gr<-rep(0,1,n)         
  hh<-trunc(n*b)         
  htm<-2*hh+1            
  ws<-matrix(0,htm,htm)  
  wk<-rep(0,1,htm)       
  xt<-matrix(0,(p+1),htm)
  xw<-xt                 
  
  #### lpf is a linear smoother. The main task is to calculate the ####
  #### weighting system at each point i. In the following only the ####
  #### weighting systems at the left boundary (i=1 to hh) and that ####
  #### at i=hh+1 will be computed. The weighting system at a point ####
  #### hh+1=<i<=n-hh is the same. The weighting system at i > n-hh ####
  #### is symmetric or asymmetric to that at the point j = n-i+1.  #### 
  
  hr<-c(hh+bb*(hh-0:hh)) #pointwise right bandwith of i,  
  
  ht<-c(1+0:hh)+hr       #pointwise total bandwith of i, 
  
  for(i in (1:(hh+1))){
    
    wk[1:ht[i]]<-(1-(((1:ht[i])-i+1)/(hr[i]+0.5))**2)**(kn-1) 
    
    for(j in (1:(p+1))){
      xt[j,1:ht[i]]<-(((1:ht[i])-i+0.5)/hh)**(j-1) 
      
      xw[j,]<-xt[j,]*t(wk)} 
    
    xa<-solve(xw%*%t(xt))%*%xw 
    ws[i,]<-xa[(v+1),] 
  } 
  ws[(hh+2):htm,]<-(-1)**(v)*ws[hh:1,htm:1]
  ws<-factorial(v)*ws*(n/hh)**(v)     
  
  ym<-matrix(0,(n-htm),htm)
  for(i in ((hh+2):(n-hh))){
    ym[(i-hh-1),]<-y[(i-hh):(i+hh)]
  }
  gr[1:(hh+1)]<-ws[1:(hh+1),]%*%cbind(y[1:htm])
  gr[(hh+2):(n-hh)]<-ym%*%cbind(ws[(hh+1),]) 
  gr[(n-hh+1):n]<-ws[(hh+2):htm,]%*%cbind(y[(n-htm+1):n])
  
  drop(gr)
}
##################################ARMA-filter#####################################SL#
arma.filt=function(x,model)
{
  if (length(model$ma)){
    x <- filter(x, c(1, -model$ma), sides = 1L)
    x[seq_along(-model$ma)]<-0
  }
  if (length(model$ar)){
    x <- filter(x, model$ar, method = "recursive")
    x
  }
}
##########################################kdf.t#############################################
kdf.t<-function(l, m, d)
{S<-0
for(i in 0:m){S1<-0
for(j in i:(l+m)){
  S1<-S1+(-1)^(j-i)*choose(l+m-i, j-i)/(2*d+j+1)*2^(2*d+j+1)
}
S<-S+2*choose(m, i)*S1/(2*d+i)
}
drop(S)
} 
#######################################cov.mle.farima###################################
cov.mle.farima<-function(eta, p, q)
{
  M<-p+q+1                    
  m<-trunc(eta[1]+0.5)
  delta<-eta[1]
  theta<-c(1,eta)
  
  # size of steps in Riemann sum: 2*pi/m.Riemann
  
  m.Riemann<- 10000
  mhalfm   <- trunc((m.Riemann-1)/2) 
  
  # size of delta for numerical calculation of derivative
  
  delta.Deriv<- 0.000000001         
  
  # partial derivatives of log f (at each Fourier frequency)
  
  lf<-matrix(1,ncol=M,nrow=mhalfm)
  f0<-fspec.farima.fourier(theta,p,q,m.Riemann)$f
  for(j in (1:M))
  {
    etaj<-eta
    etaj[j]<-etaj[j]+delta.Deriv
    thetaj<-c(1,etaj)
    fj<-fspec.farima.fourier(thetaj,p,q,m.Riemann)$f
    lf[,j]<-log(fj/f0)/delta.Deriv
  }
  
  # Calculate D
  
  Djl<-matrix(1,ncol=M,nrow=M)
  for(j in (1:M))
  {
    for(l in (1:M))
    {
      Djl[j,l]<-2*2*pi/m.Riemann*sum(lf[,j]*lf[,l])
    }   
  }
  
  # Result
  
  result<-matrix(4*pi*solve(Djl),ncol=M,nrow=M,byrow=T)
  drop(list(eta=eta,p=p,q=q,
            V=result))
}

##########################################conf.farima################################
conf.farima<-function(eta,p,q,n,alpha) 
{
  #
  z<-qnorm(1-alpha/2)
  M<-p+q+1
  
  # Covariance matrix
  V<-matrix(cov.mle.farima(eta, p, q)$V/n, ncol=M, nrow=M, byrow=T)
  
  # C.I.
  etalow<-c()
  etaup<-c()
  for(i in (1:M))
  {
    etalow<-c(etalow,eta[i]-z*sqrt(V[i,i]))
    etaup<-c(etaup,eta[i]+z*sqrt(V[i,i]))
  } 
  CI<-cbind(etalow,etaup)
  
  # result
  
  drop(list(eta=eta,p=p,q=q,n=n,alpha=alpha, 
            V=V,CI=CI))
}
########################################fspec.farima.fourier#####################
fspec.farima.fourier <- function (theta, p, q, n) 
{ 
  
  #---------parameters for the calculation of f--------
  
  sigma2<-    theta[1]
  eta   <-    theta[c(-1)]
  d     <-    eta[1]
  m     <-    trunc(d+0.5)
  delta <-    eta[1]
  phi   <-    c()
  psi   <-    c()
  
  #------   Fourier frequencies: ------------------------------
  #------   x = 2*pi*(j-1)/n (j=1,2,...,(n-1)/2) ------- 
  
  nhalf <- trunc((n-1)/2)
  x  <- 1:nhalf
  x  <- 2*pi/n*x
  
  #-----   calculation of f at Fourier frequencies   -------
  
  far   <-    rep(1,nhalf)
  fma   <-    rep(1,nhalf)
  
  if(p>0) 
  {
    phi    <-  cbind(eta[2:(p+1)])
    cosar  <- cos(cbind(x)%*%rbind(1:p))
    sinar  <- sin(cbind(x)%*%rbind(1:p))
    Rar    <- cosar%*%phi
    Iar    <- sinar%*%phi
    far    <- (1-Rar)**2 + Iar**2
  } 
  
  if(q>0) 
  {
    psi    <- cbind(eta[(p+2):(p+q+1)])
    cosma  <- cos(cbind(x)%*%rbind(1:q))
    sinma  <- sin(cbind(x)%*%rbind(1:q))
    Rma    <- cosma%*%psi
    Ima    <- sinma%*%psi
    fma    <- (1+Rma)**2 + Ima**2
  } 
  
  f.long<-sqrt((1-cos(x))**2 + sin(x)**2)**(-2*d)     
  f.short<-sigma2/(2*pi)*fma/far
  f         <- f.short*f.long
  
  drop(list(theta=theta,p=p,q=q,n=n,
            freq=x,f=f,f.long=f.long,f.short=f.short)) 
  
}
#######################################ar.coef.farima###################################SL
ar.coef.farima <- function(n, d)
{
  if(d==0)
  {
    result<-c(1,rep(0,n-1))
  }
  else
  {
    result<-c(1,gamma(1:100-d)/(gamma(1:100+1)*gamma(-d)))
    result<-c(result,( 101:max(101,(n-1)) )**(-d-1)/gamma(-d))
  }
  result<-result[1:n]
  
  drop(list(n=n,d=d,
            b=cbind(result)))
}
######################################semifar.der#####################################
semifar.der<-function(result, x, mse.RANGE, kn, bb)
{
  #########The data
  x<-x
  
  ##################### Estimating g' based on the above results
  #--- data, ti
  Cf<-result$Cf
  delta<-result$d
  
  n.DATA<-length(x)
  ti<-(1:n.DATA)/n.DATA
  b0<-result$b0
  
  #--- iteration settings
  
  nITER<-20
  icrit<-0
  
  #--- initial parameters for smoothing
  
  Kk<-9
  if(kn==1){Bk<-25/9
  Rk<-1.5}
  if(kn==2){Bk<-49/9
  Rk<-2.1429}
  if(kn==3){Bk<-9        ###81/9
  Rk<-3.1818}  ###35/11	
  
  #------- smoothing iteration
  cat("This is the smooth for estimating'.",fill=T)
  
  for(iITER in 1:nITER)
  {
    cat("iteration=",iITER,fill=T)
    
    #--- data, ti
    
    if( icrit==0 )
    {
      #---------- 3 -------- updated bandwidth for g2
      
      if(iITER==1){b2<-b0}else
      {b2<-(b0)**((7-2*delta)/(14-2*delta))}#opt. for I(g''')
      
      b2<-min(0.49,b2) #### the above bound
      
      #---------- 4 -------- kernel for g2 estimation
      
      #---------- 8 -------- estimate int(g2**2 dt)
      gk<-smooth.lpf(x, 3, 4, kn, b2, bb)###New !!!
      
      index<-max(1,trunc(mse.RANGE*n.DATA)):trunc((1-mse.RANGE)*n.DATA) 
      # avoid border effects
      gkk<-gk[index]**2
      Intgk<-sum(gkk)/n.DATA
      
      #---------- 7 -------- update b0
      
      #------- estimate g
      b0old<-b0
      
      if(delta!=0)
      {
        d<-delta
        Vc<-2*gamma(1-2*d)*sin(pi*d)
        
        
        ####For kernels of order (1, 3)
        if(kn==1){
          Rk<-(-3/2)^2*Vc*kdf.t(1,1,d)
        }  
        if(kn==2){
          Rk<-(-15/4)^2*Vc*(kdf.t(1,1,d)-2*kdf.t(3,1,d)+kdf.t(3,3,d))
        }
        if(kn==3){
          Rk<-(-105/16)^2*Vc*(kdf.t(1,1,d)+4*kdf.t(3,3,d)+kdf.t(5,5,d)
                              -4*kdf.t(3,1,d)+2*kdf.t(5,1,d)-4*kdf.t(5,3,d))
        }
        
        const1<-(Bk*Kk*Rk*(1-2*mse.RANGE)
                 *(2+1-2*delta))**(1/(7-2*delta))
      }
      else
      {
        const1<-(Bk*Rk*(1-2*mse.RANGE)*(2+1)*2*pi)**(1/7)
      }
      
      const2<-(Cf/Intgk)**(1/(7-2*delta))
      const3<-n.DATA**((2*delta-1)/(7-2*delta))
      b0<-const1*const2*const3      
      #       b0<-min(1,b0)
      b0<-min(0.49,b0) #min(1,b0)
      #       b0<-max(n.DATA**(-2/6),b0)
      b0<-max(n.DATA**(-5/7),b0)
      
      cat("Selected b0=", b0, fill=T)
      
      #--- did b0 change much? if not, stop iteration
      
      deltaITER<-0.01*max(b0,b0old)
      if( (abs(b0-b0old)<deltaITER)&(iITER>3) ){icrit<-1}
      
    }
    
    #--- end of smoothing iteration:
    
  }
  
  #--- final result
  g0<-smooth.lpf(x, 1, 2, kn, b0, bb)
  
  drop(g0)
  
}