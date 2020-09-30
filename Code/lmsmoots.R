###################################semifarima.lpf - extension of semifar.lpf with an additional MA-part################################
###################################content is heavily simplified, e.g. no more integer differencing and###############################
###################################other adjustments#################################################################################
lmsmoots = function(x, mse.RANGE,b0 = NULL, pg, kn, bb, bbi = 1, bbd = 1, IM, IF, method = "OM", A = 0.3)
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
      g0<-smooth.lpf(x, 0, pg, kn, b0, bbi)
      x.DETRENDED<-x-g0
      ############ estimation of d and filtered values ##########

      d = gph.est(x.DETRENDED, method, A)[2]
      d.all = gph.est(x.DETRENDED, method, A)
      x.filt =  dfilt(x.DETRENDED, d)


      ############## Cf ##################
      #Cf = cf0.LW.est(x.filt)$cf0.LW/(2*pi)
      Cf = d.all[1]
      ############ inflation factors ###SL
      if(convc==1){
        if(IM == "E"){
          if(IFM=="opt"){b2<-(b0)**((2 * korder + 1 -2*d)/(2 * korder + 3-2*d))}
          if(IFM=="naive"){b2<-(b0)**((2 * korder +1 -2*d)/(2 * (korder + 2) + 1-2*d))}
          if(IFM=="var"){b2<-(b0)**((1/2))}
        }
        if(IM =="M"){
          if(IFM=="opt"){b2<-(b0)*n.DATA**((2 - 4*d)/((2 * korder + 1 -2 * d) * (2* korder + 3 -2*d)))}
          if(IFM=="naive"){b2<-(b0)*n.DATA**((4 - 8*d)/((2 * korder + 1 -2*d) *(2*korder+5-2*d)))}
          if(IFM=="var"){b2<-(b0)*n.DATA**((1 - 2*d)/(4 * korder + 2 -4*d))}
        }

        b2<-min(0.49,b2)

        if(iITER<=2)
        #{gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, bb)}
        {if(pg==1){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, bbd)}
          if(pg==3){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, bbd)}}
        else
        #{gk<-smooth.lpf(x, pg+1, pg+2, kn, b2, bb)}
        {if(pg==1){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, bbd)}
          if(pg==3){gk<-smooth.lpf(x, pg+1, pg+2, 1, b2, bbd)}}

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
        #cat(const1, const2, const3,d, fill=TRUE)
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
  }else
  {
    nu<-2**(2*d)*gamma(1-2*d)*sin(pi*d)
    nu<-nu/(d*(2*d+1))
  }


  d <- d

  ####significance test trend
  g0 = smooth.lpf(x, 0, pg, kn, b0, bb)
  VARg0<-(n*b0)**(2*d-1)*nu*Cf
  CRITg0<-qnorm( (1-alpha.SEMI/2) )*sqrt(VARg0)
  SIG<-1
  MEAN<-mean(x)
  g0.limits<-c(MEAN-CRITg0,MEAN+CRITg0)

  if( max(abs(g0-MEAN)) <= CRITg0 )
  {
    SIGtext<-"Trend g0 is not significant!"
    SIG<-0
  }else
  {
    SIGtext<-"Trend g0 is significant!"
  }
  result<-list(x = x,
               mse.RANGE=mse.RANGE,
               n = n.DATA,
               d=d,
               d.all=d.all,
               g0=g0,
               b0=b0,
               Cf=Cf,
               nu=nu,
               SIG=SIG,
               SIGtext=SIGtext)
  drop(result)
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

###### App. 2: function for selecting M and estimating cf nonparametrically
###### an IPI-algorithm for estimating the spdf proposed by Buehlmann, 1996
###### For simplicity only the Bartllet window, without the c2-window, is used
cf0.LW.est=function(Xt)
{n=length(Xt)
ga=acf(Xt, lag.max=(n-1), type="covariance", plot=FALSE)$acf

nit=20
runc=1
L=(1:(nit+1))*0
L[1]=trunc(n/3)
#### the IPI-procedure for selecting L
c1=(ga[1]**2+2*sum(ga[2:n]**2))/(4*pi)
for(i in 1:nit)
{if(runc==1){
  L1=min(trunc(L[i]/n**(2/21))+2, L[1])
  x1=(0:(L1-1))/L1
  w1=1-x1 ##### the Bartlett window weights
  gai=(0:(L1-1))*ga[1:L1]*w1
  c2=3*(2*sum(gai[2:L1]**2))/(2*pi)
  L[i+1]=min(trunc(n**(1/3)*(c2/c1)**(1/3))+2, L[1])
  if(L[i+1]==L[i])
  {runc=0
  LG.opt=L[i+1]
  }
}
}
#### end of the global selection
#### the localized step
if(runc==1){LG.opt=L[nit+1]}
L1=min(trunc(LG.opt/n**(2/21))+2, L[1])
x1=(0:(L1-1))/L1
w1=1-x1 ##### the Bartlett window weights
ga1=(0:(L1-1))*ga[1:L1]*w1
c20=3*(2*sum(ga1[2:L1]))**2/(2*pi)

w0=(1+cos(pi*x1))/2 ##### Tukey-Hanning window weights
ga0=ga[1:L1]*w0
c10=(ga0[1]+2*sum(ga0[2:L1]))**2/(2*pi)

L0.opt=min(trunc(n**(1/3)*(c20/c10/2)**(1/3))+2, L[1])
#### final estimation of cf0 with L0.opt
wacf=((L0.opt+1):1)/(L0.opt+1) ##### the Bartlett window weights
acf.X=acf(Xt, lag.max=L0.opt, type="covariance", plot=FALSE)$acf
cf0.LW=2*sum(acf.X*wacf)-acf.X[1]

results=list(cf0.LW=cf0.LW, L0.opt=L0.opt, LG.opt=LG.opt)
drop(results)
}
#### end of the selection procedure for L

###############dfilt###############################
dfilt<-function(x,d){
  n=length(x)
  result<-c(1,gamma(1:100-d)/(gamma(1:100+1)*gamma(-d)))
  result<-c(result,( 101:max(101,(n-1)) )**(-d-1)/gamma(-d))
  result<-result[1:n]
  b=cbind(result)
  r<-c(x*0,x)
  r<-filter(r,b,sides=1)[(n+1):(2*n)]
  drop(r)
}
################gph.est####################################
gph.est <- function(x, method, A){
  n = length(x)
  n2 = trunc(A * n ^ (6/7))
  #browser()
  if (method == "FF") {logIj = log(per(x-mean(x))[2:(n2+1)])}
  if (method == "OM") {logIj = log(periodgm(x-mean(x))$Ilamb[2:(n2+1)])}
  lambdaj <- 2 * pi/n * 1:n2
  aj <- log(abs(1 - exp((0+1i) * lambdaj)))
  X = cbind(rep(1, n2), aj, lambdaj^2/2)
  b = ginv(t(X)%*%X)%*%t(X)
  K = sum(b[3, ] * logIj)
  C = (27 / (128 * pi^2)) ^ (.2) * (K ^ 2) ^ (-1/5)
  m = C * n ^ .8
  d.B = gph(x,m,l=1)
  M1 = lm(logIj ~ aj)
  Cf = exp(M1$coefficients[1] - digamma(1))
  d = d.B + (2 * pi ^ 2) / 9 * K * m ^ 2 / n ^ 2
  d.B.upper = d.B + 2*sqrt(pi^2/(24*m))
  d.B.lower = d.B - 2*sqrt(pi^2/(24*m))
  d.upper = d + 2*sqrt(pi^2/(24*m))
  d.lower = d - 2*sqrt(pi^2/(24*m))
  result = round(cbind(Cf, d.B, d.B.lower, d.B.upper, d, d.lower, d.upper, m),4)
  colnames(result) = c("Cf", "d.B", "d.B.L", "d.B.U", "d", "d.L", "d.U", "m")
  drop(result)
}
#############periodgm Feng###################
periodgm=function(x)
{n=length(x)
  m=trunc(n/2)
  lambda=(0:m)/n
  t=(1:n)

  Ilamb=0
  for(i in 1:m)
  {Ct=cos(lambda[i]*2*pi*t)
    St=sin(lambda[i]*2*pi*t)

    Ilambi=((sum(x*Ct))**2 + (sum(x*St))**2)/n
    Ilamb=c(Ilamb, Ilambi)
  }

  result=list(lambda=lambda, Ilamb=Ilamb)
  drop(result)
}
