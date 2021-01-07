semifarima.der_pq<-function(x, nu, mse.RANGE, kn, bb) # SL nu
{
  #########The data
  x<-x
  
  result = semifarima_pq.lpf(x, p.min = 0, p.max = 3, q.min = 0, q.max = 3, pg = 3, IF = 2, kn = 2, bb = 1)
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
  if(nu == 1){
    Kk<-9
    if(kn==1){Bk<-25/9
    Rk<-1.5}
    if(kn==2){Bk<-49/9
    Rk<-2.1429}
    if(kn==3){Bk<-9        ###81/9
    Rk<-3.1818}  ###35/11
    ### new #SL
    if(kn==4){Bk<-13.44444  ### 121/9     
    Rk<-4.405594}  
  }
  ### hier noch mal checken mit Tabelle 5.7 aus Müller 1988 #SL
  if(nu == 2){
    Kk<-144 #### Kk=(factorial(p+1))^2/(2*((p+1) - nu))
    if(kn==1){
      Bk<-(7/12)^2                  
      Rk<-22.5}
    if(kn==2){
      Bk<-(3/4)^2                
      Rk<-35}
    if(kn==3){
      Bk<-(11/12)^2             
      Rk<-59.47552}	
    if(kn==4){
      Bk<-(13/12)^2   
      Rk<-94.0724}	
  }
  
  
  
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
      {
        if(nu == 1){b2<-(b0)**((7-2*delta)/(11-2*delta))}
        if(nu == 2){b2<-(b0)**(1/2)}
      } #opt. for I(g''')
      #{b2<-(b0)*n.DATA**((2 - 4*delta)/((2 * (nu + 2) + 1 -2 * d) * (2* (nu + 2) + 3 -2*delta)))}
      
      # if(IFM=="opt"){b2<-(b0)*n.DATA**((2 - 4*d)/((2 * korder + 1 -2 * d) * (2* korder + 3 -2*d)))}
      # if(IFM=="naive"){b2<-(b0)*n.DATA**((4 - 8*d)/((2 * korder + 1 -2*d) *(2*korder+5-2*d)))}
      # if(IFM=="var"){b2<-(b0)*n.DATA**((1 - 2*d)/(4 * korder + 2 -4*d))}
      
      b2<-min(0.49,b2) #### the above bound
      
      #---------- 4 -------- kernel for g2 estimation
      
      #---------- 8 -------- estimate int(g2**2 dt)
      gk<-smooth.lpf(x, nu + 2, nu + 3, kn, b2, bb)###New !!! #SL nu
      
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
        if(nu==1 && kn==1){
          Rk<-(-3/2)^2*Vc*kdf.t(1,1,d)
        }  
        if(nu==1 && kn==2){
          Rk<-(-15/4)^2*Vc*(kdf.t(1,1,d)-2*kdf.t(3,1,d)+kdf.t(3,3,d))
        }
        if(nu==1 && kn==3){
          Rk<-(-105/16)^2*Vc*(kdf.t(1,1,d)+4*kdf.t(3,3,d)+kdf.t(5,5,d)
                              -4*kdf.t(3,1,d)+2*kdf.t(5,1,d)-4*kdf.t(5,3,d))
        }
        ### new kn = 4 #SL
        if(nu==1 && kn==4){
          Rk<-(315/32)^2*Vc*(kdf.t(1,1,d)+9*kdf.t(3,3,d)+9*kdf.t(5,5,d)+kdf.t(7,7,d)
                             -6*kdf.t(3,1,d)+6*kdf.t(5,1,d)-2*kdf.t(7,1,d)-18*kdf.t(5,3,d)
                             +6*kdf.t(7,3,d)-6*kdf.t(7,5,d))
        }
        ####For kernels of order (2, 4) # SL nachfragen
        ##### new SL
        if(nu==2 && kn==1){
          Rk<-(15/4)^2*Vc*(kdf.t(0,0,d)-
                             6*kdf.t(2,0,d)+
                             9*kdf.t(2,2,d))
        }
        if(nu==2 && kn==2){
          Rk<-(105/16)^2*Vc*(kdf.t(0,0,d)+
                               36*kdf.t(2,2,d)+
                               25*kdf.t(4,4,d)-
                               12*kdf.t(2,0,d)+
                               10*kdf.t(4,0,d)-
                               60*kdf.t(4,2,d))  
        }
        if(nu==2 && kn==3){
          Rk<-(315/32)^2*Vc*(kdf.t(0,0,d)+   
                               81*kdf.t(2,2,d)+
                               225*kdf.t(4,4,d)+
                               49*kdf.t(6,6,d)-
                               18*kdf.t(2,0,d)+
                               30*kdf.t(4,0,d)-
                               14*kdf.t(6,0,d)-
                               270*kdf.t(4,2,d)+
                               126*kdf.t(6,2,d)-
                               210*kdf.t(6,4,d))
        }  
        if(nu==2 && kn==4){
          Rk<-(3465/256)^2*Vc*(kdf.t(0,0,d)+   
                                 144*kdf.t(2,2,d)+
                                 900*kdf.t(4,4,d)+
                                 784*kdf.t(6,6,d)+
                                 81*kdf.t(8,8,d)-
                                 24*kdf.t(2,0,d)+
                                 60*kdf.t(4,0,d)-
                                 56*kdf.t(6,0,d)+
                                 18*kdf.t(8,0,d)-
                                 720*kdf.t(4,2,d)+
                                 672*kdf.t(6,2,d)-
                                 216*kdf.t(8,2,d)-
                                 1680*kdf.t(6,4,d)+
                                 540*kdf.t(8,4,d)-
                                 504*kdf.t(8,6,d))
        }  
        const1<-(Bk*Kk*Rk*(1-2*mse.RANGE)
                 *(2*nu+1-2*delta))**(1/(2*(nu+2)+1-2*delta))
      }
      else
      {
        const1<-(Bk*Rk*(1-2*mse.RANGE)*(2*nu+1)*2*pi)**(1/(2*(nu+2)+1))
      }
      
      const2<-(Cf/Intgk)**(1/(2*(nu+2)+1-2*delta))
      const3<-n.DATA**((2*delta-1)/(2*(nu+2)+1-2*delta))
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
  g0<-smooth.lpf(x, nu, nu + 1, kn, b0, bb)# SL nu
  
  drop(g0)
  
}