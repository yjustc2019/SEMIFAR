### K Methode zur Aggregation von UHF Daten mit k= N/M 
setwd("~/Arbeit/DFG/Paper4 - Partial ESEMIFAR/Data/DBKG")

name.in = c("DBKG_PV_00_1-2.txt","DBKG_PV_00_3-4.txt","DBKG_PV_01_1-2.txt","DBKG_PV_01_3-4.txt","DBKG_PV_02_1-2.txt","DBKG_PV_02_3-4.txt",
            "DBKG_PV_03_1-2.txt","DBKG_PV_03_3-4.txt","DBKG_PV_04_1-2.txt","DBKG_PV_04_3-4.txt","DBKG_PV_05_1-2.txt","DBKG_PV_05_3-4.txt",
            "DBKG_PV_06_1-2.txt","DBKG_PV_06_3-4.txt","DBKG_PV_07_1-2.txt","DBKG_PV_07_3-4.txt","DBKG_PV_08_1-2.txt","DBKG_PV_08_3-4.txt",
            "DBKG_PV_09_1-2.txt","DBKG_PV_09_3-4.txt","DBKG_PV_10_1-2.txt","DBKG_PV_10_3-4.txt","DBKG_PV_11_1-2.txt","DBKG_PV_11_3-4.txt",
            "DBKG_PV_12_1-2.txt","DBKG_PV_12_3-4.txt","DBKG_PV_13_1-2.txt","DBKG_PV_13_3-4.txt","DBKG_PV_14_1-2.txt","DBKG_PV_14_3.txt")
name.out = c("DBKG_PV_00_1-2.min.txt","DBKG_PV_00_3-4.min.txt","DBKG_PV_01_1-2.min.txt","DBKG_PV_01_3-4.min.txt","DBKG_PV_02_1-2.min.txt","DBKG_PV_02_3-4.min.txt",
             "DBKG_PV_03_1-2.min.txt","DBKG_PV_03_3-4.min.txt","DBKG_PV_04_1-2.min.txt","DBKG_PV_04_3-4.min.txt","DBKG_PV_05_1-2.min.txt","DBKG_PV_05_3-4.min.txt",
             "DBKG_PV_06_1-2.min.txt","DBKG_PV_06_3-4.min.txt","DBKG_PV_07_1-2.min.txt","DBKG_PV_07_3-4.min.txt","DBKG_PV_08_1-2.min.txt","DBKG_PV_08_3-4.min.txt",
             "DBKG_PV_09_1-2.min.txt","DBKG_PV_09_3-4.min.txt","DBKG_PV_10_1-2.min.txt","DBKG_PV_10_3-4.min.txt","DBKG_PV_11_1-2.min.txt","DBKG_PV_11_3-4.min.txt",
             "DBKG_PV_12_1-2.min.txt","DBKG_PV_12_3-4.min.txt","DBKG_PV_13_1-2.min.txt","DBKG_PV_13_3-4.min.txt","DBKG_PV_14_1-2.min.txt","DBKG_PV_14_3.min.txt")
for (j in 1:length(name.in))
{
  X=read.table(name.in[j], header=TRUE)
  
  D.Ind=unique(X$Day.N)
  L=length(D.Ind)
  M=256
  
  Results.all=matrix(0,M,L)
  
  for(l in 1:L)
  {
    
    Trading.Date=X$Day.N[X$Day.N==D.Ind[l]] 
    Trading.Hour=X$Hour[X$Day.N==D.Ind[l]]
    Trading.Minute=X$Minute[X$Day.N==D.Ind[l]]
    Trading.Second=X$Second[X$Day.N==D.Ind[l]]
    Trading.Price=X$Price[X$Day.N==D.Ind[l]]
    N.Dn=length(Trading.Price) 
    
    
    Ind.t=1:(M+1)
    Ind.t[M+1]=N.Dn
    Ind.t[2:M]=trunc((1:(M-1))/M*(N.Dn-1)+1.5) 
    #Ind.t[which(Ind.t==0)]=1
    
    
    Date=(1:M)*0
    Date[1]=Trading.Date[Ind.t[1]+1]
    Hour=(1:M)*0
    Hour[1]=Trading.Hour[Ind.t[1]+1]
    Minute=(1:M)*0
    Minute[1]=Trading.Minute[Ind.t[1]+1]
    Second=(1:M)*0
    Second[1]=Trading.Second[Ind.t[1]+1]
    Price=(1:M)*0
    Price[1]=Trading.Price[Ind.t[1]+1]
    
    for(i in 1:M)
    {
      Date[i]=Trading.Date[Ind.t[i+1]]
      Hour[i]=Trading.Hour[Ind.t[i+1]]
      Minute[i]=Trading.Minute[Ind.t[i+1]]
      Second[i]=Trading.Second[Ind.t[i+1]]
      Price[i]=Trading.Price[Ind.t[i+1]]
    }
    
    
    K=diff(Ind.t)
    
    Time=Hour+Minute/60+Second/3600
    Date=as.Date(Date, origin = '1899-12-30')
    Results.trans=data.frame(Date, Time, Price, K, stringsAsFactors = FALSE)
    assign(paste("Day",l, sep=""), Results.trans)
  }
  
  Result=Day1
  
  for (i in 2:L)
  {
    Result=rbind(Result,get(paste("Day",i,sep="")))
  }
  
  assign(paste("Res",j,sep=""),Result)
  #write.table(Result,file=name.out[j],sep="     ",row.names=FALSE,col.names=c("Date","Time","Price","K"))
}
Res.tot = Res1
for (m in 2:length(name.in))
{Res.tot = rbind(Res.tot,get(paste("Res",m,sep="")))}
write.table(Res.tot,file="DBKG00-14_3.2min.txt",sep="     ",row.names=FALSE,col.names=c("Date","Time","Price","K"))