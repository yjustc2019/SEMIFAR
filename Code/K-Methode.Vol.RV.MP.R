### K Methode zur Aggregation von UHF Daten mit k= N/M 

X=read.table("ALV_PV_00_1-2.txt", header=TRUE)

D.Ind=unique(X$Day.N)
L=length(D.Ind)
M=510



Results.all=matrix(0,L*M,4)
for(l in 1:L)
{
  
  
  Trading.Volume=X$Volume[X$Day.N==D.Ind[l]]  #### für Volumen
  Trading.Price=X$Price[X$Day.N==D.Ind[l]] #### für Preis
  N.Dn=length(Trading.Volume) #### Länge des Datensatzes
  
  
  Ind.t=1:(M+1)
  Ind.t[M+1]=N.Dn
  Ind.t[2:M]=trunc((1:(M-1))/M*N.Dn+0.5) #### (M+1) Indizes sind definiert. 
  #### Die Definition der Indizes ist für alle dieselbe.
  
  
  #######hier fängt die Schleife für Volumen an
  
  Cum.Trade.Volume=(1:M)*0	
  Mean.Trade.Volume=(1:M)*0	
  Cum.Trade.Volume[1]=sum(Trading.Volume[(Ind.t[1]+1):(Ind.t[2]-1)])
  Mean.Trade.Volume[1]=mean(Trading.Volume[(Ind.t[1]+1):(Ind.t[2]-1)])
  
  for(i in 2:(M-1))
  {
    Cum.Trade.Volume[i]=sum(Trading.Volume[Ind.t[i]:(Ind.t[i+1]-1)])  
    Mean.Trade.Volume[i]=mean(Trading.Volume[Ind.t[i]:(Ind.t[i+1]-1)])
  }
  
  Cum.Trade.Volume[M]=sum(Trading.Volume[Ind.t[M]:Ind.t[M+1]])
  Mean.Trade.Volume[M]=mean(Trading.Volume[Ind.t[M]:Ind.t[M+1]])
  
  
  
  
  
  
  Mean.Trade.Price=(1:M)*0       #### vordefinierte Matrix für den Durchschnittspreis 
  Mean.Trade.Price[1]=mean(Trading.Price[(Ind.t[1]+1):(Ind.t[2]-1)])
  
  for(i in 2:(M-1))
  {
    
    Mean.Trade.Price[i]=mean(Trading.Price[Ind.t[i]:(Ind.t[i+1]-1)])
  }
  
  
  Mean.Trade.Price[M]=mean(Trading.Price[Ind.t[M]:Ind.t[M+1]])
  
  
  #######Returns / RV
  
  Returns=diff(log(Trading.Price))
  RV=(1:M)*0
  RV[1]=sum((Returns[(Ind.t[1]):(Ind.t[2]-1)])^2)
  for(i in 2:M)
  {
    RV[i]=sum((Returns[Ind.t[i]:(Ind.t[i+1]-1)])^2)
  }
  
  K=diff(Ind.t)

  Results=cbind(Cum.Trade.Volume, Mean.Trade.Price, RV, K)
  assign(paste("Day",l, sep=""),Results)

}
Result=Day1
for (i in 2:L)
{
  Result=rbind(Result,get(paste("Day",i,sep="")))
  
}
write.table(Result,file="ALV_PV_00_1-2.min.txt",sep="     ",row.names=FALSE,col.names=c("CumVol","MeanPr","RV", "K"))