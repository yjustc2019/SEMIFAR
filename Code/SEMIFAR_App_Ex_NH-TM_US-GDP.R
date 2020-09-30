require(fracdiff)
library(arfima)
setwd("C:/Users/Letmode/Documents/Arbeit/DFG/Paper4 - Partial ESEMIFAR/Data")
#source("semifar.lpf-final-A.R")

#### US-Quarterly GDP
DAXc="N"
X = read.table("Feng-Gries-2016-D1Q-USGDP.txt", header = TRUE)
Y = X[,2]
Y0 = log(Y)

#### NH temperature, NASA data
X = read.table("NH_MeanTempChange.csv", sep=";", header = TRUE)
Y0=X[,3]

#### VIX
DAXc="N"
X=read.table("VIX_01-1990-07-2019.txt", header = TRUE)
Y0=log(X[,7])

#### US-GDP-Q-1947-Jun-2019
DAXc="N"
X = read.table("GDPC1.csv", sep=",", header = TRUE)
Y0=log(X[,2])

#### for DAX
DAXc="Y"
X0.DAX=read.table("FI-Log-GARCH-DAX-1988-2019-03.txt", header=TRUE)## matrix data with header
Year=X0.DAX[,1][X0.DAX[,7]>0&X0.DAX[,1]>=Year0]
X.DAX=as.numeric(X0.DAX[,7][X0.DAX[,7]>0&X0.DAX[,1]>=Year0])   #### re-order the observations
Ret.DAX=diff(log(X.DAX))
n=length(Ret.DAX)

n.out=250
n.in=n-n.out
mue=mean(Ret.DAX[1:(n-n.out)])
Ret.DAX.In=Ret.DAX[1:(n-n.out)]-mue
Ret.DAX.out=Ret.DAX[(n-n.out+1):n]-mue

RS=(Ret.DAX-mue)^2
Yt=log(RS)
Y0=Yt[1:n.in]


#### SEMIFAR procudrue
n = length(Y0)
pit=5
pmax=5
pe=(1:pit)*0
pe[1]=1
pg=1
kn=3
bb=1
IF=1

for(i in 1:pit)
{result = semifar.lpf(Y0, pe[i], pe[i], 0.025, pg, kn, bb, IF)
g_hat = result$g0BIC
res=Y0-g_hat
BICi=(0:pmax)*0
for(j in 0:pmax)
{            m1=fracdiff(res, nar = j, M=100)
             BICi[j+1]=-2*m1$log.likelihood+log(n)*j
}
pe[i+1]=(1:(pmax+1))[BICi==min(BICi)]-1
if(pe[i+1]==pe[i]){  
pe=pe[1:(i+1)]
break}
}

pe



par(mfrow=c(2,2))

#### only for DAX
if(DAXc=="Y")
{plot.ts(X.DAX)
plot.ts(Ret.DAX)
}


matplot(1:n, cbind(Y0, g_hat), type = "l")
result$pBIC
result$thetaBIC
result$deltaBIC
result$b0BIC
result$CI
result$TESTtext
result$BIC


acf(res, lag.max=100)



