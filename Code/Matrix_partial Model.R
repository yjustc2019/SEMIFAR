require(rgl)
require(arfima)
setwd("~/Arbeit/DFG/Paper4 - Partial ESEMIFAR/Data/DBKG")
#setwd("~/Arbeit/DFG/Paper4 - Partial ESEMIFAR/Data/1-minute-Data-examples/Examples ALV-BMW")

X = read.table("DBKG_00_14_2-3.min_last.txt", header = TRUE)
#X = read.table("ALV-1-Min-2006-09-2014.txt", header = TRUE)
Ret = diff(log(X[,1]))
mRet = Ret - mean(Ret)
lnRet = log(mRet^2)
n.col = length(X[,1])/511
Y = matrix(0, 510, n.col)
Z = matrix(0, 510, n.col)
day.ind = seq(0, length(lnRet), by = 510)
for (i in 1:n.col)
{
  Z[,i] = mRet[(day.ind[i]+1) : day.ind[i+1]]
  Y[,i] = lnRet[(day.ind[i]+1) : day.ind[i+1]]
}
y = (1:(n.col))/(n.col)*(8+9/12)+2006
x = 1:510

# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("yellow", "black") )
color = jet.colors(10000)
zfacet <- (Z[-1, -1] + Z[-1, -ncol(Z)] + Z[-nrow(Z), -1] + Z[-nrow(Z), -ncol(Z)])/4
zcol = cut(zfacet,10000)
persp3d(z = Z, x, y, col=color[zcol])
grid3d(c("x", "y+", "z"))




