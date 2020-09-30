require(fracdiff)
require(openxlsx)
setwd("~/Arbeit/DFG/Paper4 - Partial ESEMIFAR/Data/ALV")
#setwd("~/Arbeit/DFG/Paper4 - Partial ESEMIFAR/Data")

X = read.table("ALV_00_14_2-3.min_last.txt", header=TRUE)
#X = read.xlsx("test.xlsx", colNames=TRUE)
#d_name = c("2000-03-10", "2009-09-29", "2011-09-29", "2013-09-30")
d_name = c( "2005-09-30","2006-09-29","2007-09-28","2008-09-30","2009-09-30", "2010-09-30", "2011-09-30", "2012-09-28","2013-09-30","2014-09-30")
sig_test = (1:10)*0
d_par = (1:10)*0
d_time = seq(9,17.5,by = (17.5-9)/509)
for (i in 1:10)
{
Clo = X[,3][X[,1]==d_name[i] & X[,3]>0]
Ret = diff(log(Clo))
mRet = Ret - mean(Ret)
c = 1
lnRet = log(abs(mRet)+c)
#lnRet = (abs(mRet))^(1/4)
#plot.ts(cbind(Ret,lnRet))

pg = 3
kn = 2
bb = 1
IF = 2
result = semifar.lpf(lnRet, 0, 5, 0.025, pg, kn, bb, IF )
sig_test[i] = result$TESTtext
d_par[i] = result$deltaBIC
g_hat = result$g0BIC
if(i==1){par(mfrow=c(5,2),mar = c(4, 4, 4, 2))}
matplot(d_time, cbind(lnRet, g_hat),type="l", xaxt = "n", ylab = "Log Returns", xlab = "Daytime")
axis(1, at = seq(9,17,1))
title(main = paste("Log transformed data of DBKG from",d_name[i]))
nu_hat = exp(g_hat)-c
#matplot(1:510, cbind(mRet, nu_hat),type="l")
}
sig_test
d_par


