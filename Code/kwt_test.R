kwt <- function(d, kn = 1, p = 1){
  k = p + 1
  lamb = (1:(k + 2 * kn - 1))*0
  lamb[1] = (2 * factorial(2 * kn + 2)) / (factorial(kn) * factorial(kn + 1) * 2^(2*kn+3))
  for (i in 2:(k + 2 * kn - 1)){
    if((i - 1) %% 2 == 0)
      {
      if(kn - (i-1)/2 < 0) {c = 0}else{c = kn - (i-1)/2}
      lamb[i] = ((-1)^((i-1)/2) * (i + 1) * factorial(2 * kn + 2)) /
      (factorial(1+ (i-1)/2) * factorial(c) * factorial(kn +1) * 2^(2*kn+3))
    }
  }
  x = seq(-1, 1, by = 0.001)
  n = length(x)
  wt = (1:n)*0
  for(j in 1:n){
    wt[j] = lamb[1] + sum(lamb[-1] * (poly(x[j], degree = (k + 2 * kn - 2), raw = TRUE)))
  }
  zm = (1:n)*0
  z = (1:n)*0
  for(l in 1:n){
    for(m in 1:n){
      if(x[l] - x[m] != 0)
      {
        zm[m] = wt[l] * wt[m] * abs(x[l] - x[m])^((2 * d) - 1)
      }
    }
    z[l] = sum(zm)
  }
  Rk = sum(z)/(n/2)^2
  drop(Rk)
}


