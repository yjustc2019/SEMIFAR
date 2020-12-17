kdf <- function(d, kn = 1, p = 1){
  k = p + 1
  x = seq(-1, 1, by = 0.001)
  n = length(x)
  zj = (1:n)*0
  z = (1:n)*0
  for(i in 1:n){
    for(j in 1:n){
      if(x[i] - x[j] != 0)
        {
        zj[j] = (1 - x[i]^2)^(kn) * (1 - x[j]^2)^(kn) * abs(x[i] - x[j])^((2 * d) - 1)
        }
    }
    #browser()
    z[i] = sum(zj)
  }
  Rk = sum(z)/(n/2)^2
  drop(Rk)
}






