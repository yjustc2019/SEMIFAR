n= 100000000
u = (-n:n)/(n+0.5)

Bk0 = function(u){u^4*(1 -3*u^2)}
Rk0 = function(u){(1 -3*u^2)^2}

Bk1 = function(u){u^4*(-1 + 6*u^2-5*u^4)}
Rk1 = function(u){(-1 + 6*u^2-5*u^4)^2}

Bk2 = function(u){u^4*(-1 + 9*u^2-15*u^4+7*u^6)}
Rk2 = function(u){(-1 + 9*u^2-15*u^4+7*u^6)^2}

Bk3 = function(u){u^4*(-1 + 12*u^2-30*u^4+28*u^6-9*u^8)}
Rk3 = function(u){(-1 + 12*u^2-30*u^4+28*u^6-9*u^8)^2}



Bk.mu0 = 1 / (sum(Bk0(u))/n * (15/4))^2
Rk.mu0 = sum(Rk0(u))/n *(15/4)^2

Bk.mu1 = 1 / (sum(Bk1(u))/n * (105/16))^2
Rk.mu1 = sum(Rk1(u))/n *(105/16)^2

Bk.mu2 = 1 / (sum(Bk2(u))/n * (315/32))^2
Rk.mu2 = sum(Rk2(u))/n *(315/32)^2

Bk.mu3 = 1 / (sum(Bk3(u))/n * (3465/256))^2
Rk.mu3 = sum(Rk3(u))/n *(3465/256)^2

Rk = function(u){(-u+3*u^3-3*u^5+u^7)^2}
Bk = function(u){u^3*(-u+3*u^3-3*u^5+u^7)}
sum(Rk(u))/n*(315/32)^2
1/(sum(Bk(u))/n *(315/32))^2
