#install.packages("rgdal")
#install.packages("spdep")
setwd("C:\\Users\\jonel\\Bureaublad\\spatialstat\\proj3")
rm(list=ls())
library(Matrix)
library(boot)
library(matlib)
library(MASS)
source("functions.R")
load("Admin1Geography.RData")

Admin1 <- read.table('Admin1Graph.txt')
Admin2 <- read.table('Admin2Graph.txt')

Observations <- read.table('DirectEstimates.txt')

phat = Observations$V2
L = length(phat)
y = as.double(phat[2:L])


stddev = Observations$V3
V = as.double(stddev[2:L])^2

phat = inv.logit(y)

#plotAreaCol("2phat.png", 15, 15, phat, nigeriaAdm1, "observed proportion")


# b)

sigma2= 100*100
mu = y * sigma2 / (sigma2 + V)
var = sigma2 * V / (sigma2 + V)

sim_x = mvrnorm(100, mu, diag(var))
sim_p = inv.logit(sim_x) 

sim_p_median = seq(0, 0, length=37)
sim_p_cv = seq(0, 0, length=37)


for (area in 1:37) {
  array = sim_p[,area]
  sim_p_median[area] = median(array)
  sim_p_cv[area] = sd(array) / mean(array) * 100
}

#plotAreaCol("2bmedian.png", 15, 15, sim_p_median, nigeriaAdm1, "median")
#plotAreaCol("2bcv.png", 15, 15, sim_p_cv, nigeriaAdm1, "cv")
#plotAreaCol("2bdifference.png", 15, 15, phat - sim_p_median, nigeriaAdm1, "difference")



# c)

# admin1
R1 = as.matrix(Admin1)
test = seq(1, 1, length=37)
s = R1 %*% test
for (row in 1:37) {
  R1[row, row] = -s[row]
}
det(R1)
R1 = - R1
s = R1 %*% test
Q = R1

D_inv = inv(diag(V))

Qfinal = (Q + D_inv)
mu_final = inv(Q + D_inv) %*% D_inv %*% y 

#plotAreaCol("2cmean.png", 15, 15, inv.logit(mu_final), nigeriaAdm1, "difference")

sim_x = mvrnorm(100, mu_final, inv(Qfinal)  )
sim_p = inv.logit(sim_x) 

sim_p_median = seq(0, 0, length=37)
sim_p_cv = seq(0, 0, length=37)


for (area in 1:37) {
  array = sim_p[,area]
  sim_p_median[area] = median(array)
  sim_p_cv[area] = sd(array) / mean(array) * 100
}

#plotAreaCol("2cmedian.png", 15, 15, sim_p_median, nigeriaAdm1, "median")
#plotAreaCol("2ccv.png", 15, 15, sim_p_cv, nigeriaAdm1, "cv")
#plotAreaCol("2cdifference.png", 15, 15, phat - sim_p_median, nigeriaAdm1, "difference")


# d)
Sigm_inv = matrix(data = 0, nrow = 37, ncol = 37)
Sigm_inv[19, 19] = 1 / 0.1^2
Qfinal = (Q + D_inv + Sigm_inv)
mu_final = inv(Q + D_inv + Sigm_inv) %*% (D_inv + Sigm_inv) %*% y 

#plotAreaCol("2cmean.png", 15, 15, inv.logit(mu_final), nigeriaAdm1, "difference")

sim_x = mvrnorm(100, mu_final, inv(Qfinal)  )
sim_p = inv.logit(sim_x) 

sim_p_median = seq(0, 0, length=37)
sim_p_cv = seq(0, 0, length=37)


for (area in 1:37) {
  array = sim_p[,area]
  sim_p_median[area] = median(array)
  sim_p_cv[area] = sd(array) / mean(array) * 100
}

kaduna = seq(0.2, 0.5, length=37)
kaduna[19] = 1
#plotAreaCol("kaduna.png", 15, 15, kaduna, nigeriaAdm1, "kaduna")

#plotAreaCol("2dmedian.png", 15, 15, sim_p_median, nigeriaAdm1, "median")
#plotAreaCol("2dcv.png", 15, 15, sim_p_cv, nigeriaAdm1, "cv")
#plotAreaCol("2ddifference.png", 15, 15, phat - sim_p_median, nigeriaAdm1, "difference")



# e)
Q = 0.1 * R1

D_inv = inv(diag(V))

Qfinal = (Q + D_inv)
mu_final = inv(Q + D_inv) %*% D_inv %*% y 

#plotAreaCol("2cmean.png", 15, 15, inv.logit(mu_final), nigeriaAdm1, "difference")

sim_x = mvrnorm(100, mu_final, inv(Qfinal)  )
sim_p = inv.logit(sim_x) 

sim_p_median = seq(0, 0, length=37)
sim_p_cv = seq(0, 0, length=37)


for (area in 1:37) {
  array = sim_p[,area]
  sim_p_median[area] = median(array)
  sim_p_cv[area] = sd(array) / mean(array) * 100
}

#plotAreaCol("2emedian.png", 15, 15, sim_p_median, nigeriaAdm1, "median")
#plotAreaCol("2ecv.png", 15, 15, sim_p_cv, nigeriaAdm1, "cv")

Q = 10 * R1

D_inv = inv(diag(V))

Qfinal = (Q + D_inv)
mu_final = inv(Q + D_inv) %*% D_inv %*% y 

#plotAreaCol("2cmean.png", 15, 15, inv.logit(mu_final), nigeriaAdm1, "difference")

sim_x = mvrnorm(100, mu_final, inv(Qfinal)  )
sim_p = inv.logit(sim_x) 

sim_p_median = seq(0, 0, length=37)
sim_p_cv = seq(0, 0, length=37)


for (area in 1:37) {
  array = sim_p[,area]
  sim_p_median[area] = median(array)
  sim_p_cv[area] = sd(array) / mean(array) * 100
}

#plotAreaCol("2emedian2.png", 15, 15, sim_p_median, nigeriaAdm1, "median")
#plotAreaCol("2ecv2.png", 15, 15, sim_p_cv, nigeriaAdm1, "cv")


# f)
tau = 1
Q = tau * R1

D_inv = inv(diag(V))



Qfinal = (Q + D_inv)
mu_final = inv(Q + D_inv) %*% D_inv %*% y 

phat = Observations$V2
L = length(phat)
y = as.double(phat[2:L])
x = mu_final

rm(Q, mu_final, Qfinal)

f <- function (tau) {
  print(tau)
  tau = as.numeric(tau)
  
  Q = as.matrix(R1 * tau)
  det(Q)
  Qfinal = (Q + D_inv)
  mu_final = inv(Q + D_inv) %*% D_inv %*% y 
  
  #x = as.matrix(da)
  
  f1 = 18 * log(tau) - tau / 2 * (t(x) %*% R1 %*% x)
  f2 = - 1 / 2 * t(y - x) %*% D_inv %*% (y - x)
  f3 = - 1 / 2 * log(det(Qfinal))
  f4 = 1 / 2 * t(x-mu_final) %*% Qfinal %*% (x-mu_final)
  print("a")
  print(f1)
  print(f2)
  print(f3)
  print(f4)
  -(f1 + f2 + f3 + f4)
}


tau_hat = optimize(f, lower = 0.1, upper = 10)


Q = tau_hat$minimum * R1

D_inv = inv(diag(V))

Qfinal = (Q + D_inv)
mu_final = inv(Q + D_inv) %*% D_inv %*% y 

#plotAreaCol("2cmean.png", 15, 15, inv.logit(mu_final), nigeriaAdm1, "difference")

sim_x = mvrnorm(100, mu_final, inv(Qfinal)  )
sim_p = inv.logit(sim_x) 

sim_p_median = seq(0, 0, length=37)
sim_p_cv = seq(0, 0, length=37)


for (area in 1:37) {
  array = sim_p[,area]
  sim_p_median[area] = median(array)
  sim_p_cv[area] = sd(array) / mean(array) * 100
}

plotAreaCol("2fmedian.png", 15, 15, sim_p_median, nigeriaAdm1, "median")
plotAreaCol("2fcv.png", 15, 15, sim_p_cv, nigeriaAdm1, "cv")
phat = inv.logit(y)
plotAreaCol("2fdifference.png", 15, 15, phat - sim_p_median, nigeriaAdm1, "difference")



