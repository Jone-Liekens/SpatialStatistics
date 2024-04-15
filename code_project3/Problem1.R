#install.packages("rgdal")
#install.packages("spdep")
setwd("C:\\Users\\jonel\\Bureaublad\\spatialstat\\proj3")
rm(list=ls())
library(Matrix)
library(matlib)
library(MASS)
library(base)
#install.packages("ggpubr")
library("ggpubr")
source("functions.R")
load("Admin1Geography.RData")
load("Admin2Geography.RData")


Admin1 <- read.table('Admin1Graph.txt')
Admin2 <- read.table('Admin2Graph.txt')


rankMatrix(Admin1)
rankMatrix(Admin2)


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

# visualize sparsity
par(pty="s")
image(1 - (R1 != 0), col = hcl.colors(12, "BuPu", rev = TRUE),)
sum(R1 != 0)
(37*37)
sum(R1 != 0) / (37*37)


# admin2
R2 = as.matrix(Admin2)
test = seq(1, 1, length=775)
s = R2 %*% test
for (row in 1:775) {
  R2[row, row] = -s[row]
}
det(R2)
R2 = - R2
s = R2 %*% test

# visualize sparsity again
par(pty="s")
#image(1 - (R2 != 0), col = hcl.colors(12, "BuPu", rev = TRUE),)
sum(R2 != 0) 
(775*775)
sum(R2 != 0) / (775*775)


# b)

# simulate besag
eps = 1e-8
R1_eps = R1 + eps * diag(37)
inv_R1_eps = inv(R1_eps)
mu = seq(0, 0, length=37)
#ev <- eigen(inv_R1_eps)
#values <- ev$values
#values

s = mvrnorm(1, mu, inv_R1_eps)
s = s - mean(s)
#plotAreaCol("1bsim1.png", 15, 15, s, nigeriaAdm1, "x")

s = mvrnorm(1, mu, inv_R1_eps)
s = s - mean(s)
#plotAreaCol("1bsim2.png", 15, 15, s, nigeriaAdm1, "x")


# simulate gaussian
s = rnorm(37, 0, 1)
#plotAreaCol("1bsim3.png", 15, 15, s, nigeriaAdm1, "x")
s = rnorm(37, 0, 1)
#plotAreaCol("1bsim4.png", 15, 15, s, nigeriaAdm1, "x")



# c)
# simulate besag with Cholesky factorization
eps = 1e-8
N2 = 775
R2_eps = R2 + eps * diag(N2)
L = chol(R2_eps)

z = rnorm(N2, 0, 1)
s = backsolve(L, z, k = ncol(L), upper.tri = TRUE, transpose = FALSE)
s =  s - mean(s)
plotAreaCol("1csim1.png", 15, 15, s, nigeriaAdm2, "x")

z = rnorm(N2, 0, 1)
s = backsolve(L, z, k = ncol(L), upper.tri = TRUE, transpose = FALSE)
s =  s - mean(s)
plotAreaCol("1csim2.png", 15, 15, s, nigeriaAdm2, "x")


# simulate Gaussian
s = rnorm(N2, 0, 1)
plotAreaCol("1csim3.png", 15, 15, s, nigeriaAdm2, "x")
s = rnorm(N2, 0, 1)
plotAreaCol("1csim4.png", 15, 15, s, nigeriaAdm2, "x")


# d)
n = 100
results = matrix(NA, nrow = n, ncol = N2)
for (row in 1:n) {
  z = rnorm(N2, 0, 1)
  s = backsolve(L, z, k = ncol(L), upper.tri = TRUE, transpose = FALSE)
  s = s - mean(s)
  results[row,] = s
}
variances = seq(0, 0, length=N2)
correlations = seq(0, 0, length =N2)

gubio = 150
for (area in 1:N2) {
  
  variances[area] = var(results[,area])
  correlations[area] = cor(results[,150], results[,area])
}

plotAreaCol("variances.png", 15, 15, variances, nigeriaAdm2, "variance")
plotAreaCol("correlations.png", 15, 15, correlations, nigeriaAdm2, "correlation")


