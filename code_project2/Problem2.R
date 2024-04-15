#getwd()
setwd("C:\\Users\\jonel\\Bureaublad\\spatialstat\\proj2")
rm(list=ls())

library(spatial)
library(MASS)
library(fields)
library(akima)
library(geoR)
library(viridis)
library(pracma)
library(ggplot2)


probs <- read.table('obsprob.txt')
pines <- read.table('obspines.txt')

x = probs[,1]
prob_x = as.numeric(x[-1])
y = probs[,2]
prob_y = as.numeric(y[-1])
z = probs[,3]
prob_z = as.numeric(z[-1])

x = pines[,1]
pine_x = as.numeric(x[-1])
y = pines[,2]
pine_y = as.numeric(y[-1])
z = pines[,3]
pine_z = as.numeric(z[-1])

par(pty="s")
bubblePlot(prob_x, prob_y, prob_z)

par(pty="s")
alpha = prob_z
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)) + 
  geom_raster(aes(fill = alpha)) + labs(y= "y", x = "x")

par(pty="s")
pines = pine_z
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)) + 
  geom_raster(aes(fill = pines)) + labs(y= "y", x = "x")



# question c)

#sum(pine_z)
C = 1 / (sum(prob_z) / 30 / 30)
Lamb = C * sum(pine_z) / 300 / 300


xl = 0
xu = 300
yl = 0
yu = 300
space_surface = (xu - xl) * (yu - yl)


bins = matrix(data=0, nrow = 30, ncol = 30)
lambda = Lamb
dx = 10
dy = 10
lambda_0 = lambda * dx * dy
for (x_i in 1:30) {
  for (y_i in 1:30) {
    bins[x_i, y_i] <- rpois(1, lambda_0)
  }
}
num1 = sum(sum(bins))

trees = c(bins)
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)
) + 
  geom_raster(aes(fill = trees)) + labs(y= "y", x = "x") +
  scale_fill_gradientn(colours = hcl.colors(50, "YlGn"), limits = c(0, 7))
  

xs = seq(0, 0, length=num1)
ys = seq(0, 0, length=num1)
count = 1
for (x_i in 1:30) {
  for (y_i in 1:30) {
    n_trees = bins[x_i, y_i]
    for (tree in n_trees) {
      xs[count] = (x_i - 1) * 10 + 5 + runif(1, 0, 10)
      ys[count] = (y_i - 1) * 10 + 5 + runif(1, 0, 10)
      count = count + 1
    }
  }
}
par(pty="s")
plot(xs, ys, xlab='x', ylab='y')






# d)
bins = matrix(data=0, nrow = 30, ncol = 30)
lambda = Lamb
dx = 10
dy = 10
lambda_0 = lambda * dx * dy
for (x_i in 1:30) {
  for (y_i in 1:30) {
    alph = prob_z[30 * (y_i - 1) + x_i]
    obs = pine_z[30 * (y_i - 1) + x_i]
    bins[x_i, y_i] <- obs + rpois(1, lambda_0*(1-alph))
  }
}
num2 = sum(sum(bins))

trees = c(bins)
par(pty="s")
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)
) + 
  geom_raster(aes(fill = trees)) + labs(y= "y", x = "x") +
  scale_fill_gradientn(colours = hcl.colors(50, "YlGn"), limits = c(0, 7))



xs = seq(0, 0, length=num1)
ys = seq(0, 0, length=num1)
count = 1
for (x_i in 1:30) {
  for (y_i in 1:30) {
    n_trees = bins[x_i, y_i]
    for (tree in n_trees) {
      xs[count] = (x_i - 1) * 10 + 5 + runif(1, 0, 10)
      ys[count] = (y_i - 1) * 10 + 5 + runif(1, 0, 10)
      count = count + 1
    }
  }
}
par(pty="s")
plot(xs, ys, xlab='x', ylab='y')




# e)

N = 500
bins_c <- array(rep(0, N*30*30), c(N, 30, 30));  
bins_d <- array(rep(0, N*30*30), c(N, 30, 30));  

for (n in 1:N) {
  for (x_i in 1:30) {
    for (y_i in 1:30) {
      alph = prob_z[30 * (y_i - 1) + x_i]
      obs = pine_z[30 * (y_i - 1) + x_i]
      bins_c[n, x_i, y_i] <- rpois(1, lambda_0)
      bins_d[n, x_i, y_i] <- obs + rpois(1, lambda_0*(1-alph))
    }
  }
}


means_c = matrix(data=NA, nrow = 30, ncol = 30)
var_c = matrix(data=NA, nrow = 30, ncol = 30)
means_d = matrix(data=NA, nrow = 30, ncol = 30)
var_d = matrix(data=NA, nrow = 30, ncol = 30)

for (x_i in 1:30) {
  for (y_i in 1:30) {
    vals_c = bins_c[,x_i, y_i]
    means_c[x_i, y_i] = mean(vals_c)
    var_c[x_i, y_i] = var(vals_c)
    
    vals_d = bins_d[,x_i, y_i]
    means_d[x_i, y_i] = mean(vals_d)
    var_d[x_i, y_i] = var(vals_d)
  }
}


lim = max(c(means_c, means_d))
mean = c(means_c)
par(pty="s")
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)
) + 
  geom_raster(aes(fill = mean)) + labs(y= "y", x = "x") +
  scale_fill_gradientn(colours = hcl.colors(50, "YlGn"), limits = c(0, lim))



mean = c(means_d)
par(pty="s")
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)
) + 
  geom_raster(aes(fill = mean)) + labs(y= "y", x = "x") +
  scale_fill_gradientn(colours = hcl.colors(50, "YlGn"), limits = c(0, lim))



lim = max(c(var_c, var_d))
variance = c(var_c)
par(pty="s")
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)
) + 
  geom_raster(aes(fill = variance)) + labs(y= "y", x = "x") +
  scale_fill_gradientn(colours = hcl.colors(50, "YlGn"), limits = c(0, lim))



variance = c(var_d)
par(pty="s")
ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)
) + 
  geom_raster(aes(fill = variance)) + labs(y= "y", x = "x") +
  scale_fill_gradientn(colours = hcl.colors(50, "YlGn"), limits = c(0, lim))












#n_points = rpois(1, space_surface * lambda)

# binning
#x = runif(n_points, xl, xu)
#y = runif(n_points, yl, yu)

#bins = matrix(data=0, nrow = 30, ncol = 30)
#for (i in 1:n_points) {
#  x_i = x[i] %/% 10 + 1
#  y_i = y[i] %/% 10 + 1
#  s = bins[x_i, y_i]
#  bins[x_i, y_i] = s + 1
#}

#sum(sum(bins))

#alpha = c(bins)
#ggplot(data.frame(prob_x,prob_y), aes(prob_x,prob_y)) + 
#  geom_raster(aes(fill = alpha)) +labs(y= "y", x = "x", z='a')


