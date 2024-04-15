#getwd()
setwd("C:\\Users\\jonel\\Bureaublad\\spatialstat\\proj2")
rm(list=ls())
library(MASS)
library(fields)
library(akima)
library(geoR)
library(viridis)
library(pracma)
library(spatial)

xl = 0
xu = 1
yl = 0
yu = 1
area <- c(xl, xu, yl, yu)
names(area) <- c("xl", "xu", "yl", "yu")

redwood_data <- read.table('redwood.dat')
redwood <- list("x" = redwood_data[,1], "y" = redwood_data[,2], "area" = area)

par(pty="s")
ar = redwood$area
plot(redwood$x, redwood$y, xlim=c(ar[["xl"]],ar[["xu"]]), ylim=c(ar[["yl"]],ar[["yu"]]), xlab="x", ylab="y")


# neyman-scott process

#lambda_p = 10
#lambda_n = 62/lambda_p
#sigma_2 = 0.06

lambda_p = 20
lambda_n = 3.1
sigma_2 = 0.03

space_surface = 9


N = 100
l = 1000
all_L_functions = matrix(data = NA, nrow = N, ncol = l)

for (realisation in 1:N) {
  # generate number of parents:
  n_parents = rpois(1, space_surface * lambda_p)
  
  # generate locations of parents
  parent_x = runif(n_parents, xl-1, xu+1)
  parent_y = runif(n_parents, yl-1, yu+1)
  
  # for each parent, generate number of children
  n_child = rpois(n_parents, lambda_n)
  
  # for each child, generate a distance in both x and y
  total_child = sum(n_child)
  child_x = seq(0, 0, length=total_child)
  child_y = seq(0, 0, length=total_child)
  child_i = 1
  for (i in 1:n_parents){
    for (j in 1:n_child[i]) {
      dx <- rnorm(1, 0, sigma_2)
      dy <- rnorm(1, 0, sigma_2)
      child_x[child_i] <- parent_x[i] + dx
      child_y[child_i] <- parent_y[i] + dy
      child_i <- child_i + 1
    }
  }
  
  if (N - realisation < 3) {
    plot(child_x, child_y, xlim=c(xl, xu), ylim=c(yl, yu), xlab='x', ylab='y')
  }
  
  
  child_y = child_y[child_x < 1]
  child_x = child_x[child_x < 1]
  child_y = child_y[child_x > 0]
  child_x = child_x[child_x > 0]
  
  child_x = child_x[child_y < 1]
  child_y = child_y[child_y < 1]
  child_x = child_x[child_y > 0]
  child_y = child_y[child_y > 0]
  
  
  area <- c(xl, xu, yl, yu)
  names(area) <- c("xl", "xu", "yl", "yu")
  child_obj <- list("x" = child_x, "y" = child_y, "area" = area)
  
  
  par(pty="s")
  ppregion(xl, xu, yl, yu)
  temp = Kfn(child_obj, 0.7, 1000)
  
  all_L_functions[realisation,] = temp$y
}

lower_bound = seq(0, 0, length = l)
upper_bound = seq(0, 0, length = l)
for (i in 1:l) {
  sorted_ <- sort(all_L_functions[,i], decreasing = FALSE)
  lower_bound[i] <- sorted_[5]
  upper_bound[i] <- sorted_[95]
}

if (N > 10) {
  ppregion(0, 1, 0, 1)
  plot(Kfn(redwood, 0.7, 1000), type="s", xlab="distance", ylab="L(t)", col='blue', lwd=2)
  lines(temp$x, lower_bound, type="s", col='red', lwd=2)
  lines(temp$x, upper_bound, type="s", col='red', lwd=2)
  lines(c(0, 1), c(0,1), col='black')
}



