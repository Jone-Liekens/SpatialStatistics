
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

# set domain
xl = 0
xu = 1
yl = 0
yu = 1
area <- c(xl, xu, yl, yu)
names(area) <- c("xl", "xu", "yl", "yu")


# load data
cells = ppinit('cells.dat')

redwood_data <- read.table('redwood.dat')
redwood <- list("x" = redwood_data[,1], "y" = redwood_data[,2], "area" = area)

pines_data <- read.table('pines.dat', skip=3)
pines <- list("x" = pines_data[,1], "y" = pines_data[,2], "area" = area)

# part 1: plot data
par(pty="s")
ar = cells$area
plot(cells$x, cells$y, xlim=c(ar[["xl"]],ar[["xu"]]), ylim=c(ar[["yl"]],ar[["yu"]]), xlab="x", ylab="y")
par(pty="s")
plot(redwood_data[,1], redwood_data[,2], xlim=c(0,1), ylim=c(0,1), xlab="x", ylab="y")
par(pty="s")
plot(pines_data[,1], pines_data[,2], xlim=c(0,1), ylim=c(0,1), xlab="x", ylab="y")



# part 2: calculate L-function and plot it
par(pty="s")
plot(Kfn(cells, 1, 1000), type="s", xlab="distance t", ylab="L(t)", lwd=2.0, col='blue')
lines(c(0, 1), c(0, 1), col='red', lwd=2.0)

par(pty="s")
plot(Kfn(redwood, 1, 1000), type="s", xlab="distance t", ylab="L(t)", lwd=2.0, col='blue')
lines(c(0, 1), c(0, 1), col='red', lwd=2.0)

par(pty="s")
plot(Kfn(pines, 1, 1000), type="s", xlab="distance t", ylab="L(t)", lwd=2.0, col='blue')
lines(c(0, 1), c(0, 1), col='red', lwd=2.0)



# part 3: empirical 90% interval



# cells
N = 100
l = 1000
all_L_functions = matrix(data = NA, nrow = N, ncol = l)


n_cells = length(cells$x)
for (realisation in 1:N) {
  
  x = runif(n_cells, ar[["xl"]], ar[["xu"]])
  y = runif(n_cells, ar[["yl"]], ar[["yu"]])
  
  point_obj <- list("x" = x, "y" = y, "area" = area)
  
  Lfunc = Kfn(point_obj, 0.5, 1000)
  all_L_functions[realisation,] = Lfunc$y
}

lower_bound = seq(0, 0, length = l)
upper_bound = seq(0, 0, length = l)
for (i in 1:l) {
  sorted_ <- sort(all_L_functions[,i], decreasing = FALSE)
  lower_bound[i] <- sorted_[5]
  upper_bound[i] <- sorted_[95]
}

plot(Kfn(cells, 0.5, 1000), xlim=c(0,0.5), ylim=c(0, 0.5), type="s", xlab="distance t", ylab="L(t)", col='blue')
lines(Lfunc$x, lower_bound, type="s", col='red')
lines(Lfunc$x, upper_bound, type="s", col='red')





# redwood
N = 100
l = 1000
all_L_functions = matrix(data = NA, nrow = N, ncol = l)


n_cells = length(redwood$x)
for (realisation in 1:N) {
  
  x = runif(n_cells, ar[["xl"]], ar[["xu"]])
  y = runif(n_cells, ar[["yl"]], ar[["yu"]])
  
  point_obj <- list("x" = x, "y" = y, "area" = area)
  
  Lfunc = Kfn(point_obj, 0.5, 1000)
  all_L_functions[realisation,] = Lfunc$y
}

lower_bound = seq(0, 0, length = l)
upper_bound = seq(0, 0, length = l)
for (i in 1:l) {
  sorted_ <- sort(all_L_functions[,i], decreasing = FALSE)
  lower_bound[i] <- sorted_[5]
  upper_bound[i] <- sorted_[95]
}

plot(Kfn(redwood, 0.5, 1000), xlim=c(0,0.5), ylim=c(0, 0.5), type="s", xlab="distance t", ylab="L(t)", col='blue')
lines(Lfunc$x, lower_bound, type="s", col='red')
lines(Lfunc$x, upper_bound, type="s", col='red')







# pines
N = 100
l = 1000
all_L_functions = matrix(data = NA, nrow = N, ncol = l)


n_cells = length(pines$x)
for (realisation in 1:N) {
  
  x = runif(n_cells, ar[["xl"]], ar[["xu"]])
  y = runif(n_cells, ar[["yl"]], ar[["yu"]])
  
  point_obj <- list("x" = x, "y" = y, "area" = area)
  
  Lfunc = Kfn(point_obj, 0.5, 1000)
  all_L_functions[realisation,] = Lfunc$y
}

lower_bound = seq(0, 0, length = l)
upper_bound = seq(0, 0, length = l)
for (i in 1:l) {
  sorted_ <- sort(all_L_functions[,i], decreasing = FALSE)
  lower_bound[i] <- sorted_[5]
  upper_bound[i] <- sorted_[95]
}

plot(Kfn(pines, 0.5, 1000), xlim=c(0,0.5), ylim=c(0, 0.5), type="s", xlab="distance t", ylab="L(t)", col='blue')
lines(Lfunc$x, lower_bound, type="s", col='red')
lines(Lfunc$x, upper_bound, type="s", col='red')
#lines(c(0, 1), c(0, 1), col='red')







