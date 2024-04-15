rm(list=ls())
library(MASS)
library(fields)
library(akima)
library(geoR)
library(viridis)
library(pracma)

# generate grid
x1 <- rep(1:30, 30)
x2 <- rep(1:30, each = 30)  

# generate covariance matrices
cov_mx <- matrix(NA, nrow=length(x1), ncol=length(x1))
for (i in 1:length(x1))(
  for (j in 1:length(x1)){
    h = ((x1[i] - x1[j])^2 + (x2[i] - x2[j])^2)^0.5
    s2 = 2
    a = 3
    cov_mx[i,j] <- s2 * exp(-h/a)
  }
)
mean = 0 * x1

# generate sample
z = mvrnorm(1, mean, cov_mx)


# convert to correct format for plotting
data <- matrix(data=NA,nrow = length(x1), ncol = 3)
data[,1] = x1
data[,2] = x2
data[,3] = z
interp_surface <-interp(x1, x2, z)

# plot it
par(pty="s")
image.plot(interp_surface, asp = 1, xlab = "x", ylab="y")
bubblePlot(x1, x2, z)



# convert to correct format and use variog to generate variogram 
data_df <- as.data.frame(data)
temp <- as.geodata(data_df)
vario.b <- variog(temp)

# plot it
dist = seq(0, 100, length.out = 10000)
real_variogram = 2 * exp(0) - 2 * exp(-(dist / 3))
plot(vario.b, ylim = c(0, 3), type='l', col='red') 
lines(dist, real_variogram, col='blue')

q = 3;
p = 9;

# maximum likelihood estimation of parameters exponential cov function
ml = likfit(temp, ini = c(1, 2))
ml
#plot exponential semi-variogram
estimated_variogram = ml$cov.pars[1] * exp(0) - ml$cov.pars[1] * exp(-(dist / ml$cov.pars[2]))
lines(dist, estimated_variogram, col='green')

# repeat for 36 points instead of whole data set
locations36 <- sample(1:900, 36, replace=F)
new_data = data[locations36,]
new_data_df <- as.data.frame(new_data)
new_temp <- as.geodata(new_data_df)
vario.b <- variog(new_temp)
dist = seq(0, 100, length.out = 10000)
real_variogram = 2 * exp(0) - 2 * exp(-(dist / 3))

data <- matrix(data=NA,nrow = length(x1), ncol = 3)
data[,1] = x1
data[,2] = x2
data[,3] = z
interp_surface <-interp(x1, x2, z)

# plot it
par(pty="s")
image.plot(interp_surface, asp = 1, xlab = "x", ylab="y")

plot(vario.b, type='l', col='red') 
lines(dist, real_variogram, col='blue')

adjlfjlak = 17

ml = likfit(new_temp, ini = c(1, 2))
ml
estimated_variogram = ml$cov.pars[1] * exp(0) - ml$cov.pars[1] * exp(-(dist / ml$cov.pars[2]))
lines(dist, estimated_variogram, col='green')


# repeat for 9, 64 and 100 points instead of whole data set
locations9 <- sample(1:900, 9, replace=F)
new_data = data[locations9,]
new_data_df <- as.data.frame(new_data)
new_temp <- as.geodata(new_data_df)
vario.b <- variog(new_temp)
dist = seq(0, 100, length.out = 10000)
real_variogram = 2 * exp(0) - 2 * exp(-(dist / 3))
plot(vario.b, type='l', col='red') 
lines(dist, real_variogram, col='blue')

ml = likfit(new_temp, ini = c(1, 2))
ml
estimated_variogram = ml$cov.pars[1] * exp(0) - ml$cov.pars[1] * exp(-(dist / ml$cov.pars[2]))
lines(dist, estimated_variogram, col='green')



locations64 <- sample(1:900, 64, replace=F)
new_data = data[locations64,]
new_data_df <- as.data.frame(new_data)
new_temp <- as.geodata(new_data_df)
vario.b <- variog(new_temp)
dist = seq(0, 100, length.out = 10000)
real_variogram = 2 * exp(0) - 2 * exp(-(dist / 3))
plot(vario.b, type='l', col='red') 
lines(dist, real_variogram, col='blue')

ml = likfit(new_temp, ini = c(1, 2))
ml
estimated_variogram = ml$cov.pars[1] * exp(0) - ml$cov.pars[1] * exp(-(dist / ml$cov.pars[2]))
lines(dist, estimated_variogram, col='green')



locations100 <- sample(1:900, 100, replace=F)
new_data = data[locations100,]
new_data_df <- as.data.frame(new_data)
new_temp <- as.geodata(new_data_df)
vario.b <- variog(new_temp)
dist = seq(0, 100, length.out = 10000)
real_variogram = 2 * exp(0) - 2 * exp(-(dist / 3))
plot(vario.b, type='l', col='red') 
lines(dist, real_variogram, col='blue')

ml = likfit(new_temp, ini = c(1, 2))
ml
estimated_variogram = ml$cov.pars[1] * exp(0) - ml$cov.pars[1] * exp(-(dist / ml$cov.pars[2]))
lines(dist, estimated_variogram, col='green')

