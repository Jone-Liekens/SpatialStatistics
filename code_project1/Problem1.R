rm(list=ls())

library(MASS)
library(fields)
library(akima)
library(geoR)
library(viridis)
library(pracma)


# question a
h = seq(0, 100, length.out = 10000)
a = 10

fs <- matrix(NA, nrow = 8, ncol = length(h))
alfa = 1.9
s2 = 1
fs[1,] = s2 * exp(-(h / a)^alfa)
alfa = 1
s2 = 1
fs[2,] = s2 * exp(-(h / a)^alfa)
alfa = 1.9
s2 = 5
fs[3,] = s2 * exp(-(h / a)^alfa)
alfa = 1
s2 = 5
fs[4,] = s2 * exp(-(h / a)^alfa)

kapp = 1 # = nu
a = 20
phi1 = 20 / ((1 * 8)^0.5)
phi3 = 20 / ((3 * 8)^0.5)
fs[5,] <- matern(h, 1, phi1)
fs[6,] <- matern(h, 3, phi3)
fs[7,] <- 5 * matern(h, 1, phi1)
fs[8,] <- 5 * matern(h, 3, phi3)

fs[5,] <- matern(h, phi1, 1)
fs[6,] <- matern(h, phi3, 3)
fs[7,] <- 5 * matern(h, phi1, 1)
fs[8,] <- 5 * matern(h, phi3, 3)


# plot correlation functions with sigm2 = 1
plot(h, fs[5,], xlim=c(0,20), ylim=c(0,1), xlab="h", ylab="C(h)", type="l", col='blue')
lines(h, fs[1,],  type="l", col='orange')
lines(h, fs[2,], type="l", col='red')
lines(h, fs[5,], type="l", col='blue')
lines(h, fs[6,], type="l", col='green')
legend(x = c(20/2, 20), y = c(0.67, 1), col=c("orange", "red", "blue", "green"), pch = c(15,15,15,15), legend=c(
  expression(paste("Powered exponential, ", alpha, "= 1.9")),
  expression(paste("Powered exponential, ", alpha, "= 1")),
  expression(paste("Matern, ", kappa, " = 3")),
  expression(paste("Matern, ", kappa, " = 1"))
), bty="o", cex=0.8, text.font=4)


# plot correlation functions with sigm2 = 5
plot(h, fs[3,], xlim=c(0,20), ylim=c(0,5), xlab="h", ylab="C(h)", type="l", col='blue')
lines(h, fs[3,],  type="l", col='orange')
lines(h, fs[4,], type="l", col='red')
lines(h, fs[7,], type="l", col='blue')
lines(h, fs[8,], type="l", col='green')
legend(x = c(20/2, 20), y = c(0.67*5, 1*5), col=c("orange", "red", "blue", "green"), pch = c(15,15,15,15), legend=c(
  expression(paste("Powered exponential, ", alpha, "= 1.9")),
  expression(paste("Powered exponential, ", alpha, "= 1")),
  expression(paste("Matern, ", kappa, " = 3")),
  expression(paste("Matern, ", kappa, " = 1"))
), bty="o", cex=0.8, text.font=4)



# question b
# generate covariance matrices
X = seq(1, 50, 1)
cov_mx1 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx2 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx3 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx4 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx5 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx6 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx7 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_mx8 <- matrix(NA, nrow=length(X), ncol=length(X))


for (i in 1:length(X))(
  for (j in 1:length(X)){
    h <- abs(X[i] - X[j])
    
    alfa <- 1.9
    s2 <- 1
    cov_mx1[i,j] <- s2 * exp(-(h / a)^alfa)
    alfa <- 1
    s2 <- 1
    cov_mx2[i,j] <- s2 * exp(-(h / a)^alfa)
    alfa <- 1.9
    s2 <- 5
    cov_mx3[i,j] <- s2 * exp(-(h / a)^alfa)
    alfa <- 1
    s2 <- 5
    cov_mx4[i,j] <- s2 * exp(-(h / a)^alfa)
    
    phi1 = 20 / ((1 * 8)^0.5)
    phi3 = 20 / ((3 * 8)^0.5)
    cov_mx5[i,j] <- matern(h, 1, phi1)
    cov_mx6[i,j] <- matern(h, 3, phi3)
    cov_mx7[i,j] <- 5 * matern(h, 1, phi1)
    cov_mx8[i,j] <- 5 * matern(h, 3, phi3)
  }
)


# generate samples
z_value <- seq(0, 0, length.out = length(X))
sample_1 <- mvrnorm(4, z_value, cov_mx1)
sample_2 <- mvrnorm(4, z_value, cov_mx2)
sample_3 <- mvrnorm(4, z_value, cov_mx3)
sample_4 <- mvrnorm(4, z_value, cov_mx4)
sample_5 <- mvrnorm(4, z_value, cov_mx5)
sample_6 <- mvrnorm(4, z_value, cov_mx6)
sample_7 <- mvrnorm(4, z_value, cov_mx7)
sample_8 <- mvrnorm(4, z_value, cov_mx8)

# plot samples
plot(X, sample_1[1,], ylim=c(-5,5),  xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_1[2,], type="l", col='red')
lines(X, sample_1[3,], type="l", col='green')
lines(X, sample_1[4,], type="l", col='blue')

plot(X, sample_2[1,], ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_2[2,], type="l", col='red')
lines(X, sample_2[3,], type="l", col='green')
lines(X, sample_2[4,], type="l", col='blue')

plot(X, sample_3[1,],ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_3[2,], type="l", col='red')
lines(X, sample_3[3,], type="l", col='green')
lines(X, sample_3[4,], type="l", col='blue')

plot(X, sample_4[1,],ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_4[2,], type="l", col='red')
lines(X, sample_4[3,], type="l", col='green')
lines(X, sample_4[4,], type="l", col='blue')

plot(X, sample_5[1,],ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_5[2,], type="l", col='red')
lines(X, sample_5[3,], type="l", col='green')
lines(X, sample_5[4,], type="l", col='blue')

plot(X, sample_6[1,],ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_6[2,], type="l", col='red')
lines(X, sample_6[3,], type="l", col='green')
lines(X, sample_6[4,], type="l", col='blue')

plot(X, sample_7[1,],ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_7[2,], type="l", col='red')
lines(X, sample_7[3,], type="l", col='green')
lines(X, sample_7[4,], type="l", col='blue')

plot(X, sample_8[1,],ylim=c(-5,5),xlab="s", ylab="X(s)", type="l", col='orange')
lines(X, sample_8[2,], type="l", col='red')
lines(X, sample_8[3,], type="l", col='green')
lines(X, sample_8[4,], type="l", col='blue')



# question c
# observations
y1 = sample_3[1, 10]
y2 = sample_3[1, 25]
y3 = sample_3[1, 30]
Y_no_noise = c(y1, y2, y3)

X = seq(1, 50, 1)
Y = c(10, 25, 30)
cov_00_3 <- matrix(NA, nrow=length(X), ncol=length(X))
cov_01_3 <- matrix(NA, nrow=length(X), ncol=length(Y))
cov_10_3 <- matrix(NA, nrow=length(Y), ncol=length(X))
cov_11_3 <- matrix(NA, nrow=length(Y), ncol=length(Y))

for (i in 1:length(X))(
  for (j in 1:length(X)){
    h <- abs(X[i] - X[j])
    alfa <- 1.9
    s2 <- 5
    cov_00_3[i,j] <- s2 * exp(-(h / a)^alfa)
  }
)

for (i in 1:length(Y))(
  for (j in 1:length(Y)){
    h <- abs(Y[i] - Y[j])
    alfa <- 1.9
    s2 <- 5
    cov_11_3[i,j] <- s2 * exp(-(h / a)^alfa)
  }
)

for (i in 1:length(X))(
  for (j in 1:length(Y)){
    h <- abs(X[i] - Y[j])
    alfa <- 1.9
    s2 <- 5
    cov_01_3[i,j] <- s2 * exp(-(h / a)^alfa)
    cov_10_3[j,i] <- cov_01_3[i,j]
  }
)

# calculate BLUP
inv_cov_11_3 <- solve(cov_11_3)
mean = cov_01_3 %*% inv_cov_11_3 %*% Y_no_noise
cov = cov_00_3 - cov_01_3 %*% inv_cov_11_3 %*% cov_10_3
upper_ = mean + diag(cov)^0.5 * 1.645
lower_ = mean - diag(cov)^0.5 * 1.645


# plot estimations
plot(X, xlim=c(0,50), ylim=c(-7, 7), xlab="s", ylab="X(s)", upper_, type='l', col='green')
lines(X, mean,  pch = c(20,20,20), col= 'blue')
lines(X, lower_,  pch = c(20,20,20), col= 'green')
points(Y, Y_no_noise,  pch = 15, col= 'red')
lines(X, sample_3[1,])
legend(x = c(0, 18), y = c(-7, -3), col=c("black", "red", "blue", "green"), pch = c(15,15,15,15), legend=c(
  expression(paste("Real data")),
  expression(paste("Observations")),
  expression(paste("Estimated mean")),
  expression(paste("90% interval"))
), bty="o", cex=0.8, text.font=4)


# repeat estimations 100 times
samples_x <- mvrnorm(100, mean, cov)

lower_bound <- matrix(NA, 1, length(X))
upper_bound <- matrix(NA, 1, length(X))
for (i in 1:length(X)){
  sorted_x <- sort(samples_x[,i], decreasing = FALSE)
  lower_bound[i] <- sorted_x[5]
  upper_bound[i] <- sorted_x[95]
}

# plot estimations again for empiricial prediction interval
plot(X, mean, xlim=c(0,50), ylim=c(-7, 7), xlab="s", ylab="X(s)", type='l', col='blue')
points(Y, Y_no_noise,  pch = 15, col= 'red')
lines(X, lower_bound,  pch = c(20,20,20), col= 'green')
lines(X, upper_bound,  pch = c(20,20,20), col= 'green')
lines(X, sample_3[1,])
legend(x = c(0, 18), y = c(-7, -3), col=c("black", "red", "blue", "green"), pch = c(15,15,15,15), legend=c(
  expression(paste("Real data")),
  expression(paste("Observations")),
  expression(paste("Estimated mean")),
  expression(paste("90% interval"))
), bty="o", cex=0.8, text.font=4)



# repeat all of this but with noise (standard variance 0.25)
y1 = sample_3[1, 10]
y2 = sample_3[1, 25]
y3 = sample_3[1, 30]
y1 = y1 + rnorm(1, 0, 0.25)
y2 = y2 + rnorm(1, 0, 0.25)
y3 = y3 + rnorm(1, 0, 0.25)
Y_noise = c(y1, y2, y3)

sigm2_noise = 0.25 

cov_11_3 = cov_11_3 + diag(3) * sigm2_noise

inv_cov_11_3 <- solve(cov_11_3)
mean = cov_01_3 %*% inv_cov_11_3 %*% Y_noise
cov = cov_00_3 - cov_01_3 %*% inv_cov_11_3 %*% cov_10_3
cov = cov
upper_ = mean + diag(cov)^0.5 * 1.645
lower_ = mean - diag(cov)^0.5 * 1.645

plot(X, xlim=c(0,50), ylim=c(-7, 7), xlab="s", ylab="X(s)", upper_, type='l', col='green')
lines(X, mean,  pch = c(20,20,20), col= 'blue')
lines(X, lower_,  pch = c(20,20,20), col= 'green')
points(Y, Y_noise,  pch = 15, col= 'red')
lines(X, sample_3[1,])
legend(x = c(0, 18), y = c(-7, -3), col=c("black", "red", "blue", "green"), pch = c(15,15,15,15), legend=c(
  expression(paste("Real data")),
  expression(paste("Observations")),
  expression(paste("Estimated mean")),
  expression(paste("90% interval"))
), bty="o", cex=0.8, text.font=4)


samples_x <- mvrnorm(100, mean, cov)

lower_bound <- matrix(NA, 1, length(X))
upper_bound <- matrix(NA, 1, length(X))
for (i in 1:length(X)){
  sorted_x <- sort(samples_x[,i], decreasing = FALSE)
  lower_bound[i] <- sorted_x[5]
  upper_bound[i] <- sorted_x[95]
}

plot(X, mean, xlim=c(0,50), ylim=c(-7, 7), xlab="s", ylab="X(s)", type='l', col='blue')
#points(Y, Y_no_noise,  pch = c(20,20,20), col= 'red')
points(Y, Y_noise,  pch = 15, col= 'red')
lines(X, lower_bound,  pch = c(20,20,20), col= 'green')
lines(X, upper_bound,  pch = c(20,20,20), col= 'green')
lines(X, sample_3[1,])
legend(x = c(0, 18), y = c(-7, -3), col=c("black", "red", "blue", "green"), pch = c(15,15,15,15), legend=c(
  expression(paste("Real data")),
  expression(paste("Observations")),
  expression(paste("Estimated mean")),
  expression(paste("90% interval"))
), bty="o", cex=0.8, text.font=4)



# question f

samples_x <- mvrnorm(100, mean, cov)

# calculate empirically for 100 samples
sums = seq(0, 0, 100)
for (sample in 1:100) {
  sum <- 0
  for (s in 1:50) {
    if (samples_x[sample, s] > 2) {
      sum = sum + samples_x[sample, s] - 2
    }
  }
  sums[sample] = sum
}
mu = mean(sums)
sigm = var(sums)


# calculate it theoretically using integral of probablity
total = 0
probabs = seq(0, 0, 50)
for (s in 1:50) {
  var <- cov[s,s]
  mu <- mean[s]
  if (var > 1e-10) {
    
    f <- function(x) {(x-2)/sqrt(2*pi*var)*exp(-((x-mu)/var^0.5)^2/2)}
    q <- integrate(f, lower = 2, upper = 100)
    
    total = total + q$value
  } else {
    if (mu > 2) {
      total = total + mu - 2
    }
  }
}
total


xs = seq(-4, 4, length = 1000)
I = xs > 2
ys = I * (xs - 2)
plot(xs, ys, type='l', xlab='X(s)', ylab= 'A(X(s))')


