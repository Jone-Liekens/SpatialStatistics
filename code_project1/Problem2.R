rm(list=ls())
setwd("C:\\Users\\jonel\\Bureaublad\\spatialstat\\proj1")
library(MASS)
library(fields)
library(akima)
library(geoR)
library(viridis)
library(pracma)

# question a

top_data <- read.table('topo.dat')
x = unlist(top_data[,1])
y = unlist(top_data[,2])
q = top_data[,3]
z = unlist(top_data[,3])


interp_surface <-interp(x, y, z)

# visualize data
par(pty="s")
image.plot(interp_surface, asp = 1, xlab = "x", ylab="y")
bubblePlot(x, y, z)
contour(interp_surface$x, interp_surface$y, interp_surface$z)


# ordinary kriging
temp <- as.geodata(top_data)

locs <- expand.grid(1:315,1:315)

krigPred <- krige.conv(temp,locations = locs, krige = 
                         krige.control(type.krige = "ok", trend.d = "cte", trend.l = "cte",
                                       cov.model="powered.exponential", cov.pars = c(2500, 100, 1.5)))
prediction <- krigPred$predict


# plot ordinary kriging (question c)
par(pty="s")
image.plot(interp(locs$Var1,locs$Var2,prediction), xlab = "x", ylab="y")



mean = krigPred$predict[100*315 + 100]
var = krigPred$predict[100*315 + 100]
mean = krigPred$predict[(100-1)*315 + 100]
var = krigPred$predict[(100-1)*315 + 100]

# (question e)
probability = pnorm(850, mean, var^0.5)
height = qnorm(0.9, mean, var^0.5)






# universal kriging (question d)
krigPred <- krige.conv(temp,locations = locs, krige = 
                         krige.control(type.krige = "ok", trend.d = "2nd", trend.l = "2nd",
                                       cov.model="powered.exponential", cov.pars = c(2500, 100, 1.5)))
prediction <- krigPred$predict


# plot universla kriging
par(pty="s")
image.plot(interp(locs$Var1,locs$Var2,prediction), xlab = "x", ylab="y")







