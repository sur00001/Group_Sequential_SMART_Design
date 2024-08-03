#Adaptive quadrature 
#Takes way too long, need to parallelize 

library(stats)
library(mvtnorm)
library(ggplot2)
library(cubature)

# Load your data
load("C:/Users/lylae/OneDrive/Documents/Adaptive SMART paper/Group-Sequential-SMART-Research/3interim.anal.res.rda")

# Define the necessary variances and covariances (all global variables so make sure to re-run if using other scripts with similar variable names)
VAI_1 = var(int1.res$muhat.1m1) # it's essentially the same variance for every regime estimate
VC_1 = var(int1.res$muhat.0)
a = cov(int1.res$muhat.11, int1.res$muhat.1m1)                                                                   
b = cov(int1.res$muhat.m11, int1.res$muhat.m1m1)
I1 = 1 / var(int1.res$zs)
I2 = 1 / var(int2.res$zs)
I3 = 1 / var(int3.res$zs)
sim.ws1 = int1.res$ws1 = int1.res$zs*I1
sim.ws2 = int2.res$ws2 = int2.res$zs*I2
sim.ws3 = int3.res$ws3 = int3.res$zs*I3

hist(int1.res$zs,main="Interim Analysis 1 zs")
hist(int2.res$zs,main="Interim Analysis 2 zs")
hist(int3.res$zs,main="Interim Analysis 3 zs")
hist(sim.ws1)
hist(sim.ws2)
hist(sim.ws3)

theta_S = 0 


# Define the conditional density functions
f_WS2_given_WS1 = function(ws2, ws1) {
  dnorm(ws2, mean = ws1, sd = sqrt(I2 - I1))
}

f_WS3_given_WS2 = function(ws3, ws2) {
  dnorm(ws3, mean = ws2, sd = sqrt(I3 - I2))
}

# Define the marginal density function for WS1
f_WS1 = function(ws1) {
  treatments = list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
  joint_dens = matrix(NA, nrow = length(ws1), ncol = length(treatments))
  
  for (t in 1:4) {
    treat = treatments[[t]]
    a1 = treat[1]
    a2 = treat[2]
    for (j in seq_along(ws1)) {
      joint_dens[j, t] = WS1_S_joint_density(ws1[j], a1, a2, a, b, VAI=VAI_1, VC=VC_1, I=I1)
    }
  }
  
  density_WS = apply(joint_dens, 1, sum)
  return(density_WS)
}

# Compute the joint density f(W_S3, W_S2, W_S1) using cubature
joint_density_WS3_WS2_WS1 = function(ws) {
  ws1 = ws[1]
  ws2 = ws[2]
  ws3 = ws[3]
  f_WS1(ws1) * f_WS2_given_WS1(ws2, ws1) * f_WS3_given_WS2(ws3, ws2)
}

# Define the integration bounds and tolerance
lower_bounds = c(10, 10, 10)
upper_bounds = c(1000, 1000, 1000)  # Large but finite upper bounds

# Compute the joint probability P(W_S1 > lb, W_S2 > lb, W_S3 > lb)
result = adaptIntegrate(joint_density_WS3_WS2_WS1, lowerLimit = lower_bounds, upperLimit = upper_bounds, tol = 1e-10)

# Print the result
result$integral
