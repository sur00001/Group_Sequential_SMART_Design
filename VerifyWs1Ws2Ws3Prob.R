
library(stats)

# Load your data
load("C:/Users/lylae/OneDrive/Documents/Adaptive SMART paper/Group-Sequential-SMART-Research/3interim.anal.res.rda")

# Define the necessary variances and covariances (all global variables so make sure to re-run if using other scripts with similar variable names)

VAI_1 <- var(int1.res$muhat.1m1) # it's essentially the same variance for every regime estimate
VC_1 <- var(int1.res$muhat.0)
a <- cov(int1.res$muhat.11, int1.res$muhat.1m1)                                                                   
b <- cov(int1.res$muhat.m11, int1.res$muhat.m1m1)
I1 <- 1 / var(int1.res$zs)
I2 <- 1 / var(int2.res$zs)
I3 <- 1 / var(int3.res$zs)
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
lb = 10
ub = 60 #infinity takes too long and 60 is above the max ws_j simulated values anyway


#-------------------------------------------------------------------------------
##Functions to find the marginal density of W_S1
#-------------------------------------------------------------------------------

XS_S_joint_density <- function(x_s,a1,a2,a,b,VAI,VC){
  # Compute f(X_s,S)
  mean_vector <- c(0, 0) #thetas are 0 regardless of the regime bc we simulate null
  if(a1 == 1){
    cov_matrix <- matrix(c(VAI, b, b, VAI), nrow = 2)
  } else {
    cov_matrix <- matrix(c(VAI, a, a, VAI), nrow = 2)
  }
  
  # Compute the probability
  bivnorm.prob <- pmvnorm(lower = c(-Inf, -Inf), upper = c(x_s, x_s), 
                          mean = mean_vector, sigma = cov_matrix) 
  
  # Compute the mean and variance of the conditional distribution
  if(a1 == 1){
    mean_conditional <- (a / VAI) * x_s
    variance_conditional <- VAI - (a^2 / VAI)
  } else {
    mean_conditional <- (b / VAI) * x_s
    variance_conditional <- VAI - (b^2 / VAI)
  }
  
  # Compute the probability using the CDF of the normal distribution
  cond.prob <- pnorm(x_s, mean = mean_conditional, sd = sqrt(variance_conditional))
  
  f_XS_S <- dnorm(x_s, mean = 0, sd = sqrt(VAI)) * cond.prob * bivnorm.prob 
  return(as.numeric(f_XS_S))
}

WS_XS_S_joint_density <- function(x_s,w,a1,a2,a,b,VAI,VC,I) {
  f_WS_XS_S <- XS_S_joint_density(x_s,a1,a2,a,b,VAI,VC) * 
    (1 / sqrt(VC)) * dnorm((x_s - w/I) / sqrt(VC)) 
  return(as.numeric(f_WS_XS_S))
}  

WS1_S_joint_density <- function(w,a1,a2,a,b,VAI,VC,I){
  WS_XS_S_joint_density_1arg <- function(x_s) {
    sapply(x_s, function(x) WS_XS_S_joint_density(x, w, a1, a2, a, b, VAI, VC, I))
  }
  
  WS_S_joint_density <- (1/I) * integrate(WS_XS_S_joint_density_1arg, 
                                          lower = -Inf, upper = Inf)$value
  return(WS_S_joint_density)
}

# Combined function to compute the marginal density of ws1 over a range lb,ub
#Using global variables so need to ensure I have the correct components saved in the environment
f_WS1 <- function(ws1) {
  treatments <- list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
  joint_dens <- matrix(NA, nrow = length(ws1), ncol = length(treatments))
  
  for (t in 1:4) {
    treat <- treatments[[t]]
    a1 <- treat[1]
    a2 <- treat[2]
    for (j in seq_along(ws1)) {
      joint_dens[j, t] <- WS1_S_joint_density(ws1[j], a1, a2, a, b, VAI=VAI_1, VC=VC_1, I=I1)
    }
  }
  
  density_WS <- apply(joint_dens, 1, sum)
  return(density_WS)
}

#Testing that integrating f(WSI) works and matches the WS1_prob I got in the other R file
ws1.prob = integrate(f_WS1,lower = 10, upper = ub)$value 
ws1.prob #it works

#-------------------------------------------------------------------------------
# Define the density functions to do the recursive integration 
#-------------------------------------------------------------------------------
f_WS2_given_WS1 <- function(ws2, ws1) {
  dnorm(ws2, mean = ws1 + theta_S * (I2 - I1), sd = sqrt(I2 - I1))
}

f_WS3_given_WS2 <- function(ws3, ws2) {
  dnorm(ws3, mean = ws2 + theta_S * (I3 - I2), sd = sqrt(I3 - I2))
}

#-------------------------------------------------------------------------------
#Recursive integration
#-------------------------------------------------------------------------------
#Instead of Infinity, I'm doing 60 so it won't take so long
# Integrate f(W_S3 | W_S2) over W_S3
integrate_WS3 <- function(ws2_values) {
  integrate(function(ws3_values) {
    f_WS3_given_WS2(ws3_values, ws2_values)
  }, lower = lb, upper = ub)$value
}

#integrate_WS3(10) #check that function runs

# Integrate f(W_S2 | W_S1) * integrate_WS3(ws2) over W_S2
integrate_WS2 <- function(ws1_values) {
  integrate(function(ws2_values) {
    f_WS2_given_WS1(ws2_values, ws1_values) * integrate_WS3(ws2_values)
  }, lower = lb, upper = ub)$value
}

#integrate_WS2(10) #check

# Integrate f(W_S1) * integrate_WS2(ws1) over W_S1
integrate_WS1 <- function() {
  integrate(function(ws1_values) {
    f_WS1(ws1_values) * integrate_WS2(ws1_values)
  }, lower = lb, upper = ub)$value
}

# Compute the probability
probability <- integrate_WS1(); probability

#Verify with simulated ws1, ws2 and ws3 values
prop.table(table(sim.ws1>10 & sim.ws2>10 & sim.ws3>10))

#-------------------------------------------------------------------------------
# Define the marginal probability function for W_S3 > c
#-------------------------------------------------------------------------------
P_WS3_greater_than_c <- function(c) {
  integrate(function(ws2) {
    integrate(function(ws1) {
      f_WS1(ws1) * f_WS2_given_WS1(ws2, ws1)
    }, lower = -Inf, upper = Inf)$value * (1 - pnorm(c, mean = ws2 + theta_S3 * (I3 - I2), sd = sqrt(I3 - I2)))
  }, lower = -Inf, upper = Inf)$value
}

P_WS3_greater_than_c(10)
