

load("~/Adaptive SMART paper/Group-Sequential-SMART-Research/Simulated_Zsdata.rda")

library(mvtnorm)

hist(sim.zs$Zs, freq = FALSE)
# I am getting the variance and covariances from the 1000 simulations 
VAI = var(sim.zs$muhat.11)
VC = var(sim.zs$muhat.0)
a = cov(sim.zs$muhat.11,sim.zs$muhat.1m1)
b = cov(sim.zs$muhat.m11,sim.zs$muhat.m1m1)

XS_S_joint_density = function(x_s,a1,a2,a,b,VAI,VC){
  
  #------------------------------------------------------------
  # Compute f(X_s,S)
  #------------------------------------------------------------
  #-----------------------------------
  #Probability that regimes that start 
  #with different a1 are less than x_s (bivariate normal)
  #-----------------------------------
  # Define the mean vector and covariance matrix
  mean_vector <- c(0, 0) #thetas are 0 regardless of the regime bc we simulate null
  if(a1==1){cov_matrix <- matrix(c(VAI, b, b, VAI), nrow = 2) } else {cov_matrix <- matrix(c(VAI, a, a, VAI), nrow = 2)}
  
  # Compute the probability
  bivnorm.prob <- pmvnorm(lower = c(-Inf, -Inf), upper = c(x_s, x_s), 
                          mean = mean_vector, sigma = cov_matrix) #plugging in xs value from that simulation
  
  #-----------------------------------
  #P(X_a1,a2' < Xs| X_a1,a2 = x)
  #-----------------------------------
  # Compute the mean and variance of the conditional distribution
  if(a1==1){
    mean_conditional <- 0 + (a / VAI) * (x_s - 0) #this x is x_a1,a2 from a given simulation
    variance_conditional <- VAI - (a^2 / VAI)
  } else {
    mean_conditional <- 0 + (b / VAI) * (x_s - 0) #this x is x_a1,a2 
    variance_conditional <- VAI - (b^2 / VAI)}
  
  # Compute the probability using the CDF of the normal distribution
  cond.prob <- pnorm(x_s, mean = mean_conditional, sd = sqrt(variance_conditional))
  
  f_XS_S = dnorm(x_s, mean = 0, sd = sqrt(VAI))*cond.prob*bivnorm.prob #theta_a1,a2 is 0 in the null case 
  
  return(as.numeric(f_XS_S))
}

hist(sim.zs$Xs, freq = FALSE)



XS_S_joint_density(0,a1=1,a2=1,a,b,VAI=VAI,VC = VC)


treatments <- list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
density_range <- seq(from = -0.1, to = 0.20, by = 0.01) #use the range of Xs from the simulation

joint_dens <- matrix(NA, nrow = length(density_range), ncol = 4)

for (t in 1:4) {
  treat <- treatments[[t]]
  a1 <- treat[1]
  a2 <- treat[2]
  for(i in 1:length(density_range)) {
    joint_dens[i, t] <- XS_S_joint_density(density_range[i],a1,a2,a,b,VAI=VAI,VC = VC)
  }
  print(t)
}

density_XS <- apply(joint_dens, 1, sum) #marginal density of Xs

hist(sim.zs$Xs, freq = FALSE)
lines(density_range, density_XS)

ZS_XS_S_joint_density <- function(x_s,z,a1,a2,a,b,VAI,VC) {
  f_ZS_XS_S <- XS_S_joint_density(x_s,a1,a2,a,b,VAI,VC) * 
    (1 / sqrt(VC)) * dnorm((x_s - z) / sqrt(VC)) #Use xs instead of x
  return(as.numeric(f_ZS_XS_S))
}  

ZS_S_joint_density <- function(z,a1,a2,a,b,VAI,VC){
  ZS_XS_S_joint_density_1arg <- function(x_s) {
    ZS_XS_S_joint_density_1arg <- NULL
    for(j in 1:length(x_s)) {
      ZS_XS_S_joint_density_1arg <- c(ZS_XS_S_joint_density_1arg,
                                      ZS_XS_S_joint_density(x_s[j], 
                                                            z = z, a1 = a1, a2 = a2, a = a, b = b, VAI = VAI, VC = VC))
    }
    return(ZS_XS_S_joint_density_1arg)
  }
  #ZS_S_joint_density <- integrate(ZS_XS_S_joint_density_1arg ,
                                  #lower = -0.5, upper = 0.5)$value
  
  ZS_S_joint_density <- integrate(ZS_XS_S_joint_density_1arg ,
  lower = -Inf, upper = Inf)$value
  
  return(ZS_S_joint_density)
}



#ZS_XS_S_joint_density(x_s = 0,z = 0, a1 = 1,a2 = 1, a,b,VAI,VC)
ZS_S_joint_density(z = 0, a1 = 1, a2 = 1, a = a, b = b, VAI = VAI, VC = VC)

treatments <- list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
density_range <- seq(from = -0.2, to = 0.30, by = 0.01) #range of Zs from data

joint_dens <- matrix(NA, nrow = length(density_range), ncol = 4)

for (t in 1:4) {
  treat <- treatments[[t]]
  a1 <- treat[1]
  a2 <- treat[2]
  for(i in 1:length(density_range)) {
    joint_dens[i, t] <- ZS_S_joint_density(density_range[i],a1,a2,a,b,VAI=VAI,VC = VC)
  }
  print(t)
}

density_ZS <- apply(joint_dens, 1, sum)

hist(sim.zs$Zs, freq = FALSE)
lines(density_range, density_ZS)