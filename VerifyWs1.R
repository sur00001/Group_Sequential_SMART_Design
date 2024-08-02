#This code verifies f(W_S1) (derived from transforming f(Z_S1) ) and verifies Pr(W_S1)

load("~/Adaptive SMART paper/Group-Sequential-SMART-Research/Simulated_Zsdata.rda")

library(mvtnorm)

# I am getting the variance and covariances from the 1000 simulations 
sim.zs1 = sim.zs 
VAI_1 = var(sim.zs1$muhat.1m1) #it's essentially the same variance for every regime estimate
VC_1 = var(sim.zs1$muhat.0)
I_1 = 1/(var(sim.zs1$Zs))
a = cov(sim.zs1$muhat.11,sim.zs1$muhat.1m1)
b = cov(sim.zs1$muhat.m11,sim.zs1$muhat.m1m1)
sim.zs1$ws1 = sim.zs1$Zs*I_1
hist(sim.zs1$ws1, freq = FALSE)

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

WS_XS_S_joint_density <- function(x_s,w,a1,a2,a,b,VAI,VC,I) {
  #Replace z with w/I where I is the statistical information
  f_WS_XS_S <- XS_S_joint_density(x_s,a1,a2,a,b,VAI,VC) * 
    (1 / sqrt(VC)) * dnorm((x_s - w/I) / sqrt(VC)) 
  return(as.numeric(f_WS_XS_S))
}  

WS1_S_joint_density <- function(w,a1,a2,a,b,VAI,VC,I){
  WS_XS_S_joint_density_1arg <- function(x_s) {
    WS_XS_S_joint_density_1arg <- NULL
    for(j in 1:length(x_s)) {
      WS_XS_S_joint_density_1arg <- c(WS_XS_S_joint_density_1arg,
                                      WS_XS_S_joint_density(x_s[j], 
                                                            w = w, a1 = a1, a2 = a2, a = a, b = b, VAI = VAI, VC = VC,I=I))
    }
    return(WS_XS_S_joint_density_1arg)
  }
  
  WS_S_joint_density <- (1/I) * integrate(WS_XS_S_joint_density_1arg , #multiple the jacobian which is 1/I
                                  lower = -Inf, upper = Inf)$value
  
  return(WS_S_joint_density)
}

WS1_S_joint_density(w = 10, a1 = 1, a2 = 1, a = a, b = b, VAI = VAI_1, VC = VC_1, I=I_1) #check that a value is returned

#-------------------------------------------------------------------------------
# Check the density of W_S,1
#-------------------------------------------------------------------------------


treatments <- list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
density_range <- seq(from = -40, to = 62, by = 2) #range of ws1 from data

joint_dens <- matrix(NA, nrow = length(density_range), ncol = 4)

for (t in 1:4) {
  treat <- treatments[[t]]
  a1 <- treat[1]
  a2 <- treat[2]
  for(i in 1:length(density_range)) {
    joint_dens[i, t] <- WS1_S_joint_density(density_range[i],a1,a2,a,b,VAI=VAI_1,VC = VC_1,I=I_1)
  }
  print(t)
}

density_WS <- apply(joint_dens, 1, sum)

hist(sim.zs1$ws1,freq = FALSE)
lines(density_range, density_WS)


#-------------------------------------------------------------------------------
# Check the probability of W_S,1
#-------------------------------------------------------------------------------
WS1_S_joint_prob = function(a1,a2,a,b,VAI,VC,ub,lb,I){
  
  WS_S_joint_density_1arg <- function(w_s) {
    WS_S_joint_dens_1arg <- NULL
    for(j in 1:length(w_s)) {
      WS_S_joint_dens_1arg <- c(WS_S_joint_dens_1arg,
                                WS1_S_joint_density(w_s[j], 
                                                   a1 = a1, a2 = a2, a = a, b = b, VAI = VAI, VC = VC,I=I))
    }
    return(WS_S_joint_dens_1arg)
  }
  
  #Get WS1 prob by integrating over WS1 from lb to ub
  WS1_S_joint_probability <- integrate(WS_S_joint_density_1arg ,
                                       lower = lb, upper = ub)$value
  
  return(WS1_S_joint_probability)
}

WS1_S_joint_prob(a1=1,a2=1,a=a,b=b,VAI=VAI_1,VC=VC_1,ub=ub1,lb=lb1,I=I_1) #check that I get a value

#Get probability of WS1 being in (lb,ub) by summing over support of S
treatments <- list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))

WS1_prob = 0
for (t in 1:4) {
  treat <- treatments[[t]]
  WS1_prob = WS1_prob + WS1_S_joint_prob(a1= treat[1],a2= treat[2],
                                         a=a,b=b,VAI=VAI_1,VC=VC_1,ub=42.80393,lb=-Inf,I=I_1) 
  print(t)
}

#Verify results
WS1_prob
prop.table(table(sim.zs1$ws1>15 & sim.zs1$ws1<40))[2] #my prob matches the simulation

#-------------------------------------------------------------------------------
# Playing with the grid search to find the bounds
#-------------------------------------------------------------------------------

#Define the target probabilities based on alpha spending functions
alpha_star_upper <- 0.05 / 3
#alpha_star_lower <- 0.025 / 3

# Function to compute the upper prob for a given boundary
compute_upper_prob <- function(u1, a, b, VAI, VC, I) {
  upper_prob <- sum(sapply(treatments, function(treat) {
    WS1_S_joint_prob(a1 = treat[1], a2 = treat[2], a = a, b = b, VAI = VAI, VC = VC, ub = u1, lb = -Inf, I = I)
  }))
  #return(upper_prob)
  return(upper_prob - (1 - alpha_star_upper))
}
compute_upper_prob(u1=45,a=a,b=b,VAI=VAI_1,VC=VC_1,I=I_1)

compute_lower_prob <- function(l1, a, b, VAI, VC, I) {
  lower_prob <- sum(sapply(treatments, function(treat) {
    WS1_S_joint_prob(a1 = treat[1], a2 = treat[2], a = a, b = b, VAI = VAI, VC = VC, ub = Inf, lb = l1, I = I)
  }))
  return(lower_prob - alpha_star_lower)
}

# Define the range for uniroot to search for the boundaries
search_range <- c(39,45)

# Find the upper boundary u1
upper_result <- uniroot(compute_upper_prob, interval = search_range, a = a, b = b, VAI = VAI_1, VC = VC_1, I = I_1)
u1 <- upper_result$root

# Find the lower boundary l1
lower_result <- uniroot(compute_lower_prob, interval = search_range, a = a, b = b, VAI = VAI_1, VC = VC_1, I = I_1)
l1 <- lower_result$root

# Display the best l1 and u1
cat("Best l1:", l1, "\n")
cat("Best u1:", u1, "\n")


