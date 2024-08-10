
library(stats)
library(mvtnorm)
library(ggplot2)
library(dplyr)


# Load your data
load("~/Adaptive SMART paper/Group_Sequential_SMART_Design/3interim.anal.res.dmv.rda")

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
lb = 10
#ub = Inf #get "max number of subdivisions when I do infinity & takes too long 
ub=1000


#-------------------------------------------------------------------------------
##Functions to find the marginal density of W_S1
#-------------------------------------------------------------------------------

XS_S_joint_density = function(x_s,a1,a2,a,b,VAI,VC){
  # Compute f(X_s,S)
  mean_vector = c(0, 0) #thetas are 0 regardless of the regime bc we simulate null
  if(a1 == 1){
    cov_matrix = matrix(c(VAI, b, b, VAI), nrow = 2)
  } else {
    cov_matrix = matrix(c(VAI, a, a, VAI), nrow = 2)
  }
  
  # Compute the probability
  bivnorm.prob = pmvnorm(lower = c(-Inf, -Inf), upper = c(x_s, x_s), 
                          mean = mean_vector, sigma = cov_matrix) 
  
  # Compute the mean and variance of the conditional distribution
  if(a1 == 1){
    mean_conditional = (a / VAI) * x_s
    variance_conditional = VAI - (a^2 / VAI)
  } else {
    mean_conditional = (b / VAI) * x_s
    variance_conditional = VAI - (b^2 / VAI)
  }
  
  # Compute the probability using the CDF of the normal distribution
  cond.prob = pnorm(x_s, mean = mean_conditional, sd = sqrt(variance_conditional))
  
  f_XS_S = dnorm(x_s, mean = 0, sd = sqrt(VAI)) * cond.prob * bivnorm.prob 
  return(as.numeric(f_XS_S))
}

WS_XS_S_joint_density = function(x_s,w,a1,a2,a,b,VAI,VC,I) {
  f_WS_XS_S = XS_S_joint_density(x_s,a1,a2,a,b,VAI,VC) * 
    (1 / sqrt(VC)) * dnorm((x_s - w/I) / sqrt(VC)) 
  return(as.numeric(f_WS_XS_S))
}  

WS1_S_joint_density = function(w,a1,a2,a,b,VAI,VC,I){
  WS_XS_S_joint_density_1arg = function(x_s) {
    sapply(x_s, function(x) WS_XS_S_joint_density(x, w, a1, a2, a, b, VAI, VC, I))
  }
  
  WS_S_joint_density = (1/I) * integrate(WS_XS_S_joint_density_1arg, 
                                          lower = -Inf, upper = Inf)$value
  return(WS_S_joint_density)
}

# Combined function to compute the marginal density of ws1 
#Using global variables so need to ensure I have the correct components saved in the environment
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

#Testing that integrating f(WSI) works and matches the WS1_prob I got in the other R file
#ws1.prob = integrate(f_WS1,lower = 10, upper = ub)$value 
#ws1.prob #it works

#-------------------------------------------------------------------------------
# Define the conditional density functions to do the recursive integration 
#-------------------------------------------------------------------------------
f_WS2_given_WS1 = function(ws2, ws1) {
  dnorm(ws2, mean = ws1, sd = sqrt(I2 - I1)) #theta_S is 0 so I took that term out
}

f_WS3_given_WS2 = function(ws3, ws2) {
  dnorm(ws3, mean = ws2, sd = sqrt(I3 - I2))
}

#-------------------------------------------------------------------------------
#Recursive integration to find joint probability of WS1 > 10 and WS2 > 10
#-------------------------------------------------------------------------------

#First integrate over Ws2 given a fixed value of Ws1 and then integrate over Ws1

# Calculate the joint probability using Vectorize
# Vectorize applies the function to each element of the vector individually. 
#The outer integrate() integrates the vectorized function over ws1 from 10 to 10,000, 
#summing up the contributions from all relevant ws1 values to compute the overall joint probability.
#Useful when integrating over a range of values, allows the inner function to handle multiple inputs simultaneously 
#so don't need extra function() wrappers. Also slightly faster

#---------------------------------
#Joint probability of all Ws < -10
#---------------------------------
start = proc.time()

joint_prob_WS1_WS2_lm10 =integrate(Vectorize(function(ws1) {
  f_WS1(ws1) * integrate(f_WS2_given_WS1, lower = -Inf, upper = -10, ws1 = ws1,rel.tol = 1e-10)$value
}), lower = -Inf, upper = -10,rel.tol = 1e-10)$value

print(joint_prob_WS1_WS2_lm10)
print(proc.time() - start) #.008

#See if it matches the simulation 
mean(sim.ws1 < -10 & sim.ws2 < -10) #.009, close enough

#---------------------------------
#Joint probability above 10
#---------------------------------
start = proc.time()

joint_prob_WS1_WS2_g10 =integrate(Vectorize(function(ws1) {
  f_WS1(ws1) * integrate(f_WS2_given_WS1, lower = 10, upper = Inf, ws1 = ws1,rel.tol = 1e-10)$value
}), lower = 10, upper = Inf,rel.tol = 1e-10)$value

print(joint_prob_WS1_WS2_g10)
print(proc.time() - start) #.135

#See if it matches the simulation 
mean(sim.ws1>10 & sim.ws2>10) #.143, close enough

#-------------------------------------------------------------------------------
#Recursive integration to find joint probability of WS1 > 10, WS2 > 10, WS3 >10
#-------------------------------------------------------------------------------
start = proc.time()
#---------------------------------
#Less than -10
#---------------------------------
joint_prob_WS1_WS2_WS3_lm10 = integrate(Vectorize(function(ws1) {
  f_WS1(ws1) * integrate(Vectorize(function(ws2) {
    f_WS2_given_WS1(ws2, ws1) * integrate(f_WS3_given_WS2, lower = -Inf, upper = -10, ws2 = ws2, rel.tol = 1e-10)$value
  }), lower = -Inf, upper = -10, rel.tol = 1e-10)$value
}), lower = -Inf, upper = -10, rel.tol = 1e-10)$value

# Print the joint probability
print(joint_prob_WS1_WS2_WS3_lm10) #0.006028711

mean(sim.ws1< -10 & sim.ws2 < -10 &sim.ws3< -10) #.006, matches perfectly

#---------------------------------
#Greater than 10
#---------------------------------
joint_prob_WS1_WS2_WS3_g10 = integrate(Vectorize(function(ws1) {
  f_WS1(ws1) * integrate(Vectorize(function(ws2) {
    f_WS2_given_WS1(ws2, ws1) * integrate(f_WS3_given_WS2, lower = 10, upper = Inf, ws2 = ws2, rel.tol = 1e-10)$value
  }), lower = 10, upper = Inf, rel.tol = 1e-10)$value
}), lower = 10, upper = Inf, rel.tol = 1e-10)$value

# Print the joint probability
print(joint_prob_WS1_WS2_WS3_g10) #.106
print(proc.time() - start) 

mean(sim.ws1>10 & sim.ws2>10 &sim.ws3>10) #.121

#---------------------------------
#Joint Prob that WS1 in [-15,12] & WS2 in [-15,25] & WS3 > 30
#---------------------------------
joint_prob_WS1_WS2_WS3 = integrate(Vectorize(function(ws1) {
  f_WS1(ws1) * integrate(Vectorize(function(ws2) {
    f_WS2_given_WS1(ws2, ws1) * integrate(f_WS3_given_WS2, lower = 30, upper = Inf, ws2 = ws2, rel.tol = 1e-10)$value
  }), lower = -15, upper = 25, rel.tol = 1e-10)$value
}), lower = -15, upper = 12, rel.tol = 1e-10)$value 

# Print the joint probability
print(joint_prob_WS1_WS2_WS3) #.008

mean(sim.ws1> -15 &sim.ws1 <12 & sim.ws2> -15 & sim.ws2<25 & sim.ws3>30) #.005



#-------------------------------------------------------------------------------
#Check if WS1, WS2 and WS3 marginal densities match the simulated distributions
#-------------------------------------------------------------------------------
 #Only checked WS2 and WS1. WS3 is more complicated

#--------------------------
#WS1
#--------------------------

# Sequence of values for ws1 
ws1_values = seq(min(sim.ws1), max(sim.ws1), length.out = 200)

# Marginal density of ws1 for this range
ws1_density = sapply(ws1_values, f_WS1)

#Plot density over histogram
hist(sim.ws1, breaks = 30, probability = TRUE, main = "Marginal Density of ws1", xlab = "ws1", xlim = c(min(sim.ws1), max(sim.ws1)), col = "lightgray")
lines(ws1_values, ws1_density, col = "blue", lwd = 2)

#------------------------------
# Marginal density of WS2
#------------------------------

#David's code
marginal_density_WS2 = function(ws2) {
	integrate(function(ws1) {
		f_WS2_given_WS1(ws2, ws1) * f_WS1(ws1)
	}, lower = -Inf, upper = Inf, rel.tol = 10^-3)$value
}

hist(sim.ws2, breaks = 30, probability = TRUE, main = "Marginal Density of WS2", xlab = "WS2", xlim = c(min(sim.ws2), max(sim.ws2)), col = "lightgray")

ws2_values = seq(min(sim.ws2), max(sim.ws2), length.out = 100)
ws2_values = seq(-20, 30, length.out = 10)

# Compute the marginal density of WS2 for this range
ws2_density = NULL
start = proc.time()
for(i in 1:length(ws2_values)) {
	start = proc.time()
	ws2_density = c(ws2_density, marginal_density_WS2(ws2_values[i]))
	print(i)
	print(proc.time() - start)
}


lines(ws2_values, ws2_density, col = "blue", lwd = 2)


## David: not the correct form for joint density of WS2 and WS1 but you end up in the right place
## for the marginal distribution of W2
## can be quicker (only 1 integral) in the above code


#My code
WS2_joint_density = function(ws2, ws1) {
  integrate(function(ws3) {
    f_WS3_given_WS2(ws3, ws2)
  }, lower = -Inf, upper = Inf)$value * f_WS2_given_WS1(ws2, ws1)
}

marginal_density_WS2 = function(ws2) {
  integrate(function(ws1) {
    WS2_joint_density(ws2, ws1) * f_WS1(ws1)
  }, lower = -Inf, upper = Inf)$value
}

#Range of values
ws2_values = seq(min(sim.ws2), max(sim.ws2), length.out = 10)

# Compute the marginal density of WS2 for this range
start2 = proc.time()
ws2_density = sapply(ws2_values, marginal_density_WS2)
print(proc.time() - start2)

#Plot density 
hist(sim.ws2, breaks = 30, probability = TRUE, main = "Marginal Density of WS2", xlab = "WS2", xlim = c(min(sim.ws2), max(sim.ws2)), col = "lightgray")
lines(ws2_values, ws2_density, col = "blue", lwd = 2)

#------------------------------------
#Marginal probability of WS2 > 10 
#------------------------------------

#-------------------------------
## If I only use W1 and W2 ##
#-------------------------------
#Vectorize functions
f_WS2_given_WS1_vect = Vectorize(f_WS2_given_WS1)
f_WS1_vect = Vectorize(f_WS1)

marginal_density_WS2 = Vectorize(function(ws2) {
  integrate(function(ws1) {
    f_WS2_given_WS1_vect(ws2, ws1) * f_WS1_vect(ws1)
  }, lower = -100, upper = 100, rel.tol = 10^-3)$value
})

#Integrate marginal density WS2 to get the probability (integrating out WS1 from joint of WS1 and WS2)
start = proc.time()
pr.ws2.10.30 = integrate(marginal_density_WS2, lower = 10, upper = 30)$value ; print(pr.ws2.10.30)
print(proc.time() - start) #3.85 min, .28

mean(sim.ws2 >10 & sim.ws2<30) #close to simulation probability .304

#-------------------------------
## If I use joint of W1, W2 and W3 to get marginal probabilities of WS2 ##
#-------------------------------

#Vectorize 
f_WS3_given_WS2_vect = Vectorize(f_WS3_given_WS2)

WS2_joint_density_vect = Vectorize(function(ws2, ws1) {
  integrate(function(ws3) {
    f_WS3_given_WS2_vect(ws3, ws2)
  }, lower = -60, upper = 60)$value * f_WS2_given_WS1_vect(ws2, ws1)
})

marginal_density_WS2_vect = Vectorize(function(ws2) {
  integrate(function(ws1) {
    WS2_joint_density_vect(ws2, ws1) * f_WS1_vect(ws1)
  }, lower = -60, upper = 60)$value
})

#Integrate the marginal density WS2 (derived from integrating out WS1 and WS3)
start = proc.time() 
pr.ws2.fromWS3.10.15 = integrate(marginal_density_WS2_vect, lower = 10, upper = 15)$value ; print(pr.ws2.fromWS3.10.15)
print(proc.time() - start) #3.6 minutes, .145

mean(sim.ws2 >10 & sim.ws2<15) #similar to simulation probability #.155


#--------------------------
#WS3 - takes way too long
#--------------------------
#THIS IS INCORRECT becauses integrating out W1 does not simplify to f(W2)f(W3|W2) because W3 is technically conditioned on W1


#Joint density f(W_S3, W_S2) - THIS IS NOT TRUE 
WS3_joint_density = function(ws3, ws2) {
  integrate(function(ws1) {
    f_WS2_given_WS1(ws2, ws1) * f_WS1(ws1)
  }, lower = -Inf, upper = Inf, subdivisions = 100)$value * f_WS3_given_WS2(ws3, ws2)
}

#Marginal density f(W_S3)
marginal_density_WS3 = function(ws3) {
  integrate(function(ws2) {
    WS3_joint_density(ws3, ws2)
  }, lower = -Inf, upper = Inf, subdivisions = 100)$value
}

#Density for a range of ws3 values in data 
ws3_values = seq(min(sim.ws3), max(sim.ws3), length.out = 100)
ws3_density = sapply(ws3_values, marginal_density_WS3)

hist(sim.ws3, breaks = 30, probability = TRUE, main = "Marginal Density of WS3", xlab = "WS3", xlim = c(min(sim.ws3), max(sim.ws3)), col = "lightgray")
lines(ws3_values, ws3_density, col = "blue", lwd = 2)

#-------------------------------------------------------------------------------
# Define the marginal probability function for W_S3 > c

#This is incorrect. Can't just integrate out ws1 to get the joint of ws2,ws3
#because ws3 technically is conditioned on ws1 (Markov property just simplifies the conditional distribution)
#-------------------------------------------------------------------------------
P_WS3_greater_than_c = function(c) { #need to integrate out ws2 and ws1
  integrate(function(ws2) {
    integrate(function(ws1) {
      f_WS1(ws1) * f_WS2_given_WS1(ws2, ws1)
    }, lower = -Inf, upper = Inf)$value * (1 - pnorm(c, mean = ws2 + theta_S * (I3 - I2), sd = sqrt(I3 - I2)))
  }, lower = -Inf, upper = Inf)$value
}

#P_WS3_greater_than_c(10) #this takes FOREVER...


#-------------------------------------------------------------------------------
#ROUGH WORK
# 

#Rough work when trying to find marginal probabilities of WS2
# #Testing product of vector ws1 and fixed ws2 
# f_WS1(c(10,15,20)) *f_WS2_given_WS1(15, c(10,15,20))
# 
# f_WS1_vect(c(10,15,20)) *f_WS2_given_WS1_vect(c(4,5,6), c(10,15,20))
# 
# #need to debug why I get the "evaluation of function gave result of wrong length" error. Nede to check all functions return a scalar. 
# #Solution: had to Vectorize all the functions 
# # All return scalars...
# f_WS3_given_WS2(1,10) #returns scalar
# f_WS2_given_WS1(10,15) #returns scalar
# 
# WS2_joint_density(20,23) #returns scalar
# marginal_density_WS2(10) #returns scalar
# 
# #Testing f_WS1 and f_WS2_given_WS1 with vector inputs '
# ws1_test <- c(1, 2, 3)
# ws2_test <- 10
# 
# # Check the output of these functions
# print(f_WS1(ws1_test))                
# print(f_WS2_given_WS1(ws2_test, ws1_test))



# 
# #-------------------------------------------------------------------------------
# #Recursive integration (IN PARALLEL)
# 
# # Using Parallel versions of recursive integration functions
# #-------------------------------------------------------------------------------
# library(parallel)
# 
# start = proc.time()
# 
# #GETTING VERY DIFFERENT ANSWER - .02 WHEN THE SERIES WAS .60
# # Set up parallel backend
# n_cores = detectCores() - 1
# cl = makeCluster(n_cores)
# 
# 
# integrate_tight.tol_WS1_parallel = function() {
#   ws1_values = seq(lb, ub, length.out = 100)  # Divide the interval into 100 points
#   result = parSapply(cl, ws1_values, function(ws1) {
#     f_WS1(ws1) * integrate_WS2(ws1)
#   })
#   sum(result)
# }
# 
# # Export necessary functions and variables to cluster workers
# clusterExport(cl, c("f_WS1", "f_WS2_given_WS1", "f_WS3_given_WS2", "integrate_WS3", "integrate_WS2",
#                     "WS1_S_joint_density", "WS_XS_S_joint_density", "XS_S_joint_density",
#                     "lb", "ub", "I1", "I2", "I3", "VAI_1", "VC_1", "a", "b", "pmvnorm"))
# 
# #"parSapply", "integrate_tight.tol_WS1_parallel"
# 
# # Compute the probability using parallelized functions
# probability = integrate_tight.tol_WS1_parallel()
# print(probability)
# 
# # Stop the parallel backend
# stopCluster(cl)
# print(proc.time() - start)

# # 
# # #Code with cl removed from parsapply() so cluster worker does not think "cl" is a defined object
# # Parallelized integration of WS3
# integrate_WS3_parallel = function(ws2_values) {
#   parSapply(ws2_values, function(ws2) {  # = `cl` removed
#     integrate(f_WS3_given_WS2, lower = lb, upper = ub, ws2_values = ws2, rel.tol = 1e-10, subdivisions = 100)$value
#   })
# }
# 
# # Parallelized integration of WS2
# integrate_WS2_parallel = function(ws1_values) {
#   parSapply(ws1_values, function(ws1) {  # = `cl` removed
#     integrate(function(ws2) {
#       f_WS2_given_WS1(ws2, ws1) * integrate_WS3_parallel(ws2)
#     }, lower = lb, upper = ub, rel.tol = 1e-10, subdivisions = 100)$value
#   })
# }
# 
# # Parallelized integration of WS1
# integrate_tight.tol_WS1_parallel = function() {
#   result = parSapply(seq(lb, ub, length.out = 100), function(ws1) {  # = `cl` removed
#     f_WS1(ws1) * integrate_WS2_parallel(ws1)
#   })
#   sum(result)
# }
# 
