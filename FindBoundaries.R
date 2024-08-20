
#Script to find boundaries 

library(stats)
library(mvtnorm)
library(ggplot2)

#Read in functions
source("~/Adaptive SMART paper/Group_Sequential_SMART_Design/Functions4JointProbsMargDensitiesofWS.R")

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

#-------------------------------------------------------------------------------
# Grid search to find boundaries for 3 interim analyses 
#-------------------------------------------------------------------------------

# Function to compute joint probability for given boundaries for all interim analyses 
compute_joint_prob = function(u1, l1, u2, l2, u3=Inf, l3=-Inf, num.IntAnal) {

  if (num.IntAnal == 2) {
    
    joint_prob =integrate(Vectorize(function(ws1) {
      f_WS1(ws1) * integrate(f_WS2_given_WS1, lower = l2, upper = u2, ws1 = ws1,rel.tol = 1e-10)$value
    }), lower = l1, upper = u1,rel.tol = 1e-10)$value
    
  } else {
    
    joint_prob = integrate(Vectorize(function(ws1) {
      f_WS1(ws1) * integrate(Vectorize(function(ws2) {
        f_WS2_given_WS1(ws2, ws1) * integrate(f_WS3_given_WS2, lower = l3, upper = u3, ws2 = ws2, rel.tol = 1e-10)$value
      }), lower = l2, upper = u2, rel.tol = 1e-10)$value
    }), lower = l1, upper = u1, rel.tol = 1e-10)$value
  }
  return(joint_prob)
}

# How much alpha to spend for each interim analysis 
alpha_star_upper_1 = 0.05 / 3
alpha_star_upper_2 = 0.05 / 3
alpha_star_upper_3 = 0.05 / 3


#---------------------------------------------------------------
# Interim analysis 1
#---------------------------------------------------------------
# Function to compute the upper prob for interim analysis 1
find_boundary_u1 = function(u1) { #Using global variables
  
  upper_prob = integrate(f_WS1,lower = -Inf, upper=u1)$value
  #return(upper_prob)
  return(upper_prob - (1 - alpha_star_upper_1))
}

# Define the range for uniroot to search for the u1 boundaries
u1_search_range <- c(15,20)

# Grid search for u1
start = proc.time()
fixed_u1 <- uniroot(find_boundary_u1, interval = u1_search_range)$root #4 min
print(proc.time() - start)

#See if Pr(WS1 > fixed_u1) is .05/3 from simulation
mean(sim.ws1>fixed_u1)
#Probability is .01 (should be .01666) #Should try with 5000-10000 simulation reps instead of 1000

#---------------------------------------------------------------
# Interim analysis 2
#---------------------------------------------------------------

# Function to find u2 given u1 - #just use the joint distribution of Ws1, Ws2
find_boundary_u2 = function(u2) {
  upper_prob = compute_joint_prob(u2 = u2, l2 = -Inf, u1 = fixed_u1, l1 = -Inf, num.IntAnal = 2)
  return(upper_prob - (1 - alpha_star_upper_2))
}
#Need to double check this this weekend and investigate different ranges to search

# Define the range for uniroot to search for the u2 boundaries, only getting negative values
u2_search_range <- c(35,55)

# Grid search for u2, given fixed_u1
start = proc.time()
fixed_u2 = uniroot(find_boundary_u2, interval = u2_search_range)$root
print(proc.time() - start)


#---------------------------------------------------------------
# Interim analysis 3
#---------------------------------------------------------------
# Function to find u3 given u1 and u2 
find_boundary_u3 = function(u3) {
  upper_prob = compute_joint_prob(u1 = fixed_u1, l1 = -Inf, u2 = fixed_u2, l2 = -Inf, u3 = u3, l3 = -Inf)
  return(upper_prob - (1 - alpha_star_upper_3))
}

# Define the range for uniroot to search for the u2 boundaries
u3_search_range <- c(35,45)

# Grid search for u3, given fixed_u1 and fixed_u2
fixed_u3 = uniroot(find_boundary_u3, interval = u3_search_range)$root

cat("Boundary u1:", fixed_u1, "\n")
cat("Boundary u2:", fixed_u2, "\n")
cat("Boundary u3:", fixed_u3, "\n")















# 
# 
# compute_lower_prob = function(l1, a, b, VAI, VC, I) {
#   lower_prob = sum(sapply(treatments, function(treat) {
#     WS1_S_joint_prob(a1 = treat[1], a2 = treat[2], a = a, b = b, VAI = VAI, VC = VC, ub = Inf, lb = l1, I = I)
#   }))
#   return(lower_prob - alpha_star_lower)
# }
# 
# # Define the range for uniroot to search for the boundaries
# upper_search_range = c(39,45) #search range for upper boundaries
# lower_search_range = c(39,45) #search range for upper boundaries
# 
# # Find the upper boundary u1
# upper_result = uniroot(compute_upper_prob, interval = upper_search_range, a = a, b = b, VAI = VAI_1, VC = VC_1, I = I_1)
# u1 = upper_result$root
# 
# # Find the lower boundary l1
# lower_result = uniroot(compute_lower_prob, interval = search_range, a = a, b = b, VAI = VAI_1, VC = VC_1, I = I_1)
# l1 = lower_result$root
# 
# # Display the best l1 and u1
# cat("Best l1:", l1, "\n")
# cat("Best u1:", u1, "\n")
# 
# 
# 
# # #-------------------------------------------------------------------------------
# # # For interim analysis 1
# # #-------------------------------------------------------------------------------
# # 
# # #Define the target probabilities based on alpha spending functions
# # alpha_star_upper = 0.05 / 3
# # #alpha_star_lower = 0.025 / 3
# # 
# # # Function to compute the upper prob for a given boundary
# # compute_upper_prob1 = function(u1, a, b, VAI, VC, I) {
# #   upper_prob = sum(sapply(treatments, function(treat) {
# #     WS1_S_joint_prob(a1 = treat[1], a2 = treat[2], a = a, b = b, VAI = VAI, VC = VC, lb=-Inf, ub=ub, I = I)
# #   }))
# #   #return(upper_prob)
# #   return(upper_prob - (1 - alpha_star_upper))
# # }
# # #compute_upper_prob(u1=20,a=a,b=b,VAI=VAI_1,VC=VC_1,I=I1)
# # 
# # 


                            