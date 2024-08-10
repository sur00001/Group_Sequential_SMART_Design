# Functions to find joint probabilities and marginal densities of WS1, WS2 and WS3


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
# Functions to find marginal density of WS2
#-------------------------------------------------------------------------------
#Can just use WS1 and WS2 

marginal_density_WS2 = function(ws2) {
  integrate(function(ws1) {
    f_WS2_given_WS1(ws2, ws1) * f_WS1(ws1)
  }, lower = -Inf, upper = Inf, rel.tol = 10^-3)$value
}



# #If I wanted to find it with the joint of WS1, WS2 and WS3, but this is uncessary 
# WS2_joint_density = function(ws2, ws1) {
#   integrate(function(ws3) {
#     f_WS3_given_WS2(ws3, ws2)
#   }, lower = -Inf, upper = Inf)$value * f_WS2_given_WS1(ws2, ws1)
# }
# 
# marginal_density_WS2 = function(ws2) {
#   integrate(function(ws1) {
#     WS2_joint_density(ws2, ws1) * f_WS1(ws1)
#   }, lower = -Inf, upper = Inf)$value
# }