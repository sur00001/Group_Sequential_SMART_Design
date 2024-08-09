
#Script to find boundaries 

#-------------------------------------------------------------------------------
# For interim analysis 1
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


                            