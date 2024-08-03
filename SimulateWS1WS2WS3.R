#Simulate WS1, WS2, WS3 

#Load libraries
library(dplyr)
library(geepack)
library(gsDesign)
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(geepack))


#-------------------------------------------------------------------------------
# Set simulation parameters (for debugging)
#-------------------------------------------------------------------------------
# set.seed(1996)
nsim = 1000 #Number of simulations 
N.all=3000 #N is total people in the recruitment group
# int.tprop = c(.8,.1,.1) #proportion of information for each interim analysis (later need to see how this affects the alpha boundaries)
nr.prob = .3 #non-response rate in overal API group (can be different for A1 groups if A1 group outcomes are different)
control.mean = -.258
control.sd = 1.25
# delta_E = .19 #effect size (difference between E[yt|early] - E[yt|control])
delta_E = 0
delta_L = 0
delta_C = 0
delta_Em = 0
alpha = 0 #intercept when simulated y values (effect of time)
# intanal.type = "OBF"
# cont1best = 0 #indicator if you want to continue 1 best trt for interim analysis


#-------------------------------------------------------------------------------
# Function to calculate the IPW mean
#-------------------------------------------------------------------------------
ipw.est = function(df,A1.val, A2.val){
  #Randomization probabilities are .5 
  
  numerator = (df$A1 == A1.val & ((df$NR == 1 & df$A2 == A2.val) | df$NR == 0)) * df$y1
  denominator = df$A1.prob * (df$A2.prob * as.integer(df$NR==1) + as.integer(df$NR==0))
  
  #Calculate the mean
  mu_a1_a2 = mean((numerator / denominator))
  return(mu_a1_a2)
}

#-------------------------------------------------------------------------------
# Function to simulate SMART data for 4 DTRs + control for a recruitment group
#-------------------------------------------------------------------------------

sim.SMART.allDTR = function(N.all, nr.prob, control.mean, control.sd, delta_E,delta_L,delta_Em,delta_C){
  
  N.int1 = N.all/3 #how many people in each interim analysis (we will do 3 interim analyses)
  
  #Initialize data frame
  n.con = as.integer(N.int1/3)
  control = late = early = data.frame(id=1:n.con, y0=rep(0,n.con))
  late$id = control$id + n.con
  early$id = late$id + n.con
  
  control$API = "Control"; early$API = "API"; late$API = "API"
  control$A1 = "Control"; early$A1 = 1; late$A1 = -1
  
  #Only coding API,A1,A2, responder status, Y1, Y0
  
  #-------------------------
  #Simulate baseline data y0
  #-------------------------
  control$y0 = rnorm(n.con, control.mean, control.sd)
  early$y0 = rnorm(n.con,control.mean,control.sd)
  late$y0 = rnorm(n.con,control.mean,control.sd)
  api.dat = rbind(late,early)
  n.api = dim(api.dat)[1]
  
  #-------------------------
  #Simulate intermediate outcomes y.5. All deltas = 0 under null
  #E[y.5] = y0 + alpha + I(early)*delta_E*control.sd + I(late)*delta_L*control.sd + random variation (control.sd)
  #-------------------------
  #Create alternative scenario (early group has higher outcomes) and link y0 and y1
  control$y.5 = rnorm(n.con, control$y0 + alpha, control.sd)
  api.dat$y.5 = rnorm(n.api,api.dat$y0 + alpha + (api.dat$A1==1)*delta_E*control.sd +
                        (api.dat$A1==-1)*delta_L*control.sd,control.sd)
  
  #-------------------------
  # Identify non-responders
  # Code NR
  #-------------------------

    nr.threshold = quantile(api.dat$y.5,nr.prob) #identify the NR threshold based on the specified nr.prob
    api.dat$NR = ifelse(api.dat$y.5 < nr.threshold,1,0)
    control$NR = "Control"

  
  #Check non-response rates 
  table(api.dat$NR); table(api.dat$NR, api.dat$A1)
  prop.table(table(api.dat$NR, api.dat$A1))
  
  #-------------------------
  # Randomize 2nd stage trt A2 for non-responders in each A1 group
  # 1:1 randomization
  #-------------------------
  early.coach = api.dat[api.dat$NR==1 & api.dat$A1==1,] %>% sample_frac(0.5)
  early.email = anti_join(api.dat[api.dat$NR==1 & api.dat$A1==1,],early.coach)
  
  late.coach = api.dat[api.dat$NR==1 & api.dat$A1==-1,] %>% sample_frac(0.5)
  late.email = anti_join(api.dat[api.dat$NR==1 & api.dat$A1==-1,],late.coach)
  
  #Code A2
  api.dat$A2 = 0
  api.dat$A2[api.dat$id %in% c(early.coach$id,late.coach$id)] = -1
  api.dat$A2[api.dat$id %in% c(early.email$id,late.email$id)] = 1
  control$A2 = "Control"
  
  table(api.dat$A2); table(api.dat$A2,api.dat$A1)
  
  #-------------------------
  #Simulate final outcomes y1. Under null, all deltas = 0
  #-------------------------
  control$y1 = rnorm(n.con, mean = control$y.5 + alpha, sd = control.sd)
  api.dat$y1 = rnorm(n.api,
                     mean = api.dat$y.5 + alpha + (api.dat$A1==1)*delta_E*control.sd+
                       (api.dat$A1==-1)*delta_L*control.sd+
                       (api.dat$A2==1)*delta_Em*control.sd+
                       (api.dat$A2==-1)*delta_C*control.sd,
                     sd = control.sd)
 
  #-------------------------
  #Mu hats for interim analysis 1
  #------------------------- 
  api.dat$A1.prob = .5; api.dat$A2.prob = .5 #randomization probabilities (when the best trt is chosen these will be 1 for later interim analyses)
  muhat.0 = mean(control$y1)
  muhat.m11 = ipw.est(df=api.dat,A1.val=-1,A2.val = 1)
  muhat.m1m1 =  ipw.est(df=api.dat,A1.val=-1,A2.val = -1)
  muhat.11 =ipw.est(df=api.dat,A1.val=1,A2.val = 1)
  muhat.1m1 = ipw.est(df=api.dat,A1.val=1,A2.val = -1)
  mu0 = control.mean
  
  # Create the list with all variables
  int.analysis1 = data.frame(
    int.anal = 1,
    muhat.0 = muhat.0,
    muhat.m11 = muhat.m11,
    muhat.m1m1 = muhat.m1m1,
    muhat.11 = muhat.11,
    muhat.1m1 = muhat.1m1,
    mu0 = mu0,
    z11 = muhat.11 - muhat.0,
    z1m1 = muhat.1m1 - muhat.0,
    zm11 = muhat.m11 - muhat.0,
    zm1m1 = muhat.m1m1 - muhat.0,
    x11 = muhat.11 - mu0,
    xm11 = muhat.m11 - mu0,
    x1m1 = muhat.1m1 - mu0,
    xm1m1 = muhat.m1m1 - mu0,
    x0 = muhat.0 - mu0,
    xs = max(c(muhat.11 - mu0, muhat.m11 - mu0, muhat.1m1 - mu0, muhat.m1m1 - mu0)),
    zs = max(c(muhat.11 - muhat.0, muhat.1m1 - muhat.0, muhat.m11 - muhat.0, muhat.m1m1 - muhat.0)),
    S = c("11", "1-1", "-11", "-1-1")[which.max(c(muhat.11 - muhat.0, muhat.1m1 - muhat.0, muhat.m11 - muhat.0, muhat.m1m1 - muhat.0))]
  )
  # Add a1 and a2 columns based on S
  int.analysis1$Sa1 = ifelse(int.analysis1$S %in% c("11", "1-1"), 1, -1)
  int.analysis1$Sa2 = ifelse(int.analysis1$S %in% c("11", "-11"), 1, -1)
  
  #-------------------------------------------------------------------------------
  # Interim analysis 2 with best treatment (not looking at efficacy or futility here, just collecting the results to get WS2, WS3)
  #-------------------------------------------------------------------------------

  N.int2 = N.all/3
  #Initialize data frame
  n.con2 = as.integer(N.int2/2) #half in control and half gets the best AI. N is number of people in each interim analysis
  control2 = api2 = data.frame(y0=rep(0,n.con2))
  
  control2$API = "Control"; api2$API = "API"
  control2$A1 = "Control"; api2$A1 = int.analysis1$Sa1
  
  #Only coding API,A1,A2, responder status, Y1, Y0
  
  #-------------------------
  #Simulate baseline data y0
  #-------------------------
  control2$y0 = rnorm(n.con2, control.mean, control.sd)
  api2$y0 = rnorm(n.con2,control.mean,control.sd)
  n.api2 =n.con2
  
  #-------------------------
  #Simulate intermediate outcomes y.5. All deltas = 0 under null
  #E[y.5] = y0 + alpha + I(early)*delta_E*control.sd + I(late)*delta_L*control.sd + random variation (control.sd)
  #-------------------------
  #Create alternative scenario (early group has higher outcomes) and link y0 and y1
  control2$y.5 = rnorm(n.con2, control2$y0 + alpha, control.sd)
  api2$y.5 = rnorm(n.api2,api2$y0 + alpha + (api2$A1==1)*delta_E*control.sd +
                        (api2$A1==-1)*delta_L*control.sd,control.sd)
  
  #-------------------------
  # Identify non-responders
  # Code NR
  #-------------------------
  
  nr.threshold = quantile(api2$y.5,nr.prob) #identify the NR threshold based on the specified nr.prob
  api2$NR = ifelse(api2$y.5 < nr.threshold,1,0)
  control2$NR = "Control"
  
  #-------------------------
  # Assign A2 from best regime (S) for non-responders 
  #-------------------------
  #Code A2
  api2$A2 = 0
  api2$A2[api2$NR==1]= int.analysis1$Sa2
  control2$A2 = "Control"
  
 table(api2$A2,api2$A1) #check
  
  #-------------------------
  #Simulate final outcomes y1. Under null, all deltas = 0
  #-------------------------
  control2$y1 = rnorm(n.con2, mean = control2$y.5 + alpha, sd = control.sd)
  api2$y1 = rnorm(n.api2,
                     mean = api2$y.5 + alpha + (api2$A1==1)*delta_E*control.sd+
                       (api2$A1==-1)*delta_L*control.sd+
                       (api2$A2==1)*delta_Em*control.sd+
                       (api2$A2==-1)*delta_C*control.sd,
                     sd = control.sd)
  
  
  #-------------------------
  #Mu hats for interim analysis 2
  #------------------------- 
  api2$A1.prob = 1; api2$A2.prob = 1 #when the best trt is chosen, the weight are 1 for participants in recruitment periods 2 and 3
  control$id = NULL; api.dat$id = NULL
  
  #Combine data from interim analyses 1 and 2 
  api.dat2 = rbind(api.dat[api.dat$A1==int.analysis1$Sa1 & (api.dat$A2==int.analysis1$Sa2|api.dat$A2==0),],api2)
  control.dat2 = rbind(control, control2) #combine control data from 2 interim analyses
 
  muhat2.0 = mean(control.dat2$y1)
  muhat.s = ipw.est(df=api.dat2,A1.val=int.analysis1$Sa1,A2.val = int.analysis1$Sa2)
  mu0 = control.mean
  
  # Create the df with all results
  int.analysis2 = data.frame(
    int.anal = 2,
    muhat.0 = muhat2.0,
    muhat.s = muhat.s,
    mu0 = mu0,
    zs = muhat.s - muhat2.0,
    x0 = muhat2.0 - mu0,
    xs = muhat.s-mu0,
    S = int.analysis1$S,
    Sa1 = int.analysis1$Sa1,
    Sa2 = int.analysis1$Sa2
  )
  
  

  #-------------------------------------------------------------------------------
  # Interim analysis 3 with best treatment (not looking at efficacy or futility here, just collecting the results to get WS2, WS3)
  #-------------------------------------------------------------------------------
  
  N.int3 = N.all/3
  #Initialize data frame
  n.con3 = as.integer(N.int3/2) #half in control and half gets the best AI. N is number of people in each interim analysis
  control3 = api3 = data.frame(y0=rep(0,n.con3))

  control3$API = "Control"; api3$API = "API"
  control3$A1 = "Control"; api3$A1 = int.analysis1$Sa1
  
  #-------------------------
  #Simulate baseline data y0
  #-------------------------
  control3$y0 = rnorm(n.con3, control.mean, control.sd)
  api3$y0 = rnorm(n.con3,control.mean,control.sd)
  n.api3 =n.con3
  
  #-------------------------
  #Simulate intermediate outcomes y.5. All deltas = 0 under null
  #E[y.5] = y0 + alpha + I(early)*delta_E*control.sd + I(late)*delta_L*control.sd + random variation (control.sd)
  #-------------------------
  #Create alternative scenario (early group has higher outcomes) and link y0 and y1
  control3$y.5 = rnorm(n.con3, control3$y0 + alpha, control.sd)
  api3$y.5 = rnorm(n.api3,api3$y0 + alpha + (api3$A1==1)*delta_E*control.sd +
                     (api3$A1==-1)*delta_L*control.sd,control.sd)
  
  #-------------------------
  # Identify non-responders
  # Code NR
  #-------------------------
  
  nr.threshold = quantile(api3$y.5,nr.prob) #identify the NR threshold based on the specified nr.prob
  api3$NR = ifelse(api3$y.5 < nr.threshold,1,0)
  control3$NR = "Control"
  
  #-------------------------
  # Assign A2 from best regime (S) for non-responders 
  #-------------------------
  #Code A2
  api3$A2 = 0
  api3$A2[api3$NR==1]= int.analysis1$Sa2
  control3$A2 = "Control"
  
  table(api3$A2,api3$A1) #check
  
  #-------------------------
  #Simulate final outcomes y1. Under null, all deltas = 0
  #-------------------------
  control3$y1 = rnorm(n.con3, mean = control3$y.5 + alpha, sd = control.sd)
  api3$y1 = rnorm(n.api3,
                  mean = api3$y.5 + alpha + (api3$A1==1)*delta_E*control.sd+
                    (api3$A1==-1)*delta_L*control.sd+
                    (api3$A2==1)*delta_Em*control.sd+
                    (api3$A2==-1)*delta_C*control.sd,
                  sd = control.sd)
  
  
  #-------------------------
  #Mu hats for interim analysis 3 (final analysis)
  #------------------------- 
  api3$A1.prob = 1; api3$A2.prob = 1 #when the best trt is chosen, the weight are 1 for participants in recruitment periods 2 and 3
  control.dat3 = rbind(control, control2, control3) #combine all control data 
  muhat3.0 = mean(control.dat3$y1)
  
  api.dat3 = rbind(api.dat2,api3) #all the data
  muhat.s = ipw.est(df=api.dat3,A1.val=int.analysis1$Sa1,A2.val = int.analysis1$Sa2)
  mu0 = control.mean
  
  # Create the df with all results
  int.analysis3 = data.frame(
    int.anal = 3,
    muhat.0 = muhat3.0,
    muhat.s = muhat.s,
    mu0 = mu0,
    zs = muhat.s - muhat3.0,
    x0 = muhat3.0 - mu0,
    xs = muhat.s-mu0,
    S = int.analysis1$S,
    Sa1 = int.analysis1$Sa1,
    Sa2 = int.analysis2$Sa2
  )
  
  return(list(int.analysis1,int.analysis2,int.analysis3))
}

#-------------------------------------------------------------------------------
# Run simulation 
#-------------------------------------------------------------------------------
int1.res = int2.res = int3.res = data.frame()
for (i in 1:nsim){
  print(i)
  sim.dat = sim.SMART.allDTR(N.all=3000,nr.prob=nr.prob,control.mean= control.mean,control.sd=control.sd,
                             delta_E=delta_E,delta_L=delta_L,delta_Em=delta_Em,delta_C=delta_C)
  int1.res = rbind(int1.res,sim.dat[[1]])
  int2.res = rbind(int2.res,sim.dat[[2]])
  int3.res = rbind(int3.res,sim.dat[[3]])
}
int1.res$sim.no = int2.res$sim.no = int3.res$sim.no = 1:nsim

setwd("~/Adaptive SMART paper/Group-Sequential-SMART-Research")
save(int1.res,int2.res,int3.res,file="3interim.anal.res.rda")


#-------------------------------------------------------------------------------
# Is the true value theta_max = 0 under the null (all thetas = 0)? 

#this is iid though but it shouldn't be too different 
#-------------------------------------------------------------------------------

set.seed(123)

# Number of simulations - Monte Carlo to approximate the true value
n_sims = 1000000

# Simulate four standard normal variables and find their maximum
theta_hat_max = replicate(n_sims, max(rnorm(4, mean = 0, sd = 1)))

# Mean of the 1 million maximum values should approximate truth
true_theta_max = mean(theta_hat_max)

# theta_S or theta_max
true_theta_max
hist(theta_hat_max, breaks=50,freq = FALSE)

