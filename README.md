# Group_Sequential_SMART_Design

- VerifyWs1.R finds stopping boundaries, and verifies the density of Ws1 and the probability of Ws1 with the simulated Ws1 distribution (stored in Simulated_Zsdata.rda).
  
- VerifyWs1Ws2Ws3Prob.R verifies the joint probability of Ws1,Ws2 and Ws3 being within some (lb, ub) using recursive integration. Something is not right with how I calculate these probabilities, need to debug and further verify to the simulated distributions of Ws1,Ws2 and Ws3 (stored in 3interim.anal.res.rda.)

- SimulateWS1WS2WS3.R simulates 3 interim analyses to get Ws1, Ws2 and Ws3. I do not use stopping boundaries here, I'm just simulating the values across 2 interim analyses and a final analysis had the best treatment continued. 
