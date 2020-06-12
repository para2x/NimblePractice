#------------------ SETUP
library(nimble)
library(mcmcplots)
library(tidygraph)
library(ggraph)
set.seed(1)
rm(list = ls())
#TRUE model
x <- 1:100
y <- (2*x + 0.5) + rnorm(length(x),0, 0.2)

#------------------------------ NIMBLE code
# This is the function that can call the APSIM/PSIMS
predic_R_func <- function(x, slope, intercept){
 return(intercept + slope * x) 
}

# This introduces my function to nimble
predic_NIMBLE_func <- nimbleRcall(function(x = double(),
                                           slope = double(),
                                           intercept = double()
                              ){},
                     Rfun = 'predic_R_func',
                     returnType = double()
                     )

# This does the real fitting
code <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  
  for(i in 1:N) {
    predicted.y[i] <- predic_NIMBLE_func(x[i], slope, intercept) # We are estimating the mean
    y[i] ~ dnorm(predicted.y[i], sd = sigma)
  }
})

model <- nimbleModel(code, 
                     data = list(y = y),
                     inits = list(intercept = 0.5, 
                                  slope = 0.2, 
                                  sigma = 1,
                                  x = x
                                  ),
                     constants = list(N=length(x))
                     )


#--------------------------------- Visualizing the graph
# ggraph(model$graph, layout = 'auto') +
#   geom_edge_link() +
#   geom_node_point() +
#   geom_node_label(aes(label=name))+
#   theme_light()+
#   theme(legend.position = 'bottom')

#----------------------------------------------------- Compilation
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()
mcmcConf$addMonitors("x") # By default, top-level parameters are monitored. Let's get these random effects too.
DEmcmc1 <- buildMCMC(mcmcConf)
# One way to do this -----------------------------

# compiled <- compileNimble(model, DEmcmc1) # Example of compiling model and MCMC in one call, which returns a list.
# 
# DEsamples <- runMCMC(compiled$DEmcmc1,
#                      niter = 10000,
#                      nburnin = 1000,
#                      samplesAsCodaMCMC = TRUE)

#---- Alternative to the two line code for running MCMC
DEsamples <- nimbleMCMC(code=model,
                        niter = 10000,
                        nburnin = 1000,
                        samplesAsCodaMCMC = TRUE)
#summary state
nimble::samplesSummary(DEsamples)
#------------------- VIS MCMC
mcmcplots::mcmcplot(DEsamples)
