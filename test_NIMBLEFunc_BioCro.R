#------------------ SETUP
library(nimble)
library(mcmcplots)
library(tidygraph)
library(ggraph)
library(BioCro)
set.seed(1)
rm(list = ls())

output.cleaning<- function(bio.out){
  # This taken from OpBioGro {BioCro}
  ans.dat <- as.data.frame(unclass(ans)[1:11])
  sel.rows <- seq(1,nrow(ans.dat),length.out=8)
  simDat <- ans.dat[sel.rows,
                    c("ThermalT","Stem","Leaf","Root","Rhizome","Grain","LAI")]
  return(simDat$LAI)
}
#TRUE model
#------------------------------------------------------ Simulate OBS
data("weather05")
ans <- BioGro(weather05)
y <- output.cleaning(ans) 
y <- y + rnorm(length(y), 0.5) # adding errors to make obs
y <- abs(y)
#------------------------------------------------- NIMBLE code
# This is the function that can call the process-based model
predic_R_func <- function(sp){
  
  tmp.ans <- BioGro(weather05, canopyControl = canopyParms(Sp=sp))

  LAI.sim <- output.cleaning(tmp.ans)
  
 return(LAI.sim) 
}

# --------------------------------------------This introduces my function to nimble
predic_NIMBLE_func <- nimbleRcall(
  function(sp = double()){}, Rfun = 'predic_R_func',
returnType = double(1) # This means the output is a vector - double(2) means it's a matrix
)

#--------------------------------------This does the real fitting
code <- nimbleCode({

  sp ~ T(dnorm(1, 3), 0, 5) # truncated between 0-5, normal with mean=1 and sd=3
  sigma ~ dunif(0, 100)
  predicted.y[1:N] <- predic_NIMBLE_func(sp) # We are estimating the mean
  
  for(i in 1:N) { # for each obs - This is the liklihood
    y[i] ~ dnorm(predicted.y[i], sd = sigma)
  }
})
# this checks the model and registers it as an R obj
model <- nimbleModel(code, 
                     data = list(y = y),
                     inits = list(sp=1,
                                  sigma=0.5),
                     constants = list(N=length(y))
                     )

model$initializeInfo()
#--------------------------------- Visualizing the graph
ggraph(model$graph, layout = 'auto') +
  geom_edge_link() +
  geom_node_point() +
  geom_node_label(aes(label=name))+
  theme_light()+
  theme(legend.position = 'bottom')

#----------------------------------------------------- Compilation
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()
#mcmcConf$addMonitors("x") # By default, top-level parameters are monitored. Let's get these random effects too.
DEmcmc1 <- buildMCMC(mcmcConf)
# One way to do this -----------------------------

# compiled <- compileNimble(model, DEmcmc1) # Example of compiling model and MCMC in one call, which returns a list.
# 
# DEsamples <- runMCMC(compiled$DEmcmc1,
#                      niter = 10000,
#                      nburnin = 1000,
#                      samplesAsCodaMCMC = TRUE)

#---- Alternative way to the two line code for running MCMC
DEsamples <- nimbleMCMC(code=model,
                        niter = 1000,
                        nburnin = 100,
                        nchains=3,
                        samplesAsCodaMCMC = TRUE)
printErrors()
#summary state
nimble::samplesSummary(DEsamples)
#------------------- VIS MCMC
mcmcplots::mcmcplot(DEsamples)
