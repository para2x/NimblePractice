#------------------ SETUP
library(nimble)
library(mcmcplots)
library(tidygraph)
library(ggraph)
set.seed(1)

#TRUE model
x <- 1:100
y <- (2*x + 0.5) + rnorm(length(x))

plot(x, y, type='o')
#------------------------------ NIMBLE code
code <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  
  for(i in 1:N) {
    predicted.y[i] <- intercept + slope * x[i]
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
ggraph(model$graph, layout = 'auto') +
  geom_edge_link() +
  geom_node_point() +
  geom_node_label(aes(label=name))+
  theme_light()+
  theme(legend.position = 'bottom')

#----------------------------------------------------- Compilation
mcmcConf <- configureMCMC(model)
mcmcConf$printSamplers()

#mcmcConf$addMonitors("farm_effect") # By default, top-level parameters are monitored. Let's get these random effects too.
DEmcmc1 <- buildMCMC(mcmcConf)
compiled <- compileNimble(model, DEmcmc1) # Example of compiling model and MCMC in one call, which returns a list.

DEsamples <- runMCMC(compiled$DEmcmc1,
                     niter = 10000,
                     nburnin = 1000,
                     samplesAsCodaMCMC = TRUE)

nimble::samplesSummary(DEsamples)
#------------------- VIS MCMC
mcmcplots::mcmcplot(DEsamples)
