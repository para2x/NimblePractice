library(nimble)
nregions <- 10 # nrow

y <- rnorm(nregions, 5, 0.1) # this is what we try to model / Liklihood (could SSQ)


expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma   # calculate once
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    return(result)
  })


code <- nimbleCode({
  mu0 ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
  sigma2 ~ dunif(0, 100)
  rho ~ dunif(0, 5)
  #beta ~ dnorm(0, sd = 100) if there is a covariate 
  
  # latent spatial process
  mu[1:N] <- mu0*ones[1:N] # mean
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma) # cov
  x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # likelihood
  for(i in 1:N) {
    #lambda[i] <- exp(beta*z[i] + x[i]) # if there is a process model - for example, z covariate and B is now the coeff of cov
    y[i] ~ dnorm(x[i], sigma2)
  }
})

locs <-  as.matrix(data.frame(a = rnorm(nregions), # two param space
                              b = rnorm(nregions)
                              )
                   )

dists <- as.matrix(dist(locs))
dists <- dists / max(dists)  # normalize to max distance of 1

constants <- list(N = nregions,
                  dists = dists,
                  ones = rep(1, nregions)
                 )

data <- list(y = y)
inits <- list(beta = 0, mu0 = 0, sigma = 1, sigma2 = 1, rho = 0.2)

set.seed(1)

## setup initial spatially-correlated latent process values
inits$cov <- expcov(dists, inits$rho, inits$sigma)
inits$x <-  t(chol(inits$cov)) %*% rnorm(nregions)
inits$x <- inits$x[ , 1]  # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)

# 
ggraph(model$graph, layout = 'auto') +
  geom_edge_link() +
  geom_node_point() +
  geom_node_label(aes(label=name))+
  theme_light()+
  theme(legend.position = 'bottom')


#----------------
conf <- configureMCMC(model)
conf$addMonitors('x')
conf$removeSamplers('x[1:10]')
## Changing a tunable parameter in the adaptation of RW_block makes a big difference.
conf$addSampler('x[1:10]', 'RW_block', control = list(adaptFactorExponent = 0.25))

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 3000, nburnin = 500, thin = 25)

mcmcplots::mcmcplot(samples)

#----------------------------------------------------------------------
#post-hoc posterior sampling of the Gaussian process 


newlocs <- rbind(c(2.6, 6.7)
)

dist11 <- fields::rdist(locs)
dist21 <- fields::rdist(newlocs, locs)
dist22 <- fields::rdist(newlocs)


#----------------- Chris's solution

get_samples <- nimbleFunction(
  run = function(samples = double(2), z = double(1), zNew = double(1), 
                 dist11 = double(2), dist22 = double(2), dist21 = double(2)) {
    returnType(double(2))
   
    m <- dim(samples)[1]
    nstar <- dim(dist21)[2]
    output <- matrix(0, nrow = m, ncol = nstar)
    for(i in 1:m) {
      # Extract parameter values from the input matrix based on numeric column indexes.
      # Inelegant because hard-coded based on looking at column names of samples matrix.
      mu0 <- samples[i, 1]
      #beta <- samples[i, 1]
      rho <- samples[i, 2]
      sigma <- samples[i, 4]
      x <- samples[i, 5:14]
      mu <- mu0 #+ beta*z no covariate
      muNew <- mu0 #+ beta*zNew
      
      #browser()
      #--------------
      n <- length(muNew)
      sigma2 <- sigma*sigma
      C22 <- sigma2 * exp(-dist22 / rho)
      C11 <- sigma2 * exp(-dist11 / rho)
      C21 <- sigma2 * exp(-dist21 / rho)
      # Note that this could be made a bit more efficient by using the Cholesky
      # decomposition rather than solve().
      xstar <- muNew + (C21 %*% solve(C11, x - mu))[,1]
      xstar <- xstar + (t(chol(C22 - C21 %*% solve(C11, t(C21)))))[,1]
      # Get other parameters
      output[i, ] <-xstar
    }
    return(output)
  }
)


xstar_samples <- get_samples(samples, locs, newlocs, dist11, dist22, dist21)



cget_samples <- compileNimble(get_samples)
xstar_samples2 <- cget_samples(samples, z, Znew, dist11, dist22, dist21)


#---------------- My failed attempt

sample_xstar <- nimbleFunction({
  # need types for inputs
  run = function(x=double(1),
                 mu=double(0),
                 sigma=double(0),
                 rho=double(0),
                 dist11=double(2), dist22=double(2), dist21=double(2)) {
    returnType(double(2))
    
    x ~ dmnorm(mu+dist21%*%solve(dist11)%*%(x−mu),
               dist22-dist21%*%solve(dist11)%*%t(dist21))
    
    return(x)
  }
})

get_samples <- nimbleFunction(
  # need types for inputs
  run = function(samples=double(2),
                 dist11=double(2),
                 dist22=double(2),
                 dist21=double(2)
  ) {
    returnType(double(2))
    browser()
    m <- dim(samples)[1]
    nstar <- 3
    output <- matrix(0, nrow = m, ncol = nstar)
    
    for(i in 1:m) {
      # Extract parameter values from the input matrix based on numeric column indexes
      mu <- samples[i, 2]
      sigma <- samples[i, 3]
      rho <- samples[i, 4]
      x <-samples[i, 6]
      # Get other parameters
      output[i, ] <- sample_xstar(x,mu, sigma, rho, dist11, dist22, dist21)
    }
    return(output)
  }
)


xstar_samples <- get_samples(samples, dist11, dist22, dist21)

cget_samples <- compileNimble(get_samples)
xstar_samples2 <- cget_samples(samples, dist11, dist22, dist21)
