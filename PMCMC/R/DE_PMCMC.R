
#' Differential-Evolution MCMC zs with a particle filter to approximate the likelihood
#' @param prior a prior created by createPrior
#' @param iterations iterations to run
#' @param pSnooker probability of Snooker update
#' @param burnin number of iterations treated as burn-in. These iterations are not recorded in the chain.
#' @param thin thinning parameter. Determines the interval in which values are recorded.
#' @param eps small number to avoid singularity
#' @param f scaling factor gamma
#' @param parallel logical, determines weather parallel computing should be attempted (see details)
#' @param pGamma1 probability determining the frequency with which the scaling is set to 1 (allows jumps between modes)
#' @param eps.mult random term (multiplicative error)
#' @param eps.add random term
#' @param numParticles number of particles 
#' @param observations vector with data
#' @param reportIntervall intervall in which the sampler's progress is printed 
#' @param parallel either FALSE or number of nodes to be used in parallel processing
#' @param parallelOptions list with packages and objects to be exported to the cluster
#' @param zUpdateFrequency frequency with which the Z matrix is updated (= history of the chain)
#' @references ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain with snooker updater and fewer chains. Statistics and Computing http://dx.doi.org/10.1007/s11222-008-9104-9 
#' @export


DE_PMCMC <- function(iterations, model, likelihood, prior, f0, startValues, numParticles,
                     observations, names = NULL, reportIntervall = 10,
                     parallel = FALSE, parallelOptions = list(packages = NULL, objects = NULL),
                     pSnooker = 0.1, 
                     burnin = 0, 
                     thin = 1,
                     f = 2.38,
                     eps = 0,
                     pGamma1 = 0.1,
                     eps.mult =0.2,
                     eps.add = 0,
                     zUpdateFrequency = 10)
{  
  
  parLen = length(prior$sampler(1))
  X = prior$sampler(3)
  Zold = prior$sampler(parLen * 10)
  

  if (! is.matrix(X)) stop("wrong starting values")
  if (! is.matrix(Zold)) stop("wrong Z values")
  
  
  Npar <- ncol(X)
  Npar12 <- (Npar - 1)/2 # factor for Metropolis ratio DE Snooker update
  
  # M0 is initial population size of Z is the size of Z, it's the same number, only kept 2 to stay consistent with the ter Brakk & Vrugt 2008 
  M = M0 = nrow(Zold)
  Npop <- nrow(X)
  
  F2 = f/sqrt(2*Npar)
  F1 = 1.0
  rr = NULL
  r_extra = 0
  
  burnin = 0
  
  n.iter <- ceiling(iterations/Npop)
  
  lChain <- ceiling((n.iter - burnin)/thin)+1
  
  chain <- array(NA, dim=c(lChain, Npar+3, Npop))
  
  if(is.null(names)) names <- 1:Npar
  
  colnames(chain) <- c(names, "LP", "LL", "LPr")

  
  Z <- matrix(NA, nrow= M0 + floor((n.iter-1) /zUpdateFrequency) * Npop, ncol=Npar)
  
  Z[1:M,] <- Zold
  
  
  counter <- 1
  counterZ <- 0
  
  
  
  
  if(parallel){
    ## set up parallel executer 
    cl <- parallel::makeCluster(parallel)
    # Export packages and objects to cluster
    parallel::clusterCall(cl, packageFun, parallelOptions$packages)
    parallel::clusterExport(cl, varlist = parallelOptions$objects)
  }
  

  

  for(pop in 1:Npop){
    
    currentPar <- X[pop,]
    
  # Initialize start particles
  particles <-f0(numParticles)

  
  #### Evaluation of first parameter set ####
  for(i in 1:length(observations)){
    
    
    if(parallel){
      
      particles <- matrix(parallel::parRapply(cl, x = particles, FUN = model, 
                                              data = observations[i], parameters = currentPar),
                          nrow = numParticles, byrow = TRUE)
      
    }else{
      for(k in 1:numParticles){
        particles[k,] <- model(data = observations[i], parameters = currentPar, particles = particles[k,]) # update states
      }
    }
    
    
    weights[i,] <- likelihood(predicted = particles[,1],
                              observed = referenceData[i,1], error = currentPar[7]) # calculate weights for all
    # particles based on actual observations
    if(all(weights[i,] == 0)) stop("Please provide better starting values")
    
    #### Resample particles based on weights
    indX <- sample(x = 1:numParticles,size = numParticles, prob = weights[i,], replace = TRUE)
    
    particles <- particles[indX,]
    
    
  } # observations
  
  currentll <- calculateApproxLikelihood(weights) # get approximate ll values
  currentPrior <- prior$density(currentPar)
  currentlP <- currentll +  currentPrior# get posterior values
  
  
  chain[1,,pop] <- c(currentPar, currentlP, currentll, currentPrior) # write first values in chain
  
  }
  
  
  ############### Start iterations #################
  
  for (iter in 2:n.iter) {
    
    counter <- counter + 1 ## TODO this needs to be changes once 
    # thinning is included in the sampling
    
    f <- ifelse(iter%%10 == 0, 0.98, F1)
    #accept <- 0
    
      
      for (i in 1:Npop){
        # select to random different individuals (and different from i) in rr, a 2-vector
        rr <- sample.int(M, 3, replace = FALSE)
        if(runif(1) < pSnooker) {
          z <- Z[rr[3],]
          x_z <- X[i,] - z  
          D2 <- max(sum(x_z*x_z), 1.0e-300)
          projdiff <- sum((Z[rr[1],] -Z[rr[2],]) * x_z)/D2 # inner_product of difference with x_z / squared norm x_z
          gamma_snooker <- runif(1, min=1.2,max=2.2)
          x_prop <- X[i,] + gamma_snooker * projdiff * x_z
          x_z <- x_prop - z
          D2prop <- max(sum(x_z*x_z), 1.0e-300)
          r_extra <- Npar12 * (log(D2prop) - log(D2))
        } else {
          
          if ( runif(1)< pGamma1 ) { gamma_par = F1 # to be able to jump between modes
          } else {
            gamma_par = F2 * runif(Npar, min=1-eps.mult, max=1+eps.mult)    # multiplicative error to be applied to the difference
            # gamma_par = F2 
          }
          rr = sample.int(M, 2, replace = FALSE)
          if (eps.add ==0) {  # avoid generating normal random variates if possible
            x_prop = X[i,] + gamma_par * (Z[rr[1],]-Z[rr[2],]) } else {
              x_prop = X[i,] + gamma_par * (Z[rr[1],]-Z[rr[2],])  +  eps.add*rnorm(Npar,0,1)
            }
          r_extra = 0
          
        }

        
        
        # evaluate proposal 
        
        # Initialize start particles
        particles <-f0(numParticles)
        
        
        #### Make new proposal
        currentPar <- x_prop
        newPrior <- prior$density(currentPar)
        
        if(newPrior == -Inf){
          chain[iter, ,] <- chain[(iter-1), ,]
        }else{
          
          for(obs in 1:length(observations)){
            
            
            if(parallel){
              particles <- matrix(parallel::parRapply(cl, x = particles, FUN = model, 
                                                      data = observations[obs], parameters = currentPar),
                                  nrow = numParticles, byrow = TRUE)
              
            }else{
              
              for(k in 1:numParticles){
                particles[k,] <- model(data = observations[obs], parameters = currentPar, particles = particles[k,]) # update states
              }
            }
            
            weights[obs,] <- likelihood(predicted = particles[,1],
                                      observed = referenceData[obs,1], error = currentPar[7]) # calculate weights for all
            # particles based on actual observations
            
            
            if(all(weights[obs,] ==0)){
              failed = TRUE
              break
            }
            
            #### Resample particles based on weights
            indX <- sample(x = 1:numParticles,size = numParticles, prob = weights[i,], replace = TRUE)
            
            particles <- particles[indX,]
            
          } # observations
          
          
          if(failed == TRUE) {
            newll = newlP = -Inf
            failed = FALSE
          }else{
            newll <- calculateApproxLikelihood(weights) # get approximate ll values
            newlP <- newll + newPrior # get posterior values
          }
          
          if(newll == -Inf){
            chain[iter, ,] <- chain[(iter-1), ,]
          }else{
          
            
          a <- min(1, newlP / chain[(iter-1),(Npar+1),i])
          
          if(a > runif(1)){ # accept and write old values in chain
            X[i,] <- x_prop
            currentlP <- newlP
            chain[iter,,i] <- c(currentPar, currentlP, newll, newPrior)
          }else{ # reject and keep old values
            chain[iter,,i] <- chain[(iter-1), ,i]
          }
          
        } # if ll != -Inf
      } # if prior != -Inf
          
          
      } # for Npop
      

    
    
    if (iter%%zUpdateFrequency == 0) { # update history
      
      Z[( M0 + (counterZ*Npop) + 1 ):( M0 + (counterZ+1)*Npop),] <- X
      counterZ <- counterZ +1
      M <- M + Npop
    }
    # Console update
    
    if(iter %% reportIntervall == 0){
      cat("\r","Running PMCMC! Iteration " ,iter," of",n.iter,". Current logp ",
          currentlP,". Please wait!","\r")
      flush.console()
    }
    
  } # n.iter
  
  
  chain <- chain[1:counter,,]
  
  
  chain<- coda::as.mcmc.list(lapply(1:Npop,function(i) coda::as.mcmc(chain[,1:(Npar+3),i])))
  
  
  
  list(chain = chain,  X = as.matrix(X[,1:Npar]), Z = Z)
}