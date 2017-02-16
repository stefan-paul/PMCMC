
#' Function to run the particle MCMC
#' @param iterations number of iterations to run
#' @param numParticles number of particles in the particle filter
#' @param model model to make predictions about the state
#' @param likelihood likelihood function (see details)
#' @param prior log prior density
#' @param scaling scaling factor for the proposal. If not provided 1/10 of the
#' prior with is used.
#' @param f0 function that generates a set of particles for start conditions (see details)
#' @param startValues vector containing the parameter start values. Alternatively
#' a function can be provided to sample the start values.
#' @param observations vactor containing the data (e.g. PAR) to run the model.
#' @param names parameter names (otional)
#' @param parallel either FALSE number of nodes to be used in parallel processing
#' @param parallelOptions list with packages and objects to be exported to cluster
#' @export
#' @details The likelihood calculates the weigth for one set of particles. \cr\cr
#' The f0 function should accept the argument numParticles for the number of start particles.


PMCMC <- function(iterations, model, likelihood, prior, f0, startValues, numParticles,
                  scaling = NULL, observations, names = NULL, reportIntervall = 10,
                  parallel = 2, parallelOptions = list(packages = NULL, objects = NULL)){
  
  #### Initializations ###########
  failed = FALSE
  currentlP <- NA # current log posterior
  weights <- matrix(NA, nrow = length(observations), ncol = numParticles)
  
  if(is.numeric(startValues)){
    currentPar <- startValues # current parameter set
  }else{
    currentPar <- startValues()
  }
  
  if(parallel){
  ## set up parallel executer 
  cl <- parallel::makeCluster(parallel)
  # Export packages and objects to cluster
  parallel::clusterCall(cl, packageFun, parallelOptions$packages)
  parallel::clusterExport(cl, varlist = parallelOptions$objects)
  }
  
  # Initialize MCMC chain
  chain <- matrix(NA, iterations, (length(currentPar)+3))
  if(is.null(names)) names <- 1:length(currentPar)
  names <- c(names, "LP", "LL", "LP")
  colnames(chain) <- names
  
  
  # Initialize start particles
  particles <-f0(numParticles)
  
  
  ####### End initializations ####
  
  
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
    
    
    # weights[i,] <- likelihood(predicted = particles[,1],
    #                           observed = referenceData[i,1], error = currentPar[7]) # calculate weights for all
    # 
    
    weights[i,] <- likelihood(predicted = particles,
                              observed = referenceData[i,], currentPar = currentPar) # calculate weights for all
    
    
     # particles based on actual observations
    if(all(weights[i,] == 0)) stop("Please provide better starting values")
    
    #### Resample particles based on weights
    indX <- sample(x = 1:numParticles,size = numParticles, prob = weights[i,], replace = TRUE)
    
    particles <- particles[indX,]
    
    
  } # observations
  
  currentll <- calculateApproxLikelihood(weights) # get approximate ll values
  currentPrior <- prior(currentPar)
  currentlP <- currentll +  currentPrior# get posterior values
  
  
  chain[1,] <- c(currentPar, currentlP, currentll, currentPrior) # write first values in chain
  
  for(iter in 2:iterations){
    
    # Initialize start particles
    particles <-f0(numParticles)
    
    
    #### Make new proposal
    currentPar <- proposalFunction(chain[iter-1,1:(ncol(chain)-3)], scaling = scaling)
    newPrior <- prior(currentPar)
    
    if(newPrior == -Inf){
      chain[iter, ] <- chain[(iter-1), ]
    }else{
      
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
          
        # weights[i,] <- likelihood(predicted = particles[,1],
        #                           observed = referenceData[i,1], error = currentPar[7]) # calculate weights for all
        # # particles based on actual observations
        # 
        weights[i,] <- likelihood(predicted = particles,
                                  observed = referenceData[i,], currentPar = currentPar) # calculate weights for all
        
        
        if(all(weights[i,] ==0)){
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
        chain[iter, ] <- chain[(iter-1), ]
      }else{
        #### Acceptance
        
        a <- min(1, newlP / currentlP)
        
        if(a > runif(1)){ # accept and write old values in chain
          currentlP <- newlP
          chain[iter, ] <- c(currentPar, currentlP, newll, newPrior)
        }else{ # reject and keep old values
          chain[iter, ] <- chain[(iter-1), ]
        }
        
      } # if ll != -Inf
    } # if prior != -Inf
    
    if(iter %% reportIntervall == 0){
      cat("\r","Running PMCMC! Iteration " ,iter," of",iterations,". Current logp ",
          currentlP,". Please wait!","\r")
      flush.console()
    }
    
  } # iterations
  
  if(parallel) parallel::stopCluster(cl)
  
  return(coda::mcmc(chain))
  
}









##### Helper functions


#' Function to calculate the approximate likelihood based on weights
#' @param weights vector with weights
#' @return approximate likelihood
calculateApproxLikelihood <- function(weights){
  np <- length(weights[1,]) # number of particles
  
  w_hat <- numeric()
  
  for(i in 1:nrow(weights)){
    w_hat[i] <- 1/np * sum(weights[i,])
  }
  
  margLL <- sum(log(w_hat))
  return(margLL)
}

#' Function to generate proposal
#' @param parameter vector with current parameter values
#' @param scaling scaling factor for each parameter
#'
#' @return vactor with new proposal

proposalFunction <- function(parameter, scaling){
  
  proposal <- rnorm(length(parameter),mean = parameter, sd = scaling)
  return(proposal)
}


#' Function to export packages to cluster
packageFun <- function(packages = NULLL) {
  if(!is.null(packages)){
    for(i in packages) library(i, character.only = TRUE)
  }
}


