############################################################################
# all of the below is from Wright et al. 2017 (should be unchanged except for comments)
# it is including here for learning purposes, but note that the model of Wright is different
# to the NPMS data structure in a least one key way (multiple surveys of a plot within a single visit)
# another difference is the greater variability in the data structure (numbers of plots, visits per year)
############################################################################

##Functions for performing zero-augmented beta with errors (ZABE) regression
##simulations.

##See the "Read me" document for instructions on this R script and the
##related files.

##Author: Wilson Wright
##Date: 11/16/2016
##R version 3.3.1 (2016-06-21)

##Auxillary Functions
#Calculate the mean of a vector without including any zero values
mean.no.zero <- function(x){
  mean(x[which(x > 0)])
}

#Calculate the cumulative sum of a vector and replace all zeroes with ones
cumsum.no.zero <- function(x){ # op: presumably this is used as an indicator of zero value vectors
  temp <- cumsum(x)
  temp[which(temp==0)] <- 1
  return(temp)
}

#Create a data matrix for Model4 that excludes plots with no detections and
#replaces observations with NAs when the observer doesn't detect the species
u.mat.fun <- function(x.mat, u.mat){
  temp <- t(apply(x.mat*u.mat, 1, sort, decreasing=TRUE))
  temp[which(temp==0)] <- NA
  return(qlogis(temp[which(rowSums(x.mat) > 0), ]))
}

##calculate the number of divergent transitions from a model fit in Stan
count_divergences <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  sum(sapply(sampler_params, function(x) sum(x[,'divergent__']))) 
}

## slow, function compiles model as well as loading it
rstan_options(auto_write = TRUE) # avoid recompiling unchanged programs (apparently!)
model_full <- rstan::stan_model(file = "scripts/stan/model_zabe.stan")

##Function to run the simulation based on specified parameter values.
##User is able to imput sample size and observations per plot as well.
##For each iteration, function will simulate a dataset and fit three models
##based on the different approaches outlined in the paper.
##Approach 1: 'naive' - basic zero-augmented beta (ZAB) for first observer only
##Approach 2: 'ad hoc' - ZAB for averaged observer responses
##Approach 3: 'full data' based on our ZABE model using all responses
# values in comments below are OP values for testing (taken partly from Wright et al. 2017 simulation runs file)
run.sim <- function(N,        #number of plots # 100
                    J,        #number of observations per plot # 2
                    iter,     #number of iterations # 1
                    mu,       #parameter for mean of cover beta distribution # 0.25
                    phi,      #parameter for 'precision' of cover distribution # 3
                    psi,      #parameter for occupancy # 0.5
                    gamma0,   #intercept for detection logistic regression # -1.5
                    gamma1,   #slope for detection with %cover # 2
                    sigma){   #measurement variability on logit scale # 0.5
  
  #generate the cover values for every iteration
  y.mat <- matrix(rbinom(N*iter, 1, psi)*rbeta(N*iter, mu*phi, (1-mu)*phi), 
                  nrow=N, ncol=iter)
  
  #make a detection history matrix for every iteration
  x.array <- array(dim = c(N, J, iter))
  for(k in 1:iter){
    for(i in 1:N){
      for(j in 1:J){
        x.array[i, j, k] <- ifelse(y.mat[i, k] > 0,
                                   rbinom(1, 1, 
                                          plogis(gamma0 + gamma1*y.mat[i, k])),
                                   0)
      }
    }
  }
  
  #make the measurement matrix for every iteration
  u.array <- array(dim = c(N, J, iter))
  for(k in 1:iter){
    for(i in 1:N){
      for(j in 1:J){
        u.array[i, j, k] <- ifelse(y.mat[i, k] > 0,
                                   plogis(rnorm(1, qlogis(y.mat[i, k]), sigma)),
                                   0)
      }
    }
  }
  
  #create storage arrays
  output1 <- output2 <- array(NA, c(5, 4, iter))
  output3 <- array(NA, c(8, 4, iter))
  n_diverge <- matrix(NA, ncol=3, nrow=iter)
  
  #name the appropriate dimensions for all of the storage arrays
  dimnames(output1)[[1]] <- dimnames(output2)[[1]] <- c("psi", "mu", "phi", "phi_t", "lp__")
  dimnames(output1)[[2]] <- dimnames(output2)[[2]] <- c("mean", "2.5%", "97.5%", "Rhat")
  dimnames(output3)[[1]] <- c("psi", "mu", "phi", "phi_t", "sigma", "gamma0", "gamma1", "lp__")
  dimnames(output3)[[2]] <- c("mean", "2.5%", "97.5%", "Rhat")
  
  #create the progress bar to display in the R console
  pb <- txtProgressBar(min=0, max=iter, initial=0, 
                       style=3, width=40, char="+")
  
  #loop over the iterations
  for(k in 1:iter){
    #extract simulated dataset corresponding to iteration k
    u.temp <- u.array[, , k] # first iteraction, two samples consisting of recorded covers (normally distributed around true cover y.mat)
    x.temp <- x.array[, , k] # detection history for first iteration, based on underlying true cover
    u.temp[which(x.temp==0)] <- 0 # if the plant was not actually detected, then set the recorded cover to zero as well
    
    #obtain data needed for each of the three approaches investigated
    N1 <- N
    visits1 <- rep(J, N) # number of visits per plot
    pos1 <- cumsum(c(1, visits1[-N])) # vector to index the vists associated with each plot, for use with segment in Stan
    n_obs1 <- N*J # total visits
    x <- as.vector(t(x.temp))# vector of detection histories
    N2 <- length(which(rowSums(x.temp) > 0)) # number of plots with detections
    v2 <- rowSums(x.temp)[which(rowSums(x.temp)>0)] # number of dections per plot where detected
    pos2 <- cumsum(c(1, v2[-N2]))
    n_obs2 <- sum(v2)
    logit_u <- qlogis(as.vector(t(u.temp))[which(as.vector(t(u.temp))>0)])
    N1_ind <- which(rowSums(x.temp)>0)
    
    #Create Stan data lists for each of the three approaches
    data1 <- list('N1' = N1,
                  'x' = x.temp[, 1],
                  'N2' = length(which(x.temp[, 1] > 0)),
                  'y' = u.temp[which(x.temp[,1] > 0), 1])
    data2 <- list('N1' = N1,
                  'x' = apply(x.temp, 1, max),
                  'N2' = length(which(rowSums(x.temp) > 0)),
                  'y' = apply(u.temp[which(rowSums(x.temp) > 0), ], 1, mean.no.zero))
    data3 <- list('N1' = N1,
                  'visits1' = visits1,
                  'pos1' = pos1,
                  'n_obs1' = n_obs1,
                  'x' = x,
                  'N2' = N2,
                  'visits2' = v2,
                  'pos2' = pos2,
                  'n_obs2' = n_obs2,
                  'logit_u' = logit_u,
                  'N1_ind' = N1_ind)
    
    #MCMC sampling done in Stan. More iterations for warmup in ZABE model
    samps1 <- sampling(model_naive, data1, iter=1000, warmup=500,
                       pars=c("psi", "mu", "phi", "phi_t"))
    samps2 <- sampling(model_naive, data2, iter=1000, warmup=500,
                       pars=c("psi", "mu", "phi", "phi_t"))
    samps3 <- sampling(model_full, data3, iter=150, warmup=100,
                       pars=c("psi", "mu", "phi", "phi_t",
                              "sigma", "gamma0", "gamma1"))
    
    #Summarizes posterior samples, keeping the quantiles of interest
    summary1 <- summary(samps1, probs=c(0.025, 0.5, 0.975))$summary
    summary2 <- summary(samps2, probs=c(0.025, 0.5, 0.975))$summary
    summary3 <- summary(samps3, probs=c(0.025, 0.5, 0.975))$summary
    
    #save the summary matrix in the appropriate output array
    output1[, , k] <- summary1[, c("mean", "2.5%", "97.5%", "Rhat")]
    output2[, , k] <- summary2[, c("mean", "2.5%", "97.5%", "Rhat")]
    output3[, , k] <- summary3[, c("mean", "2.5%", "97.5%", "Rhat")]
    
    #save the number of divergent transitions for each model fit as well
    n_diverge[k, 1] <- count_divergences(samps1)
    n_diverge[k, 2] <- count_divergences(samps2)
    n_diverge[k, 3] <- count_divergences(samps3)
    
    setTxtProgressBar(pb, k)
    
  }
  #return a list with each array as an object
  return(list(naive=output1,
              averaged=output2,
              full.data=output3,
              n_diverge=n_diverge))
}

############################################################################
# all of the above is from Wright et al. 2017
############################################################################