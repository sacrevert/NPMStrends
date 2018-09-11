## X2. Simulation of the types of data that we are going to be modelling, using JAGS
# O.L. Pescott
# 06.09.2018
rm(list=ls())
######################################
library(R2jags)
######################################
set.seed(4323) # for reproducibility
######################################
## Define potentially useful functions
ilt<-function(x){ # inverse logit function
  exp(x)/(1+exp(x))
}
logit <- function(x){ # logit function
  log(x/(1 - x))
}
######################################

######################################
# Simulations based on those of Wright et al. 2017
######################################
N <- 100 # number of spatially unique plots
J <- 2 # number of visits to a plot within a year (assume constant for the moment)
psi <- 0.5 # true occupancy (average)
yr <- 3 # number of years
mu <- 0.25       #parameter for mean of cover beta distribution # 0.25
phi <- 3      #parameter for 'precision' of cover distribution # 3
gamma0 <- -1.5   #intercept for detection logistic regression # -1.5
gamma1 <- 2   #slope for detection with %cover # 2

# array of plot covers per visit per year
y.array <- array(dim = c(N, J, yr))
for(k in 1:yr){
      y.array[,,k] <- matrix(rbinom(N*J, 1, psi)*rbeta(N*J, mu*phi, (1-mu)*phi), 
                                 nrow=N, ncol=J)
}
# make a detection history matrix based on cover data
x.array <- array(dim = c(N, J, yr))
for(k in 1:yr){
  for(i in 1:N){
    for(j in 1:J){
      x.array[i, j, k] <- ifelse(y.array[i, j, k] > 0,
                                 rbinom(1, 1, 
                                        plogis(gamma0 + y.array[i, j, k])), # detection function 1
                                        #plogis(gamma0 + gamma1*y.array[i, j, k])), # detection function 2
                                 0)
    }
  }
}
###################################### END OF SIMS

######################################################
## Data/indicators required for running JAGS model
######################################################
N1 <- N # total number of spatially unique plots (locations)
Y <- yr # total number of years of monitoring covered by the data
 # total number of visits with positive covers
cpos <- as.vector(y.array)[which(as.vector(y.array) > 0)] # cover value for every positive cover
# indicator linking positive plot k to the visit in the total visit list; length = n.Plot.pos
plotInd <- apply(y.array, MARGIN=c(1,3), sum) # sum across visits and years within a plot (margin 1 = plot, margin 3 = year), visit information is summed
plot <- c(which(as.vector(plotInd[,1]) > 0), which(as.vector(plotInd[,2]) > 0), which(as.vector(plotInd[,3]) > 0)) # need plot indices within years
n.Plot.pos <- length(plot)
# indicator linking positive plot k to the year of the its visit in the total visit list; length = n.Plot.pos
year <- c(rep(1, length(which(as.vector(plotInd[,1]) > 0))), 
          rep(2, length(which(as.vector(plotInd[,2]) > 0))), 
          rep(3, length(which(as.vector(plotInd[,3]) > 0)))
          )
#############
### vc CALc NEEDS CHANGING WHEN THERE ARE UNEVEN VISITS NUMBERs ETC. BETWEEN YEARS
#############
V2 <- N*J*yr # total number of plot visits (# plots x # visits x # years)
# indicator linking every visit x to plot; length = V2 
plotZ <- rep(1:N, J*yr)
# indicator linking every visit x to year; length V2 
yearZ <- rep(1:3, each = N*J)
x <- as.vector(x.array) # detection indicator for every visit (length V2)

# Data list for passing to JAGS
Data <- list(N1 = N1,
            Y = Y,
            n.Plot.pos = n.Plot.pos,
            cpos = cpos,
            plot = plot,
            year = year,
            V2 = V2,
            plotZ = plotZ,
            yearZ = yearZ,
            x = x)
###################################### END OF DATA PREP

###########################################
## Initialisation of values for JAGS chains
###########################################
# Initial parameter values for JAGS
# To ensure conformity with the data all occupancies (ZI.*) are set to one
# Some other parameters are also fixed to avoid extemly small likelihoods
# but these can probably be relaxed a bit more.
zinit <- matrix(1, nrow = N1, ncol = Y)
inits.fn <- function() list(z = zinit,
                            tau.C = runif(1,1,5),
                            mu.C = 0.5,
                            gamma0 = rnorm(1,0,1)
                            )

######################################
## JAGS model
######################################
sink('scripts/JAGS/JAGS_v0.0.txt')
cat("
model{
## State model
for (i in 1:N1){ # N1 is the number of plots
  for (j in 1:Y){ # number of years
    C.S[i,j] <- z[i,j] * c.Pos[i,j] # cover including zeros
    z[i,j] ~ dbern(psi[i,j]) ## true PA state of a plot within a year depends on occupancy probability psi
    psi[i,j] ~ dunif(0,1)
    c.Pos[i,j] ~ dbeta(a.C[i,j], b.C[i,j]) T(0.00001,0.99999)
    a.C[i,j] <- mu.C * tau.C
    b.C[i,j] <- (1 - mu.C) * tau.C
    } # end of years loop
  } # end of plot loop

## Plot positive covers
for(k in 1:n.Plot.pos){ 
    ## HERE YOU CAN POTENTIALLY ACCOUNT FOR INTERVAL CENSORED NATURE OF DATA USING dinterval()
    cpos[k] ~ dbeta(a.C[plot[k], year[k]], b.C[plot[k], year[k]]) T(0.00001,0.99999) # recorded cover when present follows beta distribution
    #a.P[k] <- C.Pos[plot[k], year[k]] * tau.P # link between final state and recorded cover
    #b.P[k] <- (1 - C.Pos[plot[k], year[k]]) * tau.P # link between final state and recorded cover
  }

## Observation model for all plot visits ([within-year] detection within plots)
for (a in 1:V2){
    x[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(0.00001, p.Dec[a]), 0.99999) # trick to stop numerical problems (note that this will probably influence covars in detectability regression) -- important?
    logit(p.Dec[a]) <- gamma0 ## + covars on detectability
  }

## Priors!
gamma0 ~ dt(0, 0.01, 1)
mu.C ~ dunif(0, 1)
tau.C ~ dt(0, 0.01, 1)T(0,)

} # END MODEL
", fill = TRUE)
sink()

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_v0.0.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
# Specify parameters for which posterior samples are saved
para.names <- c('mu.C', 'tau.C', 'gamma0')
# Continue the MCMC runs with sampling
samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
## Inspect results
plot(samples)
summary(samples)
gelman.diag(samples)