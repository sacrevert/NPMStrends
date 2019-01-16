## X2b. Simulation of the types of data that we are going to be modelling, using JAGS initially
# Now with interval censoring!
#
# note that this model only uses an observation/detection model for the repeat visits within a year, not for the cover data
# whilst it is reasonable to assume that the occupancy state is stable within years (obviously this is not 100% true all the time, but is reasonable)
# it is not reasonable to assume that the repeat visits to a plot provide replicated information about a stable cover for a species (as in Wright et al. 2017,
# where there are repeat observations, i.e. Quality Assurance, within a single visit), for this reason the cover values are taken directly as accurate 
# observations and are included in the "state model" as such
#
# O.L. Pescott
# 12.09.2018
#rm(list=ls())

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
## Could put this in a simulation function (see Wright et al. 2017)
N <- 100 # number of spatially unique plots
J <- 10 # number of visits to a plot within a year (assume constant for the moment)
#psi <- 0.5 # true occupancy (average)
psi <- 1 # should be easier for model to retrieve true values of gamma0 and gamma1 if occupancy is 1 -- this appears to be true
Y <- 1 # total number of years of monitoring covered by the data
mu <- 0.5       # parameter for mean of cover beta distribution # 0.25
phi <- 10      # parameter for 'precision' of cover distribution # 3
gamma0 <- -2   # intercept for detection logistic regression # -1.5
gamma1 <- 3   # slope for detection with %cover # 2
# e.g. plogis(-2 + 3*cover) makes for greater detectability range based on covers
# if you do this but only estimate gamma0 in the model, then biases in gamma0 and mu.C estimates increase

# array of plot covers per visit per year
y.array <- array(dim = c(N, J, Y))
for(k in 1:Y){
      y.array[,,k] <- matrix(rbinom(N*J, 1, psi)*rbeta(N*J, mu*phi, (1-mu)*phi), 
                                 nrow=N, ncol=J)
}
# make a detection history matrix based on cover data
x.array <- array(dim = c(N, J, Y))
for(k in 1:Y){
  for(i in 1:N){
    for(j in 1:J){
      x.array[i, j, k] <- ifelse(y.array[i, j, k] > 0,
                                 rbinom(1, 1, 
                                        #plogis(gamma0 + y.array[i, j, k])), # detection function 1
                                        plogis(gamma0 + gamma1*y.array[i, j, k])), # detection function 2
                                 #1, # for testing: if you remove probabilistic element, then mu.C better estimated
                                 # if you make detection perfect, then gamma0 and gamma1 estiamted as ~4.5 and 2
                                 # which means that plogis(4.5 + 2*0.01) and plogis(4.5 + 2*0.99) both = ~1.00 
                                 0)
    }
  }
}
y.arrayOrig <- y.array # keep original covers
y.array[which(x.array==0)] <- 0 # if the plant was not actually detected, then set the recorded cover to zero as well
###################################### END OF SIMS

##################################################
## Data/indicators required for running JAGS model
##################################################
# total number of visits with positive covers
cpos <- matrix(NA, nrow = length(as.vector(y.array)[which(as.vector(y.array) > 0)]), ncol = 3) # empty matrix for actual covers and cover classes
cpos[,1] <- as.vector(y.array)[which(as.vector(y.array) > 0)] # actual cover value for every positive cover
cpos[,2] <- 1 # indicator stating the observation censored
cpos[,3] <- NA # NA values for latent variable
# code borrowed from Pescott et al. 2016 Power paper:
t <- c(1e-3,0.05,
       0.05,0.25,
       0.25,0.5,
       0.5,0.75,
       0.75,0.95,
       0.95,0.999)
t = matrix(t, nrow = 6, ncol = 2, byrow = TRUE)
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
intervals <- c(1,2,3,4,5,6)
tdf$int <- intervals
# library(plyr) for SQL join function (like merge but keeps order) -- although actually could use merge with sort = F
tInt <- c(1e-3,0.05,0.25,0.5,0.75,0.95,0.999)
int <- findInterval(cpos[,1], tInt) # find corresponding interval for all data points
int <- as.data.frame(int)
m1 <- plyr::join(int, tdf, by = "int", type = "left", match = "all") # change to using merge at some point
lims <- m1[,2:3] # has the intervals for all points
check <- cbind(cpos, m1); head(check); tail(check) # just for quick visual check that all is well
cpos.Cens <- rep(1, nrow(cpos))
cpos.Latent <- cpos[,3] # just NAs

n.Plot.pos <- nrow(cpos)
# indicators linking positive plot k to the correct plots and years; length = nrow(cpos)
out <- out1 <-  numeric()
for(k in 1:Y){
  for(j in 1:J){
    for(i in 1:N){
    out <- ifelse(y.array[i, j, k] > 0, i, NA)
    out1 <- c(out1, out)
    }
  }
}
plot <- out1[!is.na(out1)]

out <- out1 <-  numeric()
for(k in 1:Y){
  for(j in 1:J){
    for(i in 1:N){
    out <- ifelse(y.array[i, j, k] > 0, k, NA)
    out1 <- c(out1, out)
    }
  }
}
year <- out1[!is.na(out1)]

#############
### NEEDS CHANGING WHEN THERE ARE UNEVEN VISITS NUMBERs ETC. BETWEEN YEARS (something to consider for real data)
#############
V2 <- N*J*Y # total number of plot visits (# plots x # visits x # years)
# indicator linking every visit x to plot; length = V2 
plotZ <- rep(1:N, J*Y)
# indicator linking every visit x to year; length V2 
yearZ <- rep(1:Y, each = N*J)
x <- as.vector(x.array) # detection indicator for every visit (length V2)
# covers for all visits
y <- as.vector(y.arrayOrig) # original covers used in detectability loop (inc. zeros) for every visit (length V2)

# Data list for passing to JAGS
Data <- list(N = N,
            Y = Y,
            n.Plot.pos = n.Plot.pos,
            #cpos.Cens = cpos[,2], # indicator (is censored?)
            cpos.Cens = cpos.Cens, # indicator (is censored?)
            cpos.Latent = cpos.Latent, # NA values for latent observations
            lims = lims,
            plot = plot,
            year = year,
            V2 = V2,
            plotZ = plotZ,
            yearZ = yearZ,
            x = x,
            yOrig = y)
###################################### END OF DATA PREP

###########################################
## Initialisation of values for JAGS chains
###########################################
# Initial parameter values for JAGS
# To ensure conformity with the data all occupancies (ZI.*) are set to one
# Some other parameters are also fixed to avoid extemly small likelihoods
# but these can probably be relaxed a bit more.
zinit <- matrix(1, nrow = N, ncol = Y)
inits.fn <- function() list(z = zinit,
                            tau.C = runif(1,1,10),
                            mu.C = 0.5,
                            gamma0 = rnorm(1,0,1),
                            gamma1 = rnorm(1,0,1),
                            cpos.Latent = c(0.025,0.15,0.375,0.625,0.85,0.975)[check$int]
                            )

######################################
## JAGS model
######################################
sink('scripts/JAGS/JAGS_v0.1_cens.txt')
cat("
model{
## State model
for (i in 1:N){ # N is the number of plots
  for (j in 1:Y){ # number of years
    C.S[i,j] <- z[i,j] * c.Pos[i,j] # cover including zeros
    z[i,j] ~ dbern(psi[i,j]) ## true PA state of a plot within a year depends on occupancy probability psi
    psi[i,j] ~ dunif(0,1)
    c.Pos[i,j] ~ dbeta(a.C[i,j], b.C[i,j]) T(1e-3,0.999)
    a.C[i,j] <- mu.C * tau.C
    b.C[i,j] <- (1 - mu.C) * tau.C
    } # end of years loop
  } # end of plot loop

## Derived values from state model above (average value of c.Pos, psi and C.S per year)
for (j in 1:Y){ # number of years
    cPosAn[1,j] <- mean(c.Pos[,j]) # mean across C.Pos per year, etc.
    psiAn[1,j] <- mean(psi[,j])
    cSAn[1,j] <- mean(C.S[,j])
    }

## Plot positive covers
for(k in 1:n.Plot.pos){ 
    cpos.Cens[k] ~ dinterval(cpos.Latent[k], lims[k,])
    cpos.Latent[k] ~ dbeta(a.C[plot[k], year[k]], b.C[plot[k], year[k]]) T(1e-3,0.999) # recorded cover when present follows beta distribution
  }

## Observation model for all plot visits ([within-year] detection within plots)
for (a in 1:V2){
    x[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-3, p.Dec[a]), 0.999) # trick to stop numerical problems (note that this will probably influence covars in detectability regression) -- important?
    logit(p.Dec[a]) <- gamma0 + gamma1 * yOrig[a] ## + covars on detectability, influenced by original simulated covers (before detection)
  }

## Priors!
gamma0 ~ dt(0, 0.01, 1)
gamma1 ~ dt(0, 0.01, 1)
#gamma0 ~ dnorm(0, 1)
#gamma1 ~ dnorm(0, 1)
mu.C ~ dunif(0, 1)
tau.C ~ dt(0, 0.01, 1)T(0,)

} # END MODEL
", fill = TRUE)
sink()

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_v0.1_cens.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
# Specify parameters for which posterior samples are saved
#para.names <- c('mu.C', 'tau.C', 'gamma0')
para.names <- c('mu.C', 'tau.C', 'gamma0', 'gamma1')
#para.names <- c('psi')
#mean(summary(samples)$quantiles[1:300,3]) # mean occupancy (simulated psi value)
#mean(summary(samples)$statistics[1:300,1]) # mean occupancy (simulated psi value)
# Continue the MCMC runs with sampling
samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
## Inspect results
plot(samples)
summary(samples)
gelman.diag(samples)
