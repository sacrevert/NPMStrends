## X2c. Simulation of the types of data that we are going to be modelling, using JAGS initially, with interval censoring
# We also need to bring in the fact that some plots have no information (NAs, i.e. some plots are visited because of the habitat/level we never
# have any certain information for them)
#
# NAs across actual abundances and corresponding occupancies don't make any difference,
# nor does it matter if one is missing and not the other
# nor does it matter if one visit or year information is missing
# only important thing is if information is contradictory (e.g. detected by zero abundance)
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
# Simulations originally based on those of Wright et al. 2017
######################################
## Could put this in a simulation function (see Wright et al. 2017)
N <- 100 # number of spatially unique plots
J <- 5 # number of visits to a plot within a year (assume constant for the moment)
psi <- 0.2 # true occupancy across plots (average)
Y <- 3 # total number of years of monitoring covered by the simulated data
mu <- 0.9      # parameter for mean of cover beta distribution (i.e. mean conditional on presence) # 0.25
phi <- 10      # parameter for 'precision' of cover distribution (i.e. variance conditional on presence) # 3
gamma0 <- -1   # intercept for detection logistic regression # -1.5
gamma1 <- 1   # slope for detection with %cover # 2

# array of plot covers per visit per year
# array dimensions = "rows, columns, number of 'list' items" 
y.array <- array(dim = c(N, J, Y)) # equivalent to a Y item list, with each item having J columns and N rows
for(k in 1:Y){
      y.array[,,k] <- matrix(rbinom(N*J, 1, psi)*rbeta(N*J, mu*phi, (1-mu)*phi), 
                                 nrow=N, ncol=J) # generate 'truth'
}
# make a detection history matrix based on cover data
x.array <- array(dim = c(N, J, Y))
for(k in 1:Y){ # loop through every element of y.array
  for(i in 1:N){
    for(j in 1:J){
      x.array[i, j, k] <- ifelse(y.array[i, j, k] > 0,
                                 rbinom(1, 1, 
                                        #plogis(gamma0 + y.array[i, j, k])), # detection function 1
                                        plogis(gamma0 + gamma1*y.array[i, j, k])), # detection function 2
                                 #1, # for testing
                                 0)
    }
  }
}
y.arrayOrig <- y.array # keep original covers
y.array[which(x.array==0)] <- 0 # if the plant was not actually detected, then set the recorded cover to zero as well

## introduce some plots that are always NA (this mirrors in the real world situation for the NPMS where there are some plots that might end up being NA
# at least in the short term)
# replace the last plot (N = 100) with NAs in all reps and years
#y.array[(N-10):N, 1:J, 1:Y] <- NA
#x.array[(N-10):N, 1:J, 1:Y] <- NA
y.array[(N-10):N, J, 1:Y] <- NA
x.array[(N-10):N, J, 1:Y] <- NA

###################################### END OF DATA SIM ####

##################################################
## Data/indicators required for running JAGS model
##################################################
# total number of visits with positive covers
cpos <- matrix(NA, nrow = length(as.vector(y.array)[which(as.vector(y.array) > 0)]), ncol = 3) # empty matrix for actual covers and cover classes
cpos[,1] <- as.vector(y.array)[which(as.vector(y.array) > 0)] # actual cover value for every positive cover
cpos[,2] <- 1 # indicator stating the observation censored
cpos[,3] <- NA # NA values for latent variable
# code borrowed from Pescott et al. 2016 Power paper:
t <- c(1e-16,0.05, # boundaries of Domin scale used in NPMS
       0.05,0.25,
       0.25,0.5,
       0.5,0.75,
       0.75,0.95,
       0.95,0.9999999999999999)
t = matrix(t, nrow = 6, ncol = 2, byrow = TRUE) # turn into a matrix and then data frame
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
tdf$int <- c(1,2,3,4,5,6) # interval labels

# library(plyr) for SQL join function (like merge but keeps order) -- although actually could use merge with sort = F
tInt <- c(1e-16,0.05,0.25,0.5,0.75,0.95,0.9999999999999999)
int <- findInterval(cpos[,1], tInt) # find corresponding interval for all data points (base::findInterval)
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
    out <- ifelse(y.array[i, j, k] > 0, i, NA) # i is plot index
    out1 <- c(out1, out)
    }
  }
}
plot <- out1[!is.na(out1)]

out <- out1 <-  numeric()
for(k in 1:Y){
  for(j in 1:J){
    for(i in 1:N){
    out <- ifelse(y.array[i, j, k] > 0, k, NA) # k is year index
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
# indicator linking every visit x to year; length = V2 
yearZ <- rep(1:Y, each = N*J)
x <- as.vector(x.array) # detection indicator for every visit (length V2)
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
# Some other parameters are also fixed to avoid extremely small likelihoods
# but these can probably be relaxed a bit more.
zinit <- matrix(1, nrow = N, ncol = Y)
inits.fn <- function() list(z = zinit,
                            tau.C = runif(1,1,5),
                            mu.C = 0.5,
                            gamma0 = rnorm(1,0,1),
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
    c.Pos[i,j] ~ dbeta(a.C[i,j], b.C[i,j]) T(1e-16,0.9999999999999999)
    a.C[i,j] <- mu.C * tau.C
    b.C[i,j] <- (1 - mu.C) * tau.C
    } # end of years loop
  } # end of plot loop

## Derived values from state model above (average value of c.Pos, psi and C.S per year)
for (j in 1:Y){ # number of years
    cPosAn[j] <- mean(c.Pos[,j]) # mean across C.Pos per year, etc.
    psiAn[j] <- mean(psi[,j])
    cSAn[j] <- mean(C.S[,j])
    annOcc[j] <- (sum(z[,j]))/N
  }

## Plot positive covers
for(k in 1:n.Plot.pos){ 
    cpos.Cens[k] ~ dinterval(cpos.Latent[k], lims[k,])
    cpos.Latent[k] ~ dbeta(a.C[plot[k], year[k]], b.C[plot[k], year[k]]) T(1e-16,0.9999999999999999) # recorded cover when present follows beta distribution
  }

## Observation model for all plot visits ([within-year] detection within plots)
for (a in 1:V2){
    x[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-16, p.Dec[a]), 0.9999999999999999) # trick to stop numerical problems (note that this will probably influence covars in detectability regression) -- important?
    logit(p.Dec[a]) <- gamma0# + gamma1 * yOrig[a] ## + covars on detectability, influenced by original simulated covers (before detection)
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
para.names <- c('psi', 'gamma0', 'psiAn','annOcc')
# Continue the MCMC runs with sampling
samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
mean(summary(samples)$quantiles[1:100,3]) # mean occupancy (simulated psi value)
mean(summary(samples)$statistics[1:100,1]) # mean occupancy (simulated psi value)
## Inspect results
plot(samples)
summary(samples)
gelman.diag(samples)