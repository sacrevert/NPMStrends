## X2d. Simulation of the types of data that we are going to be modelling, using JAGS initially
# Now with interval censoring! (2c)
# Version 2d -- used for sims
# Version 2e -- fixed error in occupancy simulation and summary,and change to using y > yOrig in detection model (yOrig not available in real world)
# Version 2f -- simplified notation as per Freeman comments
#
# note that this model only uses an observation/detection model for the repeat visits within a year, not for the cover data
# whilst it is reasonable to assume that the occupancy state is stable within years (obviously this is not 100% true all the time, but is reasonable)
# it is not reasonable to assume that the repeat visits to a plot provide replicated information about a stable cover for a species (as in Wright et al. 2017,
# where there are repeat observations, i.e. Quality Assurance, within a single visit), for this reason the cover values are taken directly as accurate 
# observations and are included in the "state model" as such
#
# O.L. Pescott
# 09.01.2019
rm(list=ls())

######################################
library(R2jags)
library(mcmcplots)
#library(BayesianTools)
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
## Could put this in a simulation function (see Wright et al. 2017)
N <- 50 # number of spatially unique plots
J <- 5 # number of visits to a plot within a year (assume constant for the moment)
#psi <- 0.5 # ovreall true plot occupancy (average)
psi <- 0.25# should be easier for model to retrieve true values of gamma0 and gamma1 if occupancy is 1 -- this appears to be true
Y <- 5 # total number of years of monitoring covered by the data
mu0 <- 0.025       # parameter for mean of cover beta distribution # 0.25
tau0 <-3    # parameter for 'precision' of cover distribution # 3
gamma0 <- -2  # intercept for detection logistic regression # -1.5
gamma1 <- 3   # slope for detection with %cover # 2
# e.g. plogis(-2 + 3*cover) makes for greater detectability range based on percent covers
# if you do this but only estimate gamma0 in the model, then biases in gamma0 and mu.C estimates increase

# array of plot covers per visit per year
y.array <- array(dim = c(N, J, Y)) # plots, visits, years
for(k in 1:Y){
      y.array[,,k] <- matrix(rbeta(N*J, mu0*tau0, (1-mu0)*tau0), 
                                 nrow=N, ncol=J)
    for (i in 1:N){
      y.array[i,,k] <- rep(rbinom(1,1,psi), J) * y.array[i,,k] # plot level occupancy within years
    }
}
# make a detection history matrix based on cover data
x.array <- array(dim = c(N, J, Y))
for(k in 1:Y){
  for(i in 1:N){
    for(j in 1:J){
      x.array[i, j, k] <- ifelse(y.array[i, j, k] > 0,
                                 rbinom(1, 1, 
                                        plogis(gamma0 + gamma1*y.array[i, j, k])), # standard detection function with link between cover and detection
                                        #plogis(gamma0)) # intercept only
                                 #1), # for perfect detection
                                 0)
    }
  }
}
y.arrayOrig <- y.array # keep original covers
y.array[which(x.array==0)] <- 0 # if the plant was not actually detected, then set the recorded cover to zero as well
y.array[which(y.array<1e-16 & y.array > 0)] <- 1e-16 # put hard limit on lower cover bound
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
t <- c(1e-16,0.05, # need to tweak lower categories if starting with low percentage cover like 0.01
       0.05,0.25,
       0.25,0.5,
       0.5,0.75,
       0.75,0.95,
       0.95,0.9999999999999999)
t = matrix(t, nrow = 6, ncol = 2, byrow = TRUE)
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
intervals <- c(1,2,3,4,5,6)
tdf$int <- intervals
# library(plyr) for SQL join function (like merge but keeps order) -- although actually could use merge with sort = F
tInt <- c(1e-16,0.05,0.25,0.5,0.75,0.95,0.9999999999999999) # need to tweak lower categories if starting with low percentage cover like 0.01
int <- findInterval(cpos[,1], tInt) # find corresponding interval for all data points
int <- as.data.frame(int)
m1 <- plyr::join(int, tdf, by = "int", type = "left", match = "all") # change to using merge at some point (although merge resorts, need to add row index to resort on)
lims <- m1[,2:3] # has the intervals for all points
check <- cbind(cpos, m1); head(check); tail(check) # just for quick visual check that all is well
cposCens <- rep(1, nrow(cpos))
cposLatent <- cpos[,3] # just NAs

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
V3 <- N*Y
# indicator linking every visit x to plot; length = V2 
plotZ <- rep(1:N, J*Y)
# indicator linking every visit x to year; length V2 
yearZ <- rep(1:Y, each = N*J)
x <- as.vector(x.array) # detection indicator for every visit (length V2)
# covers for all visits
#y <- as.vector(y.arrayOrig) # original covers used in detectability loop (inc. zeros) for every visit (length V2)
y <- as.vector(y.array) # covers after non-detections -- shouldn't use previous line, as this information is not available in reality

# Data list for passing to JAGS
Data <- list(N = N,
            Y = Y,
            n.Plot.pos = n.Plot.pos,
            #cpos.Cens = cpos[,2], # indicator (is censored?)
            cposCens = cposCens, # indicator (is censored?)
            cposLatent = cposLatent, # NA values for latent observations
            lims = lims,
            plot = plot,
            year = year,
            V2 = V2,
            #V3 = V3,
            plotZ = plotZ,
            yearZ = yearZ,
            x = x,
            y = y)
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
                            #mean.p = runif(1, 0.01,0.99),
                            mean.p = runif(1,0,1),
                            tau = runif(1,1,10),
                            mu = runif(1,0,1),
                            #gamma0 = rnorm(1,0,1),
                            gamma1 = runif(1,-5,5),
                            cposLatent = c(0.025,0.15,0.375,0.625,0.85,0.975)[check$int]
                            )

######################################
## JAGS model
######################################
sink('scripts/JAGS/JAGS_vX2f_cens.txt')
cat("
model{
## State model
for (i in 1:N){ # N is the number of plots
  for (j in 1:Y){ # number of years
    C[i,j] <- z[i,j] * cPos[i,j] # cover including zeros
    z[i,j] ~ dbern(psi[i,j]) # true PA state of a plot within a year depends on occupancy probability psi
    psi[i,j] ~ dunif(0,1)
    cPos[i,j] ~ dbeta(a[i,j], b[i,j]) T(1e-16,0.9999999999999999)
    a[i,j] <- mu * tau
    b[i,j] <- (1 - mu) * tau
    } # end of years loop
} # end of plot loop


## Derived values from state model
for (j in 1:Y){ # number of years
    annOcc[j] <- (sum(z[,j]))/N # annual occupancy
}
avgOcc <- mean(annOcc[]) # average annual occupancy

## Plot positive covers
for(k in 1:n.Plot.pos){ 
    cposCens[k] ~ dinterval(cposLatent[k], lims[k,])
    cposLatent[k] ~ dbeta(a[plot[k], year[k]], b[plot[k], year[k]]) T(1e-16,0.9999999999999999) # recorded cover when present follows beta distribution
}

## Observation model for all plot visits ([within-year] detection within plots)
mY <- mean(y[]) # mean detected cover across all years for mean-centring
for (a in 1:V2){
    x[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-16, p.Dec[a]), 0.9999999999999999) # trick to stop numerical problems
    logit(p.Dec[a]) <- gamma0c + gamma1 * (y[a] - mY) # centering reduces gamma0c/gamma1 correlation
    #logit(p.Dec[a]) <- gamma0 + gamma1 * y[a]
    #logit(p.Dec[a]) <- gamma0 # without dependency of detection on abundance (requires changes to L236 and prior as well)
}
# only required if you are including detected cover in estimate of detectability
gamma0 <- gamma0c - gamma1 * mY # Recover original intercept

## Priors!
#for (j in 1:Y){
#mean.p[j] ~ dbeta(2,0.5) # regularising prior with emphasis on detection being high
#gamma0[j] <- logit(mean.p[j])
#}
mean.p ~ dbeta(1,1) # broad intercept on prob scale
gamma0c <- logit(mean.p) # transformed # note that this requires tweak to initial values
#gamma0 <- logit(mean.p) # transformed # note that this requires tweak to initial values
#gamma1 ~ dunif(-20,20) # broad uniform on logit scale
gamma1 ~ dnorm(0, 1) #tau = 1 # waekly informative prior to help with quasi-complete separation in logistic regression
mu ~ dunif(0, 1)
tau ~ dt(0, 0.01, 1)T(0,)

} ## END MODEL
", fill = TRUE)
sink()

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_vX2f_cens.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
para.names <- c('gamma1', 'gamma0', 'mu', 'tau', 'avgOcc')
samples9 <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
## Inspect results
out <- summary(samples9)
mu_sd <- out$stat[,1:2] #make columns for mean and sd
q <- out$quantile[,c(3,1,5)] #make columns for median and CI
tableOut_m9 <- as.data.frame(cbind(mu_sd,q)) #make table
tableOut_m9
#save(tableOut_m9, samples9, file = "outputs/model9.Rdata")
#tableOut[grep(rownames(tableOut), pattern = 'gamma'),]
#tableOut[grep(rownames(tableOut), pattern = 'annOcc'),]
#write.csv(tableOut, file = "outputs/tests/test1_CS.csv")

plot(samples)
gelman.diag(samples9)
## Using mcmcplots::mcmcplot
mcmcplot(samples, random = 5)
caterplot(samples, style = "plain")
samplesDF <- do.call(rbind.data.frame, samples)
cor(samplesDF)
pairs(samplesDF)
