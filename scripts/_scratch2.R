## X2d. Simulation of the types of data that we are going to be modelling, using JAGS initially
# Now with interval censoring! (2c)
# Version 2d -- investigation as to the year effects (e.g. plotting)
# currently there are no changes to the Bayesian model (as compared to versions X2b and X2c)
#
# note that this model only uses an observation/detection model for the repeat visits within a year, not for the cover data
# whilst it is reasonable to assume that the occupancy state is stable within years (obviously this is not 100% true all the time, but is reasonable)
# it is not reasonable to assume that the repeat visits to a plot provide replicated information about a stable cover for a species (as in Wright et al. 2017,
# where there are repeat observations, i.e. Quality Assurance, within a single visit), for this reason the cover values are taken directly as accurate 
# observations and are included in the "state model" as such
#
# O.L. Pescott
# 09.01.2019
#rm(list=ls())

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
# Simulations based on those of Wright et al. 2017
######################################
## Good estimates of gamma0, gamma1, mu.C and tau.C with N = 100, J = 10, psi = 1, Y = 10, mu = 0.5, phi = 10, gamma0 = -2, gamma1 = 3
#         Mean       SD  Naive SE Time-series SE
#gamma0 -1.980 0.069046 0.0017827      0.0085523
#gamma1  2.969 0.129425 0.0033417      0.0167760
#mu.C    0.536 0.002886 0.0000745      0.0001088
#tau.C  10.216 0.292774 0.0075594      0.0176602

#          2.5%     25%    50%     75%   97.5%
#gamma0 -2.1147 -2.0277 -1.982 -1.9311 -1.8434
#gamma1  2.7112  2.8822  2.974  3.0578  3.2175
#mu.C    0.5303  0.5339  0.536  0.5381  0.5414
#tau.C   9.6300 10.0310 10.221 10.3986 10.7998

## Could put this in a simulation function (see Wright et al. 2017)
N <- 100 # number of spatially unique plots
J <- 2 # number of visits to a plot within a year (assume constant for the moment)
#psi <- 0.5 # true occupancy (average)
psi <- 0.25# should be easier for model to retrieve true values of gamma0 and gamma1 if occupancy is 1 -- this appears to be true
Y <- 5 # total number of years of monitoring covered by the data
mu <- 0.025       # parameter for mean of cover beta distribution # 0.25
phi <- 3     # parameter for 'precision' of cover distribution # 3
gamma0 <- -2  # intercept for detection logistic regression # -1.5
gamma1 <- 3   # slope for detection with %cover # 2
# e.g. plogis(-2 + 3*cover) makes for greater detectability range based on percent covers
# if you do this but only estimate gamma0 in the model, then biases in gamma0 and mu.C estimates increase

# array of plot covers per visit per year
y.array <- array(dim = c(N, J, Y)) # plots, vists, years
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
                                 #1), # for testing: if you remove probabilistic element, then mu.C better estimated
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
t <- c(1e-16,0.05,
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
tInt <- c(1e-16,0.05,0.25,0.5,0.75,0.95,0.9999999999999999)
int <- findInterval(cpos[,1], tInt) # find corresponding interval for all data points
int <- as.data.frame(int)
m1 <- plyr::join(int, tdf, by = "int", type = "left", match = "all") # change to using merge at some point (although merge resorts, need to add row index to resort on)
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
plotZ <- rep(1:N, J*Y) # i, j, k =  # plots, vists, years
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
                            #mean.p = runif(1, 0.01,0.99),
                            mean.p = runif(1, 0,1),
                            tau.C = runif(1,1,10),
                            mu.C = runif(1, 0,1),
                            #gamma0 = rnorm(1,0,1),
                            gamma1 = runif(1,-5,5),
                            cpos.Latent = c(0.025,0.15,0.375,0.625,0.85,0.975)[check$int]
)

######################################
## JAGS model
######################################
sink('C:\\Users\\olipes\\Desktop\\JAGS_vX2d_cens.txt')
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
    #cPosAn[j] <- mean(c.Pos[,j]) # mean across C.Pos per year, etc.
    #psiAn[j] <- mean(psi[,j])
    #cSAn[j] <- mean(C.S[,j])
    annOcc[j] <- (sum(z[,j]))/N # avg occupancy proportion
    }
    #mCS <- mean(C.S[,]) #global C.S. mean
    #lmCS <- logit(mean(C.S[,]))
    
    ## Plot positive covers
    for(k in 1:n.Plot.pos){ 
    cpos.Cens[k] ~ dinterval(cpos.Latent[k], lims[k,])
    cpos.Latent[k] ~ dbeta(a.C[plot[k], year[k]], b.C[plot[k], year[k]]) T(1e-16,0.9999999999999999) # recorded cover when present follows beta distribution
    }
    
    ## Observation model for all plot visits ([within-year] detection within plots)
    for (a in 1:V2){
    x[a] ~ dbern(py[a]) # detectability influences detection
    py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
    p.dec[a] <- min(max(1e-16, p.Dec[a]), 0.9999999999999999) # trick to stop numerical problems (note that this *could* influence covars in detectability regression) -- important?
    #logit(p.Dec[a]) <- gamma0 + gamma1 * (logit(C.S[plotZ[a], yearZ[a]]))-lmCS ##    
    #logit(p.Dec[a]) <- gamma0c + gamma1 * ((C.S[plotZ[a], yearZ[a]])-mCS) ## centering reduces gamma0c/gamma1 correlation
    logit(p.Dec[a]) <- gamma0 + gamma1 * yOrig[a] ## + covars on detectability, influenced by original simulated covers (before detection)
    yOrig[a] ~ dunif(0,10) # Doesn't seem too make much difference to results
    }
    #gamma0 <- gamma0c - gamma1 * mCS # Recover original intercept
    
    ## Priors!
    mean.p ~ dbeta(1,1) # broad intercept on prob scale (but not going quite to +/- infinity!)
    gamma0 <- logit(mean.p) # transformed # note that this requires tweak to initial values
    gamma1 ~ dunif(-5,5) # broad uniform on logit scale (but not too broad)
    #gamma1 ~ dt(0, 0.04, 1)
    #gamma1 ~ dt(0, 0.01, 1)
    #gamma0 ~ dnorm(0, 1)
    #gamma1 ~ dnorm(0, 1)
    mu.C ~ dunif(0, 1)
    tau.C ~ dt(0, 0.01, 1)T(0,)
    #tauC.T <- pow(tau.C, -0.5)
    
    } # END MODEL
    ", fill = TRUE)
sink()

jagsModel <- jags.model(file='C:\\Users\\olipes\\Desktop\\JAGS_vX2d_cens.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
# Specify parameters for which posterior samples are saved
#para.names <- c('mu.C', 'tau.C', 'gamma0')
#para.names <- c('annOcc', 'cPosAn', 'psiAn', 'cSAn', 'gamma1', 'gamma0', 'mu.C', 'tau.C', 'tauC.T')
para.names <- c('gamma1', 'gamma0','mu.C', 'tau.C', 'annOcc')
#para.names <- c('lmCS')
#para.names <- c('psi')
#para.names <- c('z')
#para.names <- c('py')
#mean(summary(samples)$quantiles[1:300,3]) # mean occupancy (simulated psi value)
#mean(summary(samples)$statistics[1:300,1]) # mean occupancy (simulated psi value)
# Continue the MCMC runs with sampling
samples_m9 <- coda.samples(jagsModel, variable.names = para.names, n.iter = 1000)
#samples <- update(jagsModel, n.iter = 3000)
## Inspect results
#summary(samples)
## Inspect results
out <- summary(samples_m9)
mu_sd <- out$stat[,1:2] #make columns for mean and sd
q <- out$quantile[,c(3,1,5)] #make columns for median and CI
tableOut_m9 <- as.data.frame(cbind(mu_sd,q)) #make table
save(samples_m9, tableOut_m9, file ="W:\\PYWELL_SHARED\\Pywell Projects\\BRC\\_BRC_projects\\NPMS\\Analyses\\2018 08 - Per species trend analyses\\r_proj\\NPMStrends\\outputs\\model9.Rdata")
#tableOut[grep(rownames(tableOut), pattern = 'gamma'),]
#tableOut[grep(rownames(tableOut), pattern = 'annOcc'),]
#write.csv(tableOut, file = "outputs/tests/test1_CS.csv")

plot(samples)
gelman.diag(samples)

## Using mcmcplots::mcmcplot
mcmcplot(samples)
caterplot(samples, style = "plain")
samplesDF <- do.call(rbind.data.frame, samples)
cor(samplesDF)
pairs(samplesDF)
# just for gamma0 and gamma1
