################################################################################
#### 5b. Create function to make model runs across multiple species easier.#####
################################################################################
# O.L. Pescott, olipes@ceh.ac.uk
# 09.01.2019
#rm(list=ls())
######################################
list.of.packages <- c("R2jags")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(R2jags) ## obviously requires that JAGS (preferably 4.3.0) is installed


# domin things
domins <- read.csv(file = "data/dominScores.csv", header = T, stringsAsFactors = F)
#domins <- read.csv("W:/PYWELL_SHARED/Pywell Projects/BRC/_BRC_projects/NPMS/Analyses/2018 08 - Per species trend analyses/r_proj/NPMStrends/data/dominScores.csv", header = T, stringsAsFactors = F)
# t.per <- c(0.001,0.01,0.03,0.05,0.1,0.25,0.33,0.50,0.75,0.95,0.99) these are the cutpoints used in Pescott et al. 2016

## Required data 1
#N = N # the total number of plots
#Y = Y # the total number of years
#n.Plot.pos = n.Plot.pos # the total number of positive cover observations across all plot visits (= samples in NPMS database terms)
#cpos.Cens = cpos.Cens # indicator (is censored?) -- can be T/F or 1/0, doesn't matter as long as is.logical(cpos.Cens) == T
## note that for the NPMS all non-zero observations are censored; this might not be true though if we combined with CS data or similar with more "precise" cover data
#cpos.Latent = cpos.Latent # just NA values for latent observations (unknown value underlying censored observation)
#lims = lims # the limits to the intervals used for censoring (again, these should be of one type only, unless data collected under different schemes are combined) -- see Pescott et al. (2016)

##########################################################
#### Make function that can be applied across species ####
##########################################################
# Reduced list of sample data, excluding species with no data #
spForMods <- sppDatList[!names(sppDatList) %in% excludedSpp]
#
## Function for preparing a species data for JAGS model and running model
runModels_v1 <- function(i) {
  x <- spForMods[[i]]
  x$year <- format(x$date.x, "%Y")
  x <- x[order(x$year, x$plot_id), ]
  uniPlots <- unique(x$plot_id) # unique plot IDs - 31/12/2018 = 933
  # create unique plot index (useful for cross-referencing later on?)
  plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
  N <- length(uniPlots)
  # check that all years are in the range of 2015:System time
  if ( max(unique(format(x$date.x, "%Y"))) > format(Sys.time(), "%Y") ) {
    print(paste("Error. Years violation: ", names(spForMods)[i]))
  } else {
    print(paste("Years ok: ", names(spForMods)[i]))
  }
  Y <- length(unique(format(x$date.x, "%Y")))
  # All Domin data, including zeros
  yOrig <- x$dominUnify
  n.Plot.pos <- length(x$dominUnify[x$dominUnify !='0' & !is.na(x$dominUnify)]) # 300
  cpos.Cens <- rep(1, n.Plot.pos)
  cpos.Latent <- rep(NA, n.Plot.pos)
  t <- c(1e-4, 0.01,
         0.01, 0.03,
         0.03, 0.05,
         0.05, 0.1,
         0.1, 0.25,
         0.25, 0.33,
         0.33, 0.5,
         0.5, 0.75,
         0.75, 0.95, 
         0.95, 0.9999)
  tdf <- as.data.frame(matrix(t, nrow = 10, ncol = 2, byrow = TRUE))
  colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
  tdf$int <- c(1,2,3,4,5,6,7,8,9,10)
  spPos <- x[x$dominUnify !='0' & !is.na(x$dominUnify),]
  # add row number to allow resorting after merge (merge with sort = F does not actually do what we want, could also used plyer::join)
  spPos$indexPos <- 1:nrow(spPos) 
  spPos <- merge(spPos, tdf, by.x = "dominUnify", by.y = "int", all.x = T, all.y = F)
  spPos <- spPos[order(spPos$indexPos),]
  lims <- spPos[,c("L","U")] # has the lower/upper cutpoints for all intervals
  spPos <- merge(spPos, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
  spPos <- spPos[order(spPos$indexPos),]
  plot <- spPos$index
  spPos$year <- as.factor(spPos$year)
  levels(spPos$year) <- 1:length(unique(format(x$date.x, "%Y")))
  year <- spPos$year # an indicator linking a percentage cover observation to its parent year
  V2 <- nrow(x) # the total number of visits (samples), irrespective of whether there is a positive cover for a species or not
  #plotZ <- match(Achi_mill_PAN$plot_id, uniPlots)  # an indicator linking a visit to its parent (spatially unique) plot
  x <- merge(x, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
  x <- x[order(x$year, x$plot_id), ]
  plotZ <- x$index
  x$year <- as.factor(x$year)
  levels(x$year) <- 1:length(unique(format(x$date.x, "%Y")))
  yearZ <- x$year # an indicator linking a visit to its parent year
  y <- y1 <-  numeric()
  for (j in 1:nrow(x)) {
    if(is.na(x$dominUnify[j])){ 
      y <- NA
          } else if (x$dominUnify[j]==0) { 
      y <- 0
          } else {
      y <- 1
          }
    y1 <- c(y1,y)
  }
  y <- y1 # visit-level detection history (binary)
  # Prepare data for JAGS
  Data <- list(N = N,
               Y = Y,
               n.Plot.pos = n.Plot.pos,
               cpos.Cens = cpos.Cens, # indicator (is censored?)
               cpos.Latent = cpos.Latent, # NA values for latent observations
               lims = lims,
               plot = plot,
               year = year,
               V2 = V2,
               plotZ = plotZ,
               yearZ = yearZ,
               y = y,
               yOrig = yOrig)
  zinit <- matrix(1, nrow = N, ncol = Y)
  inits.fn <- function() list(z = zinit,
                              tau.C = runif(1,1,5),
                              mu.C = 0.5,
                              gamma0 = rnorm(1,0,1),
                              gamma1 = rnorm(1,0,1),
                              # cpos.Latent is approx. mid-points of the categories, used as initial values (dominUnif)
                              # intervals start from 1 (midpoint for the "zeroth" category not needed for these latent values for positive data)
                              cpos.Latent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify] )
  #for ref only
  cPos.Init <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify]
  jagsModel <- rjags::jags.model(file= 'scripts/JAGS/JAGS_v0.3_NPMS.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)
  # Specify parameters for which posterior samples are saved
  para.names <- c('mu.C', 'tau.C', 'gamma0', 'gamma1', 'cPosAn', 'psiAn', 'cSAn')
  # Continue the MCMC runs with sampling
  samples <- rjags::coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
  ## Inspect results
  out <- summary(samples)
  return(out)
}

# this will be applied across the list created in 4_extractData.R (sppDatList), minus excluded species
sppModels <- list()
#sppModels <- lapply(seq_along(spForMods[1]), function(i) runModels_v1(i)) # test
sppModels <- lapply(seq_along(spForMods[1:5]), function(i) runModels_v1(i))
names(sppModels) <- names(spForMods[1:5])


#### END OF FUNCTIONS SECTION. BELOW WAS SINGLE SPECIES DEVELOPMENT WORK ####
#######################################################################################################################################################
##############################################################
#### Single example of Achillea millefolium in grasslands ####
############################################################## #### collapse ####
#load(file = "data/Achi_mille_grassSamples_20180920.Rdata")
#load("W:/PYWELL_SHARED/Pywell Projects/BRC/_BRC_projects/NPMS/Analyses/2018 08 - Per species trend analyses/r_proj/NPMStrends/data/Achi_mille_grassSamples_20180920.Rdata")
#head(Achi_mill_PAN); tail(Achi_mill_PAN); unique(Achi_mill_PAN$dominUnify)
## 1. Achillea millefolium in grassland samples
# not sure it really matters, but sorting the data by year and then plot first may simplify downstream things
Achi_mill_PAN$year <- format(Achi_mill_PAN$date.x, "%Y")
Achi_mill_PAN <- Achi_mill_PAN[order(Achi_mill_PAN$year, Achi_mill_PAN$plot_id), ]
uniPlots <- unique(Achi_mill_PAN$plot_id) # unique plot IDs - 31/12/2018 = 933
# create unique plot index (useful for cross-referencing later on?)
plotIndex <- data.frame(plot = uniPlots, index = 1:length(uniPlots))
N <- length(uniPlots)
# check that all years are in the range of 2015:System time
if ( max(unique(format(Achi_mill_PAN$date.x, "%Y"))) > format(Sys.time(), "%Y") ) {
    print ("error")
} else {
    print("OK")
}
# 
Y <- length(unique(format(Achi_mill_PAN$date.x, "%Y")))
# All Domin data, including zeros
yOrig <- Achi_mill_PAN$dominUnify

n.Plot.pos <- length(Achi_mill_PAN$dominUnify[Achi_mill_PAN$dominUnify !='0' & !is.na(Achi_mill_PAN$dominUnify)]) # 300
## checks on n.Plot.pos
head(Achi_mill_PAN[!is.na(Achi_mill_PAN$dominUnify) & Achi_mill_PAN$dominUnify !='0',]) # data frame
head(Achi_mill_PAN$dominUnify[Achi_mill_PAN$dominUnify !='0' & !is.na(Achi_mill_PAN$dominUnify)]) # vector
## looks OK

cpos.Cens <- rep(1, n.Plot.pos)
cpos.Latent <- rep(NA, n.Plot.pos)
t <- c(1e-4, 0.01,
       0.01, 0.03,
       0.03, 0.05,
       0.05, 0.1,
       0.1, 0.25,
       0.25, 0.33,
       0.33, 0.5,
       0.5, 0.75,
       0.75, 0.95, 
       0.95, 0.9999)
t = matrix(t, nrow = 10, ncol = 2, byrow = TRUE)
tdf <- as.data.frame(t)
colnames(tdf) <- c('L','U') # 'L'ower and 'U'pper bounds of categories
tdf$int <- c(1,2,3,4,5,6,7,8,9,10)

# positive rows only
spPos <- Achi_mill_PAN[Achi_mill_PAN$dominUnify !='0' & !is.na(Achi_mill_PAN$dominUnify),]
# add row number to allow resorting after merge (merge with sort = F does not actually do what we want, could also used plyer::join)
spPos$indexPos <- 1:nrow(spPos) 
spPos <- merge(spPos, tdf, by.x = "dominUnify", by.y = "int", all.x = T, all.y = F)
spPos <- spPos[order(spPos$indexPos),]
lims <- spPos[,c("L","U")] # has the lower/upper cutpoints for all intervals

## Required data 2
#plot = plot # an indicator linking a percentage cover observation to its parent (spatially unique) plot
#year = year # an indicator linking a percentage cover observation to its parent year
#V2 = V2 # the total number if visits (samples), irrespective of whether there is a postive cover for a species or not
#plotZ = plotZ # an indicator linking a visit to its parent (spatially unique) plot
#yearZ = yearZ # an indicator linking a visit to its parent year
#x = x # visit-level detection history (binary)

#plot <- match(spPos$plot_id, Achi_mill_PAN$plot_id) # an indicator linking a percentage cover observation to its parent (spatially unique) plot
spPos <- merge(spPos, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
spPos <- spPos[order(spPos$indexPos),]
plot <- spPos$index
#yDat <- data.frame(index = 1:length(unique(format(Achi_mill_PAN$date.x, "%Y"))), year = unique(format(Achi_mill_PAN$date.x, "%Y"))) # not really needed
spPos$year <- as.factor(spPos$year)
levels(spPos$year) <- 1:length(unique(format(Achi_mill_PAN$date.x, "%Y")))
year <- spPos$year # an indicator linking a percentage cover observation to its parent year
## All visit information
V2 <- nrow(Achi_mill_PAN) # the total number of visits (samples), irrespective of whether there is a positive cover for a species or not
#plotZ <- match(Achi_mill_PAN$plot_id, uniPlots)  # an indicator linking a visit to its parent (spatially unique) plot
Achi_mill_PAN <- merge(Achi_mill_PAN, plotIndex, by.x = "plot_id", by.y = "plot", all.x = T, all.y = F)
Achi_mill_PAN <- Achi_mill_PAN[order(Achi_mill_PAN$year, Achi_mill_PAN$plot_id), ]
plotZ <- Achi_mill_PAN$index
Achi_mill_PAN$year <- as.factor(Achi_mill_PAN$year)
levels(Achi_mill_PAN$year) <- 1:length(unique(format(Achi_mill_PAN$date.x, "%Y")))
yearZ <- Achi_mill_PAN$year # an indicator linking a visit to its parent year
x <- x1 <-  numeric()
for (i in 1:nrow(Achi_mill_PAN)) {
  if(is.na(Achi_mill_PAN$dominUnify[i])){ 
    x <- NA
    #print(c(i,x)) # this was just for catching a bug
  } else if (Achi_mill_PAN$dominUnify[i]==0) { 
    x <- 0
    #print(c(i,x))
  } else {
    x <- 1
    #print(c(i,x))
  }
  x1 <- c(x1,x)
}
x <- x1 # visit-level detection history (binary)

####################
## Send data to JAGS
####################
# Data list for passing to JAGS
Data <- list(N = N,
             Y = Y,
             n.Plot.pos = n.Plot.pos,
             cpos.Cens = cpos.Cens, # indicator (is censored?)
             cpos.Latent = cpos.Latent, # NA values for latent observations
             lims = lims,
             plot = plot,
             year = year,
             V2 = V2,
             plotZ = plotZ,
             yearZ = yearZ,
             x = x,
             yOrig = yOrig)
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
                            tau.C = runif(1,1,5),
                            mu.C = 0.5,
                            gamma0 = rnorm(1,0,1),
                            gamma1 = rnorm(1,0,1),
                            # cpos.Latent is approx. mid-points of the categories, used as initial values (dominUnif)
                            # intervals start from 1 (midpoint for the "zeroth" category not needed for these latent values for positive data)
                            cpos.Latent = c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify]
)
#for ref only
cPos.Init <- c(0.001,0.025,0.04,0.075,0.175,0.29,0.375,0.625,0.85,0.975)[spPos$dominUnify]

######################################
## JAGS model
######################################
sink('scripts/JAGS/JAGS_v0.3_NPMS.txt')
cat("
    model{
    ## State model
    for (i in 1:N){ # N is the number of plots
      for (j in 1:Y){ # number of years
        C.S[i,j] <- z[i,j] * c.Pos[i,j] # cover including zeros
        z[i,j] ~ dbern(psi[i,j]) # true PA state of a plot within a year depends on occupancy probability psi
        psi[i,j] ~ dunif(0,1)
        c.Pos[i,j] ~ dbeta(a.C[i,j], b.C[i,j]) T(1e-4,0.9999)
        a.C[i,j] <- mu.C * tau.C
        b.C[i,j] <- (1 - mu.C) * tau.C
      } # end of years loop
    } # end of plot loop

    ## Derived values from state model above (average value of c.Pos, psi and C.S per year)
    for (j in 1:Y){ # number of years
      cPosAn[1,j] <- mean(c.Pos[1:N,j]) # mean across C.Pos per year, etc.
      psiAn[1,j] <- mean(psi[1:N,j])
      cSAn[1,j] <- mean(C.S[1:N,j])
    }

    ## Plot positive covers
    for(k in 1:n.Plot.pos){ 
      cpos.Cens[k] ~ dinterval(cpos.Latent[k], lims[k,])
      cpos.Latent[k] ~ dbeta(a.C[plot[k], year[k]], b.C[plot[k], year[k]]) T(1e-4,0.9999) # recorded cover when present follows beta distribution
    }
    
    ## Observation model for all plot visits ([within-year] detection within plots)
    for (a in 1:V2){
      y[a] ~ dbern(py[a]) # detectability influences detection
      py[a] <- z[plotZ[a], yearZ[a]] * p.dec[a] # true state x detectability
      p.dec[a] <- min(max(1e-4, p.Dec[a]), 0.999) # trick to stop numerical problems (note that this will probably influence covars in detectability regression) -- important?
      ## Can add observed cover (Domin scale) as covar to the following line
      logit(p.Dec[a]) <- gamma0 + gamma1 * yOrig[a] ## yOrig is reported Domin value (which can be unknown, i.e. 'NA')
      #yOrig[a] ~ dt(0, 0.01, 1) # Prior necessary to account for missing data (not sure this is the best distrbution though)
      yOrig[a] ~ dunif(0,10) # Doesn't seem to make much difference to results
    }
    
    ## Priors!
    gamma0 ~ dt(0, 0.01, 1)
    gamma1 ~ dt(0, 0.01, 1)
    mu.C ~ dunif(0, 1)
    tau.C ~ dt(0, 0.01, 1)T(0,)
    
    } # END MODEL
    ", fill = TRUE)
sink()

jagsModel <- jags.model(file= 'scripts/JAGS/JAGS_v0.3_NPMS.txt', data = Data, inits = inits.fn, n.chains = 3, n.adapt= 500)

# Specify parameters for which posterior samples are saved
para.names <- c('mu.C', 'tau.C', 'gamma0', 'gamma1', 'cPosAn', 'psiAn', 'cSAn')
# Continue the MCMC runs with sampling
samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 500)
## Inspect results
summary(samples)
gelman.diag(samples)
plot(samples)