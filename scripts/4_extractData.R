## 4. Run extraction and processing for broad habitat X and species Y
## prepare Domin data for modelling
# O.L. Pescott
# 21.08.2018, updated 10.01.2019
#rm(list=ls())

## Source file with datasets, getSamples() and spSamplePA() functions
## Might be slow, as script 3_process... itself sources script 2_getDataFrom..., which queries the Indicia database
##############################################################
################ Requires database password!!! ###############
############### or you can load existing data ################
############### from: file = "data/*"  #######################
######## See script 3_processDataFuns.R lines 10 and 11 ######
######## choose latest data files and run from there #########
##############################################################
source(file = "scripts/3_processDataFuns.R")

##############################################################
#### Single example of Achillea millefolium in grasslands ####
############################################################## #### collapse ####
#grasslands <- c("Neutral pastures and meadows", "Dry acid grassland", "Dry calcareous grassland", "Neutral damp grassland", "Lowland grassland")
#grassSamples <- getSamples(habsList = grasslands)
#Achi_mill_PAN <- spSamplePA(samples = grassSamples, species = "Achillea millefolium")
#head(Achi_mill_PAN); tail(Achi_mill_PAN)
## Save dataset (20.09.2018)
#save(Achi_mill_PAN, file = "data/Achi_mille_grassSamples_20180920.Rdata")
##

## Following can be used if database extraction not run
load(file = "data/npms1518_SamplesSpecies_2019-01-10.Rdata")
load(file = "data/npms1518_PlotsSamples_2019-01-10.Rdata")

##############################################################
############ Generalised approach using functions ############
##############################################################
## General habitats list
habs <- read.csv( file = "data/NPMShabitats_Surv_to_Broad.csv", h = T, strings = F)
# bsh = broad scale habitat
getSamps2 <- function(bsh) { getSamples(habsList = unique(habs[habs$NPMS.broad.habitat==bsh,2]))}
habSamps <- getSamps2(bsh = "Lowland grassland") # includes relevant broad and fine scale habitat names

## List of positive indicator species per fine habitat set (species subset = spp)
inds <- read.csv(file = "data/npmsIndicatorsIndicia_Aug2018.csv", header = T, stringsAsFactors = F) # repeated from script 3
selectSp <- function(bsh, sp) {tmp <- unique(inds[inds$broad.scale_habitat == bsh & inds$indicator_type == "positive",]$indiciaName)} 
# get species list for relevant broad habitat
focalSpp <- list()
focalSpp <- selectSp(bsh = "Lowland grassland")

sppDatList <- list()
# Function will print names of error-causing species to screen
# error causing species will currently just be a character element to the list, rather than a nested df
sppDatList <- lapply(focalSpp, function(x) spSamplePA(samples = habSamps, species = x))
names(sppDatList) <- focalSpp
## Species with no data for reference
excludedSpp <- unlist(lapply(seq_along(sppDatList), function(i) ifelse(which(!is.data.frame(sppDatList[[i]]))==1, names(sppDatList)[i], NULL)))
##
save(sppDatList, excludedSpp, file = "outputs/grasslandsEg_26 03 2019.Rdata")
#load(file = "outputs/grasslandsEg_26 03 2019.Rdata")
#####
# The next step is to run the model over each species data frame
# (i.e. run script 5_processData..., but over multiple species, and collecting useful results)
#####
