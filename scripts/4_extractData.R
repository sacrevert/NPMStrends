## 4. Run extraction and processing for broad habitat X and species Y
## prepare Domin data for modelling
# O.L. Pescott
# 21.08.2018
rm(list=ls())

# Source file with datasets, getSamples() and spSamplePA() functions
source(file = "scripts/3_processDataFuns.R")

## Grasslands examples
grasslands <- c("Neutral pastures and meadows", "Dry acid grassland", "Dry calcareous grassland", "Neutral damp grassland", "Lowland grassland")
grassSamples <- getSamples(habsList = grasslands)

## Get and process data for species Y
Achi_mill_PAN <- spSamplePA(samples = grassSamples, species = "Achillea millefolium")
head(Achi_mill_PAN); tail(Achi_mill_PAN)

## Save dataset (20.09.2018)
save(Achi_mill_PAN, file = "data/Achi_mille_grassSamples_20180920.Rdata")
