### X3_speciesSamples_ErrorTesting.R
## Space for testing some functions etc.

species <- "Colchicum autumnale"
samples <- grassSamples
temp <- aggregate(preferred_taxon ~ sample_id + date, data = samples, function(x) max(ifelse(x == species, 1, 0)))
presence <- temp[temp$preferred_taxon == 1,] # samples containing the taxon of interest
presSamps <- merge(presence, npms15_18plots, by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # useful for subsequent rbind
presSamps$combination <- paste(presSamps$surv_habitat,", ",presSamps$title, sep = "") # useful for subsequent rbind
absence <- temp[temp$preferred_taxon == 0,] # samples NOT containing the taxon of interest (but are these 0 or NA?)
absSamps <- merge(absence, npms15_18plots, by.x = "sample_id", by.y = "sample", all.x = T, all.y = F) # habitat/level info for every sample
absSamps$combination <- paste(absSamps$surv_habitat,", ",absSamps$title, sep = "") # create column for lookup to indsLookup
indsLookup_fil <- indsLookup[indsLookup$indiciaName == species,]
absSamps_AN <- merge(absSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# absent or NA indicator
presSamps_AN <- merge(presSamps, indsLookup_fil, by.x = "combination", by.y = "combined", all.x = T, all.y = F)# again, for rbind
absSamps_AN$preferred_taxon <- ifelse(is.na(absSamps_AN$indiciaName), "NA", 0)
# but inventory samples should always be zero if absent
absSamps_AN$preferred_taxon <- ifelse(absSamps_AN$title =="Inventory survey", 0, absSamps_AN$preferred_taxon)
samples_PAN <- rbind(presSamps_AN, absSamps_AN) # combine so we have a PAN (pres, abs, NA) view on our focal species
names(samples_PAN)[4] <- "PAN" # rename to avoid confusion
samples_PAN$PAN <- as.numeric(samples_PAN$PAN) 
# add in domins and dates etc. from original sample/species data
samples_PAN <- merge(samples_PAN, npms15_18spp[npms15_18spp$preferred_taxon == species,], by.x = "sample_id", by.y = "sample_id", all.x = T, all.y = F)
samples_PAN$domin <- ifelse(samples_PAN$PAN == 0, 0, samples_PAN$domin)
samples_PAN <- merge(samples_PAN, domins, by.x = "domin", by.y = "dominOrig", all.x = T, all.y = F)
return(samples_PAN)
