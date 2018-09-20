# samples without dates
nrow(npms15_18spp[rowSums(is.na(npms15_18spp))== 0,]) # 1718 records without associated dates (1718/101508*100 = 1.7% of all records between 2015 and August 2018)
samples_PAN


nrow(samples_PAN[rowSums(is.na(samples_PAN))== 0,]) #
