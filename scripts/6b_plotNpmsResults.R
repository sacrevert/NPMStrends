## 6b. Plot actual NPMS results from sacript 5b_processDataForTrends.R using coefplot::multiplot type appraoch
## v1 -- rough attempts for meeting on 18.01.2019
## 
# O.L. Pescott
# 17.01.2019
#rm(list=ls())
library(ggplot2)
#library(gridExtra)

load(file = "outputs/sppRuns/grasslands_16012019.Rdata") #loads list object sppModels
#load(file = "C:\\Users\\olipes\\Desktop\\grasslands_16012019.Rdata")

#1. mu.C from extract object 1
muAll <- do.call(rbind, lapply(lapply(sppModels, '[[', 1), function(x) as.data.frame(x[rownames(x) %in% c("mu.C"),])))
muAll$species <- rownames(muAll)

ggplot(muAll[c(1:43),], aes(species, `50%`)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                lwd=1.5, colour="grey", width=0) +
  #geom_errorbar(aes(ymin = quant0.25, ymax = quant0.75), 
  #              lwd=2.5, colour="gray0", width=0) +
  geom_point(size=2.0, pch=21, fill="white") +
  labs(y = "mu.C cover estimates") +
  theme_bw() +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 1145 x 973

ggplot(muAll[c(44:86),], aes(species, `50%`)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                lwd=1.5, colour="grey", width=0) +
  #geom_errorbar(aes(ymin = quant0.25, ymax = quant0.75), 
  #              lwd=2.5, colour="gray0", width=0) +
  geom_point(size=2.0, pch=21, fill="white") +
  labs(y = "mu.C cover estimates") +
  theme_bw() +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 1145 x 973