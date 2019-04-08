## 6b. Plot actual NPMS results from sacript 5b_processDataForTrends.R using coefplot::multiplot type appraoch
## v1 -- rough attempts for meeting on 18.01.2019
## 
# O.L. Pescott
# 17.01.2019
#rm(list=ls())
library(ggplot2)
#library(gridExtra)

#load(file = "outputs/sppRuns/grasslands_16012019.Rdata") #loads list object sppModels
#load(file = "C:\\Users\\olipes\\Desktop\\grasslands_16012019.Rdata")
#load(file = "outputs/sppRuns/grasslands_26032019.Rdata") 
load(file = "outputs/sppRuns/grasslands_26032019_v2_gamma0only.Rdata")

#test <- sppModels[[1]][[2]]
##################################
####1. mu.C from extract object 1
##################################
muAll <- do.call(rbind, lapply(lapply(sppModels, '[[', 2), function(x) as.data.frame(x[rownames(x) %in% c("mu.C"),])))
muAll$species <- rownames(muAll)
muAll <- muAll[order(muAll$`50%`, decreasing = T),]
muAll$species <- factor(muAll$species, levels = c(muAll$species))
muAll$plot <- rep(seq(1,2), 43)

ggplot(muAll[muAll$plot==1,], aes(species, `50%`)) + 
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

##################################
####2. annOcc years 1 and 3
##################################
aoAll1 <- do.call(rbind, lapply(lapply(sppModels, '[[', 2), function(x) as.data.frame(x[rownames(x) %in% c("annOcc[1]"),])))
aoAll1$species <- rownames(aoAll1)
aoAll1$year <- 1
aoAll1 <- aoAll1[order(aoAll1$`50%`, decreasing = T),]
aoAll1$species <- factor(aoAll1$species, levels = c(aoAll1$species))
aoAll1$plot <- rep(seq(1,2), 43)


ggplot(aoAll1[aoAll1$plot==1,], aes(species, `50%`)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                lwd=1.5, colour="red", width=0) +
  #geom_errorbar(aes(ymin = quant0.25, ymax = quant0.75), 
  #              lwd=2.5, colour="gray0", width=0) +
  scale_y_continuous(limits=c(0.4,0.65)) +
  geom_point(size=2.0, pch=21, fill="white") +
  labs(y = "Occupancy estimate: year 1") +
  theme_bw() +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 1145 x 973

aoAll3 <- do.call(rbind, lapply(lapply(sppModels, '[[', 2), function(x) as.data.frame(x[rownames(x) %in% c("annOcc[3]"),])))
aoAll3$species <- rownames(aoAll3)
aoAll3$year <- 3
ggplot(aoAll3[c(1:43),], aes(species, `50%`)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), 
                lwd=1.5, colour="grey", width=0) +
  #geom_errorbar(aes(ymin = quant0.25, ymax = quant0.75), 
  #              lwd=2.5, colour="gray0", width=0) +
  scale_y_continuous(limits=c(0.4,0.65)) +
  geom_point(size=2.0, pch=21, fill="white") +
  labs(y = "Occupancy estimate: year 3") +
  theme_bw() +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 1145 x 973

# plot years 1 and 3 together
aoAll13 <- rbind(aoAll1,aoAll3)
aoAll13 <- aoAll13[order(aoAll13$species, aoAll13$year),]
aoAll13$year <- as.factor(aoAll13$year)
ggplot(aoAll13[c(1:86),], aes(colour = year, group = year)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1, colour="grey50") +
  geom_pointrange(aes(x = species, y = `50%`, ymin = `2.5%`, ymax = `97.5%`), 
                lwd=0.5, position = position_dodge(width = 0.5)) +
  #geom_errorbar(aes(ymin = quant0.25, ymax = quant0.75), 
  #              lwd=2.5, colour="gray0", width=0) +
  scale_y_continuous(limits=c(0.35,0.75)) +
  #geom_point(size=2.0, pch=21, fill="white") +
  labs(y = "Occupancy estimate: years 1 & 3") +
  theme_bw() +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 1145 x 973


