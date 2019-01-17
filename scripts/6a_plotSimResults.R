## 6. Plot simulation results using coefplot::multiplot type appraoch
## v1 -- rough attempts for meeting on 18.01.2019
## 
# O.L. Pescott
# 16.01.2019
#rm(list=ls())
library(ggplot2)
#library(gridExtra)

load(file="outputs/model1.Rdata")
load(file="outputs/model2.Rdata")
load(file="outputs/model3.Rdata")
load(file="outputs/model4.Rdata")
load(file="outputs/model5.Rdata")
load(file="outputs/model6.Rdata")
load(file="outputs/model7.Rdata")
load(file="outputs/model8.Rdata")
load(file="outputs/model9.Rdata")

#################
### Models 1 to 3
#################
# harmonise
tableOut_m1 <- tableOut_m1[c(1,11:14),]
tableOut_m2 <- tableOut_m2[c(1,6,8:10),]
tableOut_m3 <- tableOut_m3[c(1,5,7:9),]


# Put model estimates into temporary data.frames:
model1Frame <- data.frame(Variable = rownames(tableOut_m1),
                          Median = tableOut_m1$`50%`,
                          lCI = tableOut_m1$`2.5%`,
                          uCI = tableOut_m1$`97.5%`,
                          modelName = "Model 1")
model2Frame <- data.frame(Variable = rownames(tableOut_m2),
                          Median = tableOut_m2$`50%`,
                          lCI = tableOut_m2$`2.5%`,
                          uCI = tableOut_m2$`97.5%`,
                          modelName = "Model 2")
model3Frame <- data.frame(Variable = rownames(tableOut_m3),
                          Median = tableOut_m3$`50%`,
                          lCI = tableOut_m3$`2.5%`,
                          uCI = tableOut_m3$`97.5%`,
                          modelName = "Model 3")
# Combine these data.frames
allModelFrame <- data.frame(rbind(model1Frame, model2Frame, model3Frame))  # etc.

# Plot
p1 <- ggplot(allModelFrame, aes(colour = modelName)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                              lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.key.size = unit(2,"cm")) +
  coord_flip() # export 558 x 509

#################
### Models 4 to 6
#################
# harmonise
tableOut_m4 <- tableOut_m4[c(1,6:9),]
tableOut_m5 <- tableOut_m5[c(1,6:9),]
tableOut_m6 <- tableOut_m6[c(1,6:9),]


# Put model estimates into temporary data.frames:
model4Frame <- data.frame(Variable = rownames(tableOut_m4),
                          Median = tableOut_m4$`50%`,
                          lCI = tableOut_m4$`2.5%`,
                          uCI = tableOut_m4$`97.5%`,
                          modelName = "Model 4")
model5Frame <- data.frame(Variable = rownames(tableOut_m5),
                          Median = tableOut_m5$`50%`,
                          lCI = tableOut_m5$`2.5%`,
                          uCI = tableOut_m5$`97.5%`,
                          modelName = "Model 5")
model6Frame <- data.frame(Variable = rownames(tableOut_m6),
                          Median = tableOut_m6$`50%`,
                          lCI = tableOut_m6$`2.5%`,
                          uCI = tableOut_m6$`97.5%`,
                          modelName = "Model 6")
# Combine these data.frames
allModelFrame456 <- data.frame(rbind(model4Frame, model5Frame, model6Frame))  # etc.

# Plot
p456 <- ggplot(allModelFrame456, aes(colour = modelName)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.key.size = unit(2,"cm")) +
  coord_flip() # export 558 x 509

#################
### Models 7 to 9
#################
# harmonise
tableOut_m7 <- tableOut_m7[c(1,6:9),]
tableOut_m8 <- tableOut_m8[c(1,6:9),]
tableOut_m9 <- tableOut_m9[c(1,6:9),]


# Put model estimates into temporary data.frames:
model7Frame <- data.frame(Variable = rownames(tableOut_m7),
                          Median = tableOut_m7$`50%`,
                          lCI = tableOut_m7$`2.5%`,
                          uCI = tableOut_m7$`97.5%`,
                          modelName = "Model 7")
model8Frame <- data.frame(Variable = rownames(tableOut_m8),
                          Median = tableOut_m8$`50%`,
                          lCI = tableOut_m8$`2.5%`,
                          uCI = tableOut_m8$`97.5%`,
                          modelName = "Model 8")
model9Frame <- data.frame(Variable = rownames(tableOut_m9),
                          Median = tableOut_m9$`50%`,
                          lCI = tableOut_m9$`2.5%`,
                          uCI = tableOut_m9$`97.5%`,
                          modelName = "Model 9")
# Combine these data.frames
allModelFrame789 <- data.frame(rbind(model7Frame, model8Frame, model9Frame))  # etc.
allModelFrame789$uCI[which(allModelFrame789$uCI>25)] <-25 # hard limit for plotting (not advisable!)

# Plot
p789 <- ggplot(allModelFrame789, aes(colour = modelName)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  scale_y_continuous(limits = c(-5, 25)) +
  theme(legend.key.size = unit(2,"cm")) +
  coord_flip() # export 558 x 509

