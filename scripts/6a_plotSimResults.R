## 6. Plot simulation results using coefplot::multiplot type appraoch
## v1 -- rough attempts for meeting on 18.01.2019
## 
# O.L. Pescott
# 16.01.2019
#rm(list=ls())
library(ggplot2)
library(gridExtra)

load(file="outputs/modelPD.Rdata")
load(file="outputs/model0.Rdata")
load(file="outputs/model1.Rdata")
#load(file="outputs/model1a.Rdata") # model 1a same as model 1 but with change from 0.9999999999 -> 0.99 (in simulation and JAGS code -- makes no difference)
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
tableOut <- tableOut[c(1:5),] # perfect detection
tableOut_m0 <- tableOut_m0[c(1:5),]
tableOut_m1 <- tableOut_m1[c(1:5),]
#tableOut_m1a <- tableOut_m1a[c(1,11:14),]
tableOut_m2 <- tableOut_m2[c(1:5),]
tableOut_m3 <- tableOut_m3[c(1:5),]


# Put model estimates into temporary data.frames:
modelmFrame <- data.frame(Variable = rownames(tableOut),
                          Median = tableOut$`50%`,
                          lCI = tableOut$`2.5%`,
                          uCI = tableOut$`97.5%`,
                          modelName = "Perfect detection")
model0Frame <- data.frame(Variable = rownames(tableOut_m0),
                          Median = tableOut_m0$`50%`,
                          lCI = tableOut_m0$`2.5%`,
                          uCI = tableOut_m0$`97.5%`,
                          modelName = "Model 0")
model1Frame <- data.frame(Variable = rownames(tableOut_m1),
                          Median = tableOut_m1$`50%`,
                          lCI = tableOut_m1$`2.5%`,
                          uCI = tableOut_m1$`97.5%`,
                          modelName = "Model 1")
#model1aFrame <- data.frame(Variable = rownames(tableOut_m1a),
#                          Median = tableOut_m1a$`50%`,
#                          lCI = tableOut_m1a$`2.5%`,
#                          uCI = tableOut_m1a$`97.5%`,
#                          modelName = "Model 1a")
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
allModelFrame <- data.frame(rbind(modelmFrame, model0Frame, model1Frame, model2Frame, model3Frame))  # etc.
#levels(allModelFrame$modelName) <- c("Model 3", "Model 2", "Model 1", "Model 0")
allModelFrame <- allModelFrame[order(allModelFrame$modelName, decreasing = T),]
allModelFrame[c(22:23),c(2:4)] <- NA# remove gamma0 and gamma1 rows for perfect detection model

# Plot
p <- ggplot(allModelFrame, aes(colour = modelName)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                              lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.key.size = unit(2,"cm")) +
  coord_flip() # export 558 x 509
p
### Improve display using gridExtra
p1 <- ggplot(allModelFrame[allModelFrame$Variable=="avgOcc",], aes(colour = modelName)) + 
  geom_hline(yintercept=1, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p2 <- ggplot(allModelFrame[allModelFrame$Variable=="gamma0",], aes(colour = modelName)) + 
  geom_hline(yintercept=-2, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p3 <- ggplot(allModelFrame[allModelFrame$Variable=="gamma1",], aes(colour = modelName)) + 
  geom_hline(yintercept=3, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p4 <- ggplot(allModelFrame[allModelFrame$Variable=="mu.C",], aes(colour = modelName)) + 
  geom_hline(yintercept=0.5, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p5 <- ggplot(allModelFrame[allModelFrame$Variable=="tau.C",], aes(colour = modelName)) + 
  geom_hline(yintercept=c(3,10), lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  #theme(legend.key.size = unit(0.5,"cm")) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
gridExtra::grid.arrange(p1,p2,p3,p4,p5, nrow = 2, ncol = 3)

#################
### Models 4 to 6
#################
# harmonise
tableOut_m4 <- tableOut_m4[c(1:5),]
tableOut_m5 <- tableOut_m5[c(1:5),]
tableOut_m6 <- tableOut_m6[c(1:5),]


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
p456
### Improve display using gridExtra
p6 <- ggplot(allModelFrame456[allModelFrame456$Variable=="avgOcc",], aes(colour = modelName)) + 
  geom_hline(yintercept=c(0.75,0.5,0.25), lty=2, lwd=1, color = c("#F8766D","#00BA38" ,"#619CFF")) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p7 <- ggplot(allModelFrame456[allModelFrame456$Variable=="gamma0",], aes(colour = modelName)) + 
  geom_hline(yintercept=-2, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p8 <- ggplot(allModelFrame456[allModelFrame456$Variable=="gamma1",], aes(colour = modelName)) + 
  geom_hline(yintercept=3, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p9 <- ggplot(allModelFrame456[allModelFrame456$Variable=="mu.C",], aes(colour = modelName)) + 
  geom_hline(yintercept=0.5, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p10 <- ggplot(allModelFrame456[allModelFrame456$Variable=="tau.C",], aes(colour = modelName)) + 
  geom_hline(yintercept=3, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  #theme(legend.key.size = unit(0.5,"cm")) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
gridExtra::grid.arrange(p6,p7,p8,p9,p10, nrow = 2, ncol = 3)

#################
### Models 7 to 9
#################
# harmonise
tableOut_m7 <- tableOut_m7[c(1:5),]
tableOut_m8 <- tableOut_m8[c(1:5),]
tableOut_m9 <- tableOut_m9[c(1:5),]


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
#allModelFrame789$uCI[which(allModelFrame789$uCI>25)] <-25 # hard limit for plotting (not advisable!)

# Plot
p789 <- ggplot(allModelFrame789, aes(colour = modelName)) + 
  #geom_hline(yintercept=0, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  scale_y_continuous(limits = c(-5, 25)) +
  theme(legend.key.size = unit(2,"cm")) +
  coord_flip() # export 558 x 509
p789
p11 <- ggplot(allModelFrame789[allModelFrame789$Variable=="avgOcc",], aes(colour = modelName)) + 
  geom_hline(yintercept=c(0.25), lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p12 <- ggplot(allModelFrame789[allModelFrame789$Variable=="gamma0",], aes(colour = modelName)) + 
  geom_hline(yintercept=-2, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p13 <- ggplot(allModelFrame789[allModelFrame789$Variable=="gamma1",], aes(colour = modelName)) + 
  geom_hline(yintercept=3, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p14 <- ggplot(allModelFrame789[allModelFrame789$Variable=="mu.C",], aes(colour = modelName)) + 
  geom_hline(yintercept=c(0.025,0.10,0.25), lty=2, lwd=1, color = c("#619CFF","#00BA38" ,"#F8766D")) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
p15 <- ggplot(allModelFrame789[allModelFrame789$Variable=="tau.C",], aes(colour = modelName)) + 
  geom_hline(yintercept=3, lty=2, lwd=1) +
  geom_pointrange(aes(x = Variable,y = Median,  ymin = lCI, ymax = uCI),
                  lwd = 2, position = position_dodge(width = 1/2),shape = 21, fill = "WHITE", fatten = 1) +
  #theme(legend.key.size = unit(0.5,"cm")) +
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_blank() ) +
  coord_flip() # export 558 x 509
gridExtra::grid.arrange(p11,p12,p13,p14,p15, nrow = 2, ncol = 3)

## Get colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(3)
