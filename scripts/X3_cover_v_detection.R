rm(list=ls())
## Define potentially useful functions
ilt<-function(x){ # inverse logit function
  exp(x)/(1+exp(x))
}
logit <- function(x){ # logit function
  log(x/(1 - x))
}

x <- seq(0,1,0.05)
y <- ilt(-2 + 3*x)
y1 <- -2 + 3*x
par(mfrow=c(1, 2))
#plot(x,y, xlab = "Cover", ylab = "Detection probability", ylim = c(0,1))
#plot(x,y1, xlab = "Cover", ylab = "logit(Detection probability)", ylim = c(-2,1))
lm1 <- lm(y ~ x)
#plot(lm1)
#lines(y = predict(lm1, newdata = data.frame(x=x)), x =x)

# screen out lower observations
x2 <- seq(0.5,1,0.025)
y2 <- ilt(-2 + 3*x2)
lm2 <- lm(y2 ~ x2)

#par(mfrow=c(1, 2))
plot(x,y, xlab = "Cover", ylab = "Detection probability", ylim = c(0,1), pch = 19, col = "gray")
lines(y = predict(lm1, newdata = data.frame(x=x)), x =x, col = "red", lwd = 2)
lines(y = predict(lm2, newdata = data.frame(x2=x)), x = x, col = "blue", lwd = 2)
text("(a)", x= 0.05, y =0.95)

plot(x2,y2, xlab = "Cover", ylab = "Detection probability", ylim = c(0,1), xlim = c(0,1), pch = 19, col = "gray")
lines(y = predict(lm2, newdata = data.frame(x2=x)), x = x, col = "blue", lwd = 2)
text("(b)", x= 0.05, y =0.95)

#################
##### Logit scale
#################
y1 <- -2 + 3*x
par(mfrow=c(1, 2))
#plot(x,y, xlab = "Cover", ylab = "Detection probability", ylim = c(0,1))
#plot(x,y1, xlab = "Cover", ylab = "logit(Detection probability)", ylim = c(-2,1))
lm3 <- lm(y1 ~ x)
#plot(lm1)
#lines(y = predict(lm1, newdata = data.frame(x=x)), x =x)

# screen out lower observations
x3 <- seq(0.5,1,0.025)
y3 <- -2 + 3*x3
lm4 <- lm(y3 ~ x3)

#par(mfrow=c(1, 2))
plot(x,y1, xlab = "Cover", ylab = "Logit(Detection probability)", ylim = c(-2,1), pch = 19, col = "gray")
lines(y = predict(lm3, newdata = data.frame(x=x)), x =x, col = "red", lwd = 2)
lines(y = predict(lm4, newdata = data.frame(x3=x)), x = x, col = "blue", lwd = 2)
text("(a)", x= 0.05, y =0.95)

plot(x3,y3, xlab = "Cover", ylab = "Logit(Detection probability)", ylim = c(-2,1), xlim = c(0,1), pch = 19, col = "gray")
lines(y = predict(lm4, newdata = data.frame(x3=x)), x = x, col = "blue", lwd = 2)
text("(b)", x= 0.05, y =0.95)

#plot(x,y1, xlab = "Cover", ylab = "logit(Detection probability)", ylim = c(-2,1))