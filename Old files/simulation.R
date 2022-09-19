library(MASS)
## simulation seasonal model Nicolau synchrony paper 21 march 2022
rm(list=ls())

# spring model
beta1 <- 0.5
beta2 <- 0.1
beta3 <- -0.4
beta4 <- 0.1

# fall model
gamma1 <- 0.8
gamma2 <- 0.1
gamma3 <- -0.1
gamma4 <- -0.1

sigma.eps <- 1  # winter sd
corr.s <- 0.5
cov.s <- corr.s*sigma.eps*sigma.eps
Sigma.s <- matrix(c(sigma.eps^2,cov.s,cov.s,sigma.eps^2),2,2)

sigma.eta <- 0.5 # summer sd
corr.f <- 0.3
cov.f <- corr.f*sigma.eta*sigma.eta
Sigma.f <- matrix(c(sigma.eta^2,cov.f,cov.f,sigma.eta^2),2,2)

nsim<-100000
XY <- array(0, dim=c(nsim,8))

# initial values
XY [1:2,] <- 1

for (k in 3:nsim) {
  rand.corr.s <- mvrnorm(1,mu=c(0,0),Sigma=Sigma.s)
  rand.corr.f <- mvrnorm(1,mu=c(0,0),Sigma=Sigma.f)
  XY[k,5:8]<-c(rand.corr.s, rand.corr.f)
  XY[k,1] <- beta1*XY[k-1,2] + beta2*XY[k-1,1] +
    beta3*XY[k-2,2] + beta4*XY[k-2,1] + rand.corr.s[1]
  XY[k,2] <- gamma1*XY[k,1] + gamma2*XY[k-1,2] +
    gamma3*XY[k-1,1] + gamma4*XY[k-2,2] + rand.corr.f[1]
  XY[k,3] <- beta1*XY[k-1,4] + beta2*XY[k-1,3] +
    beta3*XY[k-2,4] + beta4*XY[k-2,3] + rand.corr.s[2]
  XY[k,4] <- gamma1*XY[k,3] + gamma2*XY[k-1,4] +
    gamma3*XY[k-1,3] + gamma4*XY[k-2,4] + rand.corr.f[2]
}

ind.plot <- 800:1000
pacf(XY[ind.plot,1])
pacf(XY[ind.plot,2])
arima(XY[ind.plot,1],order=c(2,0,0))  # ca 5 years cycle
arima(XY[ind.plot,2],order=c(2,0,0))  # also...

ind.plot <- 950:1000
par(mfrow=c(1,2))
plot(ind.plot,XY[ind.plot,1], type="l")
lines(ind.plot,XY[ind.plot,3], col="red")
plot(ind.plot,XY[ind.plot,2], type="l")
lines(ind.plot,XY[ind.plot,4], col="red")

ind.plot <- 90000:100000
cor(XY[ind.plot,5],XY[ind.plot,6])  # corr spring noise
cor(XY[ind.plot,7],XY[ind.plot,8])  # corr fall noise
cor(XY[ind.plot,1],XY[ind.plot,3])  # corr spring density
cor(XY[ind.plot,2],XY[ind.plot,4])  # corr fall density

cor(arima(XY[ind.plot,1],order=c(2,0,0))$residuals,
    arima(XY[ind.plot,3],order=c(2,0,0))$residuals) # correlation residuals AR2
cor(arima(XY[ind.plot,2],order=c(2,0,0))$residuals,
    arima(XY[ind.plot,4],order=c(2,0,0))$residuals) # correlation residuals AR2

# using a function to calculate different values
corr.funct <- function(sigma.eps,corr.s,sigma.eta,corr.f, nsim) {
  cov.s <- corr.s*sigma.eps*sigma.eps
  Sigma.s <- matrix(c(sigma.eps^2,cov.s,cov.s,sigma.eps^2),2,2)
  cov.f <- corr.f*sigma.eta*sigma.eta
  Sigma.f <- matrix(c(sigma.eta^2,cov.f,cov.f,sigma.eta^2),2,2)
  XY <- array(0, dim=c(nsim,8))
  # initial values
  XY [1:2,] <- 1
  for (k in 3:nsim) {
    rand.corr.s <- mvrnorm(1,mu=c(0,0),Sigma=Sigma.s)
    rand.corr.f <- mvrnorm(1,mu=c(0,0),Sigma=Sigma.f)
    XY[k,5:8]<-c(rand.corr.s, rand.corr.f)
    XY[k,1] <- beta1*XY[k-1,2] + beta2*XY[k-1,1] +
      beta3*XY[k-2,2] + beta4*XY[k-2,1] + rand.corr.s[1]
    XY[k,2] <- gamma1*XY[k,1] + gamma2*XY[k-1,2] +
      gamma3*XY[k-1,1] + gamma4*XY[k-2,2] + rand.corr.f[1]
    XY[k,3] <- beta1*XY[k-1,4] + beta2*XY[k-1,3] +
      beta3*XY[k-2,4] + beta4*XY[k-2,3] + rand.corr.s[2]
    XY[k,4] <- gamma1*XY[k,3] + gamma2*XY[k-1,4] +
      gamma3*XY[k-1,3] + gamma4*XY[k-2,4] + rand.corr.f[2]
  }
  ind.plot<-(nsim-1000):nsim
  coeff.ars <-arima(XY[ind.plot,1],order=c(2,0,0))$coef 
  coeff.arf <-arima(XY[ind.plot,2],order=c(2,0,0))$coef 
  cor.noise.s <- cor(XY[ind.plot,5],XY[ind.plot,6])  # corr spring noise
  cor.noise.f <- cor(XY[ind.plot,7],XY[ind.plot,8])  # corr fall noise
  cor.dens.s <- cor(XY[ind.plot,1],XY[ind.plot,3])  # corr spring density
  cor.dens.f <- cor(XY[ind.plot,2],XY[ind.plot,4])  # corr fall density
  cor.ar2.s <- cor(arima(XY[ind.plot,1],order=c(2,0,0))$residuals,
                   arima(XY[ind.plot,3],order=c(2,0,0))$residuals) # correlation residuals AR2
  cor.ar2.f <- cor(arima(XY[ind.plot,2],order=c(2,0,0))$residuals,
                   arima(XY[ind.plot,4],order=c(2,0,0))$residuals)
  list(sigma.eps=sigma.eps, sigma.eta=sigma.eta,
       coef.ars=coeff.ars, coef.arf=coeff.arf,  cor.noise.s=cor.noise.s,
       cor.noise.f=cor.noise.f, cor.dens.s =cor.dens.s,
       cor.dens.f=cor.dens.f,  cor.ar2.s=cor.ar2.s,  cor.ar2.f=cor.ar2.f)
}

res1 <- corr.funct(sigma.eps=1.2,corr.s=0.5,
                   sigma.eta=0.4,corr.f=0.3, nsim=10^5)
res2 <- corr.funct(sigma.eps=0.4,corr.s=0.5,
                   sigma.eta=0.4,corr.f=0.3, nsim=10^5)
res3 <- corr.funct(sigma.eps=0.8,corr.s=0.5,
                   sigma.eta=0.8,corr.f=0.3, nsim=10^5)
res4 <- corr.funct(sigma.eps=1.2,corr.s=0.5,
                   sigma.eta=1.2,corr.f=0.3, nsim=10^5)
res5 <- corr.funct(sigma.eps=0.4,corr.s=0.5,
                   sigma.eta=1.2,corr.f=0.3, nsim=10^5)


pdf("Plots/Simulations.pdf", height=5, width=10)
par(mfrow=c(1,2))
plot(1:5,c(res1$cor.noise.s,res2$cor.noise.s,res3$cor.noise.s,res4$cor.noise.s,res5$cor.noise.s),
     ylim=c(0,0.6),xaxt="n", ylab="Correlation", pch=19, xlab="Noise Term Pairs (Spring, Fall)")
points(1:5,c(res1$cor.dens.s,res2$cor.dens.s,res3$cor.dens.s,res4$cor.dens.s,res5$cor.dens.s),
       pch=19, col=2)
points(1:5,c(res1$cor.ar2.s,res2$cor.ar2.s,res3$cor.ar2.s,res4$cor.ar2.s,res5$cor.ar2.s),
       pch=19, col=4)
abline(h=seq(0,.6,.1),lty=3, col="gray60")

legend("bottomright", # title="Correlations",
       c("True","Abundances","AR(2) Residuals"), col=c(1,2,4),pch=c(19,19,19), horiz=FALSE, cex=0.8)
axis(side=1, at=1:5, labels=c("(1.2,0.4)", "(0.4, 0.4)","(0.8, 0.8)",
                              "(1.2,1.2)","(0.4, 1.2)"))
title(main="Spring")

plot(1:5,c(res1$cor.noise.f,res2$cor.noise.f,res3$cor.noise.f,res4$cor.noise.f,res5$cor.noise.f),
     ylim=c(0,0.6),xaxt="n", ylab="Correlation", pch=19, xlab="Noise Term Pairs (Spring, Fall)")
abline(h=seq(0,.6,.1),lty=3, col="gray60")

points(1:5,c(res1$cor.dens.f,res2$cor.dens.f,res3$cor.dens.f,res4$cor.dens.f,res5$cor.dens.f),
       pch=19, col=2)
points(1:5,c(res1$cor.ar2.f,res2$cor.ar2.f,res3$cor.ar2.f,res4$cor.ar2.f,res5$cor.ar2.f), # 0.0025 as jitter
       pch=19, col=4)
legend("bottomright", # title="Correlations",
       c("True","Abundances","AR(2) Residuals"), col=c(1,2,4),pch=c(19,19,19), horiz=FALSE, cex=0.8)
axis(side=1, at=1:5, labels=c("(1.2,0.4)", "(0.4, 0.4)","(0.8, 0.8)",
                              "(1.2,1.2)","(0.4, 1.2)"))
title(main="Fall")
dev.off()
