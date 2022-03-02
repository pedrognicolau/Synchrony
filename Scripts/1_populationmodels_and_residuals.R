# Population Models -------

library(INLA)
#stackeddata <- readRDS("Data/stacked_data.rds")
stackeddata <- readRDS("Data/data4synchrony.rds")
View(stackeddata)
# expanded stacked data with extra variables
# stackeddata2 <- readRDS("~/OneDrive - UiT Office 365/Vole Synchrony/data/stacked_data2.rds")

# growth rates instead
GR_data0 <- readRDS("Vole Synchrony/data/log_centered_INLACR_growthrates.rds")

## AR(2) models =======

### Spring AR model ----
# (II) Spring

St_II <- inla(St ~ St_1 + St_2
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)

### Fall AR model ####
# (II) Fall
Ft_II <- inla(Ft ~ Ft_1 + Ft_2
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)


## AR(2) + region (III) models =======

### Spring Yearly model ----
# (III) Spring

St_III <- inla(St ~ 
                  St_1.1 + St_1.2 + St_1.3 +
                  St_2.1 + St_2.2 + St_2.3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)

### Fall Yearly model ####
# (III) Fall
Ft_III <- inla(Ft ~ 
                  Ft_1.1 + Ft_1.2 + Ft_1.3 +
                  Ft_2.1 + Ft_2.2 + Ft_2.3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)


## Season + region Models (IV) =======
### Spring seasonal model ####

# St_seasonal <- inla(S_t ~ 
#                       Ft_1R1 + Ft_1R2 + Ft_1R3 +
#                       St_1R1 + St_1R2 + St_1R3 +
#                       Ft_2R1 + Ft_2R2 + Ft_2R3 +
#                       St_2R1 + St_2R2 + St_2R3 
#                     # + f(time, mod="rw2", scale.model = TRUE,
#                     #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
#                     #, family="gaussian"
#                     ,
#                     control.predictor = list(compute = TRUE),
#                     control.compute = list(config = TRUE),
#                     data=stackeddata)
# (IV) Spring
St_IV <- inla(St ~ 
                      Ft_1.1 + Ft_1.2 + Ft_1.3 +
                      St_1.1 + St_1.2 + St_1.3 +
                      Ft_2.1 + Ft_2.2 + Ft_2.3 +
                      St_2.1 + St_2.2 + St_2.3 
                    # + f(time, mod="rw2", scale.model = TRUE,
                    #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                    #, family="gaussian"
                    ,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(config = TRUE),
                    data=stackeddata)

### Fall seasonal model ####
# Ft_seasonal <- inla(F_t ~ StR1 + StR2 + StR3 +
#                       Ft_1R1 + Ft_1R2 + Ft_1R3 +
#                       St_1R1 + St_1R2 + St_1R3 +
#                       Ft_2R1 + Ft_2R2 + Ft_2R3 
#                     # + f(time, mod="rw2", scale.model = TRUE,
#                     #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
#                     # , family="gaussian"
#                     ,
#                     control.predictor = list(compute = TRUE),
#                     control.compute = list(config = TRUE),
#                     data=stackeddata)

# (IV) Fall

Ft_IV <- inla(Ft ~ St.1 + St.2 + St.3 +
                      Ft_1.1 + Ft_1.2 + Ft_1.3 +
                      St_1.1 + St_1.2 + St_1.3 +
                      Ft_2.1 + Ft_2.2 + Ft_2.3 
                    # + f(time, mod="rw2", scale.model = TRUE,
                    #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                    # , family="gaussian"
                    ,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(config = TRUE),
                    data=stackeddata)

### List of All models ####
# popmodels <- list(St_seasonal=St_seasonal,St_year=St_year,Ft_seasonal=Ft_seasonal,Ft_year=Ft_year)
popmodels <- list(St_II=St_II,St_III=St_III,St_IV=St_IV,
                  Ft_II=Ft_II,Ft_III=Ft_III,Ft_IV=Ft_IV)



## Obtain Residuals ####

source("Scripts/0_Important_functions.R")

# obtain residuals using Bayesian sample
St_II_res <- matrix(residual_samples(stackeddata$St, St_II, ns=200, mean = TRUE),ncol=19)
Ft_II_res <- matrix(residual_samples(stackeddata$Ft, Ft_II, ns=200, mean = TRUE),ncol=19)
St_III_res <- matrix(residual_samples(stackeddata$St, St_III, ns=200, mean = TRUE),ncol=19)
Ft_III_res <- matrix(residual_samples(stackeddata$Ft, Ft_III, ns=200, mean = TRUE),ncol=19)
St_IV_res <- matrix(residual_samples(stackeddata$St, St_IV, ns=200, mean = TRUE),ncol=19)
Ft_IV_res <- matrix(residual_samples(stackeddata$Ft, Ft_IV, ns=200, mean = TRUE),ncol=19)

popres <- list(St_II=St_II_res,    Ft_II=Ft_II_res,
               St_III=St_III_res,  Ft_III=Ft_III_res,
               St_IV=St_IV_res,    Ft_IV=Ft_IV_res)

saveRDS(popres,"Data/popmodels_residuals.rds")

# Plots ####
## plot models figure 2 ####
library(plotrix)
sII <- popmodels$St_II$summary.fixed
sIII <- popmodels$St_III$summary.fixed
sIV <- popmodels$St_IV$summary.fixed

FII <- popmodels$Ft_II$summary.fixed
FIII <- popmodels$Ft_III$summary.fixed
FIV <- popmodels$Ft_IV$summary.fixed


### spring coefficients ####
datx <- sIV[-1,]

pdf("Plots/spring_seas_coefs.pdf", width=6, height = 6)
par(mar=c(2,4,1,1),mfrow=c(1,1))
plot(1, type="n", xlab="",ylab="", ylim=c(-1, 1), xlim=c(1,12), xaxt="n",
     yaxt="n", frame.plot = FALSE)

plotCI(1:12,y= datx[,4],datx[,5]-datx[,4],datx[,4]-datx[,3], col=rep(c(4,3,7),4), add=TRUE)
datx <- ss[-1,]
Axis(side=2)
abline(h=0, lty=3)
abline(v=c(3.5,6.5,9.5),lty=1)

mtext(expression(beta), side=2, line=2.5, at=0, cex=1.5)

mtext(expression(Y[t-1]), side=1, line=1, at=2, cex=1.2)
mtext(expression(X[t-1]), side=1, line=1, at=5, cex=1.2)
mtext(expression(Y[t-2]), side=1, line=1, at=8, cex=1.2)
mtext(expression(X[t-2]), side=1, line=1, at=11, cex=1.2)
dev.off()

pdf("Plots/fall_seas_coefs.pdf", width=6, height = 6)
par(mar=c(2,4,1,1),mfrow=c(1,1))
plot(1, type="n", xlab="",ylab="", ylim=c(-1, 1), xlim=c(1,12), xaxt="n",
     yaxt="n", frame.plot = FALSE)

### fall coefficients ####
datx <- FIV[-1,]

plotCI(1:12,y= datx[,4],datx[,5]-datx[,4],datx[,4]-datx[,3], col=rep(c(4,3,7),4), add=TRUE)
Axis(side=2)
abline(h=0, lty=3)
abline(v=c(3.5,6.5,9.5),lty=1)

mtext(expression(gamma), side=2, line=2.5, at=0, cex=1.5)

mtext(expression(X[t]), side=1, line=1, at=2, cex=1.2)
mtext(expression(Y[t-1]), side=1, line=1, at=5, cex=1.2)
mtext(expression(X[t-1]), side=1, line=1, at=8, cex=1.2)
mtext(expression(Y[t-2]), side=1, line=1, at=11, cex=1.2)
dev.off()

## R2 Plots


stII_r2 <- inla_R2(stackeddata$St,St_II, 1000) 
ftII_r2 <- inla_R2(stackeddata$Ft,Ft_II, 1000) 

stIII_r2 <- inla_R2(stackeddata$St,St_III, 1000) 
ftIII_r2 <- inla_R2(stackeddata$Ft,Ft_III, 1000) 

stIV_r2 <- inla_R2(stackeddata$St,St_IV, 1000) 
ftIV_r2 <- inla_R2(stackeddata$Ft,Ft_IV, 1000)

pdf("Plots/r2_4models.pdf", width=6, height = 6)
#c(bottom, left, top, right)
par(mar=c(0.2,0,0,0)+.5,mfrow=c(3,1), oma=c(5, 0,0,0))
hist(stII_r2, breaks=seq(0,.7,0.01), col=scales::alpha(4,.6), 
     xlab = "", probability = TRUE, main = "", border=F, yaxt="n",ylim = c(0,20),
     ylab="")
hist(ftII_r2, breaks=seq(0,1,0.01), add=TRUE, col=scales::alpha(2,.6), xlab = "", probability = TRUE, main = "", border=F, yaxt="n",
     ylab="")
text(x=.65,y=15," (II)", cex=1.2)


hist(stIII_r2, breaks=seq(0,0.7,0.011), col=scales::alpha(4,.7), 
     xlab = "", probability = TRUE, main = "", border=F, yaxt="n",ylim = c(0,20),
     ylab="")
hist(c, breaks=seq(0,0.7,0.011), add=TRUE, col=scales::alpha(2,.7), probability=TRUE, border=F)
text(x=.65,y=15,"(III)", cex=1.2)
mean(stIV_r2)
mean(ftII_r2)

hist(stIV_r2, breaks=seq(0,0.7,0.011), col=scales::alpha(4,.7), 
     xlab = "", probability = TRUE, main = "", border=F, yaxt="n",ylim = c(0,20),
     ylab="")
hist(ftIV_r2, breaks=seq(0,0.7,0.011), add=TRUE, col=scales::alpha(2,.7), probability=TRUE, border=F)
text(x=.65,y=15,"(IV)", cex=1.2)

mtext(expression(R^2), side=1, line=4, at=0.35, cex=1.5)
#Axis(side = 2, at=NULL)
# abline(v=median(st_r2), lwd=1, lty=2)
# hist(stII_r2, breaks=seq(0,0.7,0.011), add=TRUE, col=scales::alpha("orange1",.7), probability=TRUE, border=F)
# abline(v=median(ft_r2), lwd=1, lty=2)
dev.off()

