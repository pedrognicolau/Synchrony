# Population Models -------

library(INLA)
#stackeddata <- readRDS("Data/stacked_data.rds")
stackeddata <- readRDS("Data/data4synchrony.rds")
# expanded stacked data with extra variables
#stackeddata2 <- readRDS("~/OneDrive - UiT Office 365/Vole Synchrony/data/stacked_data2.rds")

# growth rates instead
GR_data0 <- readRDS("Vole Synchrony/data/log_centered_INLACR_growthrates.rds")
## Hansen 1999 Models =======
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

St_seasonal <- inla(St ~ 
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

Ft_seasonal <- inla(Ft ~ St.1 + St.2 + St.3 +
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

### Spring Yearly model ####
St_year <- inla(St ~ 
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
Ft_year <- inla(Ft ~ 
                  Ft_1.1 + Ft_1.2 + Ft_1.3 +
                  Ft_2.1 + Ft_2.2 + Ft_2.3 
                # + f(time, mod="rw2", scale.model = TRUE,
                #            hyper = list(theta = list(prior="pc.prec", param=c(u=1,0.01))))
                # , family="gaussian"
                ,
                control.predictor = list(compute = TRUE),
                control.compute = list(config = TRUE),
                data=stackeddata)


### List of All models ####
popmodels <- list(St_seasonal=St_seasonal,St_year=St_year,Ft_seasonal=Ft_seasonal,Ft_year=Ft_year)



## Obtain Residuals ####

source("Scripts/0_Important_functions.R")

# obtain residuals using Bayesian sample
St_seas_res <- matrix(residual_samples(stackeddata$St, St_seasonal, ns=200, mean = TRUE),ncol=19)
Ft_seas_res <- matrix(residual_samples(stackeddata$Ft, Ft_seasonal, ns=200, mean = TRUE),ncol=19)
St_year_res <- matrix(residual_samples(stackeddata$St, St_year, ns=200, mean = TRUE),ncol=19)
Ft_year_res <- matrix(residual_samples(stackeddata$Ft, Ft_year, ns=200, mean = TRUE),ncol=19)

popres <- list(St_seasonal=St_seas_res,Ft_seasonal=Ft_seas_res,St_year=St_year_res,Ft_year=Ft_year_res)

saveRDS(popres,"Data/popmodels_residuals.rds")

# Plots ####
## plot models figure 2 ####
library(plotrix)
ss <- popmodels$St_seasonal$summary.fixed
fs <- popmodels$Ft_seasonal$summary.fixed
sy <- popmodels$St_seasonal$summary.fixed
fy <- popmodels$St_seasonal$summary.fixed


### spring coefficients ####
datx <- ss[-1,]

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

datx <- fs[-1,]

### fall coefficients ####
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


st_r2 <- inla_R2(stackeddata$St,St_seasonal, 1000) 
ft_r2 <- inla_R2(stackeddata$Ft,Ft_seasonal, 1000) 

pdf("Plots/r2_seasonalmodels.pdf", width=6, height = 6)
par(mar=c(4,2,1,1),mfrow=c(1,1))
hist(st_r2, breaks=seq(0.2,0.7,0.01), ylim = c(0,20), col=scales::alpha("green4",.6), xlab = expression(R^2), probability = TRUE, main = "", border=F, yaxt="n",
     ylab="")

mtext("Number of Samples", side=2, line=0, at=10, cex=1.2)
#Axis(side = 2, at=NULL)
abline(v=median(st_r2), lwd=1, lty=2)
hist(ft_r2, breaks=seq(0,0.7,0.011), add=TRUE, col=scales::alpha("orange1",.7), probability=TRUE, border=F)
abline(v=median(ft_r2), lwd=1, lty=2)
dev.off()

