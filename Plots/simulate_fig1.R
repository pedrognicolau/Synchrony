

## generate exponential decay -----
library(scales)
seqx <- sort((c(5:40,seq(35,50,2),seq(50,80,2.5),c(51,52,65,77,81,89),seq(80,90,4))))
length(seqx)
colsss = c("orangered2","deepskyblue")
combn(12,2)
pdf("Plots/simulated_sync4.pdf", width = 7*1.2, height=2.7*1.2)
set.seed(10)
par(mfrow=c(1,4), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))

pdf("Plots/grid4sync.pdf", width = 7*1.2, height=2.7*1.2)
par(mfrow=c(1,4), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))
for(i in 1:4)
{plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
      yaxt="n", frame.plot = FALSE, main="")
  abline(h=c(0,.25,.5,.75,1), lty=3, col="gray70")
  Axis(side=1, labels=FALSE, line=1, at=c(0,90))
}

dev.off()

### Model I ####
y0 = (0.9)*exp(-0.003*seqx)+rnorm(length(seqx),0,0.03)
y01 = (0.83)*exp(-0.003*seqx)+rnorm(length(seqx),0,0.03)

syncdata <- data.frame(dist=seqx, y0=y0, y01=y01)
syncmodx1 <- INLA::inla(y0 ~ f(dist, model="rw2", scale.model = TRUE,
                              hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                 , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx2 <- INLA::inla(y01 ~ f(dist, model="rw2", scale.model = TRUE,
                        hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                 , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx <- list(syncmodx1,syncmodx2)

plot(1, type="n", xlab="",ylab="Synchrony", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="(I)")
abline(h=c(0,.25,.5,.75,1), lty=3, col="gray70")
points(seqx,y0, ylim=c(0,1), pch=19, col=alpha("red3",1))
points(seqx,y01, pch=19, col=alpha("deepskyblue4",1))
Axis(side=2, labels=TRUE, line=1, at=c(0,0.5,1),cex.axis=1.5)
Axis(side=1, labels=FALSE, line=1, at=c(0,90))

for(i in 1:2)
{  
  syncmod <- syncmodx[[i]]
  lines(seqx, syncmod$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
        ylim=c(0,.7), main="", type="l", lwd=2, col=colsss[i], lty=1)    
  min_a <- pmin(syncmod$summary.fitted.values$`0.025quant`)
  max_a <- pmax(syncmod$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      seqx,
      rev(seqx)
    ),
    c(max_a , rev(min_a)),
    col = scales::alpha(colsss[i],.5),
    border = NA
  )
}
#mtext("Population synchrony",side=2, line=3, cex=1)
#mtext("Distance",side=1, line=2, cex=1)

### Model II ###

yar1 = (0.75)*exp(-0.004*seqx)+rnorm(length(seqx),0,0.03)
yar2 = (0.65)*exp(-0.004*seqx)+rnorm(length(seqx),0,0.03)
syncdata <- data.frame(dist=seqx, y0=yar1, y01=yar2)
syncmodx1 <- INLA::inla(y0 ~ f(dist, model="rw2", scale.model = TRUE,
                               hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx2 <- INLA::inla(y01 ~ f(dist, model="rw2", scale.model = TRUE,
                                hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx <- list(syncmodx1,syncmodx2)
plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="(II)")
abline(h=c(0,.25,.5,.75,1), lty=3, col="gray70")
points(seqx, yar1, ylim=c(0,1), pch=19, col="red3")
points(seqx, yar2, ylim=c(0,1), pch=19, col="deepskyblue4")
Axis(side=1, labels=FALSE, line=1, at=c(0,90), cex=2)
for(i in 1:2)
{  
  syncmod <- syncmodx[[i]]
  lines(seqx, syncmod$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
        ylim=c(0,.7), main="", type="l", lwd=1.5, col=colsss[i])    
  min_a <- pmin(syncmod$summary.fitted.values$`0.025quant`)
  max_a <- pmax(syncmod$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      seqx,
      rev(seqx)
    ),
    c(max_a , rev(min_a)),
    col = scales::alpha(colsss[i],.5),
    border = NA
  )
}
#mtext("Distance",side=1, line=2, cex=1)

### Model III ####

yarR1 = (0.60)*exp(-0.005*seqx)+rnorm(length(seqx),0,0.03)
yarR2 = (0.40)*exp(-0.005*seqx)+rnorm(length(seqx),0,0.025)
syncdata <- data.frame(dist=seqx, y0=yarR1, y01=yarR2)
syncmodx1 <- INLA::inla(y0 ~ f(dist, model="rw2", scale.model = TRUE,
                               hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx2 <- INLA::inla(y01 ~ f(dist, model="rw2", scale.model = TRUE,
                                hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx <- list(syncmodx1,syncmodx2)

plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="(III)")
abline(h=c(0,.25,.5,.75,1), lty=3, col="gray70")
points(seqx, yarR1, ylim=c(0,1), pch=19, col="red3")
points(seqx, yarR2, ylim=c(0,1), pch=19, col="deepskyblue4")
Axis(side=1, labels=FALSE, line=1, at=c(0,90))
for(i in 1:2)
{  
  syncmod <- syncmodx[[i]]
  lines(seqx, syncmod$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
        ylim=c(0,.7), main="", type="l", lwd=1.5, col=colsss[i])    
  min_a <- pmin(syncmod$summary.fitted.values$`0.025quant`)
  max_a <- pmax(syncmod$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      seqx,
      rev(seqx)
    ),
    c(max_a , rev(min_a)),
    col = scales::alpha(colsss[i],.5),
    border = NA
  )
}
# mtext("Distance",side=1, line=2, cex=1)

### Model IV ####
yarsR1 = (0.4)*exp(-0.007*seqx)+rnorm(length(seqx),0,0.025)
yarsR2 = (0.2)*exp(-0.007*seqx)+rnorm(length(seqx),0,0.025)
syncdata <- data.frame(dist=seqx, y0=yarsR1, y01=yarsR2)
syncmodx1 <- INLA::inla(y0 ~ f(dist, model="rw2", scale.model = TRUE,
                               hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx2 <- INLA::inla(y01 ~ f(dist, model="rw2", scale.model = TRUE,
                                hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)
syncmodx <- list(syncmodx1,syncmodx2)
plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="(IV)")
abline(h=c(0,.25,.5,.75,1), lty=3, col="gray70")
points(seqx, yarsR1, ylim=c(0,1), pch=19, col="red3")
points(seqx, yarsR2, ylim=c(0,1), pch=19, col="deepskyblue4")
Axis(side=1, labels=FALSE, line=1, at=c(0,90))
for(i in 1:2)
{  
  syncmod <- syncmodx[[i]]
  lines(seqx, syncmod$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
        ylim=c(0,.7), main="", type="l", lwd=1.5, col=colsss[i])    
  min_a <- pmin(syncmod$summary.fitted.values$`0.025quant`)
  max_a <- pmax(syncmod$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      seqx,
      rev(seqx)
    ),
    c(max_a , rev(min_a)),
    col = scales::alpha(colsss[i],.5),
    border = NA
  )
}
# mtext("Distance",side=1, line=2, cex=1)
# points(seqx,y2, col=2)
# points(seqx,y0, col=2, pch=19)
dev.off()
getwd()

### WEATHER ####
pdf("Plots/simulated_weathersync2.pdf", width = 7*1.2, height=2.7*1.2)
set.seed(10)
par(mfrow=c(1,4), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))
plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="")

met3 = .1+(0.7)*exp(-0.01*seqx)+rnorm(length(seqx),0,0.07)

syncdata <- data.frame(dist=seqx, y0=met3)
syncmet <- INLA::inla(y0 ~ f(dist, model="rw2", scale.model = TRUE,
                               hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                        , family="gaussian", control.predictor = list(compute=TRUE), data=syncdata)


plot(1, type="n", xlab="",ylab="Weather Synchrony", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="")
points(seqx, met3, ylim=c(0,1), pch=19, col="darkorchid4")
# abline(lm(met3~seqx),lty=2)
Axis(side=1, labels=FALSE, line=1, at=c(0,90))
Axis(side=2, labels=TRUE, line=1, at=c(0,0.5,1),cex.axis=1.5)

lines(seqx, syncmet$summary.fitted.values$mean, xlab="", ylab=expression(Ro), 
      ylim=c(0,.7), main="", type="l", lwd=2, col="darkorchid3")    
min_a <- pmin(syncmet$summary.fitted.values$`0.025quant`)
max_a <- pmax(syncmet$summary.fitted.values$`0.975quant`)
polygon(
  c(
    seqx,
    rev(seqx)
  ),
  c(max_a , rev(min_a)),
  col = scales::alpha("darkorchid2",.5),
  border = NA
)

dev.off()


pdf("Plots/simulated_meteo.pdf", width = 8.4, height=2.76)
set.seed(10)
par(mfrow=c(1,4), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))

plot(1, type="n", xlab="",ylab="Population Synchrony", ylim=c(0, .5), xlim=c(.3,.8), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="Spring")
points(met3, yarsR2, ylim=c(0,1), pch=19, col="deepskyblue4")
abline(lm(yarsR2~met3),lty=2)
Axis(side=1, labels=TRUE, line=1, at=c(0.3,.55,.8),cex.axis=1.1)
Axis(side=2, labels=TRUE, line=1, at=c(0,.25,0.5),cex.axis=1.1)


plot(1, type="n", xlab="",ylab="",ylim=c(0, .5), xlim=c(.3,.8), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="Fall")
points(met3, yarsR1, ylim=c(0,1), pch=19, col="red3")
abline(lm(yarsR1~met3),lty=2)
Axis(side=1, labels=TRUE, line=1, at=c(0.3,.55,.8),cex.axis=1.1)

plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="")

met3 = .1+(0.7)*exp(-0.01*seqx)+rnorm(length(seqx),0,0.07)
plot(1, type="n", xlab="",ylab="", ylim=c(0, 1), xlim=c(0,90), xaxt="n",
     yaxt="n", frame.plot = FALSE, main="")
points(seqx, met3, ylim=c(0,1), pch=19, col="darkorchid4")
# abline(lm(met3~seqx),lty=2)
Axis(side=2, labels=TRUE, line=1, at=c(0,.5,1),cex.axis=1.1)

Axis(side=1, labels=FALSE, line=1, at=c(0,90))

dev.off()




# old ####
# # simulated plots 
# plot(arima.sim(n = 25, list(ar = c(-0.2, .31), ma = c(-0.2279, 0.2488)),
#           sd = sqrt(3)))
# 
# plot(1:40,exp(jitter(cos((1:40)),1000))*3.5, lty=1, pch=19, type="b")
# 
# 
# plot(arima.sim(n = 40, list(ar = c(0, .0.5)), sd = sqrt(3)))
# 
# sim <- arima.sim(n = 40, list(ar = c(0.7, -0.4858), ma = c(0, 0)),
#           sd = sqrt(0.1))
# 
# sim2 <- arima.sim(n = 40, list(ar = c(0.6, -0.4858), ma = c(0, 0)),
#                  sd = sqrt(0.05))
# 
# 
# # 
# plot(ceiling(exp(sim)^4*cfact), ylim=c(0,15))
# lines(floor(exp(sim2)*3),col=2)
# 
# 
# plot(exp(sim)*4, ylim=c(0,12), xlim=c(0,20))
# lines(exp(sim)*4*cfact, lty=2, col=2)
# 
# 
# plot(arima.sim(n = 25, list(ar = c(-0.2, .31), ma = c(-0.2279, 0.2488)),
#                sd = sqrt(3)))
# 
# set.seed(123)
# reg1<- matrix(ncol=25,nrow=5)
# for(i in 1:5)
# {
#   reg1[i,] <- arima.sim(n = 25, list(ar = c(0.5, -0.5), ma = c(0, 0)),
#                         sd = sqrt(0.5))
# }
# 
# mean1 <- apply(reg1,2,mean)
# plot(exp(mean1)*4, type="l", ylim=c(0,25))
# 
# for(i in 1:5)
# {
#   lines(smooth.spline(exp(reg1[i,])*4,spar=.4), col=2)
# }
# smooth.spline()
# cfact <-runif(40,0.3,.9)
# 
# facts <- matrix(ncol=20,nrow=5)
# for(i in 1:5)
# {
#   facts[i,] <- sample(c(0.9,0.9,0.8,0.75,1.1,1.2,0.2,0.85,1.5,0.6,1.01,0.75), 20, replace=TRUE)
# }
# 
# 
# 
# postscript("Plots/arseries_gen.eps",width = 9*1.3, height=5.5)
# par(mfrow=c(1,3), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))
# 
## ## panel 1 ----
# reg1<-c(4.2,4.4,5.8,5.7,4,4.4,8,3.5,2,2.8,4,3.5,3.7,2.7,3.3,3.5,5.4,11.2,5.3,3)
# 
# reg1ind <- matrix(ncol=20,nrow=5)
# 
# plot(1, type="n", xlab="",ylab="", ylim=c(-5, 40), xlim=c(1,20), xaxt="n",
#      yaxt="n", frame.plot = FALSE)
# 
# 
# # plot(reg1*3, type="l", ylim=c(0,12*3))
# 
# regrandom <- c(1,2,-1,-2,0.5)
# for(i in 1:5)
# {
#   reg1ind[i,] <- jitter(reg1*3+regrandom[i],15)+2
#   lines(1:20,reg1ind[i,], col="gold", lty=2)
#   points(reg1ind[i,], col="gold", pch=19)
#   
# }
# 
# 
# reg2<-reg1*cfact[1:20]
# reg2ind <- matrix(ncol=20,nrow=5)
# #lines(reg2*3, type="l", ylim=c(3,12*3))
# regrandom <- c(1,2,-1,-2,0.5)
# 
# for(i in 1:5)
# {
#   reg2ind[i,] <- jitter(reg2*2.5+regrandom[i]-2,30)-3
#   lines(reg2ind[i,], col="seagreen", lty=2)
#   points(reg2ind[i,], pch=19, col="seagreen")
# }
# 
# Axis(side=1, labels=FALSE, line=1, at=c(1,20))
# Axis(side=2, labels=FALSE, line=1, at=c(-5,40))
## ## panel 2 -----
# reg3<-c(2.2,3.4,6.2,5.0,3,6,9,1,2,2,3.5,7,3.9,2.3,3.5,3.5,6,8.2,3.3,4.3)
# 
# reg3ind <- matrix(ncol=20,nrow=5)
# plot(1, type="n", xlab="",ylab="", ylim=c(-5, 40), xlim=c(1,20), xaxt="n",
#      yaxt="n", frame.plot = FALSE)
# regrandom <- c(1,2,-1,-2,0.5)
# for(i in 1:5)
# {
#   reg3ind[i,] <- jitter(reg3*3+regrandom[i],15)+2
#   lines(reg3ind[i,], col="#E35335", lty=2)
#   points(reg3ind[i,], col="#E35335", pch=19)
#   
# }
# 
# 
# # reg4<-reg3*sample(cfact[1:20],20,TRUE)
# regrandom <- c(1,2,-1,-2,0.5)
# reg4ind <- matrix(ncol=20,nrow=5)
# 
# for(i in 1:5)
# {
#   reg4ind[i,] <- jitter(reg4*1.9+regrandom[i]-3,50)
#   lines(reg4ind[i,], col="#66CC00", lty=2)
#   points(reg4ind[i,], col="#66CC00", pch=19)
#   
# }
# dev.off()
# Axis(side=1, labels=FALSE, line=1, at=c(1,20))
# ## ----
# 
# reg1mean <- apply(reg1ind,2,mean)
# reg2mean <- apply(reg2ind,2,mean)
# reg3mean <- apply(reg3ind,2,mean)
# reg4mean <- apply(reg4ind,2,mean)
# 
# # par(mfrow=c(1,1))
# plot(1, type="n", xlab="",ylab="", ylim=c(-3, 40), xlim=c(1,20), xaxt="n",
#      yaxt="n", frame.plot = FALSE)
# lines(reg1mean+1, lwd=2, col="gold")
# lines(reg2mean, lwd=2, col="seagreen")
# lines(reg3mean+1, lwd=2,col="#E35335")
# lines(reg4mean, lwd=2, col="#66CC00")
# 
# Axis(side=1, labels=FALSE, line=1, at=c(1,20))
# Axis(side=2, labels=FALSE, line=1, at=c(-3,40))
# dev.off()
# 
# mtext("Time",side=1, line=2, cex=1.5)
# mtext("Population size",side=2, line=2, cex=1.5)
# 
# ## weather sequence ----
# pdf("Plots/weather_sim.pdf", height=4,width=8)
# plot(1, type="n", xlab="",ylab="", xlim=c(1, 20), ylim=c(-5,10), xaxt="n",
#      yaxt="n", frame.plot = FALSE)
# weathrz<-c(1,2,5,5,-3,8,5,-2,0,1,5,7,2,0,4,7,6,2,-5,-4)
# weathrx<-c(1.5,0,4.2,7.2,-2.2,8.2,5.2,2.2,0.2,0.3,2.3,1.2,2.1,0.5,2.7,4.3,6.1,2.3,0,-3)
# #1,2,5,8,9,10,13,14,19,20
# points(weathrz, pch=19, col="navy")
# lines(weathrz, pch=19, col="navy")
# lines(weathrx, pch=19, col=4)
# points(weathrx, pch=19, col=4)
# Axis(side=1, labels=FALSE, line=1, at=c(1,20))
# Axis(side=2, labels=FALSE, line=1, at=c(-5,10))
# dev.off()
