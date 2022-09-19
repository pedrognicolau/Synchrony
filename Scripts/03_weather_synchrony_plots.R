res4synch <- readRDS("Data/popmodels_residuals.rds")
icedf2 <- readRDS("Data/weather_vars.rds")
dm1 <- as.matrix(read.csv2("Data/distance_matrix.csv")[,-1])
gps <- read.csv2("Data/gps.csv")

melt_covm <- function(cor_mat, colname = ""){
  CR2 <- cor_mat
  CR2[upper.tri(CR2)] <- NA
  diag(CR2) <- NA
  nCR <- reshape2::melt(CR2, varnames = c('station1', 'station2'), na.rm = TRUE)
  colnames(nCR)[3] <- colname
  return(nCR)
}

WRain <- matrix(icedf2$winter_rain,ncol=19)
ZCw <- matrix(icedf2$zerocrosses,ncol=19)

# Melted Correlation Matrices for Weather 
ZCor <- melt_covm(cor(ZCw), colname="cor")
WRCor <- melt_covm(cor(WRain), colname="cor")

# Melted Correlation Matrices for residuals
SCorIV <- melt_covm(cor(res4synch$St_IV), colname="cor") 
FCorIV <- melt_covm(cor(res4synch$Ft_IV), colname="cor") # from model synchrony script
#SCorGR <- melt_covm(cor(res4synch$gr_spring), colname="cor") # from model synchrony script
#FCorGR <- melt_covm(cor(res4synch$gr_fall), colname="cor") # from model synchrony script


res_list <- list(SCorIV, FCorIV)
reslab <- c("Spring","Fall")

syncdata <- SCorIV
colnames(syncdata)[3] <- "SpringIV"
syncdata2 <- left_join(syncdata,FCorIV)
colnames(syncdata2)[4] <- "FallIV"
syncdata3 <- left_join(syncdata2, WRCor)
colnames(syncdata3)[5] <- "WinterRainfall"
syncdata4 <- left_join(syncdata3, ZCor)
colnames(syncdata4)[6] <- "ZeroCrosses"

syncdata9 <- arrange(syncdata4, SpringIV)
library(INLA)
mod_fall_zc <- inla(FallIV~ZeroCrosses,data=syncdata4,control.predictor = list(compute=TRUE),
                    control.compute = list(hyperpar=TRUE, config=TRUE))

mod_fall_wr <- inla(FallIV~WinterRainfall,data=syncdata4,control.predictor = list(compute=TRUE),
                    control.compute = list(hyperpar=TRUE, config=TRUE))

mod_spring_zc <- inla(SpringIV~ZeroCrosses,data=syncdata4, control.predictor = list(compute=TRUE),
                      control.compute = list(hyperpar=TRUE, config=TRUE))

mod_spring_wr <- inla(SpringIV~WinterRainfall,data=syncdata4,control.predictor = list(compute=TRUE),
                      control.compute = list(hyperpar=TRUE, config=TRUE))


modsweather <- list(spring_zc=mod_spring_zc,fall_zc=mod_fall_zc,
                    spring_wr=mod_spring_wr,fall_wr=mod_fall_wr
                    #spring_wt=mod_spring_wt,fall_st=mod_fall_st,
                    #spring_wp=mod_spring_wp,fall_sp=mod_fall_sp
                    )

# smooth -----
zcdf <- sync_df5(ZCw)
wrdf <- sync_df5(WRain)
syncmet1 <- INLA::inla(corres ~ f(dist, model="rw2", scale.model = TRUE,
                                  hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                       , family="gaussian", control.predictor = list(compute=TRUE), data=zcdf)

syncmet2 <- INLA::inla(corres ~ f(dist, model="rw2", scale.model = TRUE,
                                  hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                       , family="gaussian", control.predictor = list(compute=TRUE), data=wrdf)

par(mfrow=c(2,1))


## PLOTS ----------

varlist <- list(syncdata4$ZeroCrosses,syncdata4$ZeroCrosses,
                syncdata4$WinterRainfall,syncdata4$WinterRainfall)
seaslist <- list(syncdata4$SpringIV,syncdata4$FallIV,
                 syncdata4$SpringIV,syncdata4$FallIV)
plotlab <- c("Zero Cross Synchrony", "Zero Cross Synchrony",
             "Winter Rainfall Synchrony","Winter Rainfall Synchrony")
colsplot <- c(4,2,4,2)
mains=c("b) Spring","c) Fall","","")

# plot(sync_df5(WRain))
#plot(sync_df5(ZCw), xlab="Zero Cross Synchrony", ylab="Population synchrony",
#     pch=19, col=scales::alpha("gray30",.7))
pdf("Plots/weather_effects.pdf", width = 7.9, height= 6.2)
par(mfrow=c(2,3), mar=c(6, 2, 2, 2) + 0.1 , mai=c(.8, .6, .4, .1),oma=c(.5, 2, 2, 1.5))
for(i in 1:4)
{
  # zero cross synchrony
  if(i==1) {
    
    plot(1, type="n", xlab="",ylab="", ylim=c(-.5, 1), xlim=c(0,150), xaxt="n",
         yaxt="n", frame.plot = FALSE, main="")
    points(zcdf$dist, zcdf$corres, ylim=c(0,1), pch=19, col="gray30")
    # abline(lm(met3~seqx),lty=2)
    Axis(side=1, labels=TRUE, line=1, at=c(0,50, 100, 150))
    Axis(side=2, labels=TRUE, line=1, at=c(-.5,0,0.5,1))
    title(main = "a) Weather", cex.lab = 2,
          line = 1)
    title(ylab = "Zero Crosses Synchrony", cex.lab = 1.2,
          line = 3.5)
    title(xlab = "Distance (km)", cex.lab = 1.2, line = 3.5)
    
    
    lines(zcdf$dist, syncmet1$summary.fitted.values$mean, xlab="", ylab=expression(Ro), 
          ylim=c(0,.7), main="", type="l", lwd=2, col="gray10")    
    min_a <- pmin(syncmet1$summary.fitted.values$`0.025quant`)
    max_a0 <- pmax(syncmet1$summary.fitted.values$`0.975quant`)
    max_a <- ifelse(max_a0>1,1,max_a0)
    polygon(
      c(
        zcdf$dist,
        rev(zcdf$dist)
      ),
      c(max_a , rev(min_a)),
      col = scales::alpha("darkslategray2",.5),
      border = NA
    )
    
    # plot(sync_df5(ZCw), xlab="distance", ylab="Zero Cross Synchrony", pch=19, ylim=c(-.2,1), main="a) Weather Synchrony")
  }
  # winter rainfall synch
  if(i==3) 
    {
    
    plot(1, type="n", xlab="",ylab="", ylim=c(-.5, 1), xlim=c(0,150), xaxt="n",
         yaxt="n", frame.plot = FALSE, main="")
    points(wrdf$dist, wrdf$corres, ylim=c(0,1), pch=19, col="gray30")
    # abline(lm(met3~seqx),lty=2)
    Axis(side=1, labels=TRUE, line=1, at=c(0,50, 100, 150))
    Axis(side=2, labels=TRUE, line=1, at=c(-.5,0,0.5,1))
    title(ylab = "Winter Rainfall Synchrony", cex.lab = 1.2,
          line = 3.5)
    title(xlab = "Distance (km)", cex.lab = 1.2, line = 3.5)
    
    lines(zcdf$dist, syncmet2$summary.fitted.values$mean, xlab="", ylab=expression(Ro), 
          ylim=c(0,.7), main="", type="l", lwd=2, col="gray10")    
    min_a <- pmin(syncmet2$summary.fitted.values$`0.025quant`)
    max_a0 <- pmax(syncmet2$summary.fitted.values$`0.975quant`)
    max_a <- ifelse(max_a0>1,1,max_a0)    
    polygon(
      c(
        zcdf$dist,
        rev(zcdf$dist)
      ),
      c(max_a , rev(min_a)),
      col = scales::alpha("darkslategray2",.5),
      border = NA
    )
    
    
    
  #plot(sync_df5(WRain), xlab="distance", ylab="Winter Rainfall Synchrony", pch=19, ylim=c(-.2,1))
    }
  plot(varlist[[i]], seaslist[[i]], xlab="", ylab=" ",
       pch=19, col=scales::alpha("gray60",.7),ylim=c(-.2,.8), main=)
  title(main = mains[i], cex.lab = 3.5,
        line = 1)
  title(ylab = "Population Synchrony", cex.lab = 1.2,
        line = 2.5)
  title(xlab = plotlab[i], cex.lab = 1.2, line = 3.5)
  lines(sort(varlist[[i]]),sort(modsweather[[i]]$summary.fitted.values$mean), lty=1, col=colsplot[i], lwd=1.5)
  lines(sort(varlist[[i]]),sort(modsweather[[i]]$summary.fitted.values$`0.025quant`), lty=2, col=colsplot[i], lwd=1)
  lines(sort(varlist[[i]]),sort(modsweather[[i]]$summary.fitted.values$`0.975quant`), lty=2, col=colsplot[i], lwd=1)
}
dev.off()