# Estimate synchrony -------
library(INLA)
library(dplyr)
library(rlang)

## Read in data ---------------
source("Scripts/0_Important_functions.R")
popres <- readRDS("Data/popmodels_residuals.rds")
stackeddata <- readRDS("Data/data4synchrony.rds")

str(popres)
synclist_p = synclist_s = list()

## Compute correlation matrices ---------
for(i in 1:length(popres))
{ # St_seasonal, Ft_seasonal, St_year, Ft_year
  synclist_p[[i]] <- sync_df5(popres[[i]], method="pearson")
  # synclist_s[[i]] <- sync_df5(popres[[i]], method="spearman")
}

### Obtain raw correlations as a function of distance -----
Springcor <- sync_df5(matrix(stackeddata$St,ncol=19)) # raw
Fallcor <- sync_df5(matrix(stackeddata$Ft,ncol=19))

# append raw correlation synchrony
synclist <- append(synclist_p,list(Springcor,Fallcor))

## Model correlation in residuals as a smooth function of distance ----

modsynclist <- list()

par(mfrow=c(2,3))
i=1
for(i in 1:length(synclist)){
  dat <- synclist[[i]]
  
  colnames(dat) <- c("distance","mean.res")
  syncmod1 <- inla(mean.res ~ f(distance, model="rw2", scale.model = TRUE,
                                hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.1))))
                   , family="gaussian", control.predictor = list(compute=TRUE)
                   ,data=dat)
  msum0 <- syncmod1$summary.random$distance
  msum01 <- arrange(msum0,as.numeric(ID))
  modsynclist[[i]] <- syncmod1
  
  plot(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
       ylim=c(0,.7), main="", type="l", lwd=1.5, col=i)    
  lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.025quant`, lty=2)    
  lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.975quant`, lty=2)   
  
  
}
synclist
par(mfrow=c(2,3))
for(i in 1:6) hist(synclist[[i]][,2], breaks=seq(-1,1,0.15))
for(i in 1:6) print(shapiro.test(synclist[[i]][,2]))

plot(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
     ylim=c(0,.7), main="", type="l", lwd=1.5, col=i)    
lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.025quant`, lty=2)    
lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.975quant`, lty=2)   

### Plots of synchrony with distance --------
# spring seas, spring yr, Fall seas, Fall Yr, Spring raw, Fall Raw
pdf("Plots/scale_raw_year_seas2.pdf", width = 7.7, height= 5.4)
j=1 # parameter for the legend
mainlabel <- c("(a)", "(b)", "(c)")
par(mfrow=c(1,3), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))
i=3
for (i in c(5,3,1))
{
  
  # initiate plot
  plot(0, 0, col="white", xlim=c(0,150), ylim=c(0,.85), xlab="Distance (km)",
       ylab="",
       axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
  # if it is the first, add y lab
  if(i==5) mtext(side=2, line=3, "Synchrony", col="black", font=1, cex=1.2)
  
  #  for better comparison
  abline(h=c(seq(.1,.8,.15)), lty=2, col="gray80", cex=.6)
  
  Axis(side=1, labels=TRUE, cex.axis = .9)
  
  # add axis on the left
  if(i==5) Axis(side=2, labels=TRUE, cex.axis = .9)
  
  # plot 95% confidence intervals
  lines(
    modsynclist[[i]]$summary.random$dist$ID,
    modsynclist[[i]]$summary.fitted.values$mean,
    xlab = "Distance (km)",
    ylab = "Correlation",
    ylim = c(0, .8),
    type = "l",
    lwd = 3,
    col = 1,
    main = "",
    lty = 1
  )
  min_a <- pmin(modsynclist[[i]]$summary.fitted.values$`0.025quant`)
  max_a <- pmax(modsynclist[[i]]$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      modsynclist[[i]]$summary.random$dist$ID,
      rev(modsynclist[[i]]$summary.random$dist$ID)
    ),
    c(max_a , rev(min_a)),
    col = rgb(0, .4, 0, 0.5),
    border = NULL
  )

  
  
  min_a <- pmin(modsynclist[[i + 1]]$summary.fitted.values$`0.025quant`)
  max_a <- pmax(modsynclist[[i + 1]]$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      modsynclist[[i + 1]]$summary.random$dist$ID,
      rev(modsynclist[[1]]$summary.random$dist$ID)
    ),
    c(max_a , rev(min_a)),
    col = rgb(1, .5, 0, 0.5),
    border = NULL
  )
  
  lines(
    modsynclist[[i]]$summary.random$dist$ID,
    modsynclist[[i]]$summary.fitted.values$mean,
    col = "darkgreen",
    lty = 1,
    lwd = 3
  )
  
  lines(
    modsynclist[[i + 1]]$summary.random$dist$ID,
    modsynclist[[i + 1]]$summary.fitted.values$mean,
    col = "red4",
    lty = 1,
    lwd = 3
  )
  
  legend("top", legend=c(mainlabel[j]), bty="n", text.font = 1, cex=1.5)
  j=j+1
}

dev.off()

