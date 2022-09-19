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
str(popres)
## Compute correlation matrices ---------
for(i in 1:length(popres))
{ # St_seasonal, Ft_seasonal, St_year, Ft_year
  synclist_p[[i]] <- sync_df5(popres[[i]], method="pearson")
  # synclist_s[[i]] <- sync_df5(popres[[i]], method="spearman")
}

### Obtain raw correlations as a function of distance (I) -----
Springcor <- sync_df5(matrix(stackeddata$St,ncol=19)) # raw
Fallcor <- sync_df5(matrix(stackeddata$Ft,ncol=19))
# append raw correlation synchrony
synclist <- append(list(Springcor,Fallcor),synclist_p)

## Model correlation in residuals as a smooth function of distance ----

modsynclist <- list()

par(mfrow=c(2,4))
for(i in 1:length(synclist)){
  dat <- synclist[[i]]
  
  colnames(dat) <- c("distance","mean.res")
  syncmod1 <- inla(mean.res ~ f(distance, model="rw2", scale.model = TRUE,
                                hyper = list(theta = list(prior="pc.prec", param=c(u=.5,0.01))))
                   , family="gaussian", control.predictor = list(compute=TRUE)
                   ,data=dat)
  msum0 <- syncmod1$summary.random$distance
  msum01 <- arrange(msum0,as.numeric(ID))
  modsynclist[[i]] <- syncmod1
  
  lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
       ylim=c(0,.7), main="", type="l", lwd=1.5, col=2)    
  lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.025quant`, lty=2)    
  lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.975quant`, lty=2)   
  
  dat <- synclist[[i]]
  plot(synclist[[i]]$dist, synclist[[i]]$corres, xlab="distance", col=4, pch=19, ylim=c(-.2,1))    
  points(synclist[[i]]$dist, synclist[[i]]$corres, xlab="distance", col=2, pch=19)    
}
#synclist
# par(mfrow=c(2,4))
# for(i in 1:6) hist(synclist[[i]][,2], breaks=seq(-1,1,0.15))
# for(i in 1:6) print(shapiro.test(synclist[[i]][,2]))
# 
# plot(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$mean, xlab="distance", ylab=expression(Ro), 
#      ylim=c(0,.7), main="", type="l", lwd=1.5, col=i)    
# lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.025quant`, lty=2)    
# lines(syncmod1$summary.random$dist$ID, syncmod1$summary.fitted.values$`0.975quant`, lty=2)   

### Plots of synchrony with distance --------
# spring seas, spring yr, Fall seas, Fall Yr, Spring raw, Fall Raw
pdf("Plots/scale_sync_4models_overlayed2.pdf", width = 7.7, height= 5.4)
j=1 # parameter for the legend
mainlabel <- c("(I)", "(II)", "(III)","(IV)")
par(mfrow=c(1,4), mar=c(5, 1, 2, 2) + 0.1 , mai=c(.6, .1, .2, .1),oma=c(.5, 4, 1.5, 1.5))
for (i in c(1,3,5,7))
{
  # initiate plot
  plot(0, 0, col="white", xlim=c(0,150), ylim=c(0,1), xlab="Distance (km)",
       ylab="",
       axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
  # if it is the first, add y lab
  if(i==1) mtext(side=2, line=3, "Synchrony", col="black", font=1, cex=1.2)
  
  #  for better comparison
  abline(h=c(seq(.2,1,.2)), lty=2, col="gray50", cex=.6)
  
  Axis(side=1, labels=TRUE, cex.axis = .9, at=c(0,50,100,150))
  
  # add axis on the left
  if(i==1) Axis(side=2, labels=TRUE, cex.axis = .9)
  
  # plot dots 
  points(
    synclist[[i]]$dist, synclist[[i]]$corres, col = scales::alpha("blue3",.8), pch=19
  )
  
  points(
    synclist[[i+1]]$dist, synclist[[i+1]]$corres, col = scales::alpha("red4",.7), pch=19,
  )
  
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
    col = scales::alpha(4,.7),
    border = NA
  )

  
  
  min_a <- pmin(modsynclist[[i + 1]]$summary.fitted.values$`0.025quant`)
  max_a <- pmax(modsynclist[[i + 1]]$summary.fitted.values$`0.975quant`)
  polygon(
    c(
      modsynclist[[i + 1]]$summary.random$dist$ID,
      rev(modsynclist[[1]]$summary.random$dist$ID)
    ),
    c(max_a , rev(min_a)),
    col = scales::alpha(2,.6),
    border = NA
  )
  
  lines(
    modsynclist[[i]]$summary.random$dist$ID,
    modsynclist[[i]]$summary.fitted.values$mean,
    col = "blue4",
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



