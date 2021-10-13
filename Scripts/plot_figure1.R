crdata <- readRDS("/Users/pni026/OneDrive - UiT Office 365/Vole Synchrony/Data/Abundance_CRINLA_zeroc3.rds")





datedf <- data.frame(timepoint=1:42,year=rep(2000:2020,each=2), date=dateseq)
stationdf <- data.frame(station=factor(c(1:8,10:20)),group=c(rep(1,5),rep(2,7),rep(3,7)))

crdata2 <- left_join(crdata,datedf)
crdata3 <- left_join(crdata2,stationdf)
yrcounts <- aggregate(cbind(inlacr,counts)~year+station,data=crdata2, median)
regcounts <- aggregate(cbind(inlacr,counts)~year+group,data=crdata3, median)



# plot - median abundance per group with year ----
pdf("Plots/regional_abundances.pdf", width=12.2, height=8)
{
  par(mfrow=c(1,1), mar=c(6, 2, 0, 1) + 0.1) # B, L, T, R
  newcol=c(4,3,7)
  plot(0, 0, col="white", xlim=c(1999,2021), ylim=c(-0.25,50), ylab="",
       xlab="",
       axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
  
  for(i in 1:length(stations))
  {
    plotdata <- filter(yrcounts, station==stations[i])
    #points(plotdata$year,plotdata$inlacr, col=scales::alpha(colvec[i],.2), pch=19)
    
    #lines(plotdata$year,plotdata$inlacr, col=scales::alpha(colvec[i],.7), lwd=2)
    lines(smooth.spline(plotdata$year,plotdata$inlacr, spar=.01), col=scales::alpha(colvec[i],.7), lwd=1, lty=2)
    
    
  }
  
  for(i in 1:3)
  {
    plotdata <- filter(regcounts, group==i)
    #points(plotdata$year,plotdata$inlacr, col=scales::alpha(newcol[i],.5), pch=19)
    
    #lines(plotdata$year,plotdata$inlacr, col=scales::alpha(colvec[i],.7), lwd=2)
    lines(plotdata$year,plotdata$inlacr, col=scales::alpha(newcol[i],.9), lwd=4)
    
    
  }
  
  axis(1, cex.axis = 1.2, pos=-3)
  mtext("Abundance", side=2, line=0.5, cex=1.5)
}
dev.off()

# other related plots -----

par(mfrow=c(1,1))
stations=sort(unique(yrcounts$station))
colvec <- c(rep(4,5),rep(3,7),rep(7,7))
plot(0, 0, col="white", xlim=c(1999,2021), ylim=c(-2,100), ylab="Abundance",
     xlab="",
     axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
for(i in 1:length(stations))
{
  plotdata <- filter(yrcounts, station==stations[i])
  #points(plotdata$year,plotdata$inlacr, col=scales::alpha(colvec[i],.2), pch=19)
  
  #lines(plotdata$year,plotdata$inlacr, col=scales::alpha(colvec[i],.7), lwd=2)
  lines(smooth.spline(plotdata$year,plotdata$inlacr, spar=.3), col=scales::alpha(colvec[i],.4), lwd=1)
  
  
}

axis(1, cex.axis = 1, pos=-2)


date1 <- sort(c(seq(as.Date("2000-05-01"),as.Date("2020-05-01"),"year")))


dateseq <- sort(c(seq(as.Date("2000-05-01"),as.Date("2020-05-01"),"year"),
                  seq(as.Date("2000-08-01"),as.Date("2020-08-01"),"year")))


plot(0, 0, col="white", xlim=c(min(dateseq),max(dateseq)), ylim=c(0,40), ylab="Abundance",
     xlab="",
     axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
for(i in 1:length(stations))
{
  plotdata <- filter(crdata2, station==stations[i])
  points(plotdata$date,plotdata$counts, col=scales::alpha(colvec[i],.2), pch=19)
  # lines(plotdata$year,plotdata$inlacr, col=colvec[i], pch=19)
  
  lines(smooth.spline(plotdata$date,plotdata$counts, spar=.4), col=scales::alpha(colvec[i],.7), lwd=2)
  
  
}
axis(1, dateseq[seq(1,42,4)], format(dateseq[seq(1,42,4)], "%Y"), cex.axis = 1, pos=-2)
i=1
