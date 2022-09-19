# read in data 
library(dplyr)
library(lubridate)
# read in climate files ------
porpath <- "Data/pors klima/"
dir_pors <- sort(dir(porpath))
for(i in seq(1,length(dir_pors),2))
{

  station_weather_prec <- as.character(dir_pors[i])
  station_weather_temp <- as.character(dir_pors[i+1])
  
  dat <- read.csv(paste0(porpath,dir_pors[i]), skip=1, sep=";")
  dat2 <- read.csv(paste0(porpath,dir_pors[i+1]), skip=1, sep=";")
  
  
  stname <- strsplit(station_weather_prec, c(" "))[[1]]
  stname2 <- paste0(stname[1:(length(stname)-1)], collapse=" ")

  #variable <- strsplit(stname[length(stname)],"[.]")[[1]][1]
  #dat$var <- variable
  dat$temp_Celsius <- dat2$tm.Celcius.
  dat$station <- stname2
  if(i == 1) weather_data <- dat
  else weather_data <- rbind(weather_data,dat)
}


## minor processings -----
weather_data$datetime=strptime(weather_data$Date,"%d.%m.%Y")
weather_data$year<-lubridate::year(weather_data$datetime)
weather_data$month<-lubridate::month(weather_data$datetime)
names(weather_data) <- c("FullDate", "Prec","Temp","station","date","year","month")
weather_data$date.char <- as.character(weather_data$date)
## filter dates more than 1999 except december with is winter 2000
wdata99 <- filter(weather_data, year > 1999 | (year > 1998 & month > 11)) # we need december 1999 to compute the winter weather

# obtain total precipitation and mean daily temperature
wdata_agg <- do.call(data.frame,aggregate(cbind(Prec,Temp)~date.char+year+month+station, 
                                        FUN=function(x) c(mn = mean(x), sum = sum(x) ),data=wdata99))
wdata_agg2 <- wdata_agg[,c(1,2,3,4,6,7)]

wdata3 <- arrange(wdata_agg2, station,date.char)

# correct obviously wrong temperature values and replace with mean in adjacent days
wrongtemp <- which(wdata3$Temp.mn>50)
for (j in 1:length(wrongtemp))
{
  # replace with the mean in the adjacent days
  pos <- wrongtemp[j]
  wdata3$Temp.mn[pos] <- mean(c(wdata3$Temp.mn[pos-1],wdata3$Temp.mn[pos+1]))
}

## getting seasonal data ----
### aggregate over month and season 
seas_df <- data.frame(day=sort(unique(substr(wdata3$date.char, 6, 10))))
### label days according to season
seas_df$season <- c(rep("WINTER",81),rep("SPRING",91),rep("SUMMER",91),rep("FALL",92), rep("WINTER",11))
getseasn <- function(x) {
  daten <- substr(x, 6, 10)
  seas <- seas_df$season[seas_df$day==daten]
  return(seas)
}
### label stations
stdf <- data.frame(station=unique(wdata3$station), station2=c(1,10:19,2,20,3:8,21:24))
wdata3$season <- sapply(wdata3$date.char, FUN=getseasn)

### count number of winters ----
wdatax <- distinct(wdata3,date.char,year,month,season)

winternumber <- c()
count=0
addcounter=TRUE
for(i in 1:nrow(wdatax))
{
  if(i%%1000==0) print(i)
  if(wdatax$season[i]=="WINTER" & addcounter==TRUE) count = count + 1 ; addcounter = FALSE
  if(wdatax$season[i]!="WINTER") addcounter = TRUE
  winternumber <- c(winternumber, count)
}
wdatax$winternumber<- winternumber
wdatax$winternumber <- ifelse(wdatax$season!="WINTER",0,wdatax$winternumber)

wdata3.1 <- left_join(wdata3,stdf)
wdata3.2 <- left_join(wdata3.1, wdatax)
head(wdata3.2)

rm(wdata99,weather_data,wdata_agg,wdata_agg2, dir_pors,dat,dat2,stname,stname2)

# Plots yearly averages -----
saveRDS(wdata3.2,"Vole Synchrony/Data/temp_and_rain_1999_2020.rds")
weather_year <- wdata3.2
std <- function(x) sd(x)/sqrt(length(x))
wy1 <- filter(weather_year, year > 1999 & !(station%in%c("T1 S2","T2 S1","T3 S1","T5 S2")))
wy2_t <- aggregate(Temp.mn~station2, data=wy1,FUN=mean)
plot(wy2_t, pch=19, col=c(rep(4,5),rep(3,7),rep(7,7)))
wy2_p <- aggregate(Prec.sum~station2, data=wy1,FUN=mean)
plot(wy2_p, col=c(rep(4,5),rep(3,7),rep(7,7)), pch=19, ylab="Average Daily Precipitation")

wy2_tsd <- aggregate(Temp.mn~station2, data=wy1,FUN=std)
wy2_psd <- aggregate(Prec.sum~station2, data=wy1,FUN=std)
abs0 <- c(rep(1,14),rep(-1,5))
library(plotrix)

pdf("Vole Synchrony/Plots/weather_gradient_x.pdf", width = 8, height = 6)
par(mfrow=c(1,2), c(5, 4, 0, 1) + 0.1) # B, L, T, R
plotCI(gps$cumdist,y= wy2_t$Temp.mn,wy2_tsd$Temp.mn*1.96,wy2_tsd$Temp.mn*1.96, 
       col=c(rep(4,5),rep(3,7),rep(7,7)), pch=19, ylab="Average Daily Temperature (°C)",
       xlab="Distance (km)")
plotCI(gps$cumdist,y= wy2_p$Prec.sum,wy2_psd$Prec.sum*1.96,wy2_psd$Prec.sum*1.96, col=c(rep(4,5),rep(3,7),rep(7,7)), pch=19, ylim=c(1,2.3),
      ylab="Average Daily Precipitation (mm)", xlab="Distance (km)")
dev.off()

## temperature and precipitation per winter plot ----
winterdata <- filter(wdata3.2, winternumber>0)
temperature_per_month_station <- arrange(aggregate(Temp.mn~station+winternumber,data=winterdata, FUN=mean),station,winternumber)

stdflab <- arrange(data.frame(st=unique(wdata3$station), station=c(1,10:19,2,20,3:8,21:24)),station)
i=1

pdf("Plots/winter_patterns.pdf", height=9,width=7.7)
truecols <- c(rep(4,5),rep(3,8),rep(7,7))
par(mfrow=c(2,1), mar=c(2, 4, 2, 1) + 0.1) # B, L, T, R

for(i in stdflab$station[1:19])
{
  plot_data <- filter(temperature_per_month_station, station == stdflab$st[i])
  
  if (i == 1) plot(2000:2020,plot_data$Temp.mn, type="p", col=scales::alpha(truecols[i], .7),  ylim=c(-20,0), pch=19, ylab="Average Winter Temperature  (°C)")
  else points(2000:2020,plot_data$Temp.mn, type="b", pch=19, lty=1, col=scales::alpha(truecols[i], .7))
}

prec_per_year_station <- arrange(aggregate(Prec.sum~station+winternumber,data=winterdata, FUN=sum),station,winternumber)

stdflab <- arrange(data.frame(st=unique(wdata3$station), station=c(1,10:19,2,20,3:8,21:24)),station)

truecols <- c(rep(4,5),rep(3,8),rep(7,7))

for(i in stdflab$station[1:19])
{
  plot_data <- filter(prec_per_year_station, station == stdflab$st[i])
  
  if (i == 1) plot(2000:2020,plot_data$Prec.sum, type="p", col=scales::alpha(truecols[i], .7),  ylim=c(0,450), pch=19,
                   ylab="Total Winter Precipitation (mm)")
  else points(2000:2020,plot_data$Prec.sum, type="b", pch=19, lty=1, col=scales::alpha(truecols[i], .7))
}
dev.off()

# zero-crossings & winter rain ------

# function to compute zero cross days to use in aggregate
zcross <- function(x) sum(abs(diff(sign(x))) > 1) # sum of amount of time the difference between consecutive states is larger than abs(2)
# example -1 - (1) = -2 or 1 - (-1) = 2. It will not count -1 - 0 as a crossing

# just winter data
wdata3.3 <- filter(wdata3.2, season == "WINTER")
# how many zero crosses per winter season and station
zcrosses <- aggregate(Temp.mn~winternumber+station, FUN=zcross, data=wdata3.2)
names(zcrosses)[3] <- "zerocrosses"
# filter the zero crosses outside of winter
zcrosses2 <- filter(zcrosses, winternumber > 0)

# which days are positive and which are negative degrees celsius
wdata3.2$pos_temp <- sign(wdata3.2$Temp.mn)
# filter by data only in winter
wdataW <- filter(wdata3.2, season=="WINTER")
# remove rain on frozen days
wdataW$winter_rain <- ifelse(wdataW$pos_temp<0,0,wdataW$Prec.sum)
# get the sum of precipitation on negative days for the winter season per station
WPosPrec <- aggregate(winter_rain~winternumber+station,data=wdataW, FUN=sum)
# add that data to df with zero cross number
WRZCdf <- left_join(zcrosses2, WPosPrec)
# add trapping year corresponding to winter season
yrdf <- data.frame(year=2000:2020, winternumber=1:21)
icedf <- left_join(yrdf,WRZCdf)
icedf2 <- left_join(icedf,stdf)
icedf3 <- arrange(icedf2, station2, year)
icedf4 <- icedf3[,c(1,6,4,5,3)]
head(icedf4)
names(icedf4)[c(2,5)] <- c("station","station_name")
saveRDS(icedf4,"Data/winter_precipitation.rds")

i=1

stations <- unique(icedf4$station)
# FIGURE 2 ------
## zero  crosses -----
pdf("Plots/zerocrosses.pdf", height=4,width=6.5)
truecols <- c(rep(4,5),rep(3,8),rep(7,7))
par(mfrow=c(1,1), mar=c(4, 4, 2, 1) + 0.1) # B, L, T, R
plot(0, 0, col="white", xlim=c(1999.5,2020.5), xlab="", ylim=c(-1,30),
     ylab="", axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
abline(h=seq(5,30,5),lty=3, col="gray80")

set.seed(5)
for(i in stations[1:19])
{
  plot_data <- filter(icedf4, station == stations[i])
  
  points(2000:2020, jitter(plot_data$zerocrosses,1)+.3, type="b", col=scales::alpha(truecols[i], .5), pch=19, lty=1,
         ylab="Number of Winter Zero Crosses")
}
Axis(side=1, labels=TRUE, cex.axis = 1, pos=-2)
Axis(side=2, labels=TRUE, cex.axis = 1, at=c(0,10,20,30))
dev.off()

## winter rain ------
pdf("Plots/winterrain.pdf", height=4,width=6.5)
truecols <- c(rep(4,5),rep(3,8),rep(7,7))
par(mfrow=c(1,1), mar=c(4, 4, 2, 1) + 0.1) # B, L, T, R
summary(icedf4)
plot(0, 0, col="white", xlim=c(1999.5,2020.5), xlab="", ylim=c(-2,150),
     ylab="", axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
abline(h=seq(0,150,25),lty=3, col="gray80")
for(i in stations[1:19])
{
  plot_data <- filter(icedf4, station == stations[i])
  
  points(2000:2020, jitter(plot_data$winter_rain,1)+.5, type="b", col=scales::alpha(truecols[i], .5),  
         ylim=c(0,200), pch=19, lty=1,
                   ylab="Total Winter Rainfall (mm)")
}
Axis(side=1, labels=TRUE, cex.axis = 1, pos=-10)
Axis(side=2, labels=TRUE, cex.axis = 1, at=c(0,50,100,150))
dev.off()


######## join abundance and weather data ########

Wdata <- icedf4
Wdata2 <- filter(Wdata, !(station %in% c("T1 S2","T2 S1","T3 S1", "T5 S2")))
nrow(Wdata2)
# sstation <- function(x) as.numeric(strsplit(x,"S")[[1]][2])
saveRDS(Wdata2,"Data/weather_vars.rds")


GR_data0 <- readRDS("Vole Synchrony/data/log_centered_INLACR_growthrates.rds")
GR_data0$station <- ifelse(GR_data0$station>8,GR_data0$station+1,seasons_weather9$station)
unique(GR_data0$station)
sort(unique(seasons_weather9$station))
GRweather <- left_join(GR_data0, seasons_weather9)  

GRweather2 <- select(GRweather, -c(winternumber))
saveRDS(GRweather2, "Vole Synchrony/data/abundance_and_weather_data.rds")
s10 <- filter(seasons_weather9, station < 21)
saveRDS(s10, "Vole Synchrony/data/seasonal_weather.rds")

summary(GR_weather)
