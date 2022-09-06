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

# correct wrong temperature values and replace with mean in adjacent days
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
gps <- read.csv2("Vole synchrony/Data/gps.csv")
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


# ## abundances per season ----
# pdf("Plots/seasonal_abundances.pdf", height=4,width=7.3)
# abundanceCR <- readRDS("Data/Abundance_CRINLA_zeroc3.rds")
# abundanceCR$season <- rep(c(0,1),nrow(abundanceCR)/2)
# abundanceCR$station <- as.numeric(abundanceCR$station)
# stations <- unique(abundanceCR$station)
# 
# truecols <- c(rep(4,5),rep(3,8),rep(7,7))
# 
# par(mfrow=c(2,1), mar=c(2.2, 2, 0.8, 1.5) + 0.1) # B, L, T, R
# 
# plot(0, 0, col="white", xlim=c(1999.5,2020.5), xlab="Distance (km)", ylim=c(-2,60),
#      ylab="",
#      axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
# abline(h=c(10,20,30,40,50,60),lty=3, col="gray80")
# 
# for(i in stations)
# {
#   plot_data <- dplyr::filter(abundanceCR, station == stations[i] & season == 0)
#   
#   lines(2000:2020, plot_data$inlacr, type="b", col=scales::alpha(truecols[i], .5), pch=19, lty=1)
# }
# Axis(side=1, labels=TRUE, cex.axis = 1, pos=-5)
# Axis(side=2, labels=TRUE, cex.axis = 1, at=c(0,20,40,60))
# 
# plot(0, 0, col="white", xlim=c(1999.5,2020.5), xlab="Distance (km)", ylim=c(-2,60),
#      ylab="",
#      axes = FALSE,xaxs = "i",yaxs="i",yaxt = "n")
# abline(h=c(10,20,30,40,50,60),lty=3, col="gray80")
# 
# for(i in stations)
# {
#   plot_data <- dplyr::filter(abundanceCR, station == stations[i] & season == 1)
#   
#   lines(2000:2020, plot_data$inlacr, type="b", col=scales::alpha(truecols[i], .5), pch=19, lty=1)
# }
# Axis(side=1, labels=TRUE, cex.axis = 1, pos=-5)
# Axis(side=2, labels=TRUE, cex.axis = 1, at=c(0,20,40,60))
# dev.off()
# #### get seasonal aggregated weather data - Prec + Temp ####
# head(wdata3.2)
# 
# SDF <- distinct(wdata3.2, year,season)
# SDF$spring.no <- c(0,0,rep(1:21,each=4))
# SDF$spring.no <- ifelse(SDF$season=="SPRING",SDF$spring.no,0)
# SDF$summer.no <- c(0,0,rep(1:21,each=4))
# SDF$summer.no <- ifelse(SDF$season=="SUMMER",SDF$summer.no,0)
# SDF$fall.no <- c(0,0,rep(1:21,each=4))
# SDF$fall.no <- ifelse(SDF$season=="FALL",SDF$fall.no,0)
# 
# wdata_s2f <- left_join(wdata3.2,SDF)
# head(wdata3.2)
# seasonaltempdata <- aggregate(cbind(Temp.mn)~year+spring.no+fall.no+summer.no+station, data=wdata_s2f, FUN = mean)
# seasonalraindata <- aggregate(cbind(Temp.mn)~year+spring.no+fall.no+summer.no+station, data=wdata_s2f, FUN = sum)
# 
# #### seasonal temperature data ####
# seastd <- filter(seasonaltempdata, seasonaltempdata$spring.no+seasonaltempdata$fall.no+seasonaltempdata$summer.no>0)
# matrix(seastd$Temp.mn, ncol=3)
# 
# springdf <- filter(seastd, spring.no>0)
# colnames(springdf)[6] <- "spring_temp_mn"
# springdf2 <- springdf[,c(1,5,6)]
# summerdf <- filter(seastd, summer.no>0)
# colnames(summerdf)[6] <- "summer_temp_mn"
# summerdf2 <- summerdf[,c(1,5,6)]
# falldf <- filter(seastd, fall.no>0)
# colnames(falldf)[6] <- "fall_temp_mn"
# falldf2 <- falldf[,c(1,5,6)]
# 
# seasons_weather1 <-  left_join(springdf2,summerdf2)
# seasons_weather2 <-  left_join(seasons_weather1,falldf2)
# 
# #### seasonal rain data ####
# 
# seasfd <- filter(seasonalraindata, seasonalraindata$spring.no+seasonalraindata$fall.no+seasonalraindata$summer.no>0)
# 
# pspringdf <- filter(seasfd, spring.no>0)
# colnames(pspringdf)[6] <- "spring_prec_sum"
# pspringdf2 <- pspringdf[,c(1,5,6)]
# psummerdf <- filter(seasfd, summer.no>0)
# colnames(psummerdf)[6] <- "summer_prec_sum"
# psummerdf2 <- psummerdf[,c(1,5,6)]
# pfalldf <- filter(seasfd, fall.no>0)
# colnames(pfalldf)[6] <- "fall_prec_sum"
# pfalldf2 <- pfalldf[,c(1,5,6)]
# 
# seasons_weather3 <-  left_join(seasons_weather2,pspringdf2)
# seasons_weather4 <-  left_join(seasons_weather3,psummerdf2)
# seasons_weather5 <-  left_join(seasons_weather4,pfalldf2)
# 
# 
# #### get winter temperature / prec
# winterdata <- do.call(data.frame, aggregate(cbind(Prec.sum,Temp.mn)~winternumber+station, data=wdata3.2,
#                                     FUN = function(x) c(mn = mean(x), sum = sum(x) )))
# 
# winterdata2 <- filter(winterdata, winternumber>0)
# winterdata3 <- winterdata2[,c(1,2,4,5)]
# names(winterdata3)[3:4] <- c("winter_prec_sum","winter_temp_mn")
# 
# winterlabel <- data.frame(winternumber=1:21, year=2000:2020)
# winterdata4 <- left_join(winterlabel,winterdata3)
# seasons_weather6 <- left_join(seasons_weather5,winterdata4)
# 
# stdflab <- data.frame(st=unique(wdata3$station), station=c(1,10:19,2,20,3:8,21:24))
# colnames(seasons_weather6)[2] <- "st"
# seasons_weather7 <- left_join(seasons_weather6, stdflab)
# seasons_weather8 <- left_join(seasons_weather7, icedf4)
# seasons_weather9 <- select(seasons_weather8, -c(st))
# filter(seasons_weather9, station==20)
# 
# 
# write.csv2(seasons_weather9, "Vole Synchrony/data/seasonal_weather.csv")
# 

######## EXTRA ###########
wdata4 <- left_join(wdata3,mdata3)
wdata5 <- left_join(wdata4,seasdata)
names(wdata5)[4] <- "st"
avg_month <- mdata4
avg_seas <- distinct(wdata5, st,year,avg_seas_mntemp, avg_seas_totalprec,season)

m4p <- matrix(avg_seas$avg_seas_totalprec,ncol=23)
colss <- c("blue","red","green","gold","black")

aggtemp <- aggregate(avg_seas_mntemp ~ st+season, data=avg_seas, FUN=mean)
aggtemp$station <- c(1,10:19,2,20,3:8,21:24)
aggrain <- aggregate(avg_seas_totalprec ~ st, data=avg_seas, FUN=mean)
head(aggrain)

aggclimate <- left_join(aggtemp,aggrain)
stdflab <- data.frame(st=unique(wdata3$station), station=c(1,10:19,2,20,3:8,21:24))
wdata6 <- left_join(wdata5, stdflab)
head(wdata6)

### correlation matrices of weather ###
# wdata <- filter(mdata3, year > 1999 & station < 21)
# 
# head(wdata)
# plot(aggclimate$station,aggclimate$avg_seas_mntemp)
# summary(wdata6)
# cor(aggclimate$station,aggclimate$avg_seas_totalprec)
# 
# ###
# aggtempseas <- aggregate(avg_seas_mntemp ~ st+season, data=avg_seas, FUN=mean)
# aggrainseas <- aggregate(avg_seas_totalprec ~ st+season, data=avg_seas, FUN=mean)
# head(aggrainseas)
# 
# aggclimateseas <- left_join(aggtempseas,aggrainseas)
# head(aggclimateseas)
# 
# session :: save.session(file="syncsession.RSession")

varperyr_rain <- arrange(as.data.frame(xtabs(avg_month_totalprec~station+year+season,data=wdata6)),season,station,year)
varperyr_temp <- arrange(as.data.frame(xtabs(avg_month_mntemp~station+year+season,data=wdata6)),season,station,year)
varperyr <- data.frame(varperyr_rain,temp=varperyr_temp$Freq)
head(varperyr)
sss <- c("WINTER","SPRING","SUMMER", "FALL")
for(i in 1:4)
{
  sdf <- filter(avg_seas, season==sss[i])
  if (i==1) season_climate <- sdf
  else season_climate <- cbind(season_climate, sdf[,c(4,5)])

}
head(season_climate)
names(season_climate) <- c("year", "st", "delete","Winter_totalprec","Winter_temp",
                           "Spring_totalprec","Spring_temp",
                           "Summer_totalprec","Summer_temp",
                           "Fall_totalprec","Fall_temp")

sc_df <- left_join(season_climate, stdflab)
                   
sc_df2 <- sc_df[,-c(2:3)]
head(sc_df2)
GR_data$year <- GR_data$year+1999
head(sc_df2)
GR_weather <- left_join(GR_data,sc_df2, by=c("year","station"))
head(GR_weather)
write.csv2(GR_weather,"Vole Synchrony/Data/GR_weather.csv")
head(GR_weather)
hansenmodel_W2 <- lm(Rwt ~ F_t_1+S_t_1+F_t_2+S_t_2+Spring_totalprec+Spring_temp, data=GR_weather)
hansenmodel_S2 <- lm(F_t ~ S_t+F_t_1+S_t_1+F_t_2+Summer_totalprec+Summer_temp+Spring_totalprec+Spring_temp, data=GR_weather)

summary(hansenmodel_S2)
