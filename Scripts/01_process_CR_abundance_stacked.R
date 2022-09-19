# read and process abundance data ----
# Join abundance files computed in the CRINLA_zeroc3 file and the station data
# and process files to obtain abundance on previous years
abundanceCR <- readRDS("Data/Abundance_CRINLA_zeroc3.rds")
source("Scripts/0_Important_functions.R")
library(dplyr)
gps1 <- select(gps,station,distance,cumdist)
gps1$region <- c(rep(1,5),rep(2,7),rep(3,7))

abundanceCR$logCR <- log(abundanceCR$inlacr+1)
abundanceCR$station <- as.numeric(abundanceCR$station)
abundanceCR$station <- ifelse(abundanceCR$station<9,abundanceCR$station,abundanceCR$station+1)

abundanceCR$season <- rep(c(0,1),nrow(abundanceCR)/2)

springdata <- filter(abundanceCR, season==0)
falldata <- filter(abundanceCR, season==1)

colnames(springdata)[6] <- "St"
springdata$Ft <- falldata$logCR

springdata2 <- left_join(springdata,gps1)
seasondata <- select(springdata2,year,region,station,distance,cumdist,St,Ft)

# compute abundance on previous time points for each season ----
# spring
Stmat <- matrix(seasondata$St, ncol=19)
# St-1
St_1mat <- matrix(NA,nrow = 21, ncol=19)
St_1mat[2:21,] <- Stmat[1:20,]
# St-2
St_2mat <- matrix(NA,nrow = 21, ncol=19)
St_2mat[2:21,] <- St_1mat[1:20,]
# fall
Ftmat <- matrix(seasondata$Ft, ncol=19)
# Ft_1
Ft_1mat <- matrix(NA,nrow = 21, ncol=19)
Ft_1mat[2:21,] <- Ftmat[1:20,]
# Ft_2
Ft_2mat <- matrix(NA,nrow = 21, ncol=19)
Ft_2mat[2:21,] <- Ft_1mat[1:20,]

seasondata$St_1 <- as.vector(St_1mat)
seasondata$St_2 <- as.vector(St_2mat)
seasondata$Ft_1 <- as.vector(Ft_1mat)
seasondata$Ft_2 <- as.vector(Ft_2mat)
View(seasondata)
# df for stacked data -----
# add 0's for regional data
stackdf0 <- data.frame(seasondata,seasondata[,6:11],seasondata[,6:11],seasondata[,6:11])
head(stackdf0)
colcount <- 11
for(i in 1:3)
{
  groupx <- stackdf0$region==i
  stackdf0[!groupx,(1+colcount):(6+colcount)] <- 0
  colcount = colcount + 6
}

stackdf1 <- tibble(stackdf0)
# View(stackdf1)

saveRDS(stackdf1, "Data/data4synchrony.rds")
