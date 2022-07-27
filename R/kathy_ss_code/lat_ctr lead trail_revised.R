####################################################################################
#
# LOOK AT CHANGE IN LAT OF DIFFERENT SEGMENTS OF SPP RANGE
#
####################################################################################

library(raster)

load("C:/Users/kmills/Dropbox/NROC_data_exploration/trawl.data/Survdat_Nye2016.Rdata")
#setwd("C:/Users/jschuetz/Dropbox/NROC_data_exploration")
names(survdat)

# Years to use
year <- seq(1970, 2015)

# output df
out7015 <- data.frame()

# Species to use
trawl.spp <- c(72,141, 103,106, 143,77,503, 136,301,121, 73, 105)
trawl.spp.name <- c('SILVER HAKE', 'BLACK SEA BASS', 'SUMMER FLOUNDER', 'WINTER FLOUNDER', 'SCUP',
                  'RED HAKE', 'LONGFIN SQUID', 'ATLANTIC CROAKER', 'AMERICAN LOBSTER',
                  'ATLANTIC MACKEREL', 'ATLANTIC COD','YELLOWTAIL FLOUNDER')
 
# Pull columns
trawl.data <- survdat[, c('ID', 'EST_YEAR','SEASON', 'STRATUM', 'DECDEG_BEGLAT','DECDEG_BEGLON',
                            'SVSPP', 'CATCHSEX','BIOMASS', 'AVGDEPTH', 'ABUNDANCE')]

#STRATA FOR WHOLE OFFSHORE SURVEY
trawl.data <- trawl.data[trawl.data$STRATUM >= 01010 & trawl.data$STRATUM <= 01760,]  

# Drop specific strata
trawl.data<-trawl.data[trawl.data$STRATUM!=1310 & trawl.data$STRATUM!=1320 & trawl.data$STRATUM!=1330
                         & trawl.data$STRATUM!=1410 & trawl.data$STRATUM!=1420 & trawl.data$STRATUM!=1490,]

# check dims
dim(trawl.data)

head(trawl.data)

############### generate template of all sampling occasions 

# Unique sampling locations
template <- trawl.data[!duplicated(trawl.data[ , c("ID", "EST_YEAR", "SEASON", "STRATUM", "AVGDEPTH")]), c("ID", "EST_YEAR", "SEASON", "STRATUM", "DECDEG_BEGLAT", "DECDEG_BEGLON", "AVGDEPTH")]

############# summarize biomass for each spp in each sample


#SELECT UNIQUE ROWS BY ID, SVSPP, AND CATCHSEX TO GET RID OF LENGTH DUPLICATES
samples <- trawl.data[!duplicated(trawl.data[, c("ID", "SVSPP", "CATCHSEX")]), ]

#SUM BIOMASS OVER CATCHSEX 
samples.biomass <- aggregate(BIOMASS ~ ID + SVSPP, samples, sum)  

# Join totals back to template
spp.data <- merge(samples.biomass, template, by="ID", all.y=T)

out <- data.frame()



# Loop through species
for (s in 1:length(trawl.spp)){ # loop for species
  
  #trawl.data<- spp.data[spp.data$SVSPP==136,] 
  trawl.data<- spp.data[spp.data$SVSPP==trawl.spp[s],] # subset data for species
  
  
for (y in 1:length(year)){  # loop for years
      
      #trawl2<- trawl.data[trawl.data$EST_YEAR == 1987, ]
  
      # subset data for time
      trawl2<- trawl.data[trawl.data$EST_YEAR == year[y], ]
      
      # Drop NA ID's
      trawl2<-trawl2[!is.na(trawl2$ID),]
      
      # Order by Lat
      trawl3<-trawl2[order(trawl2$DECDEG_BEGLAT),]
      
      # Cumulative biomass along lat
      trawl3$cumsumBIOM<-cumsum(trawl3$BIOMASS)
      
      # Total annual Biomass
      biomass.total.year <- sum(trawl3$BIOMASS, na.rm=TRUE)
    
      
#SEGMENTS
     
    # Center of Biomass
    center.lat <- sum(trawl3$DECDEG_BEGLAT * (trawl3$BIOMASS/biomass.total.year))
    center.lon <- sum(trawl3$DECDEG_BEGLON*(trawl3$BIOMASS/biomass.total.year))
    
    # Trailing edge: 10%
    trail10<-trawl3[trawl3$cumsumBIOM<=.1*biomass.total.year,]
    biomass.tot10<-sum(trail10$BIOMASS, na.rm=TRUE)
    trail10.lat <- sum(trail10$DECDEG_BEGLAT*(trail10$BIOMASS/biomass.tot10))
    trail10.lon <- sum(trail10$DECDEG_BEGLON*(trail10$BIOMASS/biomass.tot10))
    
    # Trailing Edge: 25%
    trail25<-trawl3[trawl3$cumsumBIOM<=.25*biomass.total.year,]
    biomass.tot25<-sum(trail25$BIOMASS, na.rm=TRUE)
    trail25.lat <- sum(trail25$DECDEG_BEGLAT*(trail25$BIOMASS/biomass.tot25))
    trail25.lon <- sum(trail25$DECDEG_BEGLON*(trail25$BIOMASS/biomass.tot25))
    
    lead75<-trawl3[trawl3$cumsumBIOM>=.75*biomass.total.year,]
    biomass.tot75<-sum(lead75$BIOMASS, na.rm=TRUE)
    lead75.lat <- sum(lead75$DECDEG_BEGLAT*(lead75$BIOMASS/biomass.tot75))
    lead75.lon <- sum(lead75$DECDEG_BEGLON*(lead75$BIOMASS/biomass.tot75))
    
    lead90<-trawl3[trawl3$cumsumBIOM>=.9*biomass.total.year,]
    biomass.tot90<-sum(lead90$BIOMASS, na.rm=TRUE)
    lead90.lat <- sum(lead90$DECDEG_BEGLAT*(lead90$BIOMASS/biomass.tot90))
    lead90.lon <- sum(lead90$DECDEG_BEGLON*(lead90$BIOMASS/biomass.tot90))
    
    out[(s-1)*(length(year))+y, 1]<- trawl.spp[s]
    out[(s-1)*(length(year))+y, 2]<- trawl.spp.name[s]
    out[(s-1)*(length(year))+y, 3]<- year[y]
    out[(s-1)*(length(year))+y, 4] <- center.lat
    out[(s-1)*(length(year))+y, 5] <- center.lon
    out[(s-1)*(length(year))+y, 6] <- lead90.lat
    out[(s-1)*(length(year))+y, 7] <- lead90.lon
    out[(s-1)*(length(year))+y, 8] <- lead75.lat
    out[(s-1)*(length(year))+y, 9] <- lead75.lon
    out[(s-1)*(length(year))+y, 10] <- trail25.lat
    out[(s-1)*(length(year))+y, 11] <- trail25.lon
    out[(s-1)*(length(year))+y, 12] <- trail10.lat
    out[(s-1)*(length(year))+y, 13] <- trail10.lon
    
  
}
}
      
  colnames(out)<- c("TRAWL.SPP","COMNAME","YEAR","CenterLat", "CenterLon", "Lead90Lat", "Lead90Lon",
                      "Lead75Lat", "Lead75Lon", "Trail25Lat", "Trail25Lon","Trail10Lat","Trail10Lon")
 
Lead.Trail.Center<-out
write.csv(Lead.Trail.Center, file="Lead.Trail.Center.csv")

#####################
 plots
library(plotrix)
library(zoo)

 ctr.cut<-out[out$CenterLat!=0 & !is.na(out$CenterLat),]
 lead90.cut<-out[out$Lead90Lat!=0 & !is.na(out$Lead90Lat),]
 lead75.cut<-out[out$Lead75Lat!=0 & !is.na(out$Lead75Lat),]
 trail10.cut<-out[out$Trail10Lat!=0 & !is.na(out$Trail10Lat),]
 trail25.cut<-out[out$Trail25Lat!=0 & !is.na(out$Trail25Lat),]


par(mfrow=c(2,4))
s<-5
for (s in 1:length(trawl.spp)){
  data1<-trail10.cut[trail10.cut$TRAWL.SPP==trawl.spp[s],]
  min<-min(data1$Trail10Lat)
  data2<-lead90.cut[lead90.cut$TRAWL.SPP==trawl.spp[s],]
  max<-max(data2$Lead90Lat)
  plot(data2$YEAR, data2$CenterLat, type='n', pch=16, 
       ylim=c(min-0.5,max+0.5), 
       xlab='YEAR', ylab='CENTER OF LATITUDE', 
       main=trawl.spp.name[s])
  
  data3<-ctr.cut[ctr.cut$TRAWL.SPP==trawl.spp[s],]
  lines(data3$YEAR,data3$CenterLat, lwd=2)
  #lines(data3$YEAR, c(NA,NA,rollmean(data3$CenterLat,5),NA,NA),lwd=3,col="black")
  abline(lm(CenterLat~YEAR,data=data3), col="black", lwd=3)
  
  lines(data2$YEAR,data2$Lead90Lat, lwd=2, col="blue")
  #lines(data2$YEAR, c(NA,NA,rollmean(data2$Lead90Lat,5),NA,NA),lwd=3,col="blue")
  abline(lm(Lead90Lat~YEAR,data=data2), col="blue", lwd=3)
  
  data4<-lead75.cut[lead75.cut$TRAWL.SPP==trawl.spp[s],]
  lines(data4$YEAR,data4$Lead75Lat, lwd=2, col="green")
  #lines(data4$YEAR, c(NA,NA,rollmean(data4$Lead75Lat,5),NA,NA),lwd=3,col="green")
  abline(lm(Lead75Lat~YEAR,data=data4), col="green", lwd=3)
  
  lines(data1$YEAR,data1$Trail10Lat, lwd=2, col="red")
  #lines(data1$YEAR, c(NA,NA,rollmean(data1$Trail10Lat,5),NA,NA),lwd=3,col="red")
  abline(lm(Trail10Lat~YEAR,data=data1), col="red", lwd=3)
  
  data5<-trail25.cut[trail25.cut$TRAWL.SPP==trawl.spp[s],]
  lines(data5$YEAR,data5$Trail25Lat, lwd=2, col="orange")
  #lines(data5$YEAR, c(NA,NA,rollmean(data5$Trail25Lat,5),NA,NA),lwd=3,col="orange")
  abline(lm(Trail25Lat~YEAR,data=data5), col="orange", lwd=3)
  
}

trawl.spp2<-c(72,141, 103, 143,301, 105)
trawl.spp.name2<-c('SILVER HAKE', 'BLACK SEA BASS', 'SUMMER FLOUNDER', 'SCUP',
                  'AMERICAN LOBSTER', 'YELLOWTAIL FLOUNDER')

par(mfrow=c(1,1))
plot(data$YEAR, data$CENTER_OF_MASS_LATITUDE, type='n', pch=16, 
     ylim=c(min(out7015.cut$CENTER_OF_MASS_LATITUDE)-0.2,max(out7015.cut$CENTER_OF_MASS_LATITUDE)+0.2), 
     xlab='YEAR', ylab='CENTER OF LATITUDE')
for (s in 1:length(trawl.spp2)){
  data<-out7015.cut[out7015.cut$TRAWL.SPP==trawl.spp2[s],]
  lines(data$YEAR,data$CENTER_OF_MASS_LATITUDE, lwd=1, col="grey")
  lines(data$YEAR, c(NA,NA,rollmean(data$CENTER_OF_MASS_LATITUDE,5),NA,NA),lwd=3,col=s)
}




head(out.data)

for (s in 1:length(spp.code)){ # loop for species
  
  spp.data<- cleandat[cleandat$SPP_CODE == spp.code[s],]   

sort(unique(out$TRAWL.SPP))

names(out.data)

  names(out)
library(MuMIn)

centers <-  read.csv("trawl temperature centers of biomass using point locations.csv")
t <- centers$reference == "SHIFTING_DISTRIBUTION"
centers$contrast[t] <- "TRACK"
centers$contrast[!t] <- "ADAPT"

trawl.spp <- c(75, 72, 141, 401, 32, 103, 503, 77, 131, 502, 143) #federal survey codes
spp.name <- c('Pollock', 'Silver Hake', 'Black Sea Bass', 'Scallops',
              'Herring', 'Summer Flounder', 'Longfin Squid', 'Red Hake',
              'Butterfish', 'Shortfin Squid', 'Scup')

out.estimates <- data.frame()
out.models <- data.frame()

pdf("mean_T_plots_trawl_points.pdf", width = 5, height = 5.5, onefile=T)

for (s in 1:length(trawl.spp)){
  
  my.data<-centers[centers$spp==trawl.spp[s],]
  
  options(na.action = "na.fail")
  model <- lm(WEIGHTED_SST_MEAN ~  YEAR_T * contrast, my.data)
  model.table <-dredge(model, extra = list(
    "R^2", "*" = function(x) {
      s <- summary(x)
      c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
        F = s$fstatistic[[1]])
    }))
  
  set <- get.models(model.table, subset=T)
  avg <- model.avg(set)
  
  my.data$PREDS<- predict(avg, full=T)
  se <- predict(avg, se.fit=T)
  my.data$SE_2p <- se$fit+2*se$se.fit
  my.data$SE_2n <- se$fit-2*se$se.fit
  
  plot(my.data$YEAR_T, my.data$WEIGHTED_SST_MEAN, 
       pch=16, cex=0.2, 
       col="white", 
       main=paste(spp.name[s]), 
       xlab="year",
       ylab="mean SST")
  
  start <- seq(19, 335, 18)
  end <- seq(36, 342, 18)
  
  lines(my.data$YEAR_T[1:18], my.data$WEIGHTED_SST_MEAN[1:18], col=1, lwd=.2, lty=2)
  
  for(i in 1:17){
    
    lines(my.data$YEAR_T[start[i]:end[i]], my.data$WEIGHTED_SST_MEAN[start[i]:end[i]], col=2, lwd=.2, lty=2)
    
  }
  
  lines(my.data$YEAR_T[1:18], my.data$PREDS[1:18], col=1, lwd=1, lty=1)
  
  for(i in 1:17){
    
    lines(my.data$YEAR_T[start[i]:end[i]], my.data$PREDS[start[i]:end[i]], col=2, lwd=1, lty=1)
    
  }
  
  estimates <- data.frame(coefTable(avg, full=T))
  estimates$species <- spp.name[s]
  out.estimates <- rbind(out.estimates, estimates)
  
  model.selection <- data.frame(model.sel(model.table))
  model.selection$species <- spp.name[s]
  out.models <- rbind(out.models, model.selection)
  
  #my.data$RESIDS <- my.data$weighted.lat - my.data$PREDS
  #hist(my.data$RESIDS)
  
}

dev.off()

write.csv(out.estimates, "mean_T_coefficients_trawl_points.csv")
write.csv(out.models, "mean_T_model_selection_tables_trawl_points.csv", row.names=F)


##############


out.estimates <- data.frame()
out.models <- data.frame()

pdf("min_T_plots_trawl_points.pdf", width = 5, height = 5.5, onefile=T)

for (s in 1:length(trawl.spp)){
  
  my.data<-centers[centers$spp==trawl.spp[s],]
  
  options(na.action = "na.fail")
  model <- lm(WEIGHTED_SST_MIN ~  YEAR_T * contrast, my.data)
  model.table <-dredge(model, extra = list(
    "R^2", "*" = function(x) {
      s <- summary(x)
      c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
        F = s$fstatistic[[1]])
    }))
  
  set <- get.models(model.table, subset=T)
  avg <- model.avg(set)
  
  my.data$PREDS<- predict(avg, full=T)
  se <- predict(avg, se.fit=T)
  my.data$SE_2p <- se$fit+2*se$se.fit
  my.data$SE_2n <- se$fit-2*se$se.fit
  
  plot(my.data$YEAR_T, my.data$WEIGHTED_SST_MIN, 
       pch=16, cex=0.2, 
       col="white", 
       main=paste(spp.name[s]), 
       xlab="year",
       ylab="mean SST")
  
  start <- seq(19, 335, 18)
  end <- seq(36, 342, 18)
  
  lines(my.data$YEAR_T[1:18], my.data$WEIGHTED_SST_MIN[1:18], col=1, lwd=.2, lty=2)
  
  for(i in 1:17){
    
    lines(my.data$YEAR_T[start[i]:end[i]], my.data$WEIGHTED_SST_MIN[start[i]:end[i]], col=2, lwd=.2, lty=2)
    
  }
  
  lines(my.data$YEAR_T[1:18], my.data$PREDS[1:18], col=1, lwd=1, lty=1)
  
  for(i in 1:17){
    
    lines(my.data$YEAR_T[start[i]:end[i]], my.data$PREDS[start[i]:end[i]], col=2, lwd=1, lty=1)
    
  }
  
  estimates <- data.frame(coefTable(avg, full=T))
  estimates$species <- spp.name[s]
  out.estimates <- rbind(out.estimates, estimates)
  
  model.selection <- data.frame(model.sel(model.table))
  model.selection$species <- spp.name[s]
  out.models <- rbind(out.models, model.selection)
  
  #my.data$RESIDS <- my.data$weighted.lat - my.data$PREDS
  #hist(my.data$RESIDS)
  
}

dev.off()

write.csv(out.estimates, "min_T_coefficients_trawl_points.csv")
write.csv(out.models, "min_T_model_selection_tables_trawl_points.csv", row.names=F)

##########

out.estimates <- data.frame()
out.models <- data.frame()

pdf("max_T_plots_trawl_points.pdf", width = 5, height = 5.5, onefile=T)

for (s in 1:length(trawl.spp)){
  
  my.data<-centers[centers$spp==trawl.spp[s],]
  
  options(na.action = "na.fail")
  model <- lm(WEIGHTED_SST_MAX ~  YEAR_T * contrast, my.data)
  model.table <-dredge(model, extra = list(
    "R^2", "*" = function(x) {
      s <- summary(x)
      c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
        F = s$fstatistic[[1]])
    }))
  
  set <- get.models(model.table, subset=T)
  avg <- model.avg(set)
  
  my.data$PREDS<- predict(avg, full=T)
  se <- predict(avg, se.fit=T)
  my.data$SE_2p <- se$fit+2*se$se.fit
  my.data$SE_2n <- se$fit-2*se$se.fit
  
  plot(my.data$YEAR_T, my.data$WEIGHTED_SST_MAX, 
       pch=16, cex=0.2, 
       col="white", 
       main=paste(spp.name[s]), 
       xlab="year",
       ylab="mean SST")
  
  start <- seq(19, 335, 18)
  end <- seq(36, 342, 18)
  
  lines(my.data$YEAR_T[1:18], my.data$WEIGHTED_SST_MAX[1:18], col=1, lwd=.2, lty=2)
  
  for(i in 1:17){
    
    lines(my.data$YEAR_T[start[i]:end[i]], my.data$WEIGHTED_SST_MAX[start[i]:end[i]], col=2, lwd=.2, lty=2)
    
  }
  
  lines(my.data$YEAR_T[1:18], my.data$PREDS[1:18], col=1, lwd=1, lty=1)
  
  for(i in 1:17){
    
    lines(my.data$YEAR_T[start[i]:end[i]], my.data$PREDS[start[i]:end[i]], col=2, lwd=1, lty=1)
    
  }
  
  estimates <- data.frame(coefTable(avg, full=T))
  estimates$species <- spp.name[s]
  out.estimates <- rbind(out.estimates, estimates)
  
  model.selection <- data.frame(model.sel(model.table))
  model.selection$species <- spp.name[s]
  out.models <- rbind(out.models, model.selection)
  
  #my.data$RESIDS <- my.data$weighted.lat - my.data$PREDS
  #hist(my.data$RESIDS)
  
}

dev.off()

write.csv(out.estimates, "max_T_coefficients_trawl_points.csv")
write.csv(out.models, "max_T_model_selection_tables_trawl_points.csv", row.names=F)

