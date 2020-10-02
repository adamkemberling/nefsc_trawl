setwd("C:/Users/kmills/Documents/Active research/Lobsters")
load("C:/Users/kmills/Documents/Active research/CoastalSEES/Data/Survdat_Nye_allseason.Rdata")

names(survdat)
stratum.area<-read.table("C:/Users/kmills/Documents/Active research/CoastalSEES/Data/stratum.areas.csv",
                   header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE)
names(stratum.area)

head(stratum.area)

survdat$STRATUM[1:25]
survdat$AREA[1:25]

trawl.data <- survdat[, c('ID', 'EST_YEAR','SEASON', 'STRATUM', 'SVSPP', 'COMNAME', 'CATCHSEX','BIOMASS', 'ABUNDANCE')]
trawl.data <- trawl.data[trawl.data$STRATUM >= 01010 & trawl.data$STRATUM <= 01760,]  #STRATA FOR WHOLE OFFSHORE SURVEY
trawl.data <- trawl.data[trawl.data$STRATUM!=1310 & trawl.data$STRATUM!=1320 & trawl.data$STRATUM!=1330 & trawl.data$STRATUM!=1350 &
                           trawl.data$STRATUM!=1410 & trawl.data$STRATUM!=1420 & trawl.data$STRATUM!=1490,] #ELIMINATE STRATA NOT SAMPLED FOR FULL TIME SPAN
trawl.data <- trawl.data[trawl.data$SEASON=='SPRING' | trawl.data$SEASON=='FALL',]
trawl.data <- trawl.data[trawl.data$EST_YEAR >= 1970 & trawl.data$EST_YEAR < 2018,]
trawl.data <- trawl.data[!(is.na(trawl.data$BIOMASS)),]
#FOR LOBSTER SHOULD PROBABLY CUT STRATA SET FURTHER TO JUST GOM/GB

#MERGING WITH STRATUM AREA
trawl.data2<-merge(trawl.data,stratum.area,by="STRATUM")
trawl.data3<-trawl.data2[order(trawl.data2$ID, trawl.data2$SVSPP, trawl.data2$CATCHSEX),]
unique.str<-trawl.data3[!duplicated(trawl.data3[,'STRATUM']),]
totstrarea<-sum(unique.str$STRATUM_AREA)
trawl.data3$STRATIO<-trawl.data3$STRATUM_AREA/totstrarea

tows_unique<-trawl.data3[!duplicated(trawl.data3["ID"]),] 
numtows<-aggregate(ID~EST_YEAR+STRATUM,length,data=tows_unique)
colnames(numtows)<-c('EST_YEAR','STRATUM','COUNT')

#CUTTING TO JUST ONE ROW PER ID/SVSPP/CATCHSEX
bn_unique<-trawl.data3[!duplicated(trawl.data3[c("ID", "SVSPP", "CATCHSEX")]),]

#SUM BIOMASSES AND ABUNDANCES OVER CATCHSEX	
bn_unique2<-aggregate(cbind(BIOMASS, ABUNDANCE) ~ ID + SVSPP, sum, data=bn_unique)
bn_unique2<-bn_unique2[order(bn_unique2$ID),]
bn_unique3<-unique(bn_unique[,c(1:6,10,11)])
bn_unique4<-merge(bn_unique2, bn_unique3, by=c("ID", "SVSPP"))

#DOING STRATIFIED MEAN IN LOOP
years<-seq(1970,2017)
svspp<-c(103, 149, 136)
stratum<-sort(unique(trawl.data3$STRATUM))
stmean.B.loop.out<-matrix(0,length(years),length(svspp))
colnames(stmean.B.loop.out)<-as.character(svspp)
rownames(stmean.B.loop.out)<-as.character(years)
bn_setup_strata<-aggregate(BIOMASS~EST_YEAR + SVSPP + STRATUM + STRATIO, sum, data=bn_unique4)
bn_setup_strata<-merge(bn_setup_strata,numtows, by=c("EST_YEAR", "STRATUM"), all.y=TRUE)
bn_setup_strata$MBIOM.PER.TOW<-bn_setup_strata$BIOMASS/bn_setup_strata$COUNT
  #THE COUNT HERE FULLY ACCOUNTS FOR TOWS IN THE STRATA...A WORK-AROUND TO FULLY POPULATING DATAFRAME WITH ZEROS FOR
  #TOWS IN WHICH SPP WAS NOT OBSERVED
bn_setup_strata$WTMEAN.PER.TOW<-bn_setup_strata$MBIOM.PER.TOW*bn_setup_strata$STRATIO
  #WEIGHTED PROPORTIONAL SHARE FOR STRATA IN WHICH SPP WAS OBSERVED
  #WHEN SUM BELOW, STRATA WITH NO OBSERVATIONS DON'T HAVE ANYTHING TO CONTRIBUTE TO SUM
  #SO LEAVING THEIR CONTRIBUTION OUT ACCOUNTS FOR ZEROS

#THis STEP SUMS OVER THOSE WEIGHTED MEANS WHICH REALLY REFLECT PROPORTIONAL CONTRIBUTION OF EACH STRATA
#WHEN SUM, HAVE WEIGHTED MEAN STRATIFIED BIOMASS PER TOW
for(s in svspp){
  for(i in years){
    data<-bn_setup_strata[bn_setup_strata$SVSPP==s,]
    data<-data[data$EST_YEAR==i,]
    stmean.B.loop.out[which(years==i),which(svspp==s)]<-sum(data$WTMEAN.PER.TOW)
  }
}

par(mfrow=c(3,1))
plot(years, stmean.B.loop.out[,1], type='b')
plot(years, stmean.B.loop.out[,2], type='b')
plot(years, stmean.B.loop.out[,3], type='b')

#JUST AN EXAMPLE BREAKING THINGS APART FOR COD

test<-aggregate(MBIOM.PER.TOW~EST_YEAR + SVSPP, mean, data=bn_setup_strata)
test.cod<-bn_setup_strata[bn_setup_strata$SVSPP==73,]
with(test.cod, table(EST_YEAR, STRATUM)
test.cod.agg<-aggregate(MBIOM.PER.TOW ~ EST_YEAR, mean, data=test.cod)
test.cod.agg2<-aggregate(STRATUM ~ EST_YEAR, FUN="length", data=test.cod)
test.cod.agg3<-aggregate(WTMEAN.PER.TOW ~ EST_YEAR, sum, data=test.cod)
test.cod.all<-merge(test.cod.agg, test.cod.agg2, by="EST_YEAR")
test.cod.all<-merge(test.cod.all, test.cod.agg3, by="EST_YEAR")

#Mean biomass per tow without accounting for stratified survey design or fully accounting for zeros 
#(doesn't account for strata with no obs)
plot(test.cod.all$EST_YEAR, test.cod.all$MBIOM.PER.TOW, type="l") 
#Stratified mean biomass per tow accounting for stratified design and zeros
lines(test.cod.all$EST_YEAR, test.cod.all$WTMEAN.PER.TOW, col="red")
