setwd('C:/Users/kmills/Documents/Active Research/Lenfest/Ecosystem change analyses/Data')
load('all.means.out.Rdata')
load('bn.setup2.RData')
load('sppclass2.Rdata')
load('lwreg.Rdata')
load('wigleyspp2.Rdata')

#save(sppclass2, file="sppclass2.Rdata")
#save(lwreg,file="lwreg.Rdata")
#save(wigleyspp2,file="wigleyspp2.Rdata")

names(bn.setup2)
head(bn.setup2)

svspp_uniq <- sort(unique(bn.setup2$SVSPP))
svspp_uniq <- as.data.frame(svspp_uniq)
write.csv(svspp_uniq,file="svspp_uniq")

#CHECKING PORTION OF FISH REPRESENTED IN WIGLEY ET AL 2003 BASED ON TOW-LEVEL DATA
names(bn.setup2)
names(wigleyspp2)
wigley <- wigleyspp2[,c(1,4,5)]
bn <- merge(bn.setup2,wigley,by.x="SVSPP",by.y="SVSPP",all.x=TRUE)
head(bn)

#CHECKING PORTION OF FISH REPRESENTED IN WIGLEY ET AL 2003 BASED ON EXPNUMLEN 
sum(bn$NUMLEN,na.rm=TRUE)
bn.wigley <- bn[bn$wigley=="X",]
sum(bn.wigley$NUMLEN,na.rm=TRUE)
sum(bn.wigley$NUMLEN,na.rm=TRUE)/sum(bn$NUMLEN,na.rm=TRUE)
#82%

head(bn)
fish <- bn[bn$INV==0,]
sum(fish$NUMLEN,na.rm=TRUE)
fish.wigley <- fish[fish$wigley=="X",]
sum(fish.wigley$NUMLEN,na.rm=TRUE)
sum(fish.wigley$NUMLEN,na.rm=TRUE)/sum(fish$NUMLEN,na.rm=TRUE)
#94.7% of fish have L-W coefficients in Wigley

########################

#GETTING WEIGHTS FROM L-W COEFFICIENTS

fish.cut <- fish.wigley[fish.wigley$LENGTH!=0,]
fish.cut <- fish.cut[!(is.na(fish.cut$LENGTH)),]
#fish.cut[fish.cut$LENGTH==0,]

fish.cut$LLEN <- log(fish.cut$LENGTH)
class <- merge(fish.cut,sppclass2,by.x="SVSPP",by.y="SVSPP")
class.wigley <- class[class$wigley=="X",]
len.reg <- merge(class.wigley,lwreg,all.x=T)

len.reg$LWEIGHT <- len.reg$lna + (len.reg$b*len.reg$LLEN)
len.reg$WT <- exp(len.reg$LWEIGHT)
len.reg$TOTWT <- len.reg$WT * len.reg$NUMLEN

len.reg <- len.reg[,c(1:5,7,9,10,11,12,13,14,15,16,17,21,23,38,39,40)]

#SETTING UP BINS
hist(len.reg$LLEN, breaks=10, plot=TRUE)
hist(len.reg$LWEIGHT,breaks=10,plot=TRUE)
#hist(len.reg$L10WT,breaks=20,plot=TRUE)

len.reg$LENBIN2 <- 0
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>=0.5 & len.reg$LLEN<=1,1,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>1 & len.reg$LLEN<=1.5,2,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>1.5 & len.reg$LLEN<=2,3,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>2 & len.reg$LLEN<=2.5,4,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>2.5 & len.reg$LLEN<=3,5,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>3 & len.reg$LLEN<=3.5,6,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>3.5 & len.reg$LLEN<=4,7,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>4 & len.reg$LLEN<=4.5,8,len.reg$LENBIN2)
len.reg$LENBIN2 <- ifelse(len.reg$LLEN>4.5 & len.reg$LLEN<=5,9,len.reg$LENBIN2)

len.reg$WTBIN <- len.reg$LWEIGHT
len.reg$WTBIN <- ifelse(len.reg$WTBIN>4,8,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN>2 & len.reg$WTBIN<=4,7,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN>0 & len.reg$WTBIN<=2,6,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN>(-2) & len.reg$WTBIN<=0,5,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN>(-4) & len.reg$WTBIN<=(-2),4,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN>(-6) & len.reg$WTBIN<=(-4),3,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN>=(-8) & len.reg$WTBIN<=(-6),2,len.reg$WTBIN)
len.reg$WTBIN <- ifelse(len.reg$WTBIN<(-8),1,len.reg$WTBIN)

len.reg$L10WT <- log10(len.reg$INDWT)
len.reg$WT10BIN <- len.reg$L10WT
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>2.5,20,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>2 & len.reg$WT10BIN<=2.5,19,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>1.5 & len.reg$WT10BIN<=2,18,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>1 & len.reg$WT10BIN<=1.5,17,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>0.5 & len.reg$WT10BIN<=1,16,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>0 & len.reg$WT10BIN<=0.5,15,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-0.5 & len.reg$WT10BIN<=0,14,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-1 & len.reg$WT10BIN<=-0.5,13,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-1.5 & len.reg$WT10BIN<=-1,12,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-2 & len.reg$WT10BIN<=-1.5,11,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-2.5 & len.reg$WT10BIN<=-2,10,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-3 & len.reg$WT10BIN<=-2.5,9,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-3.5 & len.reg$WT10BIN<=-3,8,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-4 & len.reg$WT10BIN<=-3.5,7,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-4.5 & len.reg$WT10BIN<=-4,6,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-5 & len.reg$WT10BIN<=-4.5,5,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-5.5 & len.reg$WT10BIN<=-5,4,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-6 & len.reg$WT10BIN<=-5.5,3,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-6.5 & len.reg$WT10BIN<=-6,2,len.reg$WT10BIN)
len.reg$WT10BIN <- ifelse(len.reg$WT10BIN>-7 & len.reg$WT10BIN<=-6.5,1,len.reg$WT10BIN)

#COUNTING TOWS FOR STRATIFIED MEAN
#USES TOWS JUST WITH FISH SPECIES IN WIGLEY
tows_unique <- len.reg[!duplicated(len.reg["ID"]),] 
numtows <- aggregate(ID~EST_YEAR+STRATUM,length,data=tows_unique)
colnames(numtows) <- c('EST_YEAR','STRATUM','COUNT')

#SUM NUMLEN WITHIN EACH LENBIN2 FOR YEAR,STRATUM  

#ab.bin <- aggregate(NUMLEN ~ EST_YEAR + STRATUM + ID + LENBIN2, sum, data=len.reg)
#ab.bin$LNUMLEN <- log10(ab.bin$NUMLEN+1)
ab.bin.str.nolog <- aggregate(NUMLEN ~ EST_YEAR + STRATUM + LENBIN2, sum, data=len.reg)

#MERGE IN STRATUM RATIOS
str.wt <- len.reg[,c(4,11,12)]
str.wt.unique <- str.wt[!duplicated(str.wt["STRATUM"]),]
ab.bin.strwt.nolog <- merge(ab.bin.str.nolog,str.wt.unique,by="STRATUM")
totstrwt <- sum(str.wt.unique$STRATUM_AREA)

#CALCULATING TOTAL ABUNDANCE BY YEAR
years <- seq(1968,2013)
LENBIN2 <- c(1,2,3,4,5,6,7,8,9)
Alen.out <- matrix(0,length(years),length(LENBIN2))
colnames(Alen.out) <- as.character(LENBIN2)
rownames(Alen.out) <- as.character(years)
len_setup_strata <- aggregate(NUMLEN~EST_YEAR + LENBIN2 + STRATUM + STRATIO, sum, data=ab.bin.strwt.nolog)
len_setup_strata <- merge(len_setup_strata, numtows)
len_setup_strata$mean.ABUNDANCE <- len_setup_strata$NUMLEN / len_setup_strata$COUNT
len_setup_strata$WTMEAN <- (len_setup_strata$mean.ABUNDANCE)*len_setup_strata$STRATIO

for(x in LENBIN2){
  for(y in years){
    data <- len_setup_strata[len_setup_strata$LENBIN2==x,]
    data <- data[data$EST_YEAR==y,]
    Alen.out[which(years==y), which(LENBIN2==x)] <- log10(( sum(data$WTMEAN) / .01) * totstrwt)
  }
}

write.csv(Alen.out,"Alen.out")
save(Alen.out, file="Alen.out.Rdata")
Alen.out

par(mfrow=c(1,1))
plot(seq(1968,2013),seq(4.5,10,length=length(seq(1968,2013))),type="n")
title("Abundance by length class--all fish", cex=1.1)
lines(seq(1968,2013),Alen.out[,1],col="black") 
lines(seq(1968,2013),Alen.out[,2],col="red") #down after 2006
lines(seq(1968,2013),Alen.out[,3],col="blue") 
lines(seq(1968,2013),Alen.out[,4],col="green")  #steady upwards, jump after 2008
lines(seq(1968,2013),Alen.out[,5],col="orange") #jump in 2008
lines(seq(1968,2013),Alen.out[,6],col="purple") #steady upwards
lines(seq(1968,2013),Alen.out[,7],col="grey")
lines(seq(1968,2013),Alen.out[,8],col="pink")
lines(seq(1968,2013),Alen.out[,9],col="brown")  #come way down 1990-2000

#CLUSTERING

library(rioja)
library(vegan)
library(zoo)

#to 2013
agg.cl <- Alen.out#[2:41,1:9]
agg.cl <- rollmean(agg.cl,3)
agg.eucdist <- vegdist(agg.cl, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE) 
agg.eucdist.cl <- chclust(agg.eucdist, method = "coniss")
plot(agg.eucdist.cl,labels=seq(1970,2012),hang=-0.1,axes=F,cex=1)
axis(side=2, cex.axis=1.1)
title("Abundance by length classes--fish, 1968-2013",cex=1.1)
mtext(side=2, line=2.3, "Sum of squares", cex=1.1,las=0)
num <- bstick(agg.eucdist.cl)
num$diff <- num$dispersion-num$bstick
num.cl <- num[num$diff>=0,]
d <- dim(num.cl)
num.cl <- num.cl$nGroups[d[1]]
clus <- cutree(agg.eucdist.cl,k=num.cl)

#DEFINING PERIODS

ab.bin.strwt.nolog$period <- 0
ab.bin.strwt.nolog$period <- ifelse(ab.bin.strwt.nolog$EST_YEAR<=1989,1,ab.bin.strwt.nolog$period)
ab.bin.strwt.nolog$period <- ifelse(ab.bin.strwt.nolog$EST_YEAR>=1995 & 
    ab.bin.strwt.nolog$EST_YEAR<=2008,2,ab.bin.strwt.nolog$period)

len.reg$period <- 0
len.reg$period <- ifelse(len.reg$EST_YEAR<=1989,1,len.reg$period)
len.reg$period <- ifelse(len.reg$EST_YEAR>=1995 & len.reg$EST_YEAR<=2008,2,len.reg$period)

##ABUNDANCE BY PERIOD
#TOWS
names(ab.bin.strwt.nolog)
dim(len.reg)

fish.tows <- len.reg[len.reg$period==1 | len.reg$period==2,]
tows_unique2 <- fish.tows[!duplicated(fish.tows["ID"]),] 
numtows2 <- aggregate(ID~period+STRATUM,length,data=tows_unique2)
colnames(numtows2) <- c('period','STRATUM','COUNT')

#CALCULATING ABUND BY LEN CLASS FOR PERIODS
period <- seq(1,2)
LENBIN2 <- c(1,2,3,4,5,6,7,8,9)

stab.period <- ab.bin.strwt.nolog[ab.bin.strwt.nolog$period==1 | ab.bin.strwt.nolog$period==2,]
Alen.period.out <- matrix(0,length(period),length(LENBIN2))
colnames(Alen.period.out) <- as.character(LENBIN2)
rownames(Alen.period.out) <- as.character(period)
len_setup_strata <- aggregate(NUMLEN~period + LENBIN2 + STRATUM + STRATIO, sum, data=stab.period)
len_setup_strata <- merge(len_setup_strata,numtows)
len_setup_strata$mean.ABUNDANCE <- len_setup_strata$NUMLEN/len_setup_strata$COUNT
len_setup_strata$WTMEAN <- (len_setup_strata$mean.ABUNDANCE)*len_setup_strata$STRATIO

for(x in LENBIN2){
  for(p in period){
    data <- len_setup_strata[len_setup_strata$LENBIN2==x,]
    data <- data[data$period==p,]
    Alen.period.out[which(period==p),which(LENBIN2==x)] <- log10((sum(data$WTMEAN)/.01)*totstrwt)
  }
}
write.csv(Alen.period.out,"Alen.period.out")
save(Alen.period.out, file="Alen.period.out.Rdata")
Alen.period.out

#PLOTTING SIZE SPECTRUM--LENGTHS

plot(seq(1:9),seq(9,13,length=length(seq(1:9))),type='n')
points(Alen.period.out[1,],pch=16,col="black")
points(Alen.period.out[2,],pch=16,col="red")

##########################
#NOW DO ABOVE FOR WT

#COUNTING TOWS FOR STRATIFIED MEAN
#USES TOWS JUST WITH FISH SPECIES IN WIGLEY
tows_unique <- len.reg[!duplicated(len.reg["ID"]),] 
numtows <- aggregate(ID~EST_YEAR+STRATUM,length,data=tows_unique)
colnames(numtows) <- c('EST_YEAR','STRATUM','COUNT')

#SUM NUMLEN WITHIN EACH WTBIN FOR YEAR,STRATUM  

ab.wt10bin.str.nolog <- aggregate(NUMLEN ~ EST_YEAR + STRATUM + WT10BIN, sum, data=len.reg)

#MERGE IN STRATUM RATIOS
str.wt <- len.reg[,c(4,11,12)]
str.wt.unique <- str.wt[!duplicated(str.wt["STRATUM"]),]
ab.wt10bin.strwt.nolog <- merge(ab.wt10bin.str.nolog,str.wt.unique,by="STRATUM")
totstrwt <- sum(str.wt.unique$STRATUM_AREA)

#CALCULATING TOTAL ABUNDANCE BY YEAR
years <- seq(1968,2013)
WTBIN <- c(1:19)
Awt.out <- matrix(0,length(years),length(WTBIN))
colnames(Awt.out) <- as.character(WTBIN)
rownames(Awt.out) <- as.character(years)
wt_setup_strata <- aggregate(NUMLEN~EST_YEAR + WT10BIN + STRATUM + STRATIO, sum, data=ab.wt10bin.strwt.nolog)
wt_setup_strata <- merge(wt_setup_strata,numtows)
wt_setup_strata$mean.ABUNDANCE <- wt_setup_strata$NUMLEN/wt_setup_strata$COUNT
wt_setup_strata$WTMEAN <- (wt_setup_strata$mean.ABUNDANCE)*wt_setup_strata$STRATIO
head(wt_setup_strata)

#NOTE: 9/30/20:  In calc below, 0.01 is the standard swept area per tow in nautical miles I think.
#I don't remember the units on the stratum areas.  We should confirm these are in the same units.

for(x in WTBIN){
  for(y in years){
    data <- wt_setup_strata[wt_setup_strata$WT10BIN==x,]
    data <- data[data$EST_YEAR==y,]
    Awt.out[which(years==y),which(WTBIN==x)] <- log10((sum(data$WTMEAN)/.01)*totstrwt)
  }
}
write.csv(Awt.out,"Awt.out")
save(Awt.out, file="Awt.out.Rdata")
Awt.out

par(mfrow=c(1,1))
plot(seq(1968,2013),seq(2,10,length=length(seq(1968,2013))),type="n")
title("Abundance by weight class--all fish", cex=1.1)
lines(seq(1968,2013),Awt.out[,1],col="black") 
lines(seq(1968,2013),Awt.out[,2],col="red") 
lines(seq(1968,2013),Awt.out[,3],col="blue") 
lines(seq(1968,2013),Awt.out[,4],col="green")  
lines(seq(1968,2013),Awt.out[,5],col="orange") 
lines(seq(1968,2013),Awt.out[,6],col="purple") 
lines(seq(1968,2013),Awt.out[,7],col="grey")
lines(seq(1968,2013),Awt.out[,8],col="black") 
lines(seq(1968,2013),Awt.out[,9],col="red") 
lines(seq(1968,2013),Awt.out[,10],col="blue") 
lines(seq(1968,2013),Awt.out[,11],col="green") #above 1 g  
lines(seq(1968,2013),Awt.out[,12],col="orange") 
lines(seq(1968,2013),Awt.out[,13],col="purple") 
lines(seq(1968,2013),Awt.out[,14],col="grey") #stable
lines(seq(1968,2013),Awt.out[,15],col="black") #stable
lines(seq(1968,2013),Awt.out[,16],col="red") #declining
lines(seq(1968,2013),Awt.out[,17],col="blue") 
#lines(seq(1968,2013),Awt.out[,18],col="green")  
#lines(seq(1968,2013),Awt.out[,19],col="orange") 

#CLUSTERING

library(rioja)
library(vegan)
library(zoo)

#to 2013
agg.cl <- Awt.out[,11:17]
agg.cl <- rollmean(agg.cl,3)
agg.eucdist <- vegdist(agg.cl, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = TRUE) 
agg.eucdist.cl <- chclust(agg.eucdist, method = "coniss")
plot(agg.eucdist.cl,labels=seq(1969,2012),hang=-0.1,axes=F,cex=1)
axis(side=2, cex.axis=1.1)
title("Abundance by weight classes, 1968-2013",cex=1.1)
mtext(side=2, line=2.3, "Sum of squares", cex=1.1,las=0)
num <- bstick(agg.eucdist.cl)
num$diff <- num$dispersion-num$bstick
num.cl <- num[num$diff>=0,]
d <- dim(num.cl)
num.cl <- num.cl$nGroups[d[1]]
clus <- cutree(agg.eucdist.cl,k=num.cl)

#DEFINING PERIODS

ab.wt10bin.strwt.nolog$period <- 0
ab.wt10bin.strwt.nolog$period <- ifelse(ab.wt10bin.strwt.nolog$EST_YEAR<=1988,1,ab.wt10bin.strwt.nolog$period)
ab.wt10bin.strwt.nolog$period <- ifelse(ab.wt10bin.strwt.nolog$EST_YEAR>=1989 & 
                                    ab.wt10bin.strwt.nolog$EST_YEAR<=2009,2,ab.wt10bin.strwt.nolog$period)
ab.wt10bin.strwt.nolog$period <- ifelse(ab.wt10bin.strwt.nolog$EST_YEAR>=2010,3,ab.wt10bin.strwt.nolog$period)

len.reg$period <- 0
len.reg$period <- ifelse(len.reg$EST_YEAR<=1988,1,len.reg$period)
len.reg$period <- ifelse(len.reg$EST_YEAR>=1989 & len.reg$EST_YEAR<=2009,2,len.reg$period)
len.reg$period <- ifelse(len.reg$EST_YEAR>=2010,3,len.reg$period)

##ABUNDANCE BY PERIOD
#TOWS
fish.tows <- len.reg[len.reg$period==1 | len.reg$period==2 | len.reg$period==3,]
tows_unique2 <- fish.tows[!duplicated(fish.tows["ID"]),] 
numtows2 <- aggregate(ID~period+STRATUM,length,data=tows_unique2)
colnames(numtows2) <- c('period','STRATUM','COUNT')

#CALCULATING ABUND BY WT CLASS FOR PERIODS
period <- seq(1:3)
wt10bin <- c(11:17)

stab.period <- ab.wt10bin.strwt.nolog
Awt.period.out <- matrix(0,length(period),length(wt10bin))
colnames(Awt.period.out) <- as.character(wt10bin)
rownames(Awt.period.out) <- as.character(period)
wt_setup_strata <- aggregate(NUMLEN~period + WT10BIN + STRATUM + STRATIO, sum, data=stab.period)
wt_setup_strata <- merge(wt_setup_strata,numtows2)
wt_setup_strata$mean.ABUNDANCE <- wt_setup_strata$NUMLEN/wt_setup_strata$COUNT
wt_setup_strata$WTMEAN <- (wt_setup_strata$mean.ABUNDANCE)*wt_setup_strata$STRATIO

for(x in wt10bin){
  for(p in period){
    data <- wt_setup_strata[wt_setup_strata$WT10BIN==x,]
    data <- data[data$period==p,]
    Awt.period.out[which(period==p),which(wt10bin==x)] <- sum((data$WTMEAN) / .01) * totstrwt
  }
}
write.csv(Awt.period.out,"Awt.period.out")
save(Awt.period.out, file="Awt.period.out.Rdata")
Awt.period.out

avgAwt.period.out <- Awt.period.out
avgAwt.period.out[1,] <- Awt.period.out[1,]/20
avgAwt.period.out[2,] <- Awt.period.out[2,]/21
avgAwt.period.out[3,] <- Awt.period.out[3,]/3

LOG10avgAwt.period.out <- log10(avgAwt.period.out)

#PLOTTING SIZE SPECTRUM--WEIGHTS
par(mfrow=c(1,1))
plot(seq(11:17),seq(4,10,length=length(seq(11:17))),type='n', axes=F, xlab=" ", ylab=" ")
points(LOG10avgAwt.period.out[1,],pch=16,col="black")
points(LOG10avgAwt.period.out[2,],pch=16,col="red")
points(LOG10avgAwt.period.out[3,],pch=16,col="blue")
box()
axis(side=1,at=seq(11:17),labels=c(-1.75,1.25,0.75,0.25,-0.25,0.75,1.25),line=0,cex.axis=1.1)
axis(side=2,at=seq(4,10,length=length(seq(11:17))),labels=c(4,5,6,7,8,9,10),line=0,cex.axis=1.1)
mtext(side=1,line=2.5,"Log (mass)",cex=1.2)
mtext(side=2,line=2.5,"Log (average annual abundance)",cex=1.2)

##PLOTTING TIME SERIES BY YEARS

tows_unique <- len.reg[!duplicated(len.reg["ID"]),] 
numtows <- aggregate(ID~EST_YEAR+STRATUM,length,data=tows_unique)
colnames(numtows) <- c('EST_YEAR','STRATUM','COUNT')

#SUM NUMLEN WITHIN EACH WT10BIN FOR YEAR,STRATUM  
sz.class <- len.reg[len.reg$WT10BIN==11 | len.reg$WT10BIN==12 | len.reg$WT10BIN==13 | len.reg$WT10BIN==16 |
                len.reg$WT10BIN==17,]
sz.class$class <- 0
sz.class$class <- ifelse(sz.class$WT10BIN==11 | sz.class$WT10BIN==12 | sz.class$WT10BIN==13, 1, 
                    sz.class$class)  
sz.class$class <- ifelse(sz.class$WT10BIN==16 | sz.class$WT10BIN==17,2,sz.class$class)
ab.szclass <- aggregate(NUMLEN ~ EST_YEAR + STRATUM + class, sum, data=sz.class)

#MERGE IN STRATUM RATIOS
str.wt <- len.reg[,c(4,11,12)]
str.wt.unique <- str.wt[!duplicated(str.wt["STRATUM"]),]
ab.szclass <- merge(ab.szclass,str.wt.unique,by="STRATUM")
totstrwt <- sum(str.wt.unique$STRATUM_AREA)

#CALCULATING TOTAL ABUNDANCE BY YEAR FOR EACH CLASS
years <- seq(1968,2013)
class <- c(1,2)
Aclass.out <- matrix(0,length(years),length(class))
colnames(Aclass.out) <- as.character(class)
rownames(Aclass.out) <- as.character(years)
szclass_strata <- aggregate(NUMLEN~EST_YEAR + class + STRATUM + STRATIO, sum, data=ab.szclass)
szclass_strata <- merge(szclass_strata,numtows)
szclass_strata$mean.ABUNDANCE <- szclass_strata$NUMLEN/szclass_strata$COUNT
szclass_strata$WTMEAN <- (szclass_strata$mean.ABUNDANCE)*szclass_strata$STRATIO

for(x in class){
  for(y in years){
    data <- szclass_strata[szclass_strata$class==x,]
    data <- data[data$EST_YEAR==y,]
    Aclass.out[which(years==y),which(class==x)] <- log10((sum(data$WTMEAN)/.01)*totstrwt)
  }
}
write.csv(Aclass.out,"Aclass.out")
save(Aclass.out, file="Aclass.out.Rdata")
Aclass.out

par(mfrow=c(1,1))
plot(seq(1968,2013),seq(6,10.5,length=length(seq(1968,2013))),type="n",xlab=" ", ylab=" ")
title("Abundance by size class", cex=1.1)
lines(seq(1968,2013),Aclass.out[,1],col="black",lwd=2) 
lines(seq(1968,2013),Aclass.out[,2],col="red",lwd=2) 
mtext(side=2,line=2.5,"Log (Abundance)",cex=1.2)
mtext(side=1,line=2.5,"Year",cex=1.2)

##SIZE SPECTRUM SLOPE AND INTERCEPT

#CUTTING TO THOSE LARGER THAN 0.1 KILOGRAMS
Awt.reg <- Awt.out[,c(11:17)]
as.data.frame(Awt.reg)
dim(Awt.reg)

#rs.reg.yr.lg$labund <- log(rs.reg.yr.lg$abund)
#rs.reg.yr.lg$labund[rs.reg.yr.lg$labund=="-Inf"] <- NA

row <- seq(1,46,by=1)
yr.slopeint <- matrix(ncol=7)

for (i in row){
  data <- Awt.reg[i,]
  x <- c(-1.75,-1.25,-.75,-.25,.25,.75,1.25)
  lm <- lm(data~x,na.action=na.omit)
  coef <- coef(lm)
  se <- confint(lm)
  yr.slopeint <- rbind(yr.slopeint,c(i,coef[1],se[1,1],se[1,2],coef[2],se[2,1],se[2,2]))
}

yr.slopeint <- as.data.frame(yr.slopeint)
yr.slopeint <- yr.slopeint[2:47,2:7]
yr <- seq(1968,2013,by=1)
yr.slopeint <- cbind(yr,yr.slopeint)
colnames(yr.slopeint) <- c("year","int","l.int","h.int","slope","l.slope","h.slope")
write.csv(yr.slopeint,file="yr.slopeint")

plot(yr.slopeint$year,yr.slopeint$int,type='l')
plot(yr.slopeint$year,yr.slopeint$slope,type='l')

#PLOTTING BY YEAR
attach(yr.slopeint)
par(mfrow=c(2,2),mar=c(4,5,4,2))
#SLOPE
plot(seq(1970,2013,1),seq(min(yr.slopeint$slope,na.rm=TRUE),max(yr.slopeint$slope,na.rm=TRUE),
                          length.out=length(seq(1970,2013,1))),type="n",xlab=" ",ylab=" ",cex=1.3,cex.axis=1.3)
mtext(side=1,line=3,"Year",cex=1.3)
mtext(side=2,line=3,"Slope",cex=1.3)
points(yr.slopeint$year,yr.slopeint$slope,pch=20,col="black")
#lines(yr.slopeint$year,yr.slopeint$l.slope,lwd=1,col="grey")
#lines(yr.slopeint$year,yr.slopeint$h.slope,lwd=1,col="grey")
lines(rollmean(yr.slopeint$year,5),rollmean(yr.slopeint$slope,5),lwd=3,col="blue")

#INTERCEPT
plot(seq(1970,2013,1),seq(min(yr.slopeint$int,na.rm=TRUE),max(yr.slopeint$int,na.rm=TRUE),
                          length.out=length(seq(1970,2013,1))),type="n",xlab=" ",ylab=" ",cex=1.3,cex.axis=1.3)
mtext(side=1,line=3,"Year",cex=1.3)
mtext(side=2,line=3,"Intercept",cex=1.3)
points(yr.slopeint$year,yr.slopeint$int,pch=20,col="black")
#lines(yr.slopeint$year,yr.slopeint$l.int,lwd=1,col="grey")
#lines(yr.slopeint$year,yr.slopeint$h.int,lwd=1,col="grey")
lines(rollmean(yr.slopeint$year,5),rollmean(yr.slopeint$int,5),lwd=3,col="blue")

#PRODUCTIVITY REQUIRED

years <- seq(1968,2013)
class <- c(1,2)
Aclass.out <- matrix(0,length(years),length(class))

yr.prod <- matrix(0,length(years),6)
for(i in years){
  yr.out <- yr.slopeint[yr.slopeint$year==i,]
  a <- yr.out$int
  b <- yr.out$slope
  m0 <- 10^(-1.75)
  m1 <- 10^(1.75)
  prodreq <- (a/(7/4-b))*((m1^(7/4-b))-(m0^(7/4-b)))
  yr.prod[i-1967,] <- cbind(i,a,b,m0,m1,prodreq)
}

yr.prod <- as.data.frame(yr.prod)
colnames(yr.prod) <- c("year","a","b","m0","m1","prod")
write.csv(yr.prod,file="yr.prod")

par(mfrow=c(1,1))
plot(yr.prod$year,yr.prod$prod,type="l",xlab="Year",ylab="Productivity required")
plot(yr.prod$year,log10(yr.prod$prod),type="l",xlab="Year",ylab="Log productivity required",lwd=2,
     cex.axis=1.1)

#ADDING IN T

annT <- aggregate(sst ~ year, mean, data=GoMERSST)
yr.prod.T <- merge(yr.prod,annT,by="year",all.x=T)

names(yr.prod.T)
T1 <- mean(yr.prod.T$sst[1:21])
T2 <- mean(yr.prod.T$sst[22:42])
T3 <- mean(yr.prod.T$sst[43:46])
         
         T1
         T2
         T3

E <- 0.63
k <- 0.0008617  
 
exp(-E/k*(273+T1))-exp(-E/(k*273))

?exp


years <- seq(1968,2013)
#class <- c(1,2)
#AclassTperiod.out <- matrix(0,length(years),length(class))

yr.prod.T <- matrix(0,length(years),6)
for(i in years){
  yrT.out <- yr.slopeint[yr.slopeint$year==i,]
  a <- yrT.out$int
  b <- yrT.out$slope
  m0 <- 10^(-1.75)
  m1 <- 10^(1.75)
  prodreq <- (a/(7/4-b))*((m1^(7/4-b))-(m0^(7/4-b)))
  yr.prod[i-1967,] <- cbind(i,a,b,m0,m1,prodreq)
}



yr.prod <- as.data.frame(yr.prod)
colnames(yr.prod) <- c("year","a","b","m0","m1","prod")
write.csv(yr.prod,file="yr.prod")

par(mfrow=c(1,1))
plot(yr.prod$year,yr.prod$prod,type="l",xlab="Year",ylab="Productivity required")
plot(yr.prod$year,log10(yr.prod$prod),type="l",xlab="Year",ylab="Log productivity required",lwd=2,
     cex.axis=1.1)


plot(seq[1:9],seq(0,max(reg.out2$labund.scale,na.rm=TRUE),length.out=length(seq[1:12])),
     type="n",xlab=" ",ylab=" ",cex=1.3,cex.axis=1.3)
mtext(side=1,line=3,"log(mass)",cex=1.3)
mtext(side=2,line=3,"log(abundance)",cex=1.3)
for(i in decade){
  dec <- reg.out2[reg.out2$decade==i,]
  points(dec$bin,dec$labund.scale,pch=20,col=i)
  #slopeint <- dec.slopeint[dec.slopeint$decade==i,]
  #abline(slopeint[,2],slopeint[,5],col=i)
  lines(dec$bin,predict(lm(labund.scale~poly(bin,2),data=dec)),lwd=2,col=i)
