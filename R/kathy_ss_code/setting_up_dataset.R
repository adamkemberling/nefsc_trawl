library(tidyverse)
library(here)
library(gmRi)

# 2020 data
load(here("data/NEFSC/Survdat_Nye_Aug 2020.RData"))
mills_path <- shared.path(group = "Mills Lab", folder = "")
stratum.area <- read_csv(str_c(mills_path, "Projects/NSF_CAccel/Data/strata area.csv"), 
                       col_types = cols())
#setwd('C:/Users/Kathy/Documents/Active Research/CoastalSEES/Data')
# load('Survdat_Nye.Rdata')   #NOTE: Once loaded, file is in directory as "survdat"
# load('stratum.area.Rdata')
# setwd('C:/Users/Kathy/Documents/Active Research/Lenfest/Ecosystem change analyses/Data')
dim(survdat)



names(survdat)
bn.setup<-survdat[,c(1,2,3,8,25,26,30,31,33,34,44,45)]
bn.setup<-merge(bn.setup,stratum.area)
bn.setup<-bn.setup[order(bn.setup$ID),]

dim(bn.setup)

#SELECT YEARS, SEASON
bn.setup<-bn.setup[bn.setup$SEASON=='SPRING' | bn.setup$SEASON=='FALL',]
bn.setup<-bn.setup[bn.setup$EST_YEAR>=1968,]
bn.setup<-bn.setup[bn.setup$EST_YEAR<2014,]

dim(bn.setup)

#SELECT STRATA
#GOM and GB
bn.setup<-bn.setup[bn.setup$STRATUM>=01130 & bn.setup$STRATUM<=01400,]  
#ELIMINATE CANADIAN STRATA
bn.setup<-bn.setup[bn.setup$STRATUM!=01350 & bn.setup$STRATUM!=1310 & bn.setup$STRATUM!=1351
                   & bn.setup$STRATUM!=1320 & bn.setup$STRATUM!=1410 & bn.setup$STRATUM!=1420 
                   & bn.setup$STRATUM!=1490,] 

dim(bn.setup)
names(bn.setup)

#FIND TOTAL STRATA AREA AND CREATE "WEIGHTS" COLUMN
allstr<-bn.setup[,c(1,13)]
allstr<-allstr[!duplicated(allstr["STRATUM"]),]
totstrwt<-sum(allstr$STRATUM_AREA)
bn.setup$STRATIO<-bn.setup$STRATUM_AREA/totstrwt

#COMPUTE INDIV WT
bn.setup$BIOM.ADJ<-ifelse(bn.setup$BIOMASS==0 & bn.setup$ABUNDANCE>0,0.0001,bn.setup$BIOMASS)
bn.setup$ABUND.ADJ<-ifelse(bn.setup$ABUNDANCE==0 & bn.setup$BIOMASS>0,1,bn.setup$ABUNDANCE)
bn.setup$INDWT<-bn.setup$BIOM.ADJ/bn.setup$ABUND.ADJ

sort(unique(bn.setup$SEASON))
sort(unique(bn.setup$STRATIO))
#######################

bn.setup2<-bn.setup

#EXCLUDE THE SHRIMPS
bn.setup2<-bn.setup2[bn.setup2$SVSPP<285 | bn.setup2$SVSPP>299,]
bn.setup2<-bn.setup2[bn.setup2$SVSPP!=305 & bn.setup2$SVSPP!=306 & bn.setup2$SVSPP!=307 &
                       bn.setup2$SVSPP!=316 & bn.setup2$SVSPP!=323,]
bn.setup2<-bn.setup2[bn.setup2$SVSPP<910 | bn.setup2$SVSPP>915,]
bn.setup2<-bn.setup2[bn.setup2$SVSPP<955 | bn.setup2$SVSPP>961,]

#EXCLUDE UNIDENTIFIED FISH (OFTEN DON'T HAVE BIOMASS)
bn.setup2<-bn.setup2[bn.setup2$SVSPP!=0 & bn.setup2$SVSPP!=978 & bn.setup2$SVSPP!=979 &
                       bn.setup2$SVSPP!=980 & bn.setup2$SVSPP!=998,]

dim(bn.setup2)

#ASSIGN TAXONOMIC GROUPS
bn.setup2$SVSPP2<-bn.setup2$SVSPP
bn.setup2$tax<-NA
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=69 & bn.setup2$SVSPP2<=110,"gf",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2==773 | bn.setup2$SVSPP2==454 | bn.setup2$SVSPP2==117 | 
                        bn.setup2$SVSPP2==795 | bn.setup2$SVSPP2==866,"gf",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=003 & bn.setup2$SVSPP2<=029,"el",bn.setup2$tax)	
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2==601,"el",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=030 & bn.setup2$SVSPP2<=046,"pel",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=501 & bn.setup2$SVSPP2<=512,"pel",bn.setup2$tax)	
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2==131 |  bn.setup2$SVSPP2==135 | bn.setup2$SVSPP2==139 | 
                        bn.setup2$SVSPP2==121 | bn.setup2$SVSPP2==122 | bn.setup2$SVSPP2==123 | bn.setup2$SVSPP2==896 | 
                        bn.setup2$SVSPP2==246 | bn.setup2$SVSPP2==002 | bn.setup2$SVSPP2==067 | bn.setup2$SVSPP2==112 | 
                        bn.setup2$SVSPP2==113 | bn.setup2$SVSPP2==115 | bn.setup2$SVSPP2==205 | bn.setup2$SVSPP2==212 | 
                        bn.setup2$SVSPP2==240 | bn.setup2$SVSPP2==248 | bn.setup2$SVSPP2==251 | bn.setup2$SVSPP2==252 | 
                        bn.setup2$SVSPP2==262 | bn.setup2$SVSPP2==466 | bn.setup2$SVSPP2==856 | bn.setup2$SVSPP2==894 | 
                        bn.setup2$SVSPP2==620 | bn.setup2$SVSPP2==145 | bn.setup2$SVSPP2==054 | bn.setup2$SVSPP2==055 | 
                        bn.setup2$SVSPP2==056 | bn.setup2$SVSPP2==066 | bn.setup2$SVSPP2==120 | bn.setup2$SVSPP2==124 | 
                        bn.setup2$SVSPP2==129 | bn.setup2$SVSPP2==132 | bn.setup2$SVSPP2==133 | bn.setup2$SVSPP2==204 | 
                        bn.setup2$SVSPP2==208 | bn.setup2$SVSPP2==209 | bn.setup2$SVSPP2==228 | bn.setup2$SVSPP2==421 | 
                        bn.setup2$SVSPP2==490 | bn.setup2$SVSPP2==576 | bn.setup2$SVSPP2==582 | bn.setup2$SVSPP2==631 | 
                        bn.setup2$SVSPP2==820 | bn.setup2$SVSPP2==832 | bn.setup2$SVSPP2==851 | bn.setup2$SVSPP2==865 | 
                        bn.setup2$SVSPP2==053 | bn.setup2$SVSPP2==210 | bn.setup2$SVSPP2==876 | bn.setup2$SVSPP2==887 | 
                        bn.setup2$SVSPP2==207 | bn.setup2$SVSPP2==245 | bn.setup2$SVSPP2==750 | bn.setup2$SVSPP2==220 | 
                        bn.setup2$SVSPP2==520 | bn.setup2$SVSPP2==61 | bn.setup2$SVSPP2==127 | bn.setup2$SVSPP2==250 |
                        bn.setup2$SVSPP2==617 | bn.setup2$SVSPP2==694 | bn.setup2$SVSPP2==874, "pel", bn.setup2$tax)	
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=154 & bn.setup2$SVSPP2<=164,"dem",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=191 & bn.setup2$SVSPP2<=197,"dem",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2==001 | bn.setup2$SVSPP2==060 | bn.setup2$SVSPP2==063 | 
                        bn.setup2$SVSPP2==064 | bn.setup2$SVSPP2==384 | bn.setup2$SVSPP2==114 | bn.setup2$SVSPP2==116 | 
                        bn.setup2$SVSPP2==136 | bn.setup2$SVSPP2==141 | bn.setup2$SVSPP2==165 | bn.setup2$SVSPP2==167 | 
                        bn.setup2$SVSPP2==168 | bn.setup2$SVSPP2==170 | bn.setup2$SVSPP2==171 | bn.setup2$SVSPP2==172 | 
                        bn.setup2$SVSPP2==177 | bn.setup2$SVSPP2==180 | bn.setup2$SVSPP2==181 | bn.setup2$SVSPP2==183 | 
                        bn.setup2$SVSPP2==184 | bn.setup2$SVSPP2==185 | bn.setup2$SVSPP2==189 | bn.setup2$SVSPP2==190 | 
                        bn.setup2$SVSPP2==201 | bn.setup2$SVSPP2==213 | bn.setup2$SVSPP2==229 | bn.setup2$SVSPP2==232 | 
                        bn.setup2$SVSPP2==182 | bn.setup2$SVSPP2==166 | bn.setup2$SVSPP2==759 | bn.setup2$SVSPP2==880 | 
                        bn.setup2$SVSPP2==734 | bn.setup2$SVSPP2==144 | bn.setup2$SVSPP2==047 | bn.setup2$SVSPP2==065 | 
                        bn.setup2$SVSPP2==111 | bn.setup2$SVSPP2==169 | bn.setup2$SVSPP2==173 | bn.setup2$SVSPP2==174 | 
                        bn.setup2$SVSPP2==175 | bn.setup2$SVSPP2==199 | bn.setup2$SVSPP2==202 | bn.setup2$SVSPP2==206 | 
                        bn.setup2$SVSPP2==211 | bn.setup2$SVSPP2==264 | bn.setup2$SVSPP2==267 | bn.setup2$SVSPP2==425 | 
                        bn.setup2$SVSPP2==461 | bn.setup2$SVSPP2==557 | bn.setup2$SVSPP2==596 | bn.setup2$SVSPP2==735 | 
                        bn.setup2$SVSPP2==748 | bn.setup2$SVSPP2==797 | bn.setup2$SVSPP2==798 | bn.setup2$SVSPP2==825 | 
                        bn.setup2$SVSPP2==862 | bn.setup2$SVSPP2==134 | bn.setup2$SVSPP2==138 | bn.setup2$SVSPP2==597 | 
                        bn.setup2$SVSPP2==733 | bn.setup2$SVSPP2==852 | bn.setup2$SVSPP2==867 | bn.setup2$SVSPP2==869 | 
                        bn.setup2$SVSPP2==461 | bn.setup2$SVSPP2==557 | bn.setup2$SVSPP2==596 | bn.setup2$SVSPP2==489 | 
                        bn.setup2$SVSPP2==126 | bn.setup2$SVSPP2==142 | bn.setup2$SVSPP2==242 | bn.setup2$SVSPP2==280 | 
                        bn.setup2$SVSPP2==176 | bn.setup2$SVSPP2==143 | bn.setup2$SVSPP2==488 | bn.setup2$SVSPP2==151 | 
                        bn.setup2$SVSPP2==188 | bn.setup2$SVSPP2==231 | bn.setup2$SVSPP2==237 | bn.setup2$SVSPP2==260 | 
                        bn.setup2$SVSPP2==450 | bn.setup2$SVSPP2==452 | bn.setup2$SVSPP2==556 | bn.setup2$SVSPP2==559 | 
                        bn.setup2$SVSPP2==599 | bn.setup2$SVSPP2==230 | bn.setup2$SVSPP2==263 | bn.setup2$SVSPP2==546 | 
                        bn.setup2$SVSPP2==794 | bn.setup2$SVSPP2==221 | bn.setup2$SVSPP2==616 | bn.setup2$SVSPP2==390 |
                        bn.setup2$SVSPP2==137 | bn.setup2$SVSPP2==146 | bn.setup2$SVSPP2==239 | bn.setup2$SVSPP2==249 |
                        bn.setup2$SVSPP2==500 | bn.setup2$SVSPP2==526 | bn.setup2$SVSPP2==549 | bn.setup2$SVSPP2==595 |
                        bn.setup2$SVSPP2==838 | bn.setup2$SVSPP2==899, "dem",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=301 & bn.setup2$SVSPP2<=324,"inv",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2>=291 & bn.setup2$SVSPP2<=297,"inv",bn.setup2$tax)
bn.setup2$tax<-ifelse(bn.setup2$SVSPP2==401 | bn.setup2$SVSPP2==298 | bn.setup2$SVSPP2==348 | 
                        bn.setup2$SVSPP2==409 | bn.setup2$SVSPP2==285 | bn.setup2$SVSPP2==287 | bn.setup2$SVSPP2==299 | 
                        bn.setup2$SVSPP2==325 | bn.setup2$SVSPP2==340 | bn.setup2$SVSPP2==342 | bn.setup2$SVSPP2==517 | 
                        bn.setup2$SVSPP2==330 | bn.setup2$SVSPP2==331 | bn.setup2$SVSPP2==332 | bn.setup2$SVSPP2==335 | 
                        bn.setup2$SVSPP2==338 | bn.setup2$SVSPP2==400 | bn.setup2$SVSPP2==403 | bn.setup2$SVSPP2==604 |
                        bn.setup2$SVSPP2==402 | bn.setup2$SVSPP2==518, "inv",bn.setup2$tax)

test<-bn.setup2[is.na(bn.setup2$tax),]
dim(test)
test
sort(unique(bn.setup2$tax))

#ASSIGN ECONOMIC GROUPS
bn.setup2$econ<-"nc"
bn.setup2$econ<-ifelse(bn.setup2$SVSPP2==069 | bn.setup2$SVSPP2==072 | bn.setup2$SVSPP2==073 | 
                         bn.setup2$SVSPP2==074 | bn.setup2$SVSPP2==075 | bn.setup2$SVSPP2==076 | bn.setup2$SVSPP2==077 | 
                         bn.setup2$SVSPP2==101 | bn.setup2$SVSPP2==102 | bn.setup2$SVSPP2==105 | bn.setup2$SVSPP2==106 | 
                         bn.setup2$SVSPP2==107 | bn.setup2$SVSPP2==108 | bn.setup2$SVSPP2==155 | bn.setup2$SVSPP2==193 | 
                         bn.setup2$SVSPP2==401 | bn.setup2$SVSPP2==197 | bn.setup2$SVSPP2==032 | bn.setup2$SVSPP2==015 | 
                         bn.setup2$SVSPP2==310 | bn.setup2$SVSPP2==024 | bn.setup2$SVSPP2==025 | bn.setup2$SVSPP2==026 | 
                         bn.setup2$SVSPP2==027 | bn.setup2$SVSPP2==028 | bn.setup2$SVSPP2==023 | bn.setup2$SVSPP2==022 | 
                         bn.setup2$SVSPP2==894 | bn.setup2$SVSPP2==384 | bn.setup2$SVSPP2==301 | bn.setup2$SVSPP2==136 | 
                         bn.setup2$SVSPP2==036 | bn.setup2$SVSPP2==380 | bn.setup2$SVSPP2==141 | bn.setup2$SVSPP2==135 | 
                         bn.setup2$SVSPP2==318 | bn.setup2$SVSPP2==654 | bn.setup2$SVSPP2==143 | bn.setup2$SVSPP2==035 | 
                         bn.setup2$SVSPP2==037 | bn.setup2$SVSPP2==033 | bn.setup2$SVSPP2==034 | bn.setup2$SVSPP2==745 | 
                         bn.setup2$SVSPP2==013 | bn.setup2$SVSPP2==360 | bn.setup2$SVSPP2==928 | bn.setup2$SVSPP2==354 | 
                         bn.setup2$SVSPP2==364 | bn.setup2$SVSPP2==355 | bn.setup2$SVSPP2==358 | bn.setup2$SVSPP2==357 | 
                         bn.setup2$SVSPP2==931 | bn.setup2$SVSPP2==356 | bn.setup2$SVSPP2==359 | bn.setup2$SVSPP2==350 | 
                         bn.setup2$SVSPP2==362 | bn.setup2$SVSPP2==363 | bn.setup2$SVSPP2==366 | bn.setup2$SVSPP2==352 | 
                         bn.setup2$SVSPP2==353 | bn.setup2$SVSPP2==925 | bn.setup2$SVSPP2==017 | bn.setup2$SVSPP2==149 | 
                         bn.setup2$SVSPP2==645 | bn.setup2$SVSPP2==139 | bn.setup2$SVSPP2==103 | bn.setup2$SVSPP2==177 | 
                         bn.setup2$SVSPP2==145 | bn.setup2$SVSPP2==930 | bn.setup2$SVSPP2==121 | bn.setup2$SVSPP2==147 |
                         bn.setup2$SVSPP2==935 | bn.setup2$SVSPP2==192 | bn.setup2$SVSPP2==151 | bn.setup2$SVSPP2==403 |
                         bn.setup2$SVSPP2==409 | bn.setup2$SVSPP2==131 | bn.setup2$SVSPP2==502 | bn.setup2$SVSPP2==503, 
                       "comm", bn.setup2$econ)
sort(unique(bn.setup2$econ))

head(bn.setup2)

write.csv(bn.setup2,"bn.setup2")
save(bn.setup2, file = "bn.setup2.RData")

#UNIQUE TOWS PER YEAR
tows_unique<-bn.setup[!duplicated(bn.setup["ID"]),] 
numtows<-aggregate(ID~EST_YEAR+STRATUM,length,data=tows_unique)
colnames(numtows)<-c('EST_YEAR','STRATUM','COUNT')

write.csv(numtows,"numtows")
save(numtows, file = "numtows.RData")


