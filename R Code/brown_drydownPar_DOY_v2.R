#' Creates time series of parameters for calendar DOY
library(raster)
library(terra)

library(ncdf4)
#~~~ Create month-wise time series of drydown parameters
path_name="G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/test_regions/"
season_names=as.character(1:12)
param_names=c("m2","thetaTD","thetaWT")

#test=brick("G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/test_regions/Regional_SPL3_E_v7_NLSparBrk_R3_S4.nc")
retrieval_date=seq(as.Date("2015-03-31"),as.Date("2021-09-30"), by=2)
currRegion = 3

smapBrk = brick(paste("C:/Users/Brown/Downloads/Region_SPL3SMP_E_2dayInterp_R",currRegion,".nc", sep=""))

for (par_run in 1:3){
  eval(parse(text=paste('season_',param_names[par_run],'Brk=brick()' ,sep = "")))
  for (season_run in 1:12){
    # Convert seasonal NC files to RasterBricks    
    text_run=paste('SeasonParBrk=brick("',path_name,'Regional_SPL3_E_v7_NLSparBrk_R',currRegion,'_S',season_run,'.nc")',sep="")
    eval(parse(text=text_run))
    layerNames=c("canonical_form","thetaTD","thetaWT","thetaGW","ld","lw","m1","m2",
                 "sd_thetaTD","sd_thetaWT","sd_thetaGW","sd_ld","sd_lw","sd_m1","sd_m2",
                 "mse","CC","d","replacement","nSamp")
    names(SeasonParBrk)=layerNames
    
    # Add layer to parameter stack
    text_run2=paste('season_',param_names[par_run],'Brk=addLayer(','season_',param_names[par_run],'Brk,SeasonParBrk$',param_names[par_run],')', sep="")
    eval(parse(text=text_run2))
  }
}

#plot(rast(season_m2Brk))
#plot(rast(season_thetaTDBrk))
#plot(rast(season_thetaWTBrk))


#~~ Interpolate missing values

# Step1: Replace missing values with seasonal averages of parameters
interpM2brk=approxNA(season_m2Brk,method="linear", rule=2,NArule=1)
interpThetaTDbrk=approxNA(season_thetaTDBrk,method="linear", rule=2,NArule=1)

# interpThetaWTbrk=approxNA(season_thetaWTBrk,method="linear", rule=2,NArule=2)
# Step2: Replace missing value using seasonal max SM

#~~ a) Rolling seasonal median SM
library(parallel)
library(ClusterR)

#~~ b) Rolling seasonal max SM
rollSeasonMax=function(x){
  x=as.numeric(as.vector(x))
  if(!all(is.na(x))){
    monthTriplet=function (x){
      mseq=c(12,seq(1:12),1)
      ifelse((x !=12),return(c(mseq[which(mseq==x)[1]-1], mseq[which(mseq==x)[1]], mseq[which(mseq==x)[1]+1])),return(c(11,12,1)))
    }
    mthInd=lubridate::month(retrieval_date)
    seasonal_mean=c()
    for (runMth in c(1:12)){
      seasonal_mean[runMth]=max(x[which(mthInd %in% monthTriplet(runMth))],na.rm=TRUE)  
    }
    return(unlist(as.numeric(seasonal_mean)))
  } else{
    return(rep(NA,12))
  }
}

numCores=detectCores()  
beginCluster(n=numCores-1)    
start_time <- Sys.time() 
seasonMaxSM=clusterR(smapBrk, 
                     calc, args=list(fun=rollSeasonMax),
                     export=c("retrieval_date"),
                     progress='text');endCluster()
end_time <- Sys.time()
end_time - start_time  

#~~ c) Estimate missing values of ThetaWT
interpThetaWTbrk=season_thetaWTBrk
values(interpThetaWTbrk)[which(is.na(values(interpThetaWTbrk)))]=1.05*(values(seasonMaxSM)[which(is.na(values(interpThetaWTbrk)))])
#values(interpThetaWTbrk)[which(is.na(values(interpM2brk)))]=NA

interpThetaWTbrk=resample(interpThetaWTbrk,seasonMaxSM[[1]],method='bilinear')
interpThetaTDbrk=resample(interpThetaTDbrk,seasonMaxSM[[1]],method='bilinear')
interpM2brk=resample(interpM2brk,seasonMaxSM[[1]],method='bilinear')

#~~ Create DOY time series of the parameters
# Step1: Arrange data from seasonal stack
sampDate=seq(as.Date("2016-01-01", format="%Y-%d-%m"),as.Date("2016-31-12", format="%Y-%d-%m"),by=1)
sampSeason=lubridate::month(sampDate)

m2ts=subset(interpM2brk,sampSeason)
thetaWTts=subset(interpThetaWTbrk,sampSeason)
thetaTDts=subset(interpThetaTDbrk,sampSeason)
  

  
# Step2: Apply rolling mean
myrollfun=function(x){
  if(!all(is.na(x))){
    out=movingFun(x,31, mean,circular=TRUE,type='around')
  } else{
    out=rep(NA, length(x))
  }
  return(out)
}

library(ClusterR)
beginCluster()
start_time <- Sys.time()
m2doy=clusterR(m2ts, calc, args=list(fun=myrollfun),progress='text')
thetaWTdoy=clusterR(thetaWTts, calc, args=list(fun=myrollfun),progress='text')
thetaTDdoy=clusterR(thetaTDts, calc, args=list(fun=myrollfun),progress='text')
end_time <- Sys.time()
endCluster(); end_time - start_time #3.5 min 

# Calculate ThetaIP
thetaIPts= calc(stack(thetaWTdoy,thetaTDdoy), 
                function(x) {rowMeans(cbind(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)]))}, progress="text")


# m2out = as.numeric(extract(m2doy, SpatialPoints(cbind(-90,-55)), method = 'simple'))
# WTout = as.numeric(extract(thetaWTdoy, SpatialPoints(cbind(-100,40)), method = 'simple'))
# TDout = as.numeric(extract(thetaTDdoy, SpatialPoints(cbind(-100,40)), method = 'simple'))
# IPout = as.numeric(extract(thetaIPts, SpatialPoints(cbind(-100,40)), method = 'simple'))
# 
# y1 = 0
# y2 = 0.38
# # first plot
# plot(m2out, ylim=range(c(y1,y2)), col ='red')
# 
# # second plot  EDIT: needs to have same ylim
# par(new = TRUE)
# plot(WTout, ylim=range(c(y1,y2)), axes = FALSE, xlab = "", ylab = "", col ='blue')
# par(new = TRUE)
# plot(TDout, ylim=range(c(y1,y2)), axes = FALSE, xlab = "", ylab = "", col ='green')
# par(new = TRUE)
# plot(IPout, ylim=range(c(y1,y2)), axes = FALSE, xlab = "", ylab = "", col ='purple')
# legend(y1,y2, legend=c("m2", "WT", "TD", "IP"),
#               col=c("blue", "blue","green", "purple"))
path = "G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/paraDoy/"
writeRaster(m2doy, paste(path, "m2doy_R",currRegion,".nc", sep=""), overwrite=TRUE)
writeRaster(thetaWTdoy, paste(path, "thetaWTdoy_R",currRegion,".nc", sep=""), overwrite=TRUE)
writeRaster(thetaTDdoy, paste(path, "thetaTDdoy_R",currRegion,".nc", sep=""), overwrite=TRUE)
writeRaster(thetaIPts, paste(path, "thetaIPts_R",currRegion,".nc", sep=""), overwrite=TRUE)

