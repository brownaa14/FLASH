#Creates time series of parameters for calendar DOY

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
regionRun = currRegion

smapBrk = brick(paste("C:/Users/Brown/Documents/FLASH/9kmV2/Region_SPL3SMP_E_2dayInterp_R",currRegion,".nc", sep=""))

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
  
  remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    #1.0 is hyperparameter
    
    H <- 1.0 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  x=remove_outliers(x)
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

path = "C:/Users/Brown/Documents/FLASH/9kmV2/paraDoy/"
writeRaster(m2doy, paste(path, "m2doy_R",currRegion,".nc", sep=""), overwrite=TRUE)
writeRaster(thetaWTdoy, paste(path, "thetaWTdoy_R",currRegion,".nc", sep=""), overwrite=TRUE)
writeRaster(thetaTDdoy, paste(path, "thetaTDdoy_R",currRegion,".nc", sep=""), overwrite=TRUE)
writeRaster(thetaIPts, paste(path, "thetaIPts_R",currRegion,".nc", sep=""), overwrite=TRUE)


library(parallel)
library(ClusterR)
library(raster)
library(terra)
library(lubridate)
library(viridisLite)


###############################################
# PART 2: Calculate SMS
###############################################
regionRun=currRegion
#' For region # 3, thetaTDdoy doesnt seems to be made correctly
thetaWTdoy=brick(paste(path,"thetaWTdoy_R",regionRun,".nc",sep=""))
thetaTDdoy=brick(paste(path,"thetaTDdoy_R",regionRun,".nc",sep=""))
m2doy=brick(paste(path,"m2doy_R",regionRun,".nc",sep=""))
thetaIPdoy=brick(paste(path,"thetaIPts_R",regionRun,".nc",sep=""))

interp_brk = brick(paste("C:/Users/Brown/Documents/FLASH/9kmV2/Region_SPL3SMP_E_1dayInterp_R",regionRun,".nc", sep=""))
interp_dates=seq(as.Date("2015-03-31"),as.Date("2021-09-30"), by=1)
# thetaIPts= calc(stack(thetaWTdoy,thetaTDdoy), 
# function(x) {rowMeans(cbind(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)]))}, progress="text")

#~~ Calculate SMS 
# STEP1: Calculate (SM/thetaIP) ratio
thetaIPruntime=brick(subset(thetaIPdoy, lubridate::yday(interp_dates)))
numCores=detectCores()-1; beginCluster(n=numCores)
ratioBrk=clusterR(stack(interp_brk,thetaIPruntime), 
                  calc, args=list(fun=function(x) {x[1:(length(x)/2)]/x[(length(x)/2+1):length(x)]}),
                  #export=c(" "," "),
                  progress='text');endCluster()

# Step2: Calculate SMS
m2runtime=subset(m2doy, lubridate::yday(interp_dates))
# CHANGE    #1
nRuntime=12*sqrt(m2runtime)
smsBrk=1/(1+(ratioBrk^nRuntime))

# tsOut=as.numeric(extract(smsBrk, cbind(-115,40)))
# plot(tsOut, type="l")
# plot(smsBrk[[12]])

# 30-day Rolling mean of SMS
myrollfun=function(x){
  if(length(which(!is.na(x)))>10){
    out=movingFun(x,30, mean,circular=FALSE,type='to', na.rm=TRUE)
  } else{
    out=rep(NA, length(x))
  }
  return(out)
}
# Generate 30-day rolling average of SMS
numCores=detectCores()-1; beginCluster(n=numCores)
sms30Brk=clusterR(smsBrk, calc, args=list(fun=myrollfun),progress='text');endCluster()

# tsOut=as.numeric(extract(sms30Brk, cbind(-115,40)))
# plot(tsOut, type="l")
# 
# plot(sms30Brk[[1000]], col=viridisLite::turbo(20), zlim=c(0,1))
# points(cbind(-115,40), pch=19, col="black")

###############################################
# PART 3: Calculate RRD and FDSI
###############################################
smParStk=stack()
smParStk=addLayer(interp_brk,thetaTDdoy,thetaWTdoy)

#plot(thetaTDdoy[[200]], col=turbo(20))
#points(SpatialPoints(cbind(-115,37)), pch=19, col="white")
#x=as.numeric(extract(interp_brk[[1:366]], SpatialPoints(cbind(-115,37))))

writeRaster(smParStk,paste(path,"smParStk_",regionRun,".nc", sep=""), overwrite = TRUE)
smParStk=brick(paste(path,"smParStk_",regionRun,".nc", sep=""))

pacman::p_load(smapr, sp, raster, terra, purrr, dplyr, lubridate, zoo,twitteR,ncdf4,magick,
               ClusterR, parallel, nls2, rgdal, minpack.lm, officer, gmailr, exactextractr, sf) 

rollingm2=function(x){
  x=as.numeric(array(as.vector(x)))
  x[x=="NaN"]=NA
  
  timestep=1
  smapdim=length(x)-(366*2)
  thetaTDdim=thetaWTdim=366
  
  retrieval_date=seq(as.Date("2015-03-31"),as.Date("2021-09-30"), by=1)
  
  cellSM=x[1:smapdim]
  cellThetaTDdoy=x[(smapdim+1):(smapdim+thetaTDdim)]
  cellThetaWTdoy=x[(smapdim+thetaTDdim+1):(smapdim+thetaTDdim+thetaWTdim)]
  
  if ((length(na.omit(cellSM)>100)) & (length(na.omit(cellThetaTDdoy)>90)) & (length(na.omit(cellThetaWTdoy)>90))){
    cellThetaTDruntime=cellThetaTDdoy[lubridate::yday(retrieval_date)]
    cellThetaWTruntime=cellThetaWTdoy[lubridate::yday(retrieval_date)]
    
    ts_data=data.frame(SM=array(cellSM[-c(length(cellSM))]),
                       Loss=base::diff(cellSM, 1, na.pad = TRUE),
                       ThetaTD=cellThetaTDruntime[-c(length(cellSM))],
                       ThetaWT=cellThetaWTruntime[-c(length(cellSM))])
    ts_data[ts_data=="NaN"]=NA
    
    # Filter sensor errors
    ind=abs(ts_data$Loss)>(0.01*(max(cellSM,na.rm = T)-min(cellSM,na.rm = T)))
    ind[is.na(ind)]=FALSE
    ts_data[which(ind==FALSE),]=NA
    
    # Filter wetting data points
    ind=ts_data$Loss<0
    ind[is.na(ind)]=FALSE
    ts_data[!ind,]=NA
    ts_data$Loss=(-1/timestep)*ts_data$Loss
    
    ts_data$SMfilt=ts_data$SM
    ts_data$SMfilt[(ts_data$SM<=ts_data$ThetaTD) | (ts_data$SM>=ts_data$ThetaWT)]=NA
    
    #############
    rollNLSfit=function(df){
      if(nrow(data.frame(na.omit(df)))>=7){
        
        nlsFit=try(minpack.lm::nlsLM(Loss ~(slope*(SMfilt-ewp)), 
                                     df,
                                     start = list(slope=mean(df$Loss, na.rm=T)/mean(df$SMfilt, na.rm=T),ewp=0.1),
                                     lower=c(0.001,0.01), upper=c(20,0.65),
                                     control = list(maxiter = 1000, minFactor=1/2000, warnOnly=T)),silent = TRUE)
        
        if ("try-error" %in% class(nlsFit)){
          return(NA)
        } else {
          slope=summary(nlsFit)$parameters[1,1]
          pval=summary(nlsFit)$parameters[1,4]
          #return(slope)
          ifelse(pval<0.05,return(slope),return(NA))
        }
      } else{
        return(NA)
      }
    }
    splitRoll=function(df){return(as.data.frame(df))}
    rollapply.data.frame <- function(data, ..., fill = NULL, FUN, simplify = function(x) do.call(rbind, x)) {
      fill0 <- if (!is.null(fill)) NA
      result <- lapply(
        as.data.frame(t(rollapply(1:nrow(data), ..., fill = fill0, FUN = c))), 
        function(ix) {if (all(is.na(ix))) fill else FUN(data[ix, ])}
      )
      simplify(result)
    }
    
    # Create rolling list of 30-day data
    splitRollDF=zoo::rollapply(ts_data, 30, FUN = splitRoll,
                               by = 1,fill=NA,
                               align = c("right"), simplify = identity)
    # Map NLS fit on rolling dataset
    nlFitDf=purrr::map_dbl(splitRollDF,rollNLSfit)
    
    return(as.numeric(array(nlFitDf)))
  } else{
    return(rep(NA,smapdim-1))
  }
}


x=as.numeric(extract(smParStk, SpatialPoints(cbind(-115,37))))
plot(rollingm2(x), type="l")

numCores=detectCores()-1; beginCluster(n=numCores)
RDbrkGlobal=clusterR(smParStk,
                     calc, args=list(fun=rollingm2),
                     progress='text');endCluster()

writeRaster(RDbrkGlobal,paste("C:/Users/Brown/Documents/FLASH/9kmV2/RDbrkGlobal_",regionRun,".nc", sep=""), overwrite = TRUE)
RDbrkGlobal = brick(paste("C:/Users/Brown/Documents/FLASH/9kmV2/RDbrkGlobal_",regionRun,".nc", sep=""))
#x=as.numeric(extract(RDbrkGlobal, SpatialPoints(cbind(-115,37))))
#plot((x), type="l")

# CHANGE #2
RRD=1/(1+(m2runtime/RDbrkGlobal)^6)
values(RRD)[which((is.na(values(RRD))) & !(is.na(values(m2runtime))))]=0.5
values(RRD)[which(values(RRD)<0.5)]=0.5

#~ Calculate FDSI
FDSI30=sqrt(sms30Brk*RRD)
writeRaster(FDSI30,paste("C:/Users/Brown/Documents/FLASH/9kmV2/FDSI30_",regionRun,".nc", sep=""), overwrite = TRUE)

fdsi36 = brick("C:/Users/Brown/Documents/FLASH/FDSI.nc")


##################################
# Plotting FDSI maps
##################################
# #~ FDSI 1.0
 #fdsi36=brick("C:/Users/vinit.sehgal/Documents/GFDASbackup/FDSI.nc")
 #FDSI30=brick("G:/My Drive/TAMU/Teaching/Mentoring/AggieResearchProgram/FLASH_shared/resources/HPRC/FLASH2.0data/FDSI30_6.nc")

 # Select FDSI 1.0 raster for the plot date
 plotDate="2020-10-07" #

 # Select 36-KM raster for the specific date
 fdsi_select=fdsi36[[which(as.Date(substr(names(fdsi36),2,9), format="%Y%m%d")==plotDate)]]
 fdsi_36Ras=raster(fdsi_select)
 values(fdsi_36Ras)=values(fdsi_select)
# # Mask global raster to regional extent
 fdsi_36Ras=raster::crop(fdsi_36Ras,extent(FDSI30))

 nf=layout(mat = matrix(c(1, 1),  nrow = 1, ncol = 2))

 color_pal <- colorRampPalette(c("blue", "lightblue","white", "orange","red"))
 brk=c(0,0.2,0.4,0.5,0.6,0.71,1)

 par(mfrow=c(2,1))
 FDSI9KM=raster(FDSI30[[which(interp_dates==plotDate)]])
 values(FDSI9KM)=values(FDSI30[[which(interp_dates==plotDate)]])

 plot(rast(FDSI9KM),
      range=c(0,1),
      asp=NA,
      col=(color_pal(length(brk)-1)), # Custom colormap
      breaks=brk,
      colNA="white",
      #legend="right",
      main=plotDate)

 plot(rast(fdsi_36Ras),
      range=c(0,1),
      asp=NA,
      col=(color_pal(length(brk)-1)), # Custom colormap
      breaks=brk,
      colNA="white",
      main=plotDate
      #legend="rightoutside"
     )

 fdsi36data = as.data.frame(as(fdsi_36Ras, "SpatialPixelsDataFrame"))
 fdsi09data = as.data.frame(as(FDSI30[[which(interp_dates==plotDate)]], "SpatialPixelsDataFrame"))

# ####################
# # Plot using GGPLOT
# ####################
library(ggplot2)
   
 globalPlot=ggplot() +                           # Initialize the plot
     geom_tile(data=fdsi09data,                    # Add soil moisture data
               aes_string(x="x", y="y", fill= cut(fdsi09data$layer, breaks=brk)),alpha=1)+       # Provide X, Y and Z data
     scale_fill_manual(drop=TRUE, values=(color_pal(length(brk)-1)),
                       na.value="transparent", name="FDSI",
                       labels=levels(cut(fdsi09data$layer, breaks=brk)))+           # Use user-defined colormap                         # Set limits of colorbar
     coord_cartesian(xlim = c(-130.031, -104.9998),
                     ylim = c(33.7889, 50.0269))  +                                            #xlab("Longitude")+                            # X-axis label
     ylab("Latitude")+  xlab("Longitude")+                            # Y-axis label
     theme_bw() +
     borders("world",                              # Add global landmass boundaries
             colour="gray43",  size = 0.1,                    # Fill light-gray color to the landmass
             fill="transparent")  +                # Transparent background
     borders("state",                # Add US state borders
             colour = "gray43",      # Use light-gray color
             fill = "transparent")+
     labs(title = 'Western US Drought Outlook - FLASH 1.0 - 36km',
          subtitle = "March 30, 2021",
          caption = 'FDSI >0.71 indicates emerging/ sustained flash droughts')

   print(globalPlot)

x36=as.numeric(extract(fdsi36, SpatialPoints(cbind(-115,37))))
plot(x36, type="l")
x09=as.numeric(extract(FDSI30, SpatialPoints(cbind(-115,37))))
plot(x09, type="l")
plot(x36[1:length(x09)], x09)


