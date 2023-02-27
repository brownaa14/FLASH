library(parallel)
library(ClusterR)
library(raster)
library(terra)
library(lubridate)

###############################################
# PART 2: Calculate SMS
###############################################
thetaWTdoy=brick("G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/paraDoy/thetaWTdoy_R6.nc")
thetaTDdoy=brick("G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/paraDoy/thetaTDdoy_R6.nc")
thetaIPdoy=brick("G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/paraDoy/thetaIPts_R6.nc")
m2doy=brick("G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/paraDoy/m2doy_R6.nc")

interp_brk = brick("G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/resources/HPRC/SPL3E_1-dayInterp/Region_SPL3SMP_E_1dayInterp_R6.nc")
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
nRuntime=16*sqrt(m2runtime)
smsBrk=1/(1+(ratioBrk^nRuntime))

# 30-day Rolling mean of SMS
myrollfun=function(x){
  if(length(which(!is.na(x)))>10){
    out=movingFun(x,30, mean,circular=FALSE,type='to', na.rm=TRUE)
  } else{
    out=rep(NA, length(x))
  }
  return(out)
}
beginCluster()
sms30Brk=clusterR(smsBrk, calc, args=list(fun=myrollfun),progress='text');endCluster()


plot(smsBrk[[300]])

###############################################
# PART 3: Calculate RRD and FDSI
###############################################

dateOperate = as.Date("2021-2-15", "%y-%d-%m")

smParStk=stack()
smParStk=addLayer(interp_brk,thetaTDdoy,thetaWTdoy)
library(zoo)
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

#x=as.numeric(extract(smParStk, c(-115,37)))
#x=as.numeric(extract(smParStk, SpatialPoints(cbind(-115,37))))

rollingm2(x)

#cropped = crop(smParStk, extent(c(-130.031, -104.9998, 33.7889, 50.0269)))
numCores=detectCores()-1; beginCluster(n=numCores)
#plot(cropped[[4]])
RDbrkGlobal=clusterR(smParStk,
                     calc, args=list(fun=rollingm2),
                     #export=c("dateOperate"),
                     progress='text');endCluster()

RRD=1/(1+(m2runtime[[dim(m2runtime)[3]]]/RDbrkGlobal)^6)
values(RRD)[which((is.na(values(RRD))) & !(is.na(values(m2runtime[[dim(m2runtime)[3]]]))))]=0.5
values(RRD)[which(values(RRD)<0.5)]=0.5

#~ Calculate FDSI
FDSI30=sqrt(sms30Brk[[dim(sms30Brk)[3]]]*RRD)

plot(FDSI30[[4]])


writeRaster(FDSI30,"G:/.shortcut-targets-by-id/1_munBiOkEi1EyelAxqX1FaQ56WRqt47Z/FLASH_shared/Aidan Brown/Data/FLASH2_9km/paraDoy/FDSI30_R6.nc", overwrite = TRUE)
