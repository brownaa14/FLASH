library(raster)
library(ClusterR)
library(ncdf4) 
library(nls2)
library(minpack.lm)
library(modelr)
library(purrr)
library(broom)
library(Hmisc)
library(Rmisc)
library(data.table)
library(rgdal) 
library(parallel)
library(hydroGOF)

args <- commandArgs(trailingOnly =TRUE)
regionRun=as.numeric(args[1])
seasonRun=as.numeric(args[2])
smapNAME="SPL3SMP_E"
print(paste("regionRun=",regionRun," | smapNAME=",smapNAME," | userSeason=",seasonRun, sep = ""))

###########################
# NLS fit for SM drydowns
###########################
#~~~ SMAP raster brick
ncFileName = paste("/scratch/user/aabrown24/9km_drydownFit/regionalData_SMAPL3_E/Region_",smapNAME,"_2dayInterp_R",regionRun,".nc", sep="")
#ncFileName = paste("/scratch/user/saraant/2022/drydownFit_9km/regionalData_SMAPL3_E/Region_",smapNAME,"_2dayInterp_R",regionRun,".nc", sep="")

smap_ras = brick(ncFileName)
sm_date=seq(as.Date("2015-03-31"),as.Date("2021-09-30"), by=2)

#~~~ Model optimization function 
optimize_nls_fit=function(x){
  
  nresamp=300         # Number of resamples
  minFitSamples=50    # Minimum number of drydown samples for fitting
  testSplit=0.25      # Fraction of testing dataset
  plot_fit=FALSE      # Plot fitted curves?
  
  monthTriplet=function (x){
    mseq=c(12,seq(1:12),1)
    ifelse((x !=12),return(c(mseq[which(mseq==x)[1]-1], mseq[which(mseq==x)[1]], mseq[which(mseq==x)[1]+1])),return(c(11,12,1)))
  }
  runmonths=monthTriplet(seasonRun)
  
  #~~ Sample time series
  sm_ts=as.matrix(as.numeric(as.vector(x)))
  if(nrow(data.frame(na.omit(sm_ts)))>minFitSamples){
    ts_data=data.frame(SM=array(sm_ts[-c(length(sm_ts))]),
                       Loss=base::diff(sm_ts, 1, na.pad = TRUE), 
                       months=data.table::month(sm_date[-c(length(sm_ts))]))
    ts_data=subset(ts_data,abs(ts_data$Loss)>(0.01*(max(sm_ts,na.rm = T)-min(sm_ts,na.rm = T)))) 
    ts_data=subset(ts_data, ts_data$Loss<0)
    ts_data$Loss=(-0.5)*ts_data$Loss
    # Seasonal filtering
    ts_data=subset(ts_data, ts_data$months %in% c(runmonths))
    ts_data=subset (ts_data, select = -months)
  } else {
    ts_data=NULL
  }
  
  if (nrow(data.frame(ts_data))>100){  
    #~ Resample dataset for model fitting 
    set.seed(100); cvdata <- modelr::crossv_mc(ts_data,nresamp ,test = testSplit)
    
    #~~~ Functions
    gof=function(x){
      samp_df=as.data.frame(x)
      return(data.frame(CC=cor(samp_df$.fitted,samp_df$Loss), 
                        d=hydroGOF::d(samp_df$.fitted,samp_df$Loss), 
                        mse=hydroGOF::mse(samp_df$.fitted,samp_df$Loss)))
    }
    map_nlsFCF=function(x){
      df_in=as.data.frame(x)
      # df_in=as.data.frame(cvdata$train[[12]]); plot(df_in)
      
      limValSM=array(quantile(df_in$SM, probs = c(0.05,0.15,0.25,0.5,0.75,0.95,0.999)))
      limValLoss=array(quantile(df_in$Loss, probs = c(0.001,0.3)))
      
      RNGthetaTD=c(0.02,limValSM[4])
      RNGthetaWT=c(limValSM[3],limValSM[7])
      RNGthetaGW=c(limValSM[4],limValSM[7])
      RNGm1=RNGm2=c(0.001,1)
      
      set.seed(100);nl_mod=try(minpack.lm::nlsLM(Loss ~ ((SM <= thetaWT) * (m2*(SM-thetaTD))+
                                                           ((SM > thetaWT) & (SM < thetaGW)) * (m2*(thetaWT-thetaTD))+
                                                           (SM >= thetaGW) * (m2*(thetaWT-thetaTD)+ m1*(SM-thetaGW))),
                                                 data=df_in,
                                                 start = list(thetaTD=limValSM[1],
                                                              thetaWT=mean(RNGthetaWT),
                                                              thetaGW=mean(RNGthetaGW),
                                                              m1=mean(df_in$Loss)/mean(df_in$SM),
                                                              m2=mean(df_in$Loss)/mean(df_in$SM)),
                                                 lower=c(RNGthetaTD[1],RNGthetaWT[1],RNGthetaGW[1],RNGm1[1],RNGm2[1]),
                                                 upper=c(RNGthetaTD[2],RNGthetaWT[2],RNGthetaGW[2],RNGm1[2],RNGm2[2]),
                                                 control = nls.lm.control(maxiter = 500)),
                               silent = TRUE)
      
      # nl_mod=try(nls2::nls2(Loss ~ ((SM <= thetaWT) * (m2*(SM-thetaTD))+
      #                                 ((SM > thetaWT) & (SM < thetaGW)) * (m2*(thetaWT-thetaTD))+
      #                                 (SM >= thetaGW) * (m2*(thetaWT-thetaTD)+ m1*(SM-thetaGW))),
      #                       data=df_in,
      #                       algorithm="port",
      #                       start = list(thetaTD=limValSM[1],
      #                                    thetaWT=mean(RNGthetaWT),
      #                                    thetaGW=mean(RNGthetaGW),
      #                                    m1=mean(df_in$Loss)/mean(df_in$SM),
      #                                    m2=mean(df_in$Loss)/mean(df_in$SM)),
      #                       lower=c(RNGthetaTD[1],RNGthetaWT[1],RNGthetaGW[1],RNGm1[1],RNGm2[1]),
      #                       upper=c(RNGthetaTD[2],RNGthetaWT[2],RNGthetaGW[2],RNGm1[2],RNGm2[2]),
      #                       control = list(maxiter = 500,minFactor=1/2000, warnOnly=T)),silent = TRUE)
      
      if ("try-error" %in% class(nl_mod)){
        return(NA)
      } else {
        return(nl_mod)
      }
    }
    map_nlsCFC=function(x){
      df_in=as.data.frame(x)
      #df_in=as.data.frame(cvdata$train[[11]]); plot(df_in)
      
      # Parameter search space
      limValSM=array(quantile(df_in$SM, probs = c(0.1,0.15,0.25,0.5,0.75,0.95,0.999)))
      limValLoss=array(quantile(df_in$Loss, probs = c(0.001,0.3)))
      RNGthetaTD=c(0.02,limValSM[4])
      RNGthetaWT=c(limValSM[2],limValSM[7])
      RNGld=array(quantile(df_in$Loss, probs = c(0.001,0.3)))
      RNGm2=c(0.001,1)
      
      set.seed(100);nl_mod=try(minpack.lm::nlsLM(Loss ~(SM <= thetaTD)*ld+
                                                   ((SM > thetaTD) & (SM < thetaWT)) * (m2*(SM-thetaTD)+ld)+
                                                   (SM >= thetaWT) *(m2*(thetaWT-thetaTD)+ld),
                                                 data=df_in,
                                                 start = list(thetaTD=limValSM[1],
                                                              thetaWT=mean(RNGthetaWT),
                                                              m2=mean(df_in$Loss)/mean(df_in$SM),
                                                              ld=mean(limValLoss)),
                                                 lower=c(RNGthetaTD[1],RNGthetaWT[1],RNGm2[1],RNGld[1]),
                                                 upper=c(RNGthetaTD[2], RNGthetaWT[2],RNGm2[2],RNGld[2]),
                                                 control = nls.lm.control(maxiter = 500)),
                               silent = TRUE)
      if ("try-error" %in% class(nl_mod)){
        return(NA)
      } else {
        return(nl_mod)
      }
    }
    map_nlsFC=function(x){
      df_in=as.data.frame(x)
      #df_in=as.data.frame(cvdata$train[[11]]); plot(df_in)
      
      # Parameter search space
      limValSM=array(quantile(df_in$SM, probs = c(0.1,0.15,0.25,0.5,0.75,0.95,0.999)))
      RNGthetaTD=c(0.02,limValSM[6])
      RNGm2=c(0.001,1)
      RNGld=array(quantile(df_in$Loss, probs = c(0.001,0.25)))
      
      set.seed(100);nl_mod=try(minpack.lm::nlsLM(Loss ~((SM <= thetaTD)*ld+
                                                          (SM > thetaTD) * (ld+(m2*(SM-thetaTD)))),
                                                 data=df_in,
                                                 start = list(thetaTD=mean(RNGthetaTD),
                                                              m2=mean(df_in$Loss)/mean(df_in$SM),
                                                              ld=array(quantile(df_in$Loss, probs = c(0.02)))),
                                                 lower=c(RNGthetaTD[1],RNGm2[1],RNGld[1]),
                                                 upper=c(RNGthetaTD[2],RNGm2[2],RNGld[2]),
                                                 control = nls.lm.control(maxiter = 500)),
                               silent = TRUE)
      if ("try-error" %in% class(nl_mod)){
        return(NA)
      } else {
        return(nl_mod)
      }
    }
    map_nlsCF=function(x){
      df_in=as.data.frame(x)
      #df_in=as.data.frame(cvdata$train[[11]]); plot(df_in)
      
      # Parameter search space
      limValSM=array(quantile(df_in$SM, probs = c(0.1,0.15,0.25,0.5,0.75,0.95,0.999)))
      RNGthetaTD=c(0.02,limValSM[6])
      RNGthetaWT=c(limValSM[3],limValSM[7])
      RNGm2=c(0.001,1)
      
      set.seed(100);nl_mod=try(minpack.lm::nlsLM(Loss ~(SM < thetaWT) * (m2*(SM-thetaTD))+
                                                   (SM >= thetaWT) * (m2*(thetaWT-thetaTD)),
                                                 data=df_in,
                                                 start = list(thetaTD=mean(RNGthetaTD),
                                                              thetaWT=mean(RNGthetaWT),
                                                              m2=mean(df_in$Loss)/mean(df_in$SM)),
                                                 lower=c(RNGthetaTD[1],RNGthetaWT[1],RNGm2[1]),
                                                 upper=c(RNGthetaTD[2],RNGthetaWT[2],RNGm2[2]),
                                                 control = nls.lm.control(maxiter = 500)),
                               silent = TRUE)
      if ("try-error" %in% class(nl_mod)){
        return(NA)
      } else {
        return(nl_mod)
      }
    }
    map_nlsF=function(x){
      df_in=as.data.frame(x)
      #df_in=as.data.frame(cvdata$train[[11]]); plot(df_in)
      
      # Parameter search space
      limValSM=array(quantile(df_in$SM, probs = c(0.1,0.15,0.25,0.5,0.75,0.95,0.999)))
      RNGthetaTD=c(0.02,limValSM[4])
      RNGm2=c(0.001,1)
      
      set.seed(100);nl_mod=try(minpack.lm::nlsLM(Loss ~(m2*(SM-thetaTD)),
                                                 data=df_in,
                                                 start = list(thetaTD=mean(RNGthetaTD),
                                                              m2=mean(df_in$Loss)/mean(df_in$SM)),
                                                 lower=c(RNGthetaTD[1],RNGm2[1]),
                                                 upper=c(RNGthetaTD[2],RNGm2[2]),
                                                 control = nls.lm.control(maxiter = 500)),
                               silent = TRUE)
      if ("try-error" %in% class(nl_mod)){
        return(NA)
      } else {
        return(nl_mod)
      }
    }
    map_nlsC=function(x){
      df_in=as.data.frame(x)
      #df_in=as.data.frame(cvdata$train[[11]]); plot(df_in)
      
      #Parameter search space
      RNGld=array(quantile(df_in$Loss, probs = c(0.1,0.9)))
      
      set.seed(100);nl_mod=try(minpack.lm::nlsLM(Loss ~ld,
                                                 data=df_in,
                                                 start = list(ld=mean(RNGld)),
                                                 lower=c(RNGld[1]),
                                                 upper=c(RNGld[2]),
                                                 control = nls.lm.control(maxiter = 500)),
                               silent = TRUE)
      if ("try-error" %in% class(nl_mod)){
        return(NA)
      } else {
        return(nl_mod)
      }
    }
    
    #~ Fit NLM models
    collate_gof=list()
    collate_param=list()
    
    #~ 6:FCF 
    nls_modelsFCF <- map(cvdata$train, ~ map_nlsFCF(.))
    #~~ Apply converged models
    converged_models_FCF=nls_modelsFCF[!is.na(nls_modelsFCF)] # length(converged_models)
    test_data=cvdata$test[!is.na(nls_modelsFCF)]
    pred_ts_FCF <- map2(converged_models_FCF, test_data, ~ augment(.x, newdata = as.data.frame(.y)))
    
    #~~ Params and GOF
    collate_gof[[6]]=as.data.frame(map_df(pred_ts_FCF,gof))
    collate_param[[6]]=as.data.frame(map_dfr(converged_models_FCF,coefficients))
    
    #~~ Filter physically incorrect numbers 
    filter_in=(collate_param[[6]]$thetaGW >(1.01*collate_param[[6]]$thetaWT))
    collate_param[[6]]=collate_param[[6]][filter_in,]
    collate_gof[[6]]=collate_gof[[6]][filter_in,]
    
    filter_in=(collate_param[[6]]$thetaWT >(1.01*collate_param[[6]]$thetaTD))
    collate_param[[6]]=collate_param[[6]][filter_in,]
    collate_gof[[6]]=collate_gof[[6]][filter_in,]
    rm(converged_models_FCF,nls_modelsFCF)
    
    #~ 5:CFC
    nls_modelsCFC <- map(cvdata$train, ~ map_nlsCFC(.))
    #~~ Apply converged models
    converged_models_CFC=nls_modelsCFC[!is.na(nls_modelsCFC)] # length(converged_models)
    test_data=cvdata$test[!is.na(nls_modelsCFC)]
    pred_ts_CFC <- map2(converged_models_CFC, test_data, ~ augment(.x, newdata = as.data.frame(.y)))
    
    #~~ Params and GOF
    collate_gof[[5]]=as.data.frame(map_df(pred_ts_CFC,gof))
    collate_param[[5]]=as.data.frame(map_dfr(converged_models_CFC,coefficients))
    
    #~~ Filter physically incorrect numbers 
    filter_in=(collate_param[[5]]$thetaWT >(1.01*collate_param[[5]]$thetaTD))
    collate_param[[5]]=collate_param[[5]][filter_in,]
    collate_gof[[5]]=collate_gof[[5]][filter_in,]
    rm(converged_models_CFC,nls_modelsCFC)
    
    #~ 4:FC
    nls_modelsFC <- map(cvdata$train, ~ map_nlsFC(.))
    #~~ Apply converged models
    converged_models_FC=nls_modelsFC[!is.na(nls_modelsFC)] # length(converged_models)
    test_data=cvdata$test[!is.na(nls_modelsFC)]
    pred_ts_FC <- map2(converged_models_FC, test_data, ~ augment(.x, newdata = as.data.frame(.y)))
    
    #~~ Params and GOF
    collate_gof[[4]]=as.data.frame(map_df(pred_ts_FC,gof))
    collate_param[[4]]=as.data.frame(map_dfr(converged_models_FC,coefficients))
    rm(converged_models_FC,nls_modelsFC)
    
    #~ 3:CF
    nls_modelsCF <- map(cvdata$train, ~ map_nlsCF(.))
    #~~ Apply converged models
    converged_models_CF=nls_modelsCF[!is.na(nls_modelsCF)] # length(converged_models)
    test_data=cvdata$test[!is.na(nls_modelsCF)]
    pred_ts_CF <- map2(converged_models_CF, test_data, ~ augment(.x, newdata = as.data.frame(.y)))
    
    #~~ Params and GOF
    collate_gof[[3]]=as.data.frame(map_df(pred_ts_CF,gof))
    collate_param[[3]]=as.data.frame(map_dfr(converged_models_CF,coefficients))
    
    #~~ Filter physically incorrect numbers 
    filter_in=(collate_param[[3]]$thetaWT >(1.01*collate_param[[3]]$thetaTD))
    collate_param[[3]]=collate_param[[3]][filter_in,]
    collate_gof[[3]]=collate_gof[[3]][filter_in,]
    rm(converged_models_CF,nls_modelsCF)
    
    #~ 2:F
    nls_modelsF <- map(cvdata$train, ~ map_nlsF(.))
    #~~ Apply converged models
    converged_models_F=nls_modelsF[!is.na(nls_modelsF)] # length(converged_models)
    test_data=cvdata$test[!is.na(nls_modelsF)]
    pred_ts_F <- map2(converged_models_F, test_data, ~ augment(.x, newdata = as.data.frame(.y)))
    
    #~~ Params and GOF
    collate_gof[[2]]=as.data.frame(map_df(pred_ts_F,gof))
    collate_param[[2]]=as.data.frame(map_dfr(converged_models_F,coefficients))
    rm(converged_models_F,nls_modelsF)
    
    #~ 1:C
    nls_modelsC <- map(cvdata$train, ~ map_nlsC(.))
    #~~ Apply converged models
    converged_models_C=nls_modelsC[!is.na(nls_modelsC)] # length(converged_models)
    test_data=cvdata$test[!is.na(nls_modelsC)]
    pred_ts_C <- map2(converged_models_C, test_data, ~ augment(.x, newdata = as.data.frame(.y)))
    
    #~~ Params and GOF
    collate_gof[[1]]=as.data.frame(map_df(pred_ts_C,gof))
    collate_param[[1]]=as.data.frame(map_dfr(converged_models_C,coefficients))
    
    names(collate_gof)=rev(c("FCF","CFC","FC","CF","F","C"))
    names(collate_param)=rev(c("FCF","CFC","FC","CF","F","C"))
    rm(converged_models_C,nls_modelsC)
    
    ##########################################
    #~~ Performance comparison of fitted models
    #########################################
    if(sum(unlist(map(collate_gof, function(x){!purrr::is_empty(x)})), na.rm=TRUE)>=2){
      
      CIfun=function(x){return(CI(x$mse))}
      mseCI=as.data.frame(map_df(collate_gof, CIfun))
      bestMod=which.min(mseCI$mean)
      selectMod=min(which((mseCI$mean[1:bestMod] %between% c(mseCI$lower[bestMod],mseCI$upper[bestMod])) == TRUE))
      
      selectMod_ParMed=as.data.frame(t(apply(collate_param[[selectMod]],2, mean)))
      selectMod_GoF=as.data.frame(t(apply(collate_gof[[selectMod]],2, mean)))
      selectMod_ParSD=data.frame(t(apply(collate_param[[selectMod]],2, sd)))
      names(selectMod_ParSD)=paste("sd_",names(selectMod_ParSD),sep="")
      nSamp=nrow(collate_param[[selectMod]])
      
      #~~~~~~~~~ Import and collate parameters from all parts 
      param_names=c("canonical_form",
                    "thetaTD","thetaWT","thetaGW",
                    "ld","lw","m1","m2",
                    "sd_thetaTD","sd_thetaWT","sd_thetaGW",
                    "sd_ld","sd_lw","sd_m1","sd_m2",
                    "mse","CC","d",
                    "replacement","nSamp")
      
      collate_param_df=setNames(data.frame(matrix(ncol =20, nrow = 0)),param_names)
      mergeParFit=dplyr::bind_rows(collate_param_df,
                                   data.frame(canonical_form=selectMod,
                                              selectMod_ParMed,
                                              selectMod_ParSD, 
                                              selectMod_GoF, 
                                              replacement=ifelse(bestMod==selectMod, 0,1),
                                              nSamp=nSamp))
      
      #~ Plot fitted curved on data
      if (plot_fit ==TRUE){
        plot_fit=function(median_par,ts_data){
          
          if(ncol(as.data.frame(median_par))>0){
            
            line_plot_df=data.frame(SM=seq(0,max(ts_data$SM), by= max(ts_data$SM)/length(ts_data$SM)))
            
            if(ncol(median_par)==5){
              form="6. FCF"
              samp_out= (line_plot_df$SM <= median_par$thetaWT) * (median_par$m2*(line_plot_df$SM-median_par$thetaTD))+
                ((line_plot_df$SM > median_par$thetaWT) & (line_plot_df$SM < median_par$thetaGW)) * (median_par$m2*(median_par$thetaWT-median_par$thetaTD))+
                (line_plot_df$SM >= median_par$thetaGW) * (median_par$m2*(median_par$thetaWT-median_par$thetaTD)+ median_par$m1*(line_plot_df$SM-median_par$thetaGW))
              
            } else if (ncol(median_par)==4){
              form="5. CFC"
              samp_out=(line_plot_df$SM <= median_par$thetaTD)*median_par$ld+
                ((line_plot_df$SM > median_par$thetaTD) & (line_plot_df$SM < median_par$thetaWT)) * (median_par$m2*(line_plot_df$SM-median_par$thetaTD)+median_par$ld)+
                (line_plot_df$SM >= median_par$thetaWT) *(median_par$m2*(median_par$thetaWT-median_par$thetaTD)+median_par$ld)
              
            } else if (ncol(median_par)==3 & (all(c("thetaTD","ld","m2") %in% names(median_par)))){
              form="4. FC"
              samp_out=(line_plot_df$SM <= median_par$thetaTD)*median_par$ld+
                (line_plot_df$SM > median_par$thetaTD) * (median_par$ld+(median_par$m2*(line_plot_df$SM-median_par$thetaTD)))
              
            } else if (ncol(median_par)==3 & (all(c("thetaTD","thetaWT","m2") %in% names(median_par)))){
              form="3. CF"
              samp_out=(line_plot_df$SM < median_par$thetaWT) * (median_par$m2*(line_plot_df$SM-median_par$thetaTD))+
                (line_plot_df$SM >= median_par$thetaWT) * (median_par$m2*(median_par$thetaWT-median_par$thetaTD))
            } else if (ncol(median_par)==2){
              form="2. F"
              samp_out=(median_par$m2*(line_plot_df$SM-median_par$thetaTD))
              
            } else if (ncol(median_par)==1){
              form="1. C"
              samp_out=rep(median_par$ld, nrow(line_plot_df))
            }
            
            # Plot the fitted curve on observed data
            plot(ts_data$SM,ts_data$Loss, ylim=range(ts_data$Loss), ylab="", 
                 xlim=range(ts_data$SM), xlab="", pch=20, col="skyblue", main=form)
            par(new=TRUE)
            plot(line_plot_df$SM,samp_out, ylim=range(ts_data$Loss),ylab="Loss (v/v/t)",
                 xlim=range(ts_data$SM), xlab="SM (v/v)",col="maroon", pch=16, type="l",
                 lwd=2.4)
          } else{
            plot(NULL, xlim=c(0,1), ylim=c(0,1),ylab="Loss (v/v/t)",xlab="SM (v/v)",)
          }
        }
        par(mfrow=c(2,3))
        plot_fit(median_par=as.data.frame(t(apply(collate_param[[6]],2, median))),ts_data); grid()
        plot_fit(median_par=as.data.frame(t(apply(collate_param[[5]],2, median))),ts_data); grid()
        plot_fit(median_par=as.data.frame(t(apply(collate_param[[4]],2, median))),ts_data); grid()
        plot_fit(median_par=as.data.frame(t(apply(collate_param[[3]],2, median))),ts_data); grid()
        plot_fit(median_par=as.data.frame(t(apply(collate_param[[2]],2, median))),ts_data); grid()
        plot_fit(median_par=as.data.frame(t(apply(collate_param[[1]],2, median))),ts_data); grid()
        par(mfrow=c(1,1))
      }
      #return(mergeParFit)
      if(length(as.numeric(array(mergeParFit)))==20){
        return(as.numeric(array(mergeParFit)))
      } else{ 
        return(rep(NA,20))
      }
    } else {
      return(rep(NA,20))
    }
  } else {
    return(rep(NA,20))
  } 
}

#~~~ Apply function in parallel
numCores=detectCores()  
print(numCores)
beginCluster(n=numCores-1)    
start_time <- Sys.time() 
NLSparBrk=clusterR(smap_ras, 
                   calc, args=list(fun=optimize_nls_fit),
                   export=c("sm_date","seasonRun"),
                   progress='text');endCluster()
end_time <- Sys.time()
end_time - start_time  

#~~~ Assign names to output RasterBrick
param_names=c("canonical_form","thetaTD","thetaWT","thetaGW","ld","lw","m1","m2",
              "sd_thetaTD","sd_thetaWT","sd_thetaGW","sd_ld","sd_lw","sd_m1","sd_m2",
              "mse","CC","d","replacement","nSamp")
names(NLSparBrk)=param_names

#~~~ Export NetCDF to disk
exportNCDF=function(rb, ncName){
  r=rb[[1]] # raster taken from a first layer of a stack
  rlon=rlat=r #copy r to rlon and rlat rasters [1]][1]which will contain the longitude and latitude
  xy<-xyFromCell(r,1:length(r)) #matrix of logitudes (x) and latitudes(y)
  lon=unique(xy[,1])
  lat=unique(xy[,2])
  
  # first we write the netcdf file to disk
  writeRaster(rb,ncName,overwrite=TRUE,format="CDF",
              varname=deparse(substitute(rb)),varunit="[-]",
              longname="[-]",xname="Longitude",yname="Latitude",
              zname='Time',zunit='[-]',progress="text")
  
  # #### Adding lat, long to netCDF
  nc = nc_open(ncName, write = TRUE)
  ncdf4::ncvar_put(nc, "Longitude", lon)
  ncdf4::ncvar_put(nc, "Latitude", lat)
  nc_close(nc)
}
exportNCDF(NLSparBrk,ncName=paste("Regional_SPL3_E_v7_NLSparBrk_",seasonRun,".nc",sep=""))
