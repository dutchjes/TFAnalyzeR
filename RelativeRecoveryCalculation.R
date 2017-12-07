
## From Karin:

#' function selectISTDs() will select all ISTDs which elute within the defined RT window
#'
#' @param TFdata TraceFinder quantification file, saved as a data.frame
#' @param target name of target compound, needs to be same as name in TFdata
#' @param RT_window ISTDs are selected which elute +/- RT_window (min)
#'
selectISTDs <- function(TFdata, target, RT_window = 1.5){
  ISTDs  <- unique(subset(TFdata, TFdata$Type == "Internal Standard", select=c("Compound", "RT")))
  ISTDs$RT <- as.numeric(ISTDs$RT)
  target_RT <- as.numeric(unique(TFdata[TFdata$Compound == target, c("RT")]))
  ISTDs  <- subset(ISTDs, (ISTDs$RT< target_RT + RT_window & ISTDs$RT>target_RT - RT_window ))
  return(ISTDs$Compound)
}


#' function selects ISTDs for which rel. standard deviation < relSD.max and mean.min < mean < mean.max
#'
#' @param recoveries output of function recoveryCalc()
#' @param relSD.max
#' @param mean.min  
#' @param mean.max
#' @param sample.code
#'
selectBestISTDs <- function(recoveries, relSD.max = 25, mean.min = 80, mean.max = 120, recovery.code = recovery.code){
  
  #calculate mean & rel. standard deviation
  recoveries$mean <- round(apply(recoveries[, recovery.code], 1, mean, na.rm=TRUE ), 1)
  recoveries$relSD <- apply(recoveries[,  recovery.code], 1, sd, na.rm=TRUE )
  recoveries$relSD <- round(recoveries$relSD / recoveries$mean * 100, 1)
  recoveries <- recoveries[is.finite(recoveries$mean), ]    # delete rows, where all recoveries are NA
  
  # select best ISTDs according to criteria and order them based on relSD
  recoveries <- recoveries[recoveries$relSD < relSD.max & recoveries$mean > mean.min & recoveries$mean < mean.max, ]
  recoveries <- recoveries[ order(recoveries$relSD, decreasing = F), ]   # order ISTDs
  return(recoveries)
}


## From Jen:

### finding the best internal standard to calculate relative recovery with

### data needed for this  is a Tracefinder data file which includes:
## 1. areas of all internal standards in standards and spiked and unspiked samples
## 2. areas of all targets in standards and spiked and unspiked samples

library(reshape)

#' weighting function for the linear regression
#'
#' @param cal.levels 
#' @param calibration which type of weighting is desired? options are X, 1/X, and 1/X^2
#'
#' @return a vector is returned with the cal.levels adjusted by the desired weighting factor
#' @export
#'
#' @examples
#' 
weighting <- function(cal.levels, calibration){
  
  if(calibration == "X"){
    weights <- cal.levels
  }
  
  if(calibration == "1/X"){
    weights <- 1/cal.levels
  }
  
  if(calibration == "1/X2"){
    weights <- 1/(cal.levels^2)
  }
  
  return(weights)
}


### 
#' function for calculating response area of target vs. internal standard
#'
#' @param TFdata TraceFinder quantification file, saved as a data.frame
#' @param target name of target compound, needs to be same as name in TFdata
#' @param internal.standard name of internal standard compound, needs to be same as name in TFdata
#' @param weighting type of weighting desired, options are X, 1/X, and 1/X^2
#' @param colname.compound name of the column in TFdata with compounds names
#'
#' @return output is a data frame, which includes TFdata for target compound, TFdata for internal standard and a new column
#' named 'response.area' 
#' @export
#'
#' @examples
#' 
responseArea <- function(TFdata, target, internal.standard, weighting, colname.compound = "Compound"){
  
  data <- as.data.frame(TFdata)
  col.compound <- which(colnames(data) == colname.compound)
  target.area <- subset(data, data[,col.compound] == target)
  istd.area <- subset(data, data[,col.compound] == internal.standard)
  
  data <- merge(target.area, istd.area, by = "Filename")
  data$response.area <- as.numeric(as.character(data$Area.x)) / as.numeric(as.character(data$Area.y))

  return(data)
  
}

### 
#' function for calculating linear regression with correct calibration data and weighting
#'
#' @param TFdata TraceFinder quantification file, saved as a data.frame
#' @param target name of target compound, needs to be same as name in TFdata
#' @param internal.standard name of internal standard compound, needs to be same as name in TFdata
#' @param weighting type of weighting desired, options are X, 1/X, and 1/X^2
#' @param colname.compound name of the column in TFdata with compound names
#'
#' @return output is an lm object, based on the calibration levels and the response area, using the specified weights. if 
#' no calilbration data is input, lm = "NA"
#' @export
#'
#' @examples
#' 
calibrationCalc <- function(TFdata, target, internal.standard, weighting, colname.compound = "Compound"){
  
  data <- responseArea(TFdata, target = target, internal.standard = internal.standard, weighting = weighting, 
                       colname.compound = colname.compound)

  cal.data <- na.omit(subset(data[,c("Level.x", "response.area")], data$Sample.Type.x == "Cal Std"))
  cal.data$weights <- weighting(cal.levels = cal.data$Level.x, calibration = weighting)
  
  if(nrow(cal.data) == 0){
 #    stop("no calibration possible")
     return(lm <-"NA")
  }else{
     lm <- lm(Level.x ~ response.area, data = cal.data, weights = cal.data$weights)
     return(lm)
  }
}


#' predict sample concentrations based on new calibration model
#'
#' @param TFdata TraceFinder quantification file, saved as a data.frame
#' @param target name of target compound, needs to be same as name in TFdata
#' @param internal.standard name of internal standard compound, needs to be same as name in TFdata
#' @param weighting type of weighting desired, options are X, 1/X, 1/X^2
#' @param colname.compound name of the column in TFdata with compound names
#'
#' @return output is a data frame, which includes TFdata for target compound, TFdata for internal standard and 2 new columns
#' named 'response.area' and 'predictConc'. If prediction is not possible, then TFdata is returned
#' @export
#'
#' @examples
#' 
predictConc <- function(TFdata, target, internal.standard, weighting, colname.compound){
  
  if(is.na(TFdata)){
 #   stop("no calibration was possible")
    return(TFdata)
  }
  
  data <- responseArea(TFdata, target = target, internal.standard = internal.standard, weighting = weighting, 
                       colname.compound = colname.compound)
  

 lm <- calibrationCalc(TFdata, target = target, internal.standard = internal.standard, weighting = weighting, 
                       colname.compound = colname.compound)
  
  if(lm=="NA"){
    return(TFdata)
  }else{
    data$predictConc <- as.numeric(predict(lm, newdata = data.frame(response.area = as.numeric(as.character(data$response.area)))))
    return(data)
  }

}

#' Predicted concentrations are used to calculate relative recovery
#'
#' @param predicted.data use the output from 'predictConc' function
#' @param sample.code text code for samples 
#' @param recovery.code text code for recovery samples
#' @param standard.conc which concentration to use for the standard. options are "theoretical" or "calculated"
#' @param cal.standard text code for calibration standard to be used. use if standard.conc is set to "calculated" 
#' @param spike.level enter spike concentration. use if standard.conc is set to "theoretical"
#' @param colname.compound name of the column in predicted.data with compound names. If using the output of 'predictConc', then
#' this will need to include '.x' or '.y' at the end
#'
#' @return relative recovery, reported as a percentage
#' @export
#'
#' @examples
#' 
recoveryCalc <- function(predicted.data, sample.code, recovery.code, standard.conc = c("theoretical", "calculated"),
                         cal.standard, spike.level, colname.compound){
  
  recovery <- grep(pattern = recovery.code, x = predicted.data$Sample.Name.x)
  samples <- grep(pattern = sample.code, x = predicted.data$Sample.Name.x)
  if(standard.conc == "calculated"){
    standard <- grep(pattern = cal.standard, x = predicted.data$Sample.Name.x)
  }

  if("predictConc" %in% colnames(predicted.data) == FALSE){
    return("NA")
  }else{

  unspike.mean <- as.data.frame(cast(data = predicted.data[c(setdiff(samples, recovery)),], 
                                     formula = as.formula(paste("Sample.Name.x ~", colname.compound)), 
                                     function(x) mean(na.omit(as.numeric(as.character(x)))),
                                     value = 'predictConc'))
  
  spike.mean <- as.data.frame(cast(data = predicted.data[c(intersect(samples, recovery)),], 
                                   formula = as.formula(paste("Sample.Name.x ~", colname.compound)), 
                                   function(x) mean(na.omit(as.numeric(as.character(x)))),
                                   value = 'predictConc'))
  
  if(standard.conc == "calculated"){
    standard.mean <- as.data.frame(cast(data = predicted.data[standard,], 
                                        formula = as.formula(paste("Sample.Name.x ~", colname.compound)), 
                                        function(x) mean(na.omit(as.numeric(as.character(x)))),
                                        value = 'predictConc'))

    relRec <- ((mean(spike.mean[,2], na.rm = TRUE) - mean(unspike.mean[,2],na.rm = TRUE))/mean(standard.mean[,2], na.rm = TRUE))*100
    
    return(round(relRec, 2))
  }

  if(standard.conc == "theoretical"){
    relRec <- ((mean(spike.mean[,2], na.rm = TRUE) - mean(unspike.mean[,2],na.rm = TRUE))/as.numeric(spike.level))*100
    
    return(round(relRec, 2))
  }
  
  }
  
}


