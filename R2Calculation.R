
library(lattice)

### Evaluation of HILIC target screening data

### where data is stored and where to store
main.dir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development"
data.dir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017\\TFexports"

setwd(main.dir)

### data files available
av.data <- list.files(data.dir, full.names = TRUE)
data.name <- list.files(data.dir, full.names = FALSE)


### Functions to compute summary

### calculate R2 for every target with every ISTD

#calibration <- "X"

get.R2 <- function(x, y, cal.model="1/X"){
  
  ## x is calibration levels, that is, the expected concentrations
  ## y is the response.area that was measured for that target compound
  
  if(cal.model == "X"){
    
    regress <- lm(x~y)
    return(summary(regress)$r.squared)
    
  }
  
  if(cal.model == "1/X"){
    
    regress <- lm(x~y, weights = 1/x)
    return(summary(regress)$r.squared)
    
  }
  
  if(cal.model == "1/X2"){
    
    regress <- lm(x~y, weights = 1/(x^2))
    return(summary(regress)$r.squared)
    
  }
}

## calculate the R2 for each target with each istd

R2summary <- function(targets, istds, calibration, cal.levels){
  
  R2.summary <- data.frame(matrix("NA", nrow = length(targets), ncol = length(istds)), stringsAsFactors = FALSE)
  
  for(m in 1:length(targets)){
    
    target.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == targets[m] & 
                                                           TFdata$Sample.Type == "Cal Std" & 
                                                           TFdata$Level %in% cal.levels, 
                                                           select = c(Area)))))
    
    if(sum(!is.na(target.area)) == 0){
      R2.summary[m,] <- "ND"
      next()
    }#else{
    
    for(n in 1:length(istds)){
      
      istd.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == istds[n] & TFdata$Sample.Type == "Cal Std", select = c(Area)))))
      response.area <- target.area/istd.area
      
      ifelse(sum(!is.na(response.area))==0, R2.summary[m,n] <- "NA", R2.summary[m,n] <- get.R2(levels,response.area, cal.model = calibration))
      
    }
  }
  
  colnames(R2.summary) <- istds
  rownames(R2.summary) <- targets
  return(R2.summary)
}


## plotting of target and istd with desired calibration model

CalCurvePlot <- function(target, istd, calibration, cal.levels){
  
    target.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == target & 
                                                           TFdata$Sample.Type == "Cal Std" & 
                                                           TFdata$Level %in% cal.levels, 
                                                         select = c(Area)))))
    
    if(sum(!is.na(target.area)) == 0){
     # R2.summary[m,] <- "ND"
      print("Error: Target compound not detected")
      next()
    }
    
      istd.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == istd &
                                                           TFdata$Sample.Type == "Cal Std" &
                                                           TFdata$Level %in% cal.levels, 
                                                         select = c(Area)))))
      response.area <- target.area/istd.area
      
      R2<- c()
      
      ifelse(sum(!is.na(response.area))==0, R2 <- "NA", R2 <- get.R2(levels,response.area, cal.model = calibration))
      
      if(calibration == "X"){
        lm <- lm(response.area~cal.levels)
        plot(cal.levels, response.area, main = paste("Calibration Curve with ", calibration, sep = ""),
             xlab = "Calibration Standards", ylab = "Response Area")
        abline(lm)
        legend ("topleft", legend = paste("R2 = ", round(R2, 2), sep = ""))
      }
      
      if(calibration == "1/X"){
        lm <- lm(response.area~cal.levels, weights = 1/cal.levels)
        plot(cal.levels, response.area, main = paste("Calibration Curve with ", calibration, sep = ""),
             xlab = "1 / Calibration Standards", ylab = "Response Area")
        abline(lm)
        legend ("topleft", legend = paste("R2 = ", round(R2, 2), sep = ""))
        
      }
      
      if(calibration == "1/X2"){
        lm <- lm(response.area~cal.levels, weights = 1/(cal.levels^2))
        plot(cal.levels, response.area, main = paste("Calibration Curve with ", calibration, sep = ""),
             xlab = "1 / Calibration Standards", ylab = "Response Area")
        abline(lm)
        legend ("topleft", legend = paste("R2 = ", round(R2, 2), sep = ""))
        
      }

}



### function to find the best R2 for each target compound

best.R2 <- function(R2.summary){
  
  best.R2 <- unlist(apply(R2.summary, 1, function(x) which.max(x)))
  for(i in 1:length(best.R2)){
  best.R2[i] <- as.numeric(as.character(paste(R2.summary[i,best.R2[i]])))
  }
  return(best.R2)
}



### Calculations across all data files
all.R2s <- list()

for(f in 1:length(av.data)){
  
  ### load data
  TFdata <- read.csv(av.data[f], header = TRUE)
  TFdata <- subset(TFdata, TFdata$Peak.Label == "T1")
  
  targets <- unique(TFdata[which(TFdata$Type == "Target Compound"), 1])
  istds <- unique(TFdata[which(TFdata$Type == "Internal Standard"), 1])
  samples <- unique(TFdata[,"Sample.Name"])
  levels <- unique(TFdata[which(TFdata$Sample.Type == "Cal Std"),"Level"])
  
  
  ### calculations of all R2
  
  R2.summary.X.all <- R2summary(targets, istds, calibration = "X", cal.levels = levels)
  R2.summary.1.X.all <- R2summary(targets, istds, calibration = "1/X", cal.levels = levels)
  R2.summary.1.X.low <- R2summary(targets, istds, calibration = "1/X", cal.levels = levels[1:5])
  R2.summary.1.X.high <- R2summary(targets, istds, calibration = "1/X", cal.levels = levels[7:12])
  
  
  R2.summary.1.X.all[R2.summary.1.X.all == "1"] <- "reject"
  R2.summary.X.all[R2.summary.X.all == "1"] <- "reject"
  R2.summary.1.X.low[R2.summary.1.X.low == "1"] <- "reject"
  R2.summary.1.X.high[R2.summary.1.X.high == "1"] <- "reject"
  
  ### calculations of best R2 for each target
  
  best.R2.1.X.all <- best.R2(R2.summary.1.X.all)
  best.R2.X.all <- best.R2(R2.summary.X.all)
  best.R2.1.X.low <- best.R2(R2.summary.1.X.low)
  best.R2.1.X.high <- best.R2(R2.summary.1.X.high)
  
  all <- data.frame(R2 = c(best.R2.1.X.all, best.R2.1.X.low, best.R2.1.X.high, best.R2.X.all),
                        type = c(rep("1.X.all", length(best.R2.1.X.all)),
                                 rep("1.X.low", length(best.R2.1.X.low)),
                                 rep("1.X.high", length(best.R2.1.X.high)),
                                 rep("X", length(best.R2.X.all))
                                                )
                    )
  all$data <- paste(data.name[f])
  
  all.R2s[[f]] <- all
  
  # png(paste("BWplot_", data.name[f], ".png", sep = ""), width = 800, height = 500)
  # boxplot(R2~type, data = all.R2s, main = data.name[f])
  # dev.off()
  # 
  
}


allR2 <- do.call(rbind, all.R2s)
bwplot(R2 ~ data |type, data = allR2)

unique(allR2$data)
allR2$data <- factor(allR2$data, levels = rev(c("BEHAmide_NPStds_targets.csv", 
                                                    "BEHAmide_NPStds_neg_targets.csv",
                                                    "BEHAmide_NPBuechiStds_targets.csv",
                                                    "BEHAmide_NPBuechiStdsQEPlus_neg_targets.csv",
                                                    "BEHAmide_NPBuechiStds_pos_QE_targets.csv",
                                                    "BEHAmide_NPBuechiStds_neg_QE_targets.csv",
                                                    "BEHAmide_TWBuechiStds_pos_QE_targets.csv",
                                                    "BEHAmide_TWBuechiStds_neg_QE_targets.csv")), order = TRUE)
levels(allR2$type)
allR2$type <- factor(allR2$type, levels = rev(c("1.X.all", "X", "1.X.low", "1.X.high")), order = TRUE)

#densityplot( ~ R2 |type, groups = data, data = allR2)
stripplot(data ~ R2 |type, groups = data, data = allR2,
          par.strip.text = list(cex=1),
          strip = strip.custom(factor.levels = c("1/X Higher Calibration", 
                                                 "1/X Lower Calibration",
                                                 "X  All Calibration", 
                                                 "1/X All Calibration", "" ),
                               bg = "lightgray"),
          scales = list(y = list(cex = 0.8, labels = rev(c("NPStds positive QEPlus",
                                                     "NPStds negative QEPlus",
                                                     "NPStds Buechi positive QEPlus",
                                                     "NPStds Buechi negative QEPlus",
                                                     "NPStds Buechi positive QE",
                                                     "NPStds Buechi negative QE",
                                                     "TWStds Buechi positive QE",
                                                     "TWStds Buechi negative QE"
                                                     )),
                        x = list(cex = 1))
                        ),
          #ylab = list(label = "log10( ISTD Area in Matrix / ISTD Area in Standard )", cex = 1),
          #ylim = c(-7,7),
          xlab = list(label = expression('R'^2), cex = 1)
          )

#histogram( ~ R2 |type, data = allR2)
