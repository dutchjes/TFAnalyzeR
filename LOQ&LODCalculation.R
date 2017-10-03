


### Evaluation of HILIC target screening data

### where data is stored and where to store
main.dir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development"
data.dir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017\\TFexports"

setwd(main.dir)

### data files available
av.data <- list.files(data.dir, full.names = TRUE)
data.name <- list.files(data.dir, full.names = FALSE)

neg.data <- c(1,4,5,7) ## which files in the av.data are for negative mode data?
pos.data <- c(2,3,6,8) ## which files in the av.data are for positive mode data?

### load data
TFdata <- read.csv(av.data[1], header = TRUE)
TFdata <- subset(TFdata, TFdata$Peak.Label == "T1")

targets <- unique(TFdata[which(TFdata$Type == "Target Compound"), 1])
istds <- unique(TFdata[which(TFdata$Type == "Internal Standard"), 1])
samples <- unique(TFdata[,"Sample.Name"])
levels <- unique(TFdata[which(TFdata$Sample.Type == "Cal Std"),"Level"])

target.rts <- c(unlist(subset(TFdata, TFdata$Type == "Target Compound" & 
                      # TFdata$Sample.Type == "Cal Std" & 
                       TFdata$Level == 500,
                     select = c(Actual.RT))))



#### calculation of lowest cal standard detected

LOQsummary <- function(TFdata, targets, istds, cal.levels){
  
  loq.summary <- data.frame(matrix("NA", nrow = length(targets), ncol = length(istds)), stringsAsFactors = FALSE)
  
  for(m in 1:length(targets)){
    
    target.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == targets[m] & 
                                                           TFdata$Sample.Type == "Cal Std" & 
                                                           TFdata$Level %in% cal.levels, 
                                                         select = c(Area)))))
    
    if(sum(!is.na(target.area)) == 0){
      loq.summary[m,] <- "ND"
      next()
    }#else{
    
    for(n in 1:length(istds)){
      
      istd.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == istds[n] & TFdata$Sample.Type == "Cal Std", select = c(Area)))))
      response.area <- target.area/istd.area
      
      ifelse(sum(!is.na(response.area))==0, 
             loq.summary[m,n] <- "NA", loq.summary[m,n] <- min(cal.levels[!is.na(response.area)]))
      
    }
  }
  
  colnames(loq.summary) <- istds
  rownames(loq.summary) <- targets
  return(loq.summary)
}




### function to find the lowest LOQ for each target compound

lowest.loq <- function(loq.summary){
  
  lowestloq <- apply(loq.summary, 1, function(x) min(x))
  lowestloq[unlist(lapply(lowestloq, length))==0] <- NA
  lowestloq <- unlist(lowestloq)
  
  # for(i in 1:length(lowestloq)){
  #   ifelse(is.na(lowestloq[i]), next(),
  #   lowestloq[i] <- as.numeric(as.character(paste(loq.summary[i,lowestloq[i]]))))
  # }
  
  return(as.numeric(lowestloq))
}


### function to find the lowest LOD (detection of target compound)

lowest.lod <- function(TFdata, targets, cal.levels){
  
  lod.summary <- c()
  
  for(m in 1:length(targets)){
    
    target.area <- as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == targets[m] & 
                                                           TFdata$Sample.Type == "Cal Std" & 
                                                           TFdata$Level %in% cal.levels, 
                                                         select = c(Area)))))
    
    ifelse(sum(!is.na(target.area)) == 0, 
           lod.summary[m] <- "ND", lod.summary[m] <- as.numeric(as.character(min(cal.levels[!is.na(target.area)]))))
  }
    
  names(lod.summary) <- targets
  return(as.numeric(lod.summary))
}


#### Calculations across data files

best.loq <- list()
best.lod <- list()

for(f in 1:length(av.data)){
  
  ### load data
  TFdata <- read.csv(av.data[f], header = TRUE)
  TFdata <- subset(TFdata, TFdata$Peak.Label == "T1")
  
  targets <- unique(TFdata[which(TFdata$Type == "Target Compound"), 1])
  istds <- unique(TFdata[which(TFdata$Type == "Internal Standard"), 1])
  samples <- unique(TFdata[,"Sample.Name"])
  levels <- unique(TFdata[which(TFdata$Sample.Type == "Cal Std"),"Level"])
  
  target.rts <- as.numeric(as.character(unlist(subset(TFdata, TFdata$Type == "Target Compound" & 
                                TFdata$Level == 500,
                                select = c(Actual.RT)))))
  
  loqsummary <- LOQsummary(TFdata = TFdata, targets = targets, istds = istds, cal.levels = levels)
  best.loq[[f]] <- list(data = lowest.loq(loq.summary = loqsummary),
                        rts = target.rts)
  best.lod[[f]] <- list(data = lowest.lod(TFdata = TFdata, targets = targets, cal.levels = levels),
                        rts = target.rts)

}


library(lattice)

#### Negative mode comparison

neg <- neg.data
n.neg <- length(best.loq[[neg[1]]]$data)
neg.loq <- data.frame(loqs = unlist(lapply(best.loq[neg], function(l) l$data)),
                      data = c(rep("NPBuechi_QE", n.neg),
                               rep("NPBuechi_QEPlus", n.neg),
                               rep("NPStds_QEPlus", n.neg),
                               rep("TWBuechi_QE", n.neg)
                      ),
                      rts = unlist(lapply(best.loq[neg], function(s) s$rts))
                      
)
boxplot(log10(loqs) ~ data, data = neg.loq)
histogram(~loqs | data, data = neg.loq, nint = 100)
xyplot(loqs ~ rts, data = neg.loq, groups = data)

n.neg <- length(best.lod[[neg[1]]])
neg.lod <- data.frame(lods = unlist(lapply(best.lod[neg], function(l) l$data)),
                      data = c(rep("NPBuechi_QE", n.neg),
                               rep("NPBuechi_QEPlus", n.neg),
                               rep("NPStds_QEPlus", n.neg),
                               rep("TWBuechi_QE", n.neg)
                      ),
                      rts = unlist(lapply(best.lod[neg], function(s) s$rts))
                      
)
boxplot(log10(lods) ~ data, data = neg.lod)
histogram(~lods | data, data = neg.lod, nint = 100)
xyplot(lods ~ rts, data = neg.lod, groups = data, ylim = c(-5,100), auto.key = TRUE)

##### Positive mode comparison

pos <- pos.data
n.pos <- unlist(lapply(best.loq[pos], function(x) length(x$data)))
pos.loq <- data.frame(loqs = unlist(lapply(best.loq[pos], function(l) l$data)),
                      data = c(rep("NPBuechi_QE", n.pos[1]),
                               rep("NPBuechi_QEPlus", n.pos[2]),
                               rep("NPStds_QEPlus", n.pos[3]),
                               rep("TWBuechi_QE", n.pos[4])
                      ),
                      rts = unlist(lapply(best.loq[pos], function(s) s$rts))
)
boxplot(log10(loqs) ~ data, data = pos.loq)
histogram(~loqs | data, data = pos.loq, nint = 100)
xyplot(loqs ~ rts, data = pos.loq, groups = data, auto.key = TRUE, pch = 16)

n.pos <- unlist(lapply(best.lod[pos], function(x) length(x$data)))
pos.lod <- data.frame(lods = unlist(lapply(best.lod[pos], function(l) l$data)),
                      data = c(rep("NPBuechi_QE", n.pos[1]),
                               rep("NPBuechi_QEPlus", n.pos[2]),
                               rep("NPStds_QEPlus", n.pos[3]),
                               rep("TWBuechi_QE", n.pos[4])
                      ),
                      rts = unlist(lapply(best.lod[pos], function(s) s$rts))
)
boxplot(log10(lods) ~ data, data = pos.lod)
histogram(~lods | data, data = pos.lod, nint = 100)
xyplot(lods ~ rts, data = pos.lod, groups = data, auto.key = TRUE, pch = 16)


### summary tables

sum.lod <- lapply(best.lod, function(x) table(x$data))
sum.loq <- lapply(best.loq, function(x) table(x$data))

summary.lod <- data.frame(matrix(NA, nrow = length(av.data), ncol = length(levels)+1))
summary.loq <- data.frame(matrix(NA, nrow = length(av.data), ncol = length(levels)+1))

colnames(summary.lod) <- c(levels, "ND")
colnames(summary.loq) <- c(levels, "ND")
rownames(summary.lod) <- data.name
rownames(summary.loq) <- data.name

for(f in 1:length(av.data)){
  locs <- sapply(names(sum.loq[[f]]), function(x) which(x == colnames(summary.loq)))
  summary.loq[f,locs] <- sum.loq[[f]]
  
  locs <- as.integer(sapply(names(sum.lod[[f]]), function(x) which(x == colnames(summary.lod))))
  summary.lod[f, locs] <- sum.lod[[f]]
}

rownames(summary.lod)

new.rows <- c(6,5,3,4,2,1,8,7)
summary.lod <- summary.lod[new.rows,]
summary.loq <- summary.loq[new.rows,]

write.csv(summary.lod, "LODsummary.csv")
write.csv(summary.loq, "LOQsummary.csv")

png("LODComparison2.png", width = 3300, height = 1500, res = 300)
barplot(t(summary.lod), names.arg = c(rep(c("pos", "neg"),4)), xlim = c(0,12), legend = rownames(t(summary.lod)), # xaxt = "n",
        ylab = "Number of Compounds", col = topo.colors(length(rownames(t(summary.lod)))+1), main = "LOD Comparison",
        args.legend = list(title = "Concentration (ug/L)", bty = "n")
        )
mtext(side = 1, text = c("NPStds \n (QEPlus)", "NPBuechi \n (QEPlus)", "NPBuechi \n (QE)", "TWBuechi \n (QE)"), line = 3,
      at = c(1.3, 3.7, 6.2, 8.6))
dev.off()

png("LOQComparison.png", width = 1000, height = 600, res = 100)
barplot(t(summary.loq), names.arg = c(rep(c("pos", "neg"),4)), xlim = c(0,12), legend = rownames(t(summary.loq)), # xaxt = "n",
        ylab = "Number of Compounds", col = topo.colors(length(rownames(t(summary.loq)))+1), main = "LOQ Comparison",
        args.legend = list(title = "Concentration (ug/L)", bty = "n")
        )
mtext(side = 1, text = c("NPStds \n (QEPlus)", "NPBuechi \n (QEPlus)", "NPBuechi \n (QE)", "TWBuechi \n (QE)"), line = 3,
     at = c(1.3, 3.7, 6.2, 8.6))
dev.off()



