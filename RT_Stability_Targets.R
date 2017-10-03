
library(reshape)
library(lattice)

### Evaluation of HILIC target screening data

### where data is stored and where to store
main.dir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017"
data.dir <- "Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017\\TFexports"

setwd(main.dir)

### data files available
av.data <- list.files(data.dir, full.names = TRUE)
data.name <- list.files(data.dir, full.names = FALSE)

neg.data <- c(1,4,5,7) ## which files in av.data are for negative mode data?
pos.data <- c(2,3,6,8) ## which files in av.data are for positive mode data?

#### compiling all data
target.rts.all <- list()

for(f in 1:length(av.data)){

## load data
TFdata <- read.csv(av.data[f], header = TRUE)
TFdata <- subset(TFdata, TFdata$Peak.Label == "T1")

targets <- unique(TFdata[which(TFdata$Type == "Target Compound"), 1])
istds <- unique(TFdata[which(TFdata$Type == "Internal Standard"), 1])
samples <- unique(TFdata[,"Sample.Name"])
levels <- unique(TFdata[which(TFdata$Sample.Type == "Cal Std"),"Level"])


target.rts <- data.frame(matrix(NA, nrow = length(targets), ncol = 4))
colnames(target.rts) <- c("TargetCompound", "MeanRT", "StdRT", "RSDRT")

for(m in 1:length(targets)){
  
  target.rts[m,1] <- paste(targets[m])
  
  target.rts[m,2] <- mean(as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == targets[m] & 
                                                              TFdata$Sample.Type == "Cal Std" & 
                                                              TFdata$Level %in% levels, 
                                                            select = c(Actual.RT))))))
  
  target.rts[m,3] <- sd(as.numeric(as.character(unlist(subset(TFdata, TFdata[,1] == targets[m] & 
                                                               TFdata$Sample.Type == "Cal Std" & 
                                                               TFdata$Level %in% levels, 
                                                             select = c(Actual.RT))))))
  
  target.rts[m,4] <- target.rts[m,3] / target.rts[m,2]
}
  
target.rts.all[[f]] <- target.rts

#densityplot(target.rts[,2])
#densityplot(target.rts[,3])
# #densityplot(target.rts[,4])
# png(paste("Histogram_RTRSD_", data.name[f], ".png", sep = ""), width = 1100, height = 600)
# hist.now <- hist(target.rts[,4])
# hist(target.rts[,4], labels = TRUE, ylim = c(0,max(hist.now$counts)+10),
#      col = "lightblue", main = "Retention Time RSD of Target Compounds",
#      xlab = "RSD", n=20, xlim = c(0,round(max(na.omit(target.rts[,4])),1)), sub = paste(data.name[f]))
# legend("topright", legend = paste("n = ", length(na.omit(target.rts[,2])), sep = ""))
#dev.off()

}


## negative ESI
mode <- unlist(lapply(target.rts.all, function(x) nrow(x)))

neg.rts <- data.frame()
for(f in 1:length(neg.data)){
  
  take <- neg.data[f]
  
  if(f == 1){
    neg.rts <- data.frame(target.rts.all[[take]][,c(1,2)])
  }else{
      neg.rts <- data.frame(cbind(neg.rts, target.rts.all[[take]][,2]))
    }
}

##here you need to give names to the columns of your table
## column 1 always corresponds to the Target Compound so that should remain the same
## all subsequent columns are the different data files processed as defined by av.data[[neg.data]]

colnames(neg.rts) <- c("TargetCompound", "NPBuechi_QE", "NPBuechi_QEPlus", "NPStds_QEPlus", "TWBuechi_QE")
neg.rts <- neg.rts[,c(1,4,3,2,5)] ## i reordered the columns for plotting later

##decide which column of RTs should be the reference RT. Here I used neg.rts$NPStds_QEPlus...
neg.rts.shift <- melt(apply(neg.rts[,3:5], 2, function(x) neg.rts$NPStds_QEPlus - x))


png("RTShifts_neg.png", width = 1100, height = 500)
bwplot(value ~ X2, data = neg.rts.shift,
       ylab = list(cex = 1.3, label = "RT Shift of Target Compounds (minutes)"),
       xlab = list(cex = 1.3, label = NULL),
       scales = list(x = list(cex = 1.2),
                     y = list(cex = 1.2)),
       panel = function(x, y, ...) { 
         panel.abline(h=0, lwd = 2)
         panel.bwplot(x, y, fill = "lightblue",...) 
         nn <- table(x) 
         panel.text(paste("n =", nn), x = seq_along(nn), 
                    y = current.panel.limits()$y[1], pos = 3)
       }
       )
dev.off()


png("RTShiftsHistogram_neg.png", width = 1100, height = 500)
histogram( ~ value| X2, data = neg.rts.shift, labels = TRUE, type = "count", nint = 50,
           ylab = list(label = "Counts", cex = 1.3),
           xlab = list(label = "RT Shift (minutes)", cex = 1.3),
           scales = list(x = list(cex = 1.2),
                         y = list(cex = 1.2)),
           panel=function(x, ...){
             panel.histogram(x=x,...)
             m <- length(na.omit(x))
             panel.text(x=10,y=20,labels=paste("n=",m, sep = ""))
           }
           )
dev.off()


### positive ESI
pos.rts <- data.frame()
for(f in 1:length(pos.data)){

  take <- pos.data[f]

  if(f == 1){
    pos.rts <- data.frame(target.rts.all[[take]][,c(1,2)])
  }else{
    pos.rts <- data.frame(merge(pos.rts, target.rts.all[[take]], by = "TargetCompound"))[,1:(f+1)]
  }
}

##here you need to give names to the columns of your table
## column 1 always corresponds to the Target Compound so that should remain the same
## all subsequent columns are the different data files processed as defined by av.data[[pos.data]]

colnames(pos.rts) <- c("TargetCompound", "NPBuechi_QE", "NPBuechi_QEPlus", "NPStds_QEPlus", "TWBuechi_QE")
pos.rts <- pos.rts[,c(1,4,3,2,5)] ## i reordered the columns for plotting later

##decide which column of RTs should be the reference RT. Here I used pos.rts$NPStds_QEPlus...
pos.rts.shift <- melt(apply(pos.rts[,3:5], 2, function(x) pos.rts$NPStds_QEPlus - x))

png("RTShifts_pos.png", width = 1100, height = 500)
bwplot(value ~ X2, data = pos.rts.shift,
       ylab = list(cex = 1.3, label = "RT Shift of Target Compounds (minutes)"),
       xlab = list(cex = 1.3, label = NULL),
       scales = list(x = list(cex = 1.2),
                     y = list(cex = 1.2)),
       panel = function(x, y, ...) { 
         panel.abline(h=0, lwd = 2)
         panel.bwplot(x, y, fill = "lightblue",...) 
         nn <- table(x) 
         panel.text(paste("n =", nn), x = seq_along(nn), 
                    y = current.panel.limits()$y[1], pos = 3)
       }
)
dev.off()

png("RTShiftsHistogram_pos.png", width = 1100, height = 500)
histogram( ~ value| X2, data = pos.rts.shift, labels = TRUE, type = "count", nint = 50,
           ylab = list(label = "Counts", cex = 1.3),
           xlab = list(label = "RT Shift (minutes)", cex = 1.3),
           scales = list(x = list(cex = 1.2),
                         y = list(cex = 1.2)),
           panel=function(x, ...){
             panel.histogram(x=x,...)
             m<-length(na.omit(x))
             panel.text(x=3,y=20,labels=paste("n=",m, sep = ""))
           }
           )
dev.off()

### combining data from both positive and negative ESI modes
all.rt.shifts <- rbind(neg.rts.shift, pos.rts.shift)

png("RTShiftsHistogram_all.png", width = 1100, height = 500)
histogram( ~ value| X2, data = all.rt.shifts, labels = TRUE, type = "count", nint = 50,
           xlim = c(-15,10),
           ylab = list(label = "Counts", cex = 1.3),
           xlab = list(label = "RT Shift (minutes)", cex = 1.3),
           scales = list(x = list(cex = 1.2, tick.number = 10),
                         y = list(cex = 1.2)),
           panel=function(x, ...){
             panel.histogram(x=x,...)
             m<-length(na.omit(x))
             panel.text(x=-10,y=40,labels=paste("n = ",m, sep = ""))
           }
)
dev.off()

histogram.computations <-
  function(x, breaks, equal.widths = TRUE,
           type = "density", nint, ...)
  {
    if (is.null(breaks))
    {
      breaks <-
        if (is.factor(x)) seq_len(1 + nlevels(x)) - 0.5
      else if (equal.widths) do.breaks(range(x, finite = TRUE), nint)
      else quantile(x, 0:nint/nint, na.rm = TRUE)
    }
    hist(x, breaks = breaks, plot = FALSE)
  }

haha <- histogram( ~ value| X2, data = all.rt.shifts, labels = TRUE, type = "count", nint = 50,
                  xlim = c(-15,10),
                  ylab = list(label = "Counts", cex = 1.3),
                  xlab = list(label = "RT Shift (minutes)", cex = 1.3),
                  scales = list(x = list(cex = 1.2, tick.number = 10),
                                y = list(cex = 1.2)),
                  panel=function(x, ...){
                    panel.histogram(x=x,...)
                    m<-length(na.omit(x))
                    panel.text(x=-10,y=40,labels=paste("n = ",m, sep = ""))
                  }
)


# colfunc <- colorRampPalette(c("white", "blue"))
# colRT <- colfunc(length(seq(-8,8,0.5)))
# my.legend <- packGrob(
#   draw.colorkey(key=list(col = colRT, at = do.breaks(range(-8,8),20))), 
#   textGrob("My title", x = 0, y = 0.5, just = c("left", "centre")), 
#   height=unit(2, "lines"),side="top", dynamic=T)

## what are the names of the different scenarios?
levels(all.rt.shifts$X2) 

##reorder the names in the way you want for plotting
all.rt.shifts$X2 <- factor(all.rt.shifts$X2, levels = c("NPBuechi_QEPlus", "NPBuechi_QE", "TWBuechi_QE"), order = TRUE)

## adding some formatting to the names...
levels(all.rt.shifts$X2) <- c("NP Buechi\n QEPlus", "NP Buechi\n QE", "TW Buechi\n QE")


png("RTShiftsHistogram2_all.png", width = 3300, height = 1500, res = 300)
histogram( ~ value| X2, data = all.rt.shifts, labels = TRUE, type = "count", breaks = seq(-25,15,0.5), #nint = 50,
                   xlim = c(-8,8),
                   ylab = list(label = "Counts", cex = 1.3),
                   xlab = list(label = "RT Shift (minutes)", cex = 1.3),
                   scales = list(x = list(cex = 1.2, tick.number = 10),
                                 y = list(cex = 1.2)),
           strip = strip.custom(bg = "lightgray", par.strip.text = list(cex = 1)),
           par.settings = list(layout.heights = list(strip = 2.45)), ## this makes the strip taller because I used two line names
           #legend=list(right=list(fun=my.legend)),
           #par.strip.text = list(label = c("NP Buechi", "NP Buechi", "TW Buechi")),
                   panel=function(x, ...){
                     panel.histogram(x=x,# col = colRT,
                                     ...)
                     h <- histogram.computations(x=x,...)
                   #  par.settings = list(box.rectangle = list(fill = colRT))
                     m<-length(na.omit(x))
                     panel.text(x=-10,y=40,labels=paste("n = ",m, sep = ""))
                     panel.text(x = h$breaks+0.35, y = h$counts+1.5, labels = ifelse(h$counts==0, "", h$counts), cex = 0.7)
                   }
           
)
dev.off()

# png("RTShiftsHistogram3_all.png", width = 1100, height = 500, res = 300)
# histogram( ~ value| X2, data = all.rt.shifts, labels = TRUE, type = "count", breaks = seq(-25,15,0.5), #nint = 50,
#            xlim = c(-8,8),
#            ylab = list(label = "Counts", cex = 1.3),
#            xlab = list(label = "RT Shift (minutes)", cex = 1.3),
#            scales = list(x = list(cex = 1.2, tick.number = 10),
#                          y = list(cex = 1.2)),
#            strip = strip.custom(bg = "lightgray"),
#            
#            panel=function(x, ...){
#              panel.histogram(x=x,...)
#              h <- histogram.computations(x=x,...)
#              m<-length(na.omit(x))
#              panel.text(x=-10,y=40,labels=paste("n = ",m, sep = ""))
#              panel.text(x = h$breaks+0.35, y = h$counts+1.5, 
#                         labels = ifelse(h$counts==0, "", round(h$breaks, 1)), cex = 0.75)
#            }
# )
# dev.off()
