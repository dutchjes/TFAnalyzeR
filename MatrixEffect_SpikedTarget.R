

library(lattice)
library(reshape)
library(gridExtra)
library(plotrix)
library(datapasta)

#### Calculation of Recovery of Target Compound

## where to save everything?
setwd("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\1st test - May 2017\\Targets")

## where is the tracefinder file?
sample.data <- read.csv("HILIC_pos_everything_samples_BEHAmide_targets.csv", header = TRUE)

sample.data <- subset(sample.data, sample.data$Peak.Label == "T1") ## select only monoisotopic peaks

istd.data <- subset(sample.data, sample.data$Type == "Internal Standard")
target.data <- subset(sample.data, sample.data$Type == "Target Compound")

## for each relevant sample type, define a vector
blanks <- grep(pattern = "ACN-", x = target.data$Sample.Name)
STDA <- grep(pattern = "STDA_", x = target.data$Sample.Name)
STDB <- grep(pattern = "STDB_", x = target.data$Sample.Name)
STDC <- grep(pattern = "STDC_", x = target.data$Sample.Name)
VKB <- grep(pattern = "VKB_", x = target.data$Sample.Name)
NachO3 <- grep(pattern = "NachO3_", x = target.data$Sample.Name)
GAK5b <- grep(pattern = "GAK5b_", x = target.data$Sample.Name)
Dry <- grep(pattern = "Dry_", x = target.data$Sample.Name)
Wet <- grep(pattern = "Wet_", x = target.data$Sample.Name)
Recovery <- grep(pattern = "_A_", x = target.data$Sample.Name)


### function for calculating the average area of spiked target compound
### NOTE: in the following analysis I used the Area output from TF. BUT normally it would be correct the use the
### concentration. I just didn't do a quantification and so used the areas instead. But in the analysis below,
### define peak.area = 'Concentration' and then it will be calculated correctly

targetSpike.mean <- function(data, formula, peak.area, ...){
  
  colMeans(cast(formula, data = data, 
                function(x) as.numeric(as.character(x)), value = peak.area),
           na.rm = TRUE)
  
}

### Selecting which pure standard to compare to
## in the end you need a data frame named 'NPstd' which contains only the target compund in the standard you selected

### which standard to compare to? CC9 is the closest to the spike amt, absolute ng on the column
# here for the pure NP standard from Test 1 (was dried)
# NPstd <- target.data[c(STDC),]
# # NPstd <- targetSpike.mean(data = NPstd, formula = Sample.Name ~ Compound,
# #                           peak.area = 'Response.Ratio', na.rm = TRUE)
# NPstd <- targetSpike.mean(data = NPstd, formula = Sample.Name ~ Compound,
#                           peak.area = 'Area', na.rm = TRUE)

## here for the pure NP standard from Test 2 (was not dried)
# std.data <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017\\TFexports\\BEHAmide_NPStds_targets.csv",
#                      header = TRUE)
# levels(std.data$Sample.Name)
# NPstd <- subset(std.data, std.data$Sample.Name == "CC9_01" & std.data$Type == "Target Compound" & std.data$Peak.Label == "T1")

# here for the Buechi NP standard
std.data <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017\\TFexports\\BEHAmide_NPBuechiStds_targets.csv",
                     header = TRUE)
levels(std.data$Sample.Name)
NPstd <- subset(std.data, std.data$Sample.Name == "CC9_NP_01" & std.data$Type == "Target Compound" & std.data$Peak.Label == "T1")



### calculate area of each compound in spiked samples
### selected the data which was a recovery sample from the right sample type and buechi treatment
spike.VKB.wet <- target.data[c(intersect(VKB, intersect(Wet, Recovery))),]
spike.NachO3.wet <- target.data[c(intersect(NachO3, intersect(Wet, Recovery))),]
spike.GAK5b.wet <- target.data[c(intersect(GAK5b, intersect(Wet, Recovery))),]

spike.VKB.dry <- target.data[c(intersect(VKB, intersect(Dry, Recovery))),]
spike.NachO3.dry <- target.data[c(intersect(NachO3, intersect(Dry, Recovery))),]
spike.GAK5b.dry <- target.data[c(intersect(GAK5b, intersect(Dry, Recovery))),]

### Test 2 has more target compounds, subset to keep the same set
# NPstd <- NPstd[which(c(NPstd$ï..Compound %in% spike.VKB.wet$Compound)==TRUE),]
# NPstd <- as.vector(NPstd[,c("Area")])

### calculate average of each compound in spiked samples (since there are multiple injections)
## separate for each sample type and treatment you're interested in
## see NOTE in lines 36-38 about defining peak.area
spike.VKB.wet <- targetSpike.mean(data = spike.VKB.wet, formula = Sample.Name ~ Compound,
                                  peak.area = 'Area', na.rm = TRUE)

spike.NachO3.wet <- targetSpike.mean(data = spike.NachO3.wet, formula = Sample.Name ~ Compound,
                                     peak.area = 'Area', na.rm = TRUE)

spike.GAK5b.wet <- targetSpike.mean(data = spike.GAK5b.wet, formula = Sample.Name ~ Compound,
                                    peak.area = 'Area', na.rm = TRUE)

spike.VKB.dry <- targetSpike.mean(data = spike.VKB.dry, formula = Sample.Name ~ Compound,
                                  peak.area = 'Area', na.rm = TRUE)

spike.NachO3.dry <- targetSpike.mean(data = spike.NachO3.dry, formula = Sample.Name ~ Compound,
                                     peak.area = 'Area', na.rm = TRUE)

spike.GAK5b.dry <- targetSpike.mean(data = spike.GAK5b.dry, formula = Sample.Name ~ Compound,
                                    peak.area = 'Area', na.rm = TRUE)


### calculate area of each compound in unspiked samples
## separate for each sample type and treatment you're interested in
nospike.VKB.wet <- target.data[c(setdiff(intersect(Wet, VKB), Recovery)),]
nospike.NachO3.wet <- target.data[c(setdiff(intersect(Wet, Recovery), NachO3)),]
nospike.GAK5b.wet <- target.data[c(setdiff(intersect(Wet, Recovery), GAK5b)),]

nospike.VKB.dry <- target.data[c(setdiff(intersect(Dry, VKB), Recovery)),]
nospike.NachO3.dry <- target.data[c(setdiff(intersect(Dry, Recovery), NachO3)),]
nospike.GAK5b.dry <- target.data[c(setdiff(intersect(Dry, Recovery), GAK5b)),]

### calculate average of each compound in unspiked samples (since there are multiple injections)
## separate for each sample type and treatment you're interested in
## see NOTE in lines 36-38 about defining peak.area
nospike.VKB.wet <- targetSpike.mean(data = nospike.VKB.wet, formula = Sample.Name ~ Compound,
                                    peak.area = 'Area', na.rm = TRUE)

nospike.NachO3.wet <- targetSpike.mean(data = nospike.NachO3.wet, formula = Sample.Name ~ Compound,
                                       peak.area = 'Area', na.rm = TRUE)

nospike.GAK5b.wet <- targetSpike.mean(data = nospike.GAK5b.wet, formula = Sample.Name ~ Compound,
                                      peak.area = 'Area', na.rm = TRUE)

nospike.VKB.dry <- targetSpike.mean(data = nospike.VKB.dry, formula = Sample.Name ~ Compound,
                                    peak.area = 'Area', na.rm = TRUE)

nospike.NachO3.dry <- targetSpike.mean(data = nospike.NachO3.dry, formula = Sample.Name ~ Compound,
                                       peak.area = 'Area', na.rm = TRUE)

nospike.GAK5b.dry <- targetSpike.mean(data = nospike.GAK5b.dry, formula = Sample.Name ~ Compound,
                                      peak.area = 'Area', na.rm = TRUE)


### remove values where spiked area is not minimum 1.5 times higher than the unspike area (check for correct spike level)
error.margin <- 1.5 ## is the spike level correct? set ratio here

take.spike.VKB.wet <- which(spike.VKB.wet > error.margin * nospike.VKB.wet)
take.spike.NachO3.wet <- which(spike.NachO3.wet > error.margin * nospike.NachO3.wet)
take.spike.GAK5b.wet <- which(spike.GAK5b.wet > error.margin * nospike.GAK5b.wet)

take.spike.VKB.dry <- which(spike.VKB.dry > error.margin * nospike.VKB.dry)
take.spike.NachO3.dry <- which(spike.NachO3.dry > error.margin * nospike.NachO3.dry)
take.spike.GAK5b.dry <- which(spike.GAK5b.dry > error.margin * nospike.GAK5b.dry)

good.spikes <- c(length(take.spike.VKB.wet), length(take.spike.NachO3.wet), length(take.spike.GAK5b.wet),
                 length(take.spike.VKB.dry), length(take.spike.NachO3.dry), length(take.spike.GAK5b.dry))

#### difference of spiked area and unspiked area
spike.VKB.wet <- spike.VKB.wet - nospike.VKB.wet
spike.NachO3.wet <- spike.NachO3.wet - nospike.NachO3.wet
spike.GAK5b.wet <- spike.GAK5b.wet - nospike.GAK5b.wet

spike.VKB.dry <- spike.VKB.dry - nospike.VKB.dry
spike.NachO3.dry <- spike.NachO3.dry - nospike.NachO3.dry
spike.GAK5b.dry <- spike.GAK5b.dry - nospike.GAK5b.dry


### Matrix Effect

### remove values where the spiked area - unspiked area is negative
spike.VKB.wet[spike.VKB.wet < 0] <- NA
spike.NachO3.wet[spike.NachO3.wet < 0] <- NA
spike.GAK5b.wet[spike.GAK5b.wet < 0] <- NA

spike.VKB.dry[spike.VKB.dry < 0] <- NA
spike.NachO3.dry[spike.NachO3.dry < 0] <- NA
spike.GAK5b.dry[spike.GAK5b.dry < 0] <- NA


#### calculate spiked area in comparison to a standard
recovery.VKB.wet <- spike.VKB.wet / (as.numeric(NPstd) )
recovery.NachO3.wet <- spike.NachO3.wet / (as.numeric(NPstd) )
recovery.GAK5b.wet <- spike.GAK5b.wet / (as.numeric(NPstd))

recovery.VKB.dry <- spike.VKB.dry / (as.numeric(NPstd) )
recovery.NachO3.dry <- spike.NachO3.dry / (as.numeric(NPstd) )
recovery.GAK5b.dry <- spike.GAK5b.dry / (as.numeric(NPstd) )


### combine data
spike.amts <- data.frame(areas = c(spike.VKB.wet[take.spike.VKB.wet], 
                                   spike.NachO3.wet[take.spike.NachO3.wet],
                                   spike.GAK5b.wet[take.spike.GAK5b.wet],
                                   spike.VKB.dry[take.spike.VKB.dry],
                                   spike.NachO3.dry[take.spike.NachO3.dry],
                                   spike.GAK5b.dry[take.spike.GAK5b.dry]
                                   ),
                         sampletype = c(rep("VKB", good.spikes[1]),
                                        rep("NachO3", good.spikes[2]),
                                        rep("GAK5b", good.spikes[3]),
                                        rep("VKB", good.spikes[4]),
                                        rep("NachO3", good.spikes[5]),
                                        rep("GAK5b", good.spikes[6])
                                        ),
                         buechi = c(rep("wet", sum(good.spikes[1:3])),
                                    rep("dry", sum(good.spikes[4:6]))
                                    ),
                         recovery = c(recovery.VKB.wet[take.spike.VKB.wet],
                                      recovery.NachO3.wet[take.spike.NachO3.wet],
                                      recovery.GAK5b.wet[take.spike.GAK5b.wet],
                                      recovery.VKB.dry[take.spike.VKB.dry],
                                      recovery.NachO3.dry[take.spike.NachO3.dry],
                                      recovery.GAK5b.dry[take.spike.GAK5b.dry])
                         )


#### look at data
bwplot(log10(areas) ~ sampletype | buechi, data = spike.amts)
# not terribly informative since all targets have their own area



bwplot(recovery ~ sampletype | buechi, data = spike.amts
       #,ylim = c(-1,1)
       )

### better bwplot

## how to order the bwplot?
levels(spike.amts$sampletype)
spike.amts$sampletype <- factor(spike.amts$sampletype, levels = c("VKB", "NachO3", "GAK5b"), order = TRUE)

### depending on which standard is used, the output of this analysis is either the matrix effect
### or the absolute recovery

# png("MEinBuechi_CC9_Test2_SMmethod.png", width = 3300, height = 1500, res = 300)
# bwplot(log10(recovery) ~ sampletype | buechi, data = spike.amts, #ylim = c(-2,5),
#        par.strip.text = list(cex=1),
#        strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
#                             bg = "lightgray"),
#        ylab = list(label = "log10( Target Spiked Area / Target Standard Area )", cex = 1),
#        ylim = c(-7,7),
#        xlab = list(label = NULL, cex = 1),
#        scales = list(x = list(cex = 1, labels = c(expression("Before CAS"),
#                                                   expression('After O'[3]), expression("After GAK5b"))),
#                      y = list(cex = 1)),
#        #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
#        panel = function(x = x, y = y,...){
#          panel.abline(h=c(log10(1))
#                       #, col = "red", lty = 3, lwd = 3
#          )
#          panel.bwplot(x = x,y = y, ..., fill = rev(c(blues9[3:6])),
#                       # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
#                       par.settings = list(y = grid::unit(110, "mm")),
#                       panel.axis = list(outside = TRUE))
#           nn <- table(x)
#           panel.text(paste("n=", nn), x = seq_along(nn),
#                      y = current.panel.limits()$y[1], pos = 3)
#        }
# )
# dev.off()

# png("AbsRecoveryinBuechi_CC9_Test2_SMmethod.png", width = 3300, height = 1500, res = 300)
# bwplot(recovery*100 ~ sampletype | buechi, data = spike.amts, #ylim = c(-2,5),
#        par.strip.text = list(cex=1),
#        strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
#                             bg = "lightgray"),
#        ylab = list(label = "log10( Target Spiked Area / Target Standard Area )", cex = 1),
#        ylim = c(-10,1000),
#        xlab = list(label = NULL, cex = 1),
#        scales = list(x = list(cex = 1, labels = c(expression("Before CAS"),
#                                                   expression('After O'[3]), expression("After GAK5b"))),
#                      y = list(cex = 1)),
#        #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
#        panel = function(x = x, y = y,...){
#          panel.abline(h=c(log10(1))
#                       #, col = "red", lty = 3, lwd = 3
#          )
#          panel.bwplot(x = x,y = y, ..., fill = rev(c(blues9[3:6])),
#                       # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
#                       par.settings = list(y = grid::unit(110, "mm")),
#                       panel.axis = list(outside = TRUE))
#          nn <- table(x)
#          panel.text(paste("n=", nn), x = seq_along(nn),
#                     y = current.panel.limits()$y[1], pos = 3)
#        }
# )
# dev.off()


# png("MEinBuechi_vsNPBuechi_Test2_SMmethod.png", width = 3300, height = 1500, res = 300)
# bwplot(log10(recovery) ~ sampletype | buechi, data = spike.amts, #ylim = c(-2,5),
#        par.strip.text = list(cex=1),
#        strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
#                             bg = "lightgray"),
#        ylab = list(label = "log10( Target Spiked Area / Target Standard Area )", cex = 1),
#        ylim = c(-7,7),
#        xlab = list(label = NULL, cex = 1),
#        scales = list(x = list(cex = 1, labels = c(expression("Before CAS"),
#                                                   expression('After O'[3]), expression("After GAK5b"))),
#                      y = list(cex = 1)),
#        #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
#        panel = function(x= x, y= y, ...){
#          panel.abline(h=c(log10(1))
#                       #, col = "red", lty = 3, lwd = 3
#                       )
#          panel.bwplot(x = x, y = y, ..., fill = rev(c(blues9[3:6])),
#                       # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
#                       par.settings = list(y = grid::unit(110, "mm")),
#                       panel.axis = list(outside = TRUE))
#                    nn <- table(x)
#                    panel.text(paste("n=", nn), x = seq_along(nn),
#                               y = current.panel.limits()$y[1], pos = 3)
#        }
# )
# dev.off()

# png("MEinBuechi_spikedTargetvsNPBuechi.png", width = 3300, height = 1500, res = 300)
# bwplot(recovery ~ sampletype | buechi, data = spike.amts, #ylim = c(-2,5),
#        par.strip.text = list(cex=1),
#        strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
#                             bg = "lightgray"),
#        ylab = list(label = "log10( Target Spiked Area / Target Standard Area )", cex = 1),
#         ylim = c(-0.05,1),
#        xlab = list(label = NULL, cex = 1),
#        scales = list(x = list(cex = 1, labels = c(expression("Before CAS"),
#                                                   expression('After O'[3]), expression("After GAK5b"))),
#                      y = list(cex = 1)),
#        #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
#        panel = function(...,box.ratio){
#          panel.abline(h=c(log10(1))
#                       #, col = "red", lty = 3, lwd = 3
#          )
#          panel.bwplot(..., fill = rev(c(blues9[3:6])),
#                       # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
#                       par.settings = list(y = grid::unit(110, "mm")),
#                       panel.axis = list(outside = TRUE))
#        }
# )
# dev.off()


upper.rec <- 120
lower.rec <- 80
cast(sampletype ~ buechi, data = spike.amts, value = 'recovery', 
     function(x) length(which(x < upper.rec & x > lower.rec)))


# 
# densityplot(~ recovery | buechi, groups = sampletype, data = spike.amts 
#        #ylim = c(1000,-1000)
# )
# 
# 
# histogram(~ recovery | buechi, groups = sampletype, data = spike.amts, 
#           type = "p",
#           panel = function(...)panel.superpose(...,panel.groups = panel.histogram,
#                                                col = c("cyan", "magenta", "yellow"), alpha = 0.4),
#           auto.key = list(columns = 3, rectangles = FALSE,
#                           col = c("cyan", "magenta", "yellow3"))
#           )


#### comparison to RP
### to compare get RPLC ME values
RP.MEs <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\RPLC_MEs.csv",
                   header = TRUE)

spike.amts <- data.frame(areas = c(spike.VKB.wet[take.spike.VKB.wet], 
                                   spike.NachO3.wet[take.spike.NachO3.wet],
                                   spike.GAK5b.wet[take.spike.GAK5b.wet],
                                   spike.VKB.dry[take.spike.VKB.dry],
                                   spike.NachO3.dry[take.spike.NachO3.dry],
                                   spike.GAK5b.dry[take.spike.GAK5b.dry],
                                   rep("NA", nrow(RP.MEs))
                                   ),
sampletype = c(rep("VKB", good.spikes[1]),
               rep("NachO3", good.spikes[2]),
               rep("GAK5b", good.spikes[3]),
               rep("VKB", good.spikes[4]),
               rep("NachO3", good.spikes[5]),
               rep("GAK5b", good.spikes[6]),
               as.character(RP.MEs$sampletype)),
buechi = c(rep("wet", sum(good.spikes[1:3])),
           rep("dry", sum(good.spikes[4:6])),
           as.character(RP.MEs$ME_method)),
recovery = c(recovery.VKB.wet[take.spike.VKB.wet],
             recovery.NachO3.wet[take.spike.NachO3.wet],
             recovery.GAK5b.wet[take.spike.GAK5b.wet],
             recovery.VKB.dry[take.spike.VKB.dry],
             recovery.NachO3.dry[take.spike.NachO3.dry],
             recovery.GAK5b.dry[take.spike.GAK5b.dry],
             RP.MEs$value)
)

png("MEinBuechi_vsNPBuechi_Test2_SMmethod_wRPLCMEs.png", width = 3300, height = 1500, res = 300)
bwplot(log10(recovery) ~ sampletype | buechi, data = spike.amts, #ylim = c(-2,5),
       par.strip.text = list(cex=1),
       # strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
       #                      bg = "lightgray"),
       ylab = list(label = "log10( Target Spiked Area / Target Standard Area )", cex = 1),
       ylim = c(-7,7),
       xlab = list(label = NULL, cex = 1),
       # scales = list(x = list(cex = 1, labels = c(expression("Before CAS"),
       #                                            expression('After O'[3]), expression("After GAK5b"))),
       #               y = list(cex = 1)),
       # #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
       panel = function(x= x, y= y, ...){
         panel.abline(h=c(log10(1))
                      #, col = "red", lty = 3, lwd = 3
         )
         panel.bwplot(x = x, y = y, ..., fill = rev(c(blues9[3:6])),
                      # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
                      par.settings = list(y = grid::unit(110, "mm")),
                      panel.axis = list(outside = TRUE))
         nn <- table(x)
         panel.text(paste("n=", nn), x = seq_along(nn),
                    y = current.panel.limits()$y[1], pos = 3)
       }
)
dev.off()
