

library(lattice)
library(reshape)
library(gridExtra)
library(plotrix)

### Calculating Matrix Effect based on ISTD Area

### Where to save everything?
setwd("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\1st test - May 2017\\Targets")

## where is tracefinder data file?
sample.data <- read.csv("HILIC_pos_everything_samples_BEHAmide_targets.csv", header = TRUE)

sample.data <- subset(sample.data, sample.data$Peak.Label == "T1")  ## selecting only monoisotopic peaks

istd.data <- subset(sample.data, sample.data$Type == "Internal Standard")

## for each relevant sample type, define a vector
blanks <- grep(pattern = "ACN-", x = istd.data$Sample.Name)
STDA <- grep(pattern = "STDA_", x = istd.data$Sample.Name)
STDB <- grep(pattern = "STDB_", x = istd.data$Sample.Name)
STDC <- grep(pattern = "STDC_", x = istd.data$Sample.Name)
VKB <- grep(pattern = "VKB_", x = istd.data$Sample.Name)
NachO3 <- grep(pattern = "NachO3_", x = istd.data$Sample.Name)
GAK5b <- grep(pattern = "GAK5b_", x = istd.data$Sample.Name)
Dry <- grep(pattern = "Dry_", x = istd.data$Sample.Name)
Wet <- grep(pattern = "Wet_", x = istd.data$Sample.Name)




#boxplot(as.numeric(as.character(Area)) ~ Sample.Name, data = istd.data)

## ISTD areas in different sample types
istds.wet <- istd.data[c(intersect(VKB, Wet), intersect(NachO3, Wet), intersect(GAK5b, Wet)),]
istds.wet$matrix <- c(rep("VKB", length(intersect(VKB, Wet))), 
                      rep("NachO3", length(intersect(NachO3, Wet))),
                      rep("GAK5b", length(intersect(GAK5b, Wet))))
istds.wet$enrich <- "wet"
boxplot(log10(as.numeric(as.character(Area))) ~ matrix, data = istds.wet)


istds.dry <- istd.data[c(intersect(VKB, Dry), intersect(NachO3, Dry), intersect(GAK5b, Dry)),]
istds.dry$matrix <- c(rep("VKB", length(intersect(VKB, Dry))), 
                      rep("NachO3", length(intersect(NachO3, Dry))),
                      rep("GAK5b", length(intersect(GAK5b, Dry))))
istds.dry$enrich <- "dry"
boxplot(log10(as.numeric(as.character(Area))) ~ matrix, data = istds.dry)

## combine all data and plot
istds.samples <- rbind(istds.wet, istds.dry)
bwplot(log10(as.numeric(as.character(Area))) ~ matrix | enrich, data = istds.samples,
       par.strip.text = list(cex=1.3, label = c("Dry", "Wet")),
       ylab = list(label = "log10(Internal Standard Area)", cex = 1.3),
       xlab = list(label = "Sample Type", cex = 1.3),
       scales = list(x = list(cex = 1.2),
                     y = list(cex = 1.2)),
       par.settings = list(box.rectangle = list(fill = c("lightblue")))
)



## ISTD areas in seletected standard
# istds.std <- istd.data[c(STDC),]
# istds.std$matrix <- "STD"
# istds.std$enrich <- "NP"

# # here for the Buechi NP standard
std.data <- read.csv("Q:\\Abteilungsprojekte\\uchem\\Projekte\\Aktuelle Projekte\\Advanced_Treatment_wastewater\\Non-target BAFU\\Sampling\\HILIC Method Development\\2nd test - July 2017\\TFexports\\BEHAmide_NPBuechiStds_targets.csv",
                     header = TRUE)
std.data <- subset(std.data, std.data$Peak.Label == "T1" & std.data$Type == "Internal Standard")
colnames(std.data)[1] <- "Compound"
levels(std.data$Sample.Name)
STDs <- grep(pattern = "CC", x = std.data$Sample.Name)
istds.std <- std.data[STDs,]

### Test 2 has more target compounds, subset to keep the same set
istds.std <- istds.samples[which(c(istds.std$Compound %in% istds.wet$Compound)==TRUE),]
#istds.std <- as.vector(istds.std[,c("Area")])

istds.std$matrix <- "STD"
istds.std$enrich <- "Buechi"

istds.std <- istds.std[,c(paste(colnames(istds.wet)))]


## combine all data and plot
istds.all <- rbind(istds.std, istds.wet, istds.dry)

bwplot(log10(as.numeric(as.character(Area))) ~ matrix | enrich, data = istds.all,
       par.strip.text = list(cex=1.3, label = c("Dry", "Wet")),
       ylab = list(label = "log10(Internal Standard Area)", cex = 1.3),
       xlab = list(label = "Sample Type", cex = 1.3),
       scales = list(x = list(cex = 1.2),
                     y = list(cex = 1.2)),
       par.settings = list(box.rectangle = list(fill = c("lightblue")))
)


### standardizing ISTD areas
std.means <- as.data.frame(cast(istds.std, Compound ~ . , function(x) mean(na.omit(as.numeric(as.character(x)))), value = "Area"))
dry.means <- as.data.frame(cast(istds.dry, Compound ~ matrix, function(x) mean(na.omit(as.numeric(as.character(x)))), value = "Area"))
wet.means <- as.data.frame(cast(istds.wet, Compound ~ matrix, function(x) mean(na.omit(as.numeric(as.character(x)))), value = "Area"))

dry.means <- apply(dry.means[,2:4], 2, function(x) x/std.means[,2])
rownames(dry.means) <- std.means[,1]
dry.means <- melt(dry.means)
dry.means$type <- "dry"

wet.means <- apply(wet.means[,2:4], 2, function(x) x/std.means[,2])
rownames(wet.means) <- std.means[,1]
wet.means <- melt(wet.means)
wet.means$type <- "wet"

bwplot(value ~ X2, data = wet.means, ylim = c(-5,20))
bwplot(value ~ X2, data = dry.means, ylim = c(-5,20))

all.means <- rbind(dry.means, wet.means)

## how to order the bwplot?
levels(all.means$X2)
all.means$X2 <- factor(all.means$X2, levels = c("VKB", "NachO3", "GAK5b"), order = TRUE)

## remove those values where ISTD was not detected
all.means <- na.omit(all.means)


lattice.options(layout.heights = list(bottom.padding = list(x = 1)))

png("MatrixEffect_log10ISTD_NPBuechi.png", width = 3300, height = 1500, res = 300)
bwplot(log10(value) ~ X2 | type, data = all.means, #ylim = c(-2,5), 
       par.strip.text = list(cex=1),
       strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
                            bg = "lightgray"),
       ylab = list(label = "log10( ISTD Area in Matrix / ISTD Area in Standard )", cex = 1),
       ylim = c(-7,7),
       xlab = list(label = NULL, cex = 1),
       scales = list(x = list(cex = 1, labels = c(expression("Before CAS"), 
                                                  expression('After O'[3]), expression("After GAK5b"))),
                     y = list(cex = 1)),
       #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
       panel = function(x = x, y = y,...,box.ratio){
         panel.abline(h=log10(1)
                      #, col = "red", lty = 3, lwd = 3
                      )
         panel.bwplot(x=x, y=y,..., fill = rev(c(blues9[3:6])),
                      # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
                      par.settings = list(y = grid::unit(110, "mm")),
                      panel.axis = list(outside = TRUE))
          nn <- table(x)
          panel.text(paste("n=", nn), x = seq_along(nn),
                     y = current.panel.limits()$y[1], pos = 3)
       }
)
dev.off()

png("MatrixEffect2_log10ISTD.png", width = 3300, height = 1500, res = 300)
bwplot(log10(1/value) ~ X2 | type, data = all.means, #ylim = c(-2,5), 
       par.strip.text = list(cex=1),
       strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
                            bg = "lightgray"),
       ylab = list(label = "log10( ISTD Area in Standard / ISTD Area in Matrix )", cex = 1),
       ylim = c(-6,6),
       xlab = list(label = NULL, cex = 1),
       scales = list(x = list(cex = 1, labels = c(expression("Before CAS"), 
                                                  expression('After O'[3]), expression("After GAK5b"))),
                     y = list(cex = 1)),
       #trellis.par.set((layout.heights = list(bottom.padding = list(x = 3)))),
       panel = function(...,box.ratio){
         panel.abline(h=log10(1), col = "red", lty = 3, lwd = 3)
         panel.bwplot(..., fill = rev(c(blues9[3:6])),
                      # par.settings = list(layout.heights = list(bottom.padding = list(x=3))),
                      par.settings = list(y = grid::unit(110, "mm")),
                      panel.axis = list(outside = TRUE))
       }
)
dev.off()

# png("MatrixEffect_ISTD_Zoom.png", width = 1100, height = 500)
# bwplot(log10(value) ~ X2 | type, data = all.means, ylim = c(-5,5), 
#        par.strip.text = list(cex=1.3),
#        strip = strip.custom(factor.levels = c("Dry Buechi Enrichment", "Wet Buechi Enrichment"),
#                             bg = "lightgray"),
#        ylab = list(label = "ISTD Area in Standard / ISTD Area in Matrix", cex = 1.3),
#        xlab = list(label = NULL, cex = 1.3),
#        scales = list(x = list(cex = 1.2),
#                      y = list(cex = 1.2)),
#        panel = function(...,box.ratio){
#          panel.abline(h=log10(1))
#          panel.bwplot(..., fill = c(blues9[3:6]))}
# )
# dev.off()

