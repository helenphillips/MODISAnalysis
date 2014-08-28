# Analysis of EVI and NDVI data for a single pixel

#setwd()

library(lme4)
library(MuMIn) ## for RsquareGLMM function
library(car) ## for Anova function

## library(roquefort) - for plotting. Not publically available
source("/Users/hp1111/PhD/AdriannasFunctions/model_plot.R")


file <- list.files(path = "MOD13Q1_IndividualPixel", pattern = "^MODIS_Data") ## The file that begins MODIS_Data is the original data with MODIS appended
file
# file2 ## Need to ensure there is only 1 file, or that we know which file we want to load

PREDICTS <- read.csv(file.path("MOD13Q1_IndividualPixel", file), header = TRUE)
names(PREDICTS)
nrow(PREDICTS)
str(PREDICTS)



###### Appending additional summary information ### 
file2 <- list.files(path = "MOD13Q1_IndividualPixel", pattern = "^MODIS_Summary")
file2
file2 <- read.csv(file.path("MOD13Q1_IndividualPixel", file2), header = TRUE)

names(file2)
nrow(file2) ## Twice the number of rows as downloaded subsets - EVI and NDVI are on separate rows 

evi <- file2[file2$data.band == "250m_16_days_EVI",]
ndvi <- file2[file2$data.band == "250m_16_days_NDVI",]

## Add all additional columns to the dataset (SD, minimum, maximum etc.), rename, and remove unwanted rows
PREDICTS <- merge(PREDICTS, evi, by.x = "ID", by.y = "ID", all.x = TRUE)
PREDICTS$lat.y <- NULL
PREDICTS$long.y <- NULL
PREDICTS$end.date.y <- NULL
PREDICTS$data.band <- NULL
newnames <- names(PREDICTS)[28:length(names(PREDICTS))]
newnames <- paste("EVI.", newnames, sep = "")
names(PREDICTS)[28:length(names(PREDICTS))] <- newnames

PREDICTS <- merge(PREDICTS, ndvi, by.x = "ID", by.y = "ID", all.x = TRUE)
PREDICTS$lat <- NULL
PREDICTS$long <- NULL
PREDICTS$end.date <- NULL
PREDICTS$start.date.x <- NULL
PREDICTS$data.band <- NULL
newnames <- names(PREDICTS)[35:length(names(PREDICTS))]
newnames <- paste("NDVI.", newnames, sep = "")
names(PREDICTS)[35:length(names(PREDICTS))] <- newnames

summary(PREDICTS$EVI.mean.band)
summary(PREDICTS$NDVI.mean.band)

### To check that everything is correct, PREDICTS$X250m_16_days_EVI_pixel1 and PREDICTS$EVI.mean.band should the same numbers (both the mean)
sum(PREDICTS$X250m_16_days_EVI_pixel1 == PREDICTS$EVI.mean.band, na.rm = TRUE)
## same with NDVI
sum(PREDICTS$X250m_16_days_NDVI_pixel1 == PREDICTS$NDVI.mean.band, na.rm = TRUE)
## Good.

## Data exploration
summary(PREDICTS$EVI.mean.band)
# NAs present, so let's remove them
PREDICTS <-PREDICTS[complete.cases(PREDICTS$EVI.mean.band),]

summary(PREDICTS$NDVI.mean.band)
## No NAs present


summary(PREDICTS$EVI.band.yield)
summary(PREDICTS$NDVI.band.yield)
## another 8 data points missing
PREDICTS <-PREDICTS[complete.cases(PREDICTS$EVI.band.yield),]

par(mfrow = c(2,2))
hist(PREDICTS$EVI.mean.band, xlab = "EVI", main = "Histogram of mean EVI")
hist(PREDICTS$NDVI.mean.band, xlab = "NDVI", main = "Histogram of mean NDVI")
hist(PREDICTS$EVI.band.yield, xlab = "EVI", main = "Histogram of EVI yield")
hist(PREDICTS$NDVI.band.yield, xlab = "NDVI", main = "Histogram of NDVI yield") 


par(mfrow = c(2,2))
plot(PREDICTS$EVI.mean.band, col = PREDICTS$Source_ID, ylab = "EVI", main = "Scatterplot of EVI by source")
plot(PREDICTS$NDVI.mean.band, col = PREDICTS$Source_ID, ylab = "NDVI", main = "Scatterplot of NDVI by source")
plot(PREDICTS$EVI.band.yield, col = PREDICTS$Source_ID, ylab = "EVI", main = "Scatterplot of EVI by source")
plot(PREDICTS$NDVI.band.yield, col = PREDICTS$Source_ID, ylab = "NDVI", main = "Scatterplot of NDVI by source")
## Expect to see the horizontal bands in colour, but not nessecarily that each source was consistently higher or lower (bar the few in the middle). Homogeneity within sources.
## A random effects structure will deal with this homogeneity
## As technically the data points are nested within each source, it will also deal with that

## How many publications are in the dataset
length(unique(PREDICTS$Source_ID))
length(unique(paste(PREDICTS$Source_ID, PREDICTS$Study_number))) ## Same number of studies and sources, so random effects will only need the one

plot(PREDICTS$EVI.mean.band, PREDICTS$NDVI.mean.band, xlab = "EVI", ylab = "NDVI")
abline(a=0, b=1)
## As expected, strong relationship between the two (EVI is similar to NDVI but with a small correction factor), so both will be analysed independently

plot(PREDICTS$Higher_taxon, PREDICTS$Species_richness, ylab = "Species Richness")
## Clear distinction between higher taxon of species richness. Again, higher taxon can be a random or fixed effect to deal with that

table(PREDICTS$Higher_taxon)
## Too few Arachnida, so remove them from the dataset
PREDICTS <- droplevels(PREDICTS[PREDICTS$Higher_taxon != "Arachnida",])

## Create folder for figures ##
if(file.exists("Figures") == FALSE)
	dir.create("Figures")
	
dir.create(file.path("Figures", "EVISinglePixel"))
dir.create(file.path("Figures", "NDVISinglePixel"))
dir.create(file.path("Figures", "SupplementaryFigures"))

####### Univariate models (+higher taxon)######



## EVI mean
cols_higher_taxon <- c("#000000", "#9932CC", "#FFD700", "#006400")
VIs.e1 <- glmer(Species_richness ~ EVI.mean.band * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.e1b <- glmer(Species_richness ~ EVI.mean.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e1, VIs.e1b) ## p value < 0.05 need interaction

summary(VIs.e1)
model_plots(VIs.e1)

pdf(file.path("Figures", "EVISinglePixel", "Mean.pdf"))
PlotContEffects(VIs.e1, data = PREDICTS,effects = "EVI.mean.band", xlab = "Mean EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()
r.squaredGLMM(VIs.e1) # 
      # R2m       R2c 
# 0.5871584 0.8097819


## EVI max

VIs.e2 <- glmer(Species_richness ~  Higher_taxon * EVI.max.band + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.e2b <- glmer(Species_richness ~ EVI.max.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e2, VIs.e2b) ## need interaction

pdf(file.path("Figures", "EVISinglePixel", "Max.pdf"))
PlotContEffects(VIs.e2, data = PREDICTS,effects = "EVI.max.band", xlab = "Maximum EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()
summary(VIs.e2)
r.squaredGLMM(VIs.e2)
      # R2m       R2c 
# 0.6022174 0.7996292 


### EVI yield

VIs.e3 <- glmer(Species_richness ~ EVI.band.yield * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.e3b <- glmer(Species_richness ~ EVI.band.yield + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e3, VIs.e3b) ## not significant ( p > 0.05)

VIs.e3c <- glmer(Species_richness ~ EVI.band.yield + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e3b, VIs.e3c) ## Need higher taxa
VIs.e3d <- glmer(Species_richness ~ Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e3b, VIs.e3d) ## do need yield
# Stick with the additive model

summary(VIs.e3b)

pdf(file.path("Figures", "EVISinglePixel", "Yield.pdf"))
PlotContEffects(VIs.e3b, data = PREDICTS,effects = "EVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "EVI yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
dev.off()
r.squaredGLMM(VIs.e2b)
      # R2m       R2c 
# 0.5774434 0.8018797 





### NDVI
#Mean

VIs.n1 <- glmer(Species_richness ~ NDVI.mean.band * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.n1b <- glmer(Species_richness ~ NDVI.mean.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n1, VIs.n1b) ## Significantly different

summary(VIs.n1)
r.squaredGLMM(VIs.n1)
      # R2m       R2c 
# 0.5696208 0.8253525

pdf(file.path("Figures", "NDVISinglePixel", "Mean.pdf"))
PlotContEffects(VIs.n1, data = PREDICTS,effects = "NDVI.mean.band", xlab = "mean NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


# Max


VIs.n2 <- glmer(Species_richness ~ NDVI.max.band * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.n2b <- glmer(Species_richness ~ NDVI.max.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n2, VIs.n2b) ## significantly different

summary(VIs.n2)
r.squaredGLMM(VIs.n2)
      # R2m       R2c 
# 0.5954573 0.8041141 

pdf(file.path("Figures", "NDVISinglePixel", "Max.pdf"))
PlotContEffects(VIs.n2, data = PREDICTS,effects = "NDVI.max.band", xlab = "Maximum NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()

# Yield

VIs.n3 <- glmer(Species_richness ~ NDVI.band.yield * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.n3b <- glmer(Species_richness ~ NDVI.band.yield + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n3, VIs.n3b) ## don't need interaction

VIs.n3c <- glmer(Species_richness ~ NDVI.band.yield + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
anova(VIs.n3b, VIs.n3c) # can't remove
VIs.n3d <- glmer(Species_richness ~ Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n3b, VIs.n3d) # can't remove
# Stick with additive model

summary(VIs.n3b)
r.squaredGLMM(VIs.n3b)

pdf(file.path("Figures", "NDVISinglePixel", "Yield.pdf"))
PlotContEffects(VIs.n3b, data = PREDICTS,effects = "NDVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "NDVI yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
dev.off()




## Creating supplementary figures
pdf(file.path("Figures", "SupplementaryFigures", "SinglePixels.pdf"))
par(mfrow = c(3, 2))
PlotContEffects(VIs.e1, data = PREDICTS,effects = "EVI.mean.band", xlab = "Mean EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
text(0.56, 70, "(a)", cex = 1.5)

PlotContEffects(VIs.n1, data = PREDICTS,effects = "NDVI.mean.band", xlab = "Mean NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
text(0.83, 80, "(b)", cex = 1.5)


PlotContEffects(VIs.e2, data = PREDICTS,effects = "EVI.max.band", xlab = "Maximum EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
text(0.83, 80, "(c)", cex = 1.5)


PlotContEffects(VIs.n2, data = PREDICTS,effects = "NDVI.max.band", xlab = "Maximum NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
text(0.9, 70, "(d)", cex = 1.5)


PlotContEffects(VIs.e3b, data = PREDICTS,effects = "EVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "EVI Yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
text(0.28, 18, "(e)", cex = 1.5)


PlotContEffects(VIs.n3b, data = PREDICTS,effects = "NDVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "NDVI Yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
text(0.41, 21, "(f)", cex = 1.5)

dev.off()