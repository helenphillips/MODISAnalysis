# Analysis of EVI and NDVI data for a single pixel

#setwd()

library(lme4)
library(MuMIn) ## for RsquareGLMM function
library(car) ## for Anova function

## library(roquefort) - for plotting. Not publically available
model_plot <- function(mod.for.plot){
	par(mfrow = c(1,3))
	qqnorm(resid(mod.for.plot))
	qqline(resid(mod.for.plot), col = 2)
  	plot(fitted(mod.for.plot), resid(mod.for.plot),xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs fitted")
  	abline(h=0, lty=2)
  	lines(smooth.spline(fitted(mod.for.plot), resid(mod.for.plot)), col = "red")
	hist(resid(mod.for.plot))}


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
PREDICTS.e1 <- glmer(Species_richness ~ EVI.mean.band * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PREDICTS.e1b <- glmer(Species_richness ~ EVI.mean.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.e1, PREDICTS.e1b) ## p value < 0.05 need interaction

summary(PREDICTS.e1)
model_plot(PREDICTS.e1)

pdf(file.path("Figures", "EVISinglePixel", "Mean.pdf"))
PlotContEffects(PREDICTS.e1, data = PREDICTS,effects = "EVI.mean.band", xlab = "Mean EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()
r.squaredGLMM(PREDICTS.e1) # 
      # R2m       R2c 
# 0.5871584 0.8097819


## EVI max

PREDICTS.e2 <- glmer(Species_richness ~  Higher_taxon * EVI.max.band + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PREDICTS.e2b <- glmer(Species_richness ~ EVI.max.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.e2, PREDICTS.e2b) ## need interaction

pdf(file.path("Figures", "EVISinglePixel", "Max.pdf"))
PlotContEffects(PREDICTS.e2, data = PREDICTS,effects = "EVI.max.band", xlab = "Maximum EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()
summary(PREDICTS.e2)
model_plot(PREDICTS.e2)
r.squaredGLMM(PREDICTS.e2)
      # R2m       R2c 
# 0.6022174 0.7996292 


### EVI yield

PREDICTS.e3 <- glmer(Species_richness ~ EVI.band.yield * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PREDICTS.e3b <- glmer(Species_richness ~ EVI.band.yield + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.e3, PREDICTS.e3b) ## not significant ( p > 0.05)

PREDICTS.e3c <- glmer(Species_richness ~ EVI.band.yield + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.e3b, PREDICTS.e3c) ## Need higher taxa
PREDICTS.e3d <- glmer(Species_richness ~ Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.e3b, PREDICTS.e3d) ## do need yield
# Stick with the additive model

summary(PREDICTS.e3b)
model_plot(PREDICTS.e3b)

pdf(file.path("Figures", "EVISinglePixel", "Yield.pdf"))
PlotContEffects(PREDICTS.e3b, data = PREDICTS,effects = "EVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "EVI yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
dev.off()
r.squaredGLMM(PREDICTS.e2b)
      # R2m       R2c 
# 0.5774434 0.8018797 





### NDVI
#Mean

PREDICTS.n1 <- glmer(Species_richness ~ NDVI.mean.band * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PREDICTS.n1b <- glmer(Species_richness ~ NDVI.mean.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.n1, PREDICTS.n1b) ## Significantly different

summary(PREDICTS.n1)
model_plot(PREDICTS.n1)

r.squaredGLMM(PREDICTS.n1)
      # R2m       R2c 
# 0.5696208 0.8253525

pdf(file.path("Figures", "NDVISinglePixel", "Mean.pdf"))
PlotContEffects(PREDICTS.n1, data = PREDICTS,effects = "NDVI.mean.band", xlab = "mean NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


# Max


PREDICTS.n2 <- glmer(Species_richness ~ NDVI.max.band * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PREDICTS.n2b <- glmer(Species_richness ~ NDVI.max.band + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.n2, PREDICTS.n2b) ## significantly different

summary(PREDICTS.n2)
model_plot(PREDICTS.n2)

r.squaredGLMM(PREDICTS.n2)
      # R2m       R2c 
# 0.5954573 0.8041141 

pdf(file.path("Figures", "NDVISinglePixel", "Max.pdf"))
PlotContEffects(PREDICTS.n2, data = PREDICTS,effects = "NDVI.max.band", xlab = "Maximum NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()

# Yield

PREDICTS.n3 <- glmer(Species_richness ~ NDVI.band.yield * Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PREDICTS.n3b <- glmer(Species_richness ~ NDVI.band.yield + Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.n3, PREDICTS.n3b) ## don't need interaction

PREDICTS.n3c <- glmer(Species_richness ~ NDVI.band.yield + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
anova(PREDICTS.n3b, PREDICTS.n3c) # can't remove
PREDICTS.n3d <- glmer(Species_richness ~ Higher_taxon + (1|SS), data = PREDICTS, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PREDICTS.n3b, PREDICTS.n3d) # can't remove
# Stick with additive model

summary(PREDICTS.n3b)
model_plot(PREDICTS.n3b)

r.squaredGLMM(PREDICTS.n3b)

pdf(file.path("Figures", "NDVISinglePixel", "Yield.pdf"))
PlotContEffects(PREDICTS.n3b, data = PREDICTS,effects = "NDVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "NDVI yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
dev.off()




## Creating supplementary figures
pdf(file.path("Figures", "SupplementaryFigures", "SinglePixels.pdf"), height = 11.7, width = 8.3)
par(mfrow = c(3, 2))
PlotContEffects(PREDICTS.e1, data = PREDICTS,effects = "EVI.mean.band", xlab = "Mean EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(a) ", side=3, line=-1.5, adj=1.0, cex=1.2)

PlotContEffects(PREDICTS.n1, data = PREDICTS,effects = "NDVI.mean.band", xlab = "Mean NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(b) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PREDICTS.e2, data = PREDICTS,effects = "EVI.max.band", xlab = "Maximum EVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(c) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PREDICTS.n2, data = PREDICTS,effects = "NDVI.max.band", xlab = "Maximum NDVI", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PREDICTS$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(d) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PREDICTS.e3b, data = PREDICTS,effects = "EVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "EVI Yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
mtext("(e) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PREDICTS.n3b, data = PREDICTS,effects = "NDVI.band.yield", otherFactors = list(Higher_taxon = levels(PREDICTS$Higher_taxon)[1]), xlab = "NDVI Yield", ylab = "Species Richness", byFactor=NULL, logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,line.cols=cols_higher_taxon)
mtext("(f) ", side=3, line=-1.5, adj=1.0, cex=1.2)

dev.off()