# Analysis of EVI and NDVI data for a single pixel

#setwd()

library(lme4)
library(MuMIn) ## for RsquareGLMM function
library(car) ## for Anova function

source("/Users/hp1111/PhD/Functions/NEWplot_lmer_means.R")
source("/Users/hp1111/PhD/AdriannasFunctions/model_plot.R")


file <- list.files(path = "MOD13Q1_IndividualPixel", pattern = "^MODIS_Data") ## The file that begins MODIS_Data is the original data with MODIS appended
file
# file2 ## Need to ensure there is only 1 file, or that we know which file we want to load

VIs <- read.csv(file.path("MOD13Q1_IndividualPixel", file), header = TRUE)
names(VIs)
nrow(VIs)
str(VIs)



###### Appending additional summary information ### 
file2 <- list.files(path = "MOD13Q1_IndividualPixel", pattern = "^MODIS_Summary")
file2
file2 <- read.csv(file.path("MOD13Q1_IndividualPixel", file2), header = TRUE)

names(file2)
nrow(file2) ## Twice the number of rows as downloaded subsets - EVI and NDVI are on separate rows 

evi <- file2[file2$data.band == "250m_16_days_EVI",]
ndvi <- file2[file2$data.band == "250m_16_days_NDVI",]

## Add all additional columns to the dataset (SD, minimum, maximum etc.), rename, and remove unwanted rows
VIs <- merge(VIs, evi, by.x = "ID", by.y = "ID", all.x = TRUE)
VIs$lat.y <- NULL
VIs$long.y <- NULL
VIs$end.date.y <- NULL
VIs$data.band <- NULL
newnames <- names(VIs)[28:length(names(VIs))]
newnames <- paste("EVI.", newnames, sep = "")
names(VIs)[28:length(names(VIs))] <- newnames

VIs <- merge(VIs, ndvi, by.x = "ID", by.y = "ID", all.x = TRUE)
VIs$lat <- NULL
VIs$long <- NULL
VIs$end.date <- NULL
VIs$start.date.x <- NULL
VIs$data.band <- NULL
newnames <- names(VIs)[35:length(names(VIs))]
newnames <- paste("NDVI.", newnames, sep = "")
names(VIs)[35:length(names(VIs))] <- newnames

summary(VIs$EVI.mean.band)
summary(VIs$NDVI.mean.band)

### To check that everything is correct, VIs$X250m_16_days_EVI_pixel1 and VIs$EVI.mean.band should the same numbers (both the mean)
sum(VIs$X250m_16_days_EVI_pixel1 == VIs$EVI.mean.band, na.rm = TRUE)
## same with NDVI
sum(VIs$X250m_16_days_NDVI_pixel1 == VIs$NDVI.mean.band, na.rm = TRUE)
## Good.

## Data exploration
summary(VIs$EVI.mean.band)
# NAs present, so let's remove them
VIs <-VIs[complete.cases(VIs$EVI.mean.band),]

summary(VIs$NDVI.mean.band)
## No NAs present


summary(VIs$EVI.band.yield)
summary(VIs$NDVI.band.yield)
## another 8 data points missing
VIs <-VIs[complete.cases(VIs$EVI.band.yield),]

par(mfrow = c(2,2))
hist(VIs$EVI.mean.band, xlab = "EVI", main = "Histogram of mean EVI")
hist(VIs$NDVI.mean.band, xlab = "NDVI", main = "Histogram of mean NDVI")
hist(VIs$EVI.band.yield, xlab = "EVI", main = "Histogram of EVI yield")
hist(VIs$NDVI.band.yield, xlab = "NDVI", main = "Histogram of NDVI yield") ## tiny bit of skew, but nothing to worry about


par(mfrow = c(2,2))
plot(VIs$EVI.mean.band, col = VIs$Source_ID, ylab = "EVI", main = "Scatterplot of EVI by source")
plot(VIs$NDVI.mean.band, col = VIs$Source_ID, ylab = "NDVI", main = "Scatterplot of NDVI by source")
plot(VIs$EVI.band.yield, col = VIs$Source_ID, ylab = "EVI", main = "Scatterplot of EVI by source")
plot(VIs$NDVI.band.yield, col = VIs$Source_ID, ylab = "NDVI", main = "Scatterplot of NDVI by source")
## Expect to see the horizontal bands in colour, but not nessecarily that each source was consistently higher or lower (bar the few in the middle). Homogeneity within sources.
## A random effects structure will deal with this homogeneity
## As technically the data points are nested within each source, it will also deal with that

## How many publications are in the dataset
length(unique(VIs$Source_ID))
length(unique(paste(VIs$Source_ID, VIs$Study_number))) ## Same number of studies and sources, so random effects will only need the one

plot(VIs$EVI.mean.band, VIs$NDVI.mean.band, xlab = "EVI", ylab = "NDVI")
abline(a=0, b=1)
## As expected, strong relationship between the two (EVI is similar to NDVI but with a small correction factor), so both will be analysed independently

plot(VIs$Higher_taxon, VIs$Species_richness, ylab = "Species Richness")
## Clear distinction between higher taxon of species richness. Again, higher taxon can be a random or fixed effect to deal with that

table(VIs$Higher_taxon)
## Too few Arachnida, so remove them from the dataset
VIs <- droplevels(VIs[VIs$Higher_taxon != "Arachnida",])

## Create folder for figures ##
if(file.exists("Figures") == FALSE)
	dir.create("Figures")
	
dir.create(file.path("Figures", "EVISinglePixel"))
dir.create(file.path("Figures", "NDVISinglePixel"))

####### Univariate models (+higher taxon)######

## EVI mean
cols_higher_taxon <- c("blue", "darkorchid", "gold", "darkgreen")
VIs.e1 <- glmer(Species_richness ~ EVI.mean.band * Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.e1b <- glmer(Species_richness ~ EVI.mean.band + Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e1, VIs.e1b) ## p value < 0.05 need interaction

summary(VIs.e1)
model_plots(VIs.e1)

pdf(file.path("Figures", "EVISinglePixel", "Mean.pdf"))
NEWplot_lmer_means(VIs.e1, var.of.interest = "EVI.mean.band", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Mean EVI", CIcol = cols_higher_taxon, meancol = cols_higher_taxon, line.width = 1)
legend("topright", legend = levels(VIs$Higher_taxon), col = cols_higher_taxon, lwd = 2)
dev.off()
r.squaredGLMM(VIs.e1) # 
      # R2m       R2c 
# 0.5871584 0.8097819


## EVI max

VIs.e2 <- glmer(Species_richness ~  Higher_taxon * EVI.max.band + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.e2b <- glmer(Species_richness ~ EVI.max.band + Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e2, VIs.e2b) ## need interaction

pdf(file.path("Figures", "EVISinglePixel", "Max.pdf"))
NEWplot_lmer_means(VIs.e2, var.of.interest = "EVI.max.band", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Max EVI", CIcol = cols_higher_taxon, meancol = cols_higher_taxon, line.width = 2)
legend("topright", legend = levels(VIs$Higher_taxon), col = cols_higher_taxon, lwd = 2)
dev.off()
summary(VIs.e2)
r.squaredGLMM(VIs.e2)
      # R2m       R2c 
# 0.6022174 0.7996292 


### EVI yield

VIs.e3 <- glmer(Species_richness ~ EVI.band.yield * Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.e3b <- glmer(Species_richness ~ EVI.band.yield + Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e3, VIs.e3b) ## not significant ( p > 0.05)

VIs.e3c <- glmer(Species_richness ~ EVI.band.yield + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e3b, VIs.e3c) ## Need higher taxa
VIs.e3d <- glmer(Species_richness ~ Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.e3b, VIs.e3d) ## do need yield
# Stick with the additive model

summary(VIs.e3b)

pdf(file.path("Figures", "EVISinglePixel", "Yield.pdf"))
NEWplot_lmer_means(VIs.e3b, var.of.interest = "EVI.band.yield", length = 50, ylabel = "Species Richness", xlabel = "EVI yield", CIcol = "black", meancol = "black", line.width = 1)
dev.off()
r.squaredGLMM(VIs.e2b)
      # R2m       R2c 
# 0.5774434 0.8018797 





### Without arachnids


VIs.n1 <- glmer(Species_richness ~ NDVI.mean.band * Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.n1b <- glmer(Species_richness ~ NDVI.mean.band + Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n1, VIs.n1b) ## Significantly different

summary(VIs.n1)
r.squaredGLMM(VIs.n1)
      # R2m       R2c 
# 0.5696208 0.8253525

pdf(file.path("Figures", "EVISinglePixel", "Mean.pdf"))
NEWplot_lmer_means(VIs.n1, var.of.interest = "NDVI.mean.band", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Mean NDVI", CIcol = cols_higher_taxon, meancol = cols_higher_taxon, line.width = 2)
legend("topright", legend = levels(VIs$Higher_taxon), col = cols_higher_taxon, lwd = 2)
dev.off()






VIs.n2 <- glmer(Species_richness ~ NDVI.max.band * Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.n2b <- glmer(Species_richness ~ NDVI.max.band + Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n2, VIs.n2b) ## significantly different

summary(VIs.n2)
r.squaredGLMM(VIs.n2)
      # R2m       R2c 
# 0.5954573 0.8041141 

pdf(file.path("Figures", "EVISinglePixel", "Max.pdf"))
NEWplot_lmer_means(VIs.n2, var.of.interest = "NDVI.max.band", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Max NDVI", CIcol = cols_higher_taxon, meancol = cols_higher_taxon, line.width = 2)
legend("topright", legend = levels(VIs$Higher_taxon), col = cols_higher_taxon, lwd = 2)
dev.off()



VIs.n3 <- glmer(Species_richness ~ NDVI.band.yield * Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VIs.n3b <- glmer(Species_richness ~ NDVI.band.yield + Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n3, VIs.n3b) ## don't need interaction

VIs.n3c <- glmer(Species_richness ~ NDVI.band.yield + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
anova(VIs.n3b, VIs.n3c) # can't remove
VIs.n3d <- glmer(Species_richness ~ Higher_taxon + (1|SS), data = VIs, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VIs.n3b, VIs.n3d) # can't remove
# Stick with additive model

summary(VIs.n3b)
r.squaredGLMM(VIs.n3b)

pdf(file.path("Figures", "EVISinglePixel", "Yield.pdf"))
NEWplot_lmer_means(VIs.n3b, var.of.interest = "NDVI.band.yield", length = 50, ylabel = "Species Richness", xlabel = "NDVI yield", CIcol = "black", meancol = "black", line.width = 1)
dev.off()
