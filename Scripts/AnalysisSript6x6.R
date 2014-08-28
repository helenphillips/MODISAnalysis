library(lme4)
library(nlme)
library(MuMIn)
#library(roquefort) # private library used for plotting

# setwd()

source("/Users/hp1111/PhD/AdriannasFunctions/model_plot.R")
source("~/Scripts/creating_weight_matrix.R")



file <- list.files(path = "MOD13Q1_6x6km", pattern = "^MODIS_Data") ## The file that begins MODIS_Data is the original data with MODIS appended
file2 <- list.files(path = "MOD13Q1_6x6km", pattern = "^MODIS_Summary") ## The file that begins MODIS_Summary contains additional summarised data, such as Yield, standard deviation, min and max (if any were selected). This will need to be appeneded onto the original data

file
file2 ## Need to ensure there is only 1 file, or that we know which file we want to load




PREDICTS <- read.csv(file.path("MOD13Q1_6x6km", file), header = TRUE) ## a lot of columns
length(names(PREDICTS)) ## Let's strip it down...

PRED <- PREDICTS[,1:25]

file2 <- read.csv(file.path("MOD13Q1_6x6km", file2), header = TRUE)
names(file2)

evi <- file2[file2$data.band == "250m_16_days_EVI",]
ndvi <- file2[file2$data.band == "250m_16_days_NDVI",]

## Calculate the average of all the pixels at each site
evimean <- data.frame(tapply(evi$mean.band, evi$ID, mean, na.rm = TRUE))
names(evimean) <- "Arith.EVI.mean"
evimean$ID <- row.names(evimean)
PRED <- merge(PRED, evimean, by.x= "ID", by.y = "ID", all.x = TRUE)

ndvimean <- data.frame(tapply(ndvi$mean.band, ndvi$ID, mean, na.rm = TRUE))
names(ndvimean) <- "Arith.NDVI.mean"
ndvimean$ID <- row.names(ndvimean)
PRED <- merge(PRED, ndvimean, by.x= "ID", by.y = "ID", all.x = TRUE)






## Can calculate the average EVI yield as well 
eviyield <- data.frame(tapply(evi$band.yield, evi$ID, mean, na.rm = TRUE))
names(eviyield) <- "Arith.EVI.yield"
eviyield$ID <- row.names(eviyield)
PRED <- merge(PRED, eviyield, by.x= "ID", by.y = "ID", all.x = TRUE)

ndviyield <- data.frame(tapply(ndvi$band.yield, ndvi$ID, mean, na.rm = TRUE))
names(ndviyield) <- "Arith.NDVI.yield"
ndviyield$ID <- row.names(ndviyield)
PRED <- merge(PRED, ndviyield, by.x= "ID", by.y = "ID", all.x = TRUE)

## minimums
evimin <- data.frame(tapply(evi$min.band, evi$ID, mean, na.rm = TRUE))
names(evimin) <- "Arith.min.EVI"
evimin$ID <- row.names(evimin)
PRED <- merge(PRED, evimin, by.x= "ID", by.y = "ID", all.x = TRUE)

ndvimin <- data.frame(tapply(ndvi$min.band, ndvi$ID, mean, na.rm = TRUE))
names(ndvimin) <- "Arith.min.NDVI"
ndvimin$ID <- row.names(ndvimin)
PRED <- merge(PRED, ndvimin, by.x= "ID", by.y = "ID", all.x = TRUE)

## maximums
evimax <- data.frame(tapply(evi$max.band, evi$ID, mean, na.rm = TRUE))
names(evimax) <- "Arith.max.EVI"
evimax$ID <- row.names(evimax)
PRED <- merge(PRED, evimax, by.x= "ID", by.y = "ID", all.x = TRUE)

ndvimax <- data.frame(tapply(ndvi$max.band, ndvi$ID, mean, na.rm = TRUE))
names(ndvimax) <- "Arith.max.NDVI"
ndvimax$ID <- row.names(ndvimax)
PRED <- merge(PRED, ndvimax, by.x= "ID", by.y = "ID", all.x = TRUE)



## Now to create the weighting vector, and where this might get complicated.
## The "tile" of pixels is now a single row of numbers, where the middle number is the central pixel

evi.tile <- PREDICTS[,c(24, 651:ncol(PREDICTS))] ## these are the means (and ID)
ndvi.tile <- PREDICTS[,c(24, 26:650)]

## This functions creates a weight vector based on exponential distance decay

	weights <- weight_matrix(evi.tile[1, 2:ncol(evi.tile)])
	
	# Long way of going through each row and calculating the spatial mean based on the weighted mean
	
	for(i in 1:nrow(evi.tile)){
		weightedmean <- weighted.mean(evi.tile[i, 2:ncol(evi.tile)], weights, na.rm = TRUE)
		PRED$Weighted.mean.EVI[i] <- weightedmean
	}
		
	for(i in 1:nrow(ndvi.tile)){
		weightedmean <- weighted.mean(ndvi.tile[i, 2:ncol(ndvi.tile)], weights, na.rm = TRUE)
		PRED$Weighted.mean.NDVI[i] <- weightedmean
	}


## weighted yields
weightedeviyields <- as.data.frame(tapply(evi$band.yield, evi$ID, weighted.mean, w = weights, na.rm = TRUE))
names(weightedeviyields) <- "Weighted.EVI.yield"
weightedeviyields$ID <- row.names(weightedeviyields)
PRED <- merge(PRED, weightedeviyields, by.x= "ID", by.y = "ID", all.x = TRUE)

weightedndviyields <- as.data.frame(tapply(ndvi$band.yield, ndvi$ID, weighted.mean, w = weights, na.rm = TRUE))
names(weightedndviyields) <- "Weighted.NDVI.yield"
weightedndviyields$ID <- row.names(weightedndviyields)
PRED <- merge(PRED, weightedndviyields, by.x= "ID", by.y = "ID", all.x = TRUE)


## weighted minimums
weightedevimin <- as.data.frame(tapply(evi$min.band, evi$ID, weighted.mean, w = weights, na.rm = TRUE))
names(weightedevimin) <- "Weighted.min.EVI"
weightedevimin$ID <- row.names(weightedevimin)
PRED <- merge(PRED, weightedevimin, by.x= "ID", by.y = "ID", all.x = TRUE)

weightedndvimin <- as.data.frame(tapply(ndvi$min.band, ndvi$ID, weighted.mean, w = weights, na.rm = TRUE))
names(weightedndvimin) <- "Weighted.min.NDVI"
weightedndvimin$ID <- row.names(weightedndvimin)
PRED <- merge(PRED, weightedndvimin, by.x= "ID", by.y = "ID", all.x = TRUE)

## weighted maximum

weightedevimax <- as.data.frame(tapply(evi$max.band, evi$ID, weighted.mean, w = weights, na.rm = TRUE))
names(weightedevimax) <- "Weighted.max.EVI"
weightedevimax$ID <- row.names(weightedevimax)
PRED <- merge(PRED, weightedevimax, by.x= "ID", by.y = "ID", all.x = TRUE)

weightedndvimax <- as.data.frame(tapply(ndvi$max.band, ndvi$ID, weighted.mean, w = weights, na.rm = TRUE))
names(weightedndvimax) <- "Weighted.max.NDVI"
weightedndvimax$ID <- row.names(weightedndvimax)
PRED <- merge(PRED, weightedndvimax, by.x= "ID", by.y = "ID", all.x = TRUE)


##########

summary(PRED$Arith.EVI.mean)
summary(PRED$Arith.EVI.yield) ## 33 NAs in the yield (will be the same for all)
PRED<-PRED[complete.cases(PRED$Arith.EVI.yield),]

par(mfrow=c(4,2))
hist(PRED$Arith.EVI.mean)
hist(PRED$Arith.NDVI.mean)
hist(PRED$Arith.EVI.yield)
hist(PRED$Arith.NDVI.yield)
hist(PRED$Weighted.mean.EVI)
hist(PRED$Weighted.mean.NDVI)
hist(PRED$Weighted.EVI.yield)
hist(PRED$Weighted.NDVI.yield)

## We assume that we seem the same structure in terms of homogeneity from sources

table(PRED$Higher_taxon)
## again, not enough Arachnida for the models, so remove those data
PRED <- droplevels(PRED[PRED$Higher_taxon != "Arachnida",])

## Create folders for figures
if(file.exists("Figures") == FALSE)
	dir.create("Figures")
	
dir.create(file.path("Figures", "EVITile"))
dir.create(file.path("Figures", "NDVITile"))


#### EVI MODELS #####
# Arithmetic

cols_higher_taxon <- c("#000000", "#9932CC", "#FFD700", "#006400")

#mean
PRED.e1 <- glmer(Species_richness ~ Arith.EVI.mean * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.e1b <- glmer(Species_richness ~ Arith.EVI.mean + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e1, PRED.e1b) # significant
summary(PRED.e1)
model_plot(PRED.e1)
r.squaredGLMM(PRED.e1)

pdf(file.path("Figures", "EVITile", "ARITH_mean.pdf"))
PlotContEffects(PRED.e1, data = PRED,effects = "Arith.EVI.mean", xlab = "Mean EVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()



#yield

PRED.e2 <- glmer(Species_richness ~ Arith.EVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.e2b <- glmer(Species_richness ~ Arith.EVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e2, PRED.e2b) ## significant

summary(PRED.e2)
r.squaredGLMM(PRED.e2)

pdf(file.path("Figures", "EVITile", "ARITH_yield.pdf"))
PlotContEffects(PRED.e2, data = PRED,effects = "Arith.EVI.yield", xlab = "EVI Yield (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


#max
PRED.e3 <- glmer(Species_richness ~ Arith.max.EVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PRED.e3b <- glmer(Species_richness ~ Arith.max.EVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e3, PRED.e3b) ## significant

summary(PRED.e3)
r.squaredGLMM(PRED.e3)

pdf(file.path("Figures", "EVITile", "ARITH_Max.pdf"))
PlotContEffects(PRED.e3, data = PRED,effects = "Arith.max.EVI", xlab = "Max EVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()






## Weighted (exponential decay model)
#mean
PRED.e4 <- glmer(Species_richness ~ Weighted.mean.EVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.e4b <- glmer(Species_richness ~ Weighted.mean.EVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e4, PRED.e4b) # significant
summary(PRED.e4)
r.squaredGLMM(PRED.e4)

pdf(file.path("Figures", "EVITile", "WEIGHTED_Mean.pdf"))
PlotContEffects(PRED.e4, data = PRED,effects = "Weighted.mean.EVI", xlab = "Mean EVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


#yield

PRED.e5 <- glmer(Species_richness ~ Weighted.EVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) # didn't converge

PRED.e5b <- glmer(Species_richness ~ Weighted.EVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e5, PRED.e5b) ## significant
summary(PRED.e5)
r.squaredGLMM(PRED.e5)

pdf(file.path("Figures", "EVITile", "WEIGHTED_yield.pdf"))
PlotContEffects(PRED.e5, data = PRED,effects = "Weighted.EVI.yield", xlab = "EVI yield (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


#max
PRED.e6 <- glmer(Species_richness ~ Weighted.max.EVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
PRED.e6b <- glmer(Species_richness ~ Weighted.max.EVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e6, PRED.e6b) ## significant
summary(PRED.e6)
r.squaredGLMM(PRED.e6)


pdf(file.path("Figures", "EVITile", "WEIGHTED_Max.pdf"))
PlotContEffects(PRED.e6, data = PRED,effects = "Weighted.max.EVI", xlab = "Maximum EVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()

## Supplementary figures
pdf(file.path("Figures", "SupplementaryFigures", "EVITile.pdf"), height = 11.7, width = 8.3)

par(mfrow = c(3,2))
PlotContEffects(PRED.e3, data = PRED,effects = "Arith.max.EVI", xlab = "Maximum EVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(a) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PRED.e6, data = PRED,effects = "Weighted.max.EVI", xlab = "Maximum EVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(b) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PRED.e1, data = PRED,effects = "Arith.EVI.mean", xlab = "Mean EVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(c) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PRED.e4, data = PRED,effects = "Weighted.mean.EVI", xlab = "Mean EVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(d) ", side=3, line=-1.5, adj=1.0, cex=1.2)

PlotContEffects(PRED.e2, data = PRED,effects = "Arith.EVI.yield", xlab = "EVI Yield (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(e) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PRED.e5, data = PRED,effects = "Weighted.EVI.yield", xlab = "EVI yield (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(f) ", side=3, line=-1.5, adj=1.0, cex=1.2)

dev.off()



#### NDVI MODELS ####
# Averaged
#mean
PRED.n1 <- glmer(Species_richness ~ Arith.NDVI.mean * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.n1b <- glmer(Species_richness ~ Arith.NDVI.mean + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.n1, PRED.n1b) # significant

summary(PRED.n1)
r.squaredGLMM(PRED.n1)

pdf(file.path("Figures", "NDVITile", "ARITH_Mean.pdf"))
PlotContEffects(PRED.n1, data = PRED,effects = "Arith.NDVI.mean", xlab = "Mean NDVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()



#yield

PRED.n2 <- glmer(Species_richness ~ Arith.NDVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.n2b <- glmer(Species_richness ~ Arith.NDVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
anova(PRED.n2, PRED.n2b) ## significant

summary(PRED.n2)
r.squaredGLMM(PRED.n2)

pdf(file.path("Figures", "NDVITile", "ARITH_yield.pdf"))
PlotContEffects(PRED.n2, data = PRED,effects = "Arith.NDVI.yield", xlab = "NDVI yield (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()

#max
PRED.n3 <- glmer(Species_richness ~ Arith.max.NDVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
PRED.n3b <- glmer(Species_richness ~ Arith.max.NDVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.n3, PRED.n3b) ## significant
summary(PRED.n3)
r.squaredGLMM(PRED.n3)

pdf(file.path("Figures", "NDVITile", "ARITH_max.pdf"))
PlotContEffects(PRED.n3, data = PRED,effects = "Arith.max.NDVI", xlab = "Maximum NDVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()



# Weighted (exponential decay model)
#mean
PRED.n4 <- glmer(Species_richness ~ Weighted.mean.NDVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.n4b <- glmer(Species_richness ~ Weighted.mean.NDVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.n4, PRED.n4b) # significant
summary(PRED.n4)
r.squaredGLMM(PRED.n4)

pdf(file.path("Figures", "NDVITile", "WEIGHTED_mean.pdf"))
PlotContEffects(PRED.n4, data = PRED,effects = "Weighted.mean.NDVI", xlab = "Mean NDVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()

#yield

PRED.n5 <- glmer(Species_richness ~ Weighted.NDVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

PRED.n5b <- glmer(Species_richness ~ Weighted.NDVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.n5, PRED.n5b) ## significant
summary(PRED.n5)
r.squaredGLMM(PRED.n5)


pdf(file.path("Figures", "NDVITile", "WEIGHTED_yield.pdf"))
PlotContEffects(PRED.n5, data = PRED,effects = "Weighted.NDVI.yield", xlab = "NDVI yield (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


#max
PRED.n6 <- glmer(Species_richness ~ Weighted.max.NDVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) # 
PRED.n6b <- glmer(Species_richness ~ Weighted.max.NDVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(PRED.e6, PRED.e6b) ## significant
summary(PRED.n6)
r.squaredGLMM(PRED.n6)


pdf(file.path("Figures", "NDVITile", "WEIGHTED_max.pdf"))
PlotContEffects(PRED.n6, data = PRED,effects = "Weighted.max.NDVI", xlab = "Maximum NDVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
dev.off()


## Supplementary Figures
pdf(file.path("Figures", "SupplementaryFigures", "NDVITile.pdf"), height = 11.7, width = 8.3)

par(mfrow=c(3,2))
PlotContEffects(PRED.n3, data = PRED,effects = "Arith.max.NDVI", xlab = "Maximum NDVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(a) ", side=3, line=-1.5, adj=1.0, cex=1.2)

PlotContEffects(PRED.n6, data = PRED,effects = "Weighted.max.NDVI", xlab = "Maximum NDVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(b) ", side=3, line=-1.5, adj=1.0, cex=1.2)

PlotContEffects(PRED.n1, data = PRED,effects = "Arith.NDVI.mean", xlab = "Mean NDVI (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(c) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PRED.n4, data = PRED,effects = "Weighted.mean.NDVI", xlab = "Mean NDVI (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(d) ", side=3, line=-1.5, adj=1.0, cex=1.2)

PlotContEffects(PRED.n2, data = PRED,effects = "Arith.NDVI.yield", xlab = "NDVI yield (Arithmetic mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)
mtext("(e) ", side=3, line=-1.5, adj=1.0, cex=1.2)


PlotContEffects(PRED.n5, data = PRED,effects = "Weighted.NDVI.yield", xlab = "NDVI yield (Weighted mean)", ylab = "Species Richness", 
	byFactor="Higher_taxon", logLink="e",plotRug=FALSE,seMultiplier=1.96, params=list(),axis.log="y",ylim=NULL,
	line.cols=cols_higher_taxon)
legend("topleft", legend = levels(PRED$Higher_taxon), col = cols_higher_taxon, lwd = 4)mtext("(a) ", side=3, line=-1.5, adj=1.0, cex=1.5)
mtext("(f) ", side=3, line=-1.5, adj=1.0, cex=1.2)

dev.off()