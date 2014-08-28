PREDICTSlibrary(lme4)
library(nlme)
library(MuMIn)

source("/Users/hp1111/PhD/Functions/NEWplot_lmer_means.R")
source("/Users/hp1111/PhD/AdriannasFunctions/model_plot.R")
source("/Users/hp1111/PhD/MODISTools Paper/Analysis/MyAnalysis/creating_weight_matrix.R")



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

higher_taxon_cols <- c("blue", "darkorchid", "gold", "darkgreen")

#mean
VI.e1 <- glmer(Species_richness ~ Arith.EVI.mean * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.e1b <- glmer(Species_richness ~ Arith.EVI.mean + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e1, VI.e1b) # significant
summary(VI.e1)
model_plot(VI.e1)
r.squaredGLMM(VI.e1)

pdf(file.path("Figures", "EVITile", "ARITH_mean.pdf"))
NEWplot_lmer_means(VI.e1, var.of.interest = "Arith.EVI.mean", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Mean EVI (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()



#yield

VI.e2 <- glmer(Species_richness ~ Arith.EVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.e2b <- glmer(Species_richness ~ Arith.EVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e2, VI.e2b) ## significant

summary(VI.e2)
r.squaredGLMM(VI.e2)

pdf(file.path("Figures", "EVITile", "ARITH_yield.pdf"))
NEWplot_lmer_means(VI.e2, var.of.interest = "Arith.EVI.yield", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "EVI yield (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()


#max
VI.e3 <- glmer(Species_richness ~ Arith.max.EVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VI.e3b <- glmer(Species_richness ~ Arith.max.EVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e3, VI.e3b) ## significant

summary(VI.e3)
r.squaredGLMM(VI.e3)

pdf(file.path("Figures", "EVITile", "ARITH_Max.pdf"))
NEWplot_lmer_means(VI.e3, var.of.interest = "Arith.max.EVI", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Arith.max.EVI (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()






## Weighted (exponential decay model)
#mean
VI.e4 <- glmer(Species_richness ~ Weighted.mean.EVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.e4b <- glmer(Species_richness ~ Weighted.mean.EVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e4, VI.e4b) # significant
summary(VI.e4)
r.squaredGLMM(VI.e4)

pdf(file.path("Figures", "EVITile", "WEIGHTED_Mean.pdf"))
NEWplot_lmer_means(VI.e4, var.of.interest = "Weighted.mean.EVI", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Weighted mean EVI (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()


#yield

VI.e5 <- glmer(Species_richness ~ Weighted.EVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) # didn't converge

VI.e5b <- glmer(Species_richness ~ Weighted.EVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e5, VI.e5b) ## significant
summary(VI.e5)
r.squaredGLMM(VI.e5)

pdf(file.path("Figures", "EVITile", "WEIGHTED_yield.pdf"))
NEWplot_lmer_means(VI.e5, var.of.interest = "Weighted.EVI.yield", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "EVI yield (weighted mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()


#max
VI.e6 <- glmer(Species_richness ~ Weighted.max.EVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
VI.e6b <- glmer(Species_richness ~ Weighted.max.EVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e6, VI.e6b) ## significant
summary(VI.e6)
r.squaredGLMM(VI.e6)


pdf(file.path("Figures", "EVITile", "WEIGHTED_Max.pdf"))
NEWplot_lmer_means(VI.e6, var.of.interest = "Weighted.max.EVI", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Max EVI (weighted mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()



#### NDVI MODELS ####
# Averaged
#mean
VI.n1 <- glmer(Species_richness ~ Arith.NDVI.mean * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.n1b <- glmer(Species_richness ~ Arith.NDVI.mean + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.n1, VI.n1b) # significant

summary(VI.n1)
r.squaredGLMM(VI.n1)

pdf(file.path("Figures", "EVITile", "ARITH_Mean.pdf"))
NEWplot_lmer_means(VI.n1, var.of.interest = "Arith.NDVI.mean", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Mean NDVI (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()



#yield

VI.n2 <- glmer(Species_richness ~ Arith.NDVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.n2b <- glmer(Species_richness ~ Arith.NDVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
anova(VI.n2, VI.n2b) ## significant

summary(VI.n2)
r.squaredGLMM(VI.n2)

pdf(file.path("Figures", "EVITile", "ARITH_yield.pdf"))
NEWplot_lmer_means(VI.n2, var.of.interest = "Arith.NDVI.yield", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "NDVI yield (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()

#max
VI.n3 <- glmer(Species_richness ~ Arith.max.NDVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) 
VI.n3b <- glmer(Species_richness ~ Arith.max.NDVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.n3, VI.n3b) ## significant
summary(VI.n3)
r.squaredGLMM(VI.n3)

pdf(file.path("Figures", "EVITile", "ARITH_max.pdf"))
NEWplot_lmer_means(VI.n3, var.of.interest = "Arith.max.NDVI", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Max NDVI (arithmetic mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()



# Weighted (exponential decay model)
#mean
VI.n4 <- glmer(Species_richness ~ Weighted.mean.NDVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.n4b <- glmer(Species_richness ~ Weighted.mean.NDVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.n4, VI.n4b) # significant
summary(VI.n4)
r.squaredGLMM(VI.n4)

pdf(file.path("Figures", "EVITile", "WEIGHTED_mean.pdf"))
NEWplot_lmer_means(VI.n4, var.of.interest = "Weighted.mean.NDVI", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Mean NDVI (weighted mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()

#yield

VI.n5 <- glmer(Species_richness ~ Weighted.NDVI.yield * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

VI.n5b <- glmer(Species_richness ~ Weighted.NDVI.yield + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.n5, VI.n5b) ## significant
summary(VI.n5)
r.squaredGLMM(VI.n5)


pdf(file.path("Figures", "EVITile", "WEIGHTED_yield.pdf"))
NEWplot_lmer_means(VI.n5, var.of.interest = "Weighted.NDVI.yield", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "NDVI yield (weighted mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()


#max
VI.n6 <- glmer(Species_richness ~ Weighted.max.NDVI * Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))) # 
VI.n6b <- glmer(Species_richness ~ Weighted.max.NDVI + Higher_taxon + (1|SS), data = PRED, family = poisson, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
anova(VI.e6, VI.e6b) ## significant
summary(VI.n6)
r.squaredGLMM(VI.n6)


pdf(file.path("Figures", "EVITile", "WEIGHTED_max.pdf"))
NEWplot_lmer_means(VI.n6, var.of.interest = "Weighted.max.NDVI", var.interaction = "Higher_taxon", length = 50, ylabel = "Species Richness", xlabel = "Max NDVI (weighted mean)", CIcol = higher_taxon_cols, meancol = higher_taxon_cols, line.width = 2)
legend("topright", legend = levels(PRED$Higher_taxon), col = higher_taxon_cols, lwd = 2)
dev.off()

