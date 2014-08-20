## Analysis done on a MAC computer (10.8.5).

# Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/usr/texbin",sep=":")) 

library(devtools)
install_github("MODISTools", "seantuck12", ref = "master") # Has author ID and UpdateSubsets fix
library(MODISTools)

sessionInfo()
# MODISTools_0.94.2

## setwd()


dat <- read.csv("Data/TemperateSites_2014-08-14.csv")


## Making column names match what is needed by MODISTools
names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date <- as.Date(dat$end.date)

## Creating my unique ID for each unique site
dat$ID <- paste(dat$SS, dat$lat, dat$long, dat$end.date, sep=" ")

###############################
#### EVI and NDVI DOWNLOAD ####
###############################

## Specifying my product, and associated bands
products <- "MOD13Q1"
bands <- c("250m_16_days_pixel_reliability", "250m_16_days_NDVI", "250m_16_days_EVI")

if(file.exists("MOD13Q1_IndividualPixel") == FALSE)
	dir.create("MOD13Q1_IndividualPixel")

## Call to download the subsets
MODISSubsets(dat, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD13Q1_IndividualPixel", TimeSeriesLength = 2)


# To make sure that all time-series were downloaded, run the data frame through update subsets
# to create a new data frame of all the site that might have failed to download first time through
take2 <- UpdateSubsets(dat,StartDate = FALSE, Dir = "MOD13Q1_IndividualPixel")


# Now to use MODISSubsets again
MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD13Q1_IndividualPixel", TimeSeriesLength = 2)


############################################
## If all bands use the same QA band, they can be summarised together
MODISSummaries(dat, Dir = "MOD13Q1_IndividualPixel", Product = products, Bands = bands[2:3], ValidRange = c(-2000, 10000), NoDataFill = -3000, ScaleFactor = 0.0001, StartDate = FALSE, Interpolate = TRUE, Yield = TRUE, QualityScreen = TRUE, QualityBand = "250m_16_days_pixel_reliability", QualityThreshold = 0, DiagnosticPlot = TRUE)

## Should now be a file with the MODIS data appended to our original data in the folder "MOD13Q1_IndividualPixel", ready for analysis

############
########## Larger Tiles ############

if(file.exists("MOD13Q1_6x6km") == FALSE)
	dir.create("MOD13Q1_6x6km")
	
## Call to download the subsets
MODISSubsets(dat, Products = products, Bands = bands, Size = c(3,3), SaveDir = "MOD13Q1_6x6km", TimeSeriesLength = 2)


# To make sure that all time-series were downloaded, run the data frame through update subsets
# to create a new data frame of all the site that might have failed to download first time through
take2 <- UpdateSubsets(dat,StartDate = FALSE, Dir = "MOD13Q1_6x6km")

# Now to use MODISSubsets again
MODISSubsets(take2, Products = products, Bands = bands, Size = c(3,3), SaveDir = "MOD13Q1_6x6km", TimeSeriesLength = 2)



############################################
## If all bands use the same QA band, they can be summarised together
MODISSummaries(dat, Dir = "MOD13Q1_6x6km", Product = products, Bands = bands[2:3], ValidRange = c(-2000, 10000), NoDataFill = -3000, ScaleFactor = 0.0001, StartDate = FALSE, Interpolate = TRUE, Yield = TRUE, QualityScreen = TRUE, QualityBand = "250m_16_days_pixel_reliability", QualityThreshold = 1)

## Should now be a file with the MODIS data appended to our original data in the folder "MOD13Q1_6x6km", ready for analysis
