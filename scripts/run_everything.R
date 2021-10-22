# Backbone script to reproduce analysis for phenology variance paper

setwd("/home/michael/Documents/Grad School/Research Projects/pheno_variance")

### load libraries
library(lubridate)
library(raster)
library(visreg)
library(dplyr)
library(lmtest)
library(data.table)
library(mgsub)
library(Hmisc)
library(arm)
library(mclust)
library(lme4)
library(lmerTest)
library(RColorBrewer)
library(car)
library(MuMIn)
library(beepr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(sf)
library(cowplot)
library(autoimage)
library(png)
library(DescTools)
library(flextable)
library(animation)
library(textclean)
library(readxl)
library(pbapply)
library(ncdf4)
library(Taxonstand)
library(stringr)
library(BIEN)
library(taxize)
library(curl)
library(devtools)
install_github("willpearse/phest")
library(phest)
library(ape)
library(caper)
library(phytools)
library(quantreg)
print("Done loading libraries")

'%!in%' <- Negate('%in%')

### data cleaning
source("scripts/japan_data_cleaning.R")
source("scripts/korea_data_cleaning.R")
source("scripts/manomet_data_cleaning_all.R")
source("scripts/nectar_data_cleaning.R") 
source("scripts/rmbl_data_cleaning.R")
source("scripts/rothamsted_data_cleaning.R")
source("scripts/ussr_data_cleaning.R")
source("scripts/pep_data_cleaning.R")
gc()
print("Done cleaning data")

### retrieving climate data and assign coords
source("scripts/microclimate.R")
print("Done getting climate data")

### standardizing phenophase names, subsetting to spring, and combining into one file
source("scripts/standardize_data.R")
print("Done standardizing phenophase names")

### retrieving trait data
source("scripts/plant_traits.R") 
source("scripts/bird_traits.R")
print("Done getting trait data")

### phylogeny
source("scripts/phylogeny.R") # come back to this one you've run analysis... or just cut phylogeny
print("Done getting phylogeny")

### analysis
source("script/data_analysis.R")
print("Done running analysis")

### figures
source("script/figures.R")
print("Done making figures")


