###
# Cleaning up NECTAR (Lizzie Wolkovich) data for analysis
# This dataset is already so clean
# Thanks Lizzie W! :)
#
# Outcomes:
# Formatted into columns: station, year, species, stage, phenology, DOY
# 
#
###

#setwd("/Users/Michael/Documents/Grad School/Research Projects/Phenology variability")
#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

nectar <- read.csv("raw_data/NECTAR_raw.csv")
sp_list_all <- read.csv("raw_data/NECTAR_all_species_list.csv")
sp_harvard <- read.csv("raw_data/NECTAR_harvard_species_list.csv")
sp_luquillo <- read.csv("raw_data/NECTAR_luquillo_species_list.csv") # I have no idea what's going on here

full_length <- nrow(nectar)
#station, year, species, stage, phenology, DOY
nectar_formatted <- data.frame(station = rep(NA,full_length),
                               year = rep(NA,full_length),
                               species = rep(NA,full_length),
                               stage = rep(NA,full_length),
                               phenology = rep(NA,full_length),
                               DOY = rep(NA,full_length))

nectar_formatted$station <- as.factor(nectar$site)
nectar_formatted$year <- as.numeric(nectar$year)
nectar_formatted$species <- as.factor(apply(nectar,1,function(x) paste(x[7],x[8])))
nectar_formatted$stage <- as.factor(nectar$event)
nectar_formatted$phenology <- as.character(nectar$date)
nectar_formatted$DOY <- as.numeric(nectar$doy)

write.csv(nectar_formatted,"clean_data/nectar_ready_for_analysis.csv")

# cleanup
rm(sp_luquillo, nectar, nectar_formatted)

