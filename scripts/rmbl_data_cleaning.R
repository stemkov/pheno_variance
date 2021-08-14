###
# Cleaning up RMBL flower data for analysis
#
#
# Outcomes:
#
# Got rid of first column, which was just redundant sp./event code
# Formatted into columns: station, year, species, stage, phenology, DOY
# Assuming the phenology dates are represented as "MMDD". Need to check on this
# Also assuming that "0" phenology dates are actually "NA" ß
# I need to figure out what the "event" codes mean
#
###

#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

library(plyr)

load("raw_data/rmbl.RData")
rmbl <- all.data
rm(all.data)

rmbl_2016 <- rmbl[which(rmbl$year == 2016),]
rmbl_na <- apply(rmbl_2016, 1, function(x) ifelse(sum(is.na(x) > 0), TRUE, FALSE))
just_counts <- rmbl[which(rmbl$floralcount > 0),]
all_counts <- rep(just_counts$doy, just_counts$floralcount)
non_zeros <- rmbl[which(rmbl$floralcount > 0),]
unique_combs <- unique(non_zeros[,1:3])

rmbl_formatted <- data.frame(station = rep(NA,nrow(unique_combs)),
                             year = rep(NA,nrow(unique_combs)),
                             species = rep(NA,nrow(unique_combs)),
                             estimate = rep(NA,nrow(unique_combs)),
                             lower = rep(NA,nrow(unique_combs)),
                             upper = rep(NA,nrow(unique_combs)))

for(i in 1:nrow(unique_combs)){
  #i <- 1
  y <- unique_combs[i,1]
  site <- unique_combs[i,2]
  sp <- unique_combs[i,3]
  
  print(paste(i,site,sp,y)) 
  
  rmbl_sub <- rmbl[which(rmbl$species == sp & rmbl$plot == site & rmbl$year == y),]
  if(sum(rmbl_sub$floralcount != 0) == 0) next() #making sure not to try to estimate when there are not non-zero records... though I should have filtered them all out already
  doys <- rep(rmbl_sub$doy, rmbl_sub$floralcount)
  if(length(doys) < 4) next() #making sure phest has enough to work with
  phest <- weib.limit(head(doys,50))
  
  rmbl_formatted$station[i] <- as.character(site)
  rmbl_formatted$year[i] <- y
  rmbl_formatted$species[i] <- as.character(sp)
  rmbl_formatted$estimate[i] <- phest[[1]]
  rmbl_formatted$lower[i] <- phest[[2]]
  rmbl_formatted$upper[i] <- phest[[3]]
  
}

rmbl_formatted_cleaned <- rmbl_formatted[-which(is.na(rmbl_formatted$estimate) | is.na(rmbl_formatted$lower) | is.na(rmbl_formatted$upper)),]

# removing species/site combinations where there are fewer than 10 estimates, following Pearse et al 2017
# this is painfully convoluted programming, but it works ¯\_(ツ)_/¯
n_estimates <- plyr::count(rmbl_formatted_cleaned, c("station", "species"))
estimates_over_10 <- n_estimates[which(n_estimates$freq >= 10),]
sites_species <- paste(estimates_over_10$station, estimates_over_10$species)
rmbl_formatted_cleaned$sites_species <- paste(rmbl_formatted_cleaned$station, rmbl_formatted_cleaned$species)
rmbl_formatted_cleaned <- rmbl_formatted_cleaned[which(rmbl_formatted_cleaned$sites_species %in% sites_species),]
rmbl_formatted_cleaned <- rmbl_formatted_cleaned[,colnames(rmbl_formatted_cleaned) != "sites_species"]

write.csv(rmbl_formatted_cleaned, "clean_data/rmbl_ready_for_analysis.csv")

# cleanup
rm(rmbl, all_counts, just_counts, non_zeros, rmbl_na, rmbl_2016, rmbl_formatted, rmbl_formatted_cleaned)

