###
# Cleaning up PEP data for analysis
#
# BBCH is the phenology event code
# PEP_ID is station information
#
# Outcomes:
#
# subsetted to just early season stages (first leaves and first flowers)
#
# there are 2 commonly recorded stages for first leaves and three for first flowers, so...
# picked the most commonly used leaf and flower stage for each species
# These are mostly easy to call because of huge unevenness in the stage data almost all species, but there is one possibly problematic species:
# Secale cereale. flowering 60 and 61 are close in coverage. It's rye, so it's just closely monitored
# overall, this results in the loss of about 1% of data, so I'm not too worried
# 
#
###

#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

# bbch <- read.table("Raw Data/PEP725_DE_108_010/PEP725_BBCH.csv", sep=";", header=T)
# pep_id <- read.table("Raw Data/PEP725_DE_108_010/PEP725_DE_stations.csv", sep=";", header=T)
bbch <- fread("raw_data/PEP725_DE_108_010/PEP725_BBCH.csv", stringsAsFactors = FALSE)
pep_id <- fread("raw_data/PEP725_DE_108_010/PEP725_DE_stations.csv", stringsAsFactors = FALSE)

#setwd("/home/michael/Documents/Grad School/Research Projects/pheno_variance/raw_data/PEP725_records_Stemkovski")
root_dir <- getwd()
setwd(paste0(root_dir,"/raw_data/PEP725_records_Stemkovski"))

data <- fread(paste0(getwd(), "/pep725_Stemkovski_0.csv"), stringsAsFactors = FALSE)
files <- list.files()[-1]
file_numbers <- gsub("\\D","",files)
files <- files[order(as.numeric(file_numbers))]

for (i in files){
  to_add <- fread(i, stringsAsFactors = FALSE)
  data <- rbindlist(list(data, to_add))
}
rm(to_add)
setwd(root_dir)

data <- data[phase_id %in% c(10,11,60),] # same BBCH's as in the USSR data
data <- data[year >= 1958,]
n_sites_full <- length(unique(data$s_id))
#data[, location := do.call(paste,round(.SD, 1)), .SDcols=c("lat", "lon")] # coordinates are rounded down to 1 decimal point (~11km) for unique location - this actually only cuts out 23 out of 14423 sites
data[, location := do.call(paste,round(.SD, 0)), .SDcols=c("lat", "lon")]
n_sites_subset <- length(unique(data$location))

sites <- data[, .(records = .N, location = unique(location), lat = unique(lat), lon = unique(lon)), by=s_id]
sites <- sites[, .SD[which.max(records)], by=location] # some sites are at (or close to) the same coordinates, so I'm selecting the site that has the most records at each location - then I can really call them sites
data <- data[s_id %in% sites$s_id,]

`%!in%` <- Negate(`%in%`)

# doing a bit of pruning so that I have to rely less on the cleaning function later on
data <- data[, if (.N >= 10 & max(year) >= 2010) .SD, by = .(s_id, species, phase_id)] # 10 or more years and at least up to 2010
data <- data[species %!in% unique(data[cult_season != 0, species]),] # getting rid of agricultural cereals

data <- data[, phase_id:=as.character(phase_id)]
data[phase_id == 60, "phase_id"] <- "first_flower"
data[phase_id %in% c(10, 11), "phase_id"] <- "first_leaf"

# getting temp data
#data[, med_month := get.median.month(date), .(s_id, species, phase_id)]
n_combs <- uniqueN(data[,c("s_id", "species", "phase_id")])
data[, med_month := {cat(.GRP/n_combs*100,"%\n"); get.median.month(date)}, .(s_id, species, phase_id)]
n_combs <- uniqueN(data[,c("s_id", "med_month", "year")])
#data$tmax <- with(data, pbmapply(get.terraclim, lat=latitude, lon=longitude, year=year, month=med_month, check_water=F)) # you know you could make this faster if you do it by location/month combinations
data[, tmax := {cat(.GRP/n_combs*100,"%\n"); get.terraclim(lat, lon, year, med_month, check_water = F)}, .(s_id, med_month, year)]
n_combs <- uniqueN(data[,c("s_id")])
data[, seasonality := {cat(.GRP/n_combs*100,"% seasonality\n"); seasonality(mean(lat), mean(lon))}, by=s_id]
data[, mean_temp := {cat(.GRP/n_combs*100,"% mean_temp\n"); seasonality(mean(lat), mean(lon), tmax_mean=T)}, by=s_id]

#seasonality(trends[1,lat], trends[1,lon], both=T)

##### may be over water - check it (40.0666, 12.916)

write.csv(data, "clean_data/data_100km.csv", row.names = FALSE)

formatted_data <- data.table(site = data$s_id,
                             species = data$species,
                             phenophase = data$phase_id,
                             doy = data$day,
                             year = data$year,
                             tmax = data$tmax,
                             dataset = "pep",
                             lat = data$lat,
                             lon = data$lon)

pep_sites <- data[, .(records = .N, location = unique(location), lat = unique(lat), lon = unique(lon), seasonality = mean(seasonality), mean_temp = mean(mean_temp)), by=s_id]
pep_sites <- data.table(site = pep_sites$s_id,
                        seasonality = pep_sites$seasonality,
                        mean_temp = pep_sites$mean_temp)

write.csv(formatted_data, "clean_data/pep_clean.csv", row.names = FALSE)
write.csv(pep_sites, "clean_data/pep_sites.csv", row.names = FALSE) # doing this so that I don't have to recalculate seasonality and mean_temp later

rm(data, formatted_data, pep_id)

