### calculating microclimate a month of the phenophase for all datasets 

sort(sapply(ls(),function(x){object.size(get(x))}))

global.counter <- function(reset = FALSE){
  if(!exists("global_counter")){
    global_counter <<- 0
    print(global_counter)
  } 
  if(reset){
    global_counter <<- 0
  } 
  if(global_counter >= 0){
    global_counter <<- global_counter + 1
    print(global_counter, replace=T)
  }
}

### load in all datasets - I already calculated PEP and USSR microclimate in data cleaning scripts
korea_data <- fread("clean_data/korea_ready_for_analysis.csv", stringsAsFactors = FALSE)
japan_data <- fread("clean_data/japan_ready_for_analysis.csv", stringsAsFactors = FALSE)
manomet_data <- fread("clean_data/manomet_ready_for_analysis.csv", stringsAsFactors = FALSE) #it's "almost" because it doesn't account for sampling effort... but we may never get this data and sampling effort is probably kept pretty constant at manomet
rmbl_data <- fread("clean_data/rmbl_ready_for_analysis.csv", stringsAsFactors = FALSE)
nectar_data <- fread("clean_data/nectar_ready_for_analysis.csv", stringsAsFactors = FALSE)
rothamsted_data <- fread("clean_data/rothamsted_clean.csv", stringsAsFactors = FALSE)

# separate coordinate info, if needed
korea_station_info <- read.csv("raw_data/StationInfo_removed_special.csv", stringsAsFactors = FALSE)
rothamsted_info <- as.data.frame(read_xlsx("raw_data/rothamsted/SitesYearsLocations.xlsx"))
rothamsted_info <- data.frame(site = rothamsted_info$`Site Name`,
                              lat = rothamsted_info$Latitude,
                              lon = rothamsted_info$Longitude)

### pulling out spatial data for all datasets
match.coords <- function(site_id, dataset="korea"){
  if(dataset == "korea"){
    this_station <- korea_station_info[which(korea_station_info$country == "Korea" & korea_station_info$KoreanID == site_id),c("latdec","londec")]
    lat <- as.numeric(this_station["latdec"])
    lon <- as.numeric(this_station["londec"])
  }
  if(dataset == "japan"){
    this_station <- korea_station_info[which(korea_station_info$country == "Japan" & korea_station_info$ID == site_id),c("latdec","londec")]
    lat <- as.numeric(this_station["latdec"])
    lon <- as.numeric(this_station["londec"])
  }
  if(dataset == "nectar"){
    this_station <- nectar_coords[site_id,]
    lat <- as.numeric(this_station["lat"])
    lon <- as.numeric(this_station["lon"])
  }
  if(dataset == "rothamsted"){
    this_station <- rothamsted_info[which(rothamsted_info$site == site_id),]
    lat <- as.numeric(this_station["lat"])
    lon <- as.numeric(this_station["lon"])
  }
  results <- c(lat = lat, lon = lon)
  return(results)
}

# inserting and/or standardizing coordinates
korea_coords <- mapply(match.coords, site_id = korea_data$station, dataset = "korea")
korea_coords <- t(korea_coords)
korea_data$lat <- korea_coords[,1]
korea_data$lon <- korea_coords[,2]
korea_data$dataset <- "korea"
korea_data[which(korea_data$Station == 102),"lat"] <- 37.948798 # the point for site 102 is slightly off of the island - I manually picked the point in the middle of the small island
korea_data[which(korea_data$Station == 102),"lon"] <- 124.674469
apply(korea_data, 2, function(x) sum(is.na(x)))
korea_data <- korea_data[-which(korea_data$station == 164),] # station 164 doesn't have coordinates associated with it, so I'm just removing it
apply(korea_data, 2, function(x) sum(is.na(x))) # no more NAs

japan_codes <- as.data.frame(read_excel("raw_data/Japan Weather/StationNos.xls")) # the station numbers don't match the numbers in the station info spreadsheet. Luckily there's this link
japan_data$pheno_station <- sapply(japan_data$station, function(x) japan_codes[which(japan_codes$`Weather #` == x),1]) # matching up station number station number in coordinates spreadsheet
japan_coords <- mapply(match.coords, site_id = japan_data$pheno_station, dataset = "japan")
japan_coords <- t(japan_coords)
japan_data$lat <- japan_coords[,1]
japan_data$lon <- japan_coords[,2]
japan_data$dataset <- "japan"
japan_data[which(japan_data$pheno_station == 12),"lat"] <- 41.865184 # this point was very close to the coast and needed to be moved slightly inland
japan_data[which(japan_data$pheno_station == 12),"lon"] <- 140.130728
japan_data <- japan_data[-which(is.na(japan_data$DOY)),]
apply(japan_data, 2, function(x) sum(is.na(x))) # no NAs

# manomet & RMBL sites are just one location
manomet_data$lat <- 41.919934 # taken from google maps
manomet_data$lon <- -70.542633 # taken from google maps
manomet_data$dataset <- "manomet"

rmbl_data$lat <- 38.955942 # taken from google maps
rmbl_data$lon <- -106.988305 # taken from google maps
rmbl_data$dataset <- "rmbl"

# kivach_data$lat <- 62.276376 # taken from google maps
# kivach_data$lon <- 33.978323 # taken from google maps
# kivach_data$dataset <- "kivach"

rothamsted_coords <- mapply(match.coords, site_id = rothamsted_data$site, dataset = "rothamsted")
table(rothamsted_data$site)
rothamsted_coords <- t(rothamsted_coords)
rothamsted_data$lat <- rothamsted_coords[,1]
rothamsted_data$lon <- rothamsted_coords[,2]

#colnames(pep_data) <- gsub("long", "lon", colnames(pep_data)) # standardizing name
#pep_data$dataset <- "pep"

# coordinated obtained from the supplement for Cook et al. 2012 (Ecosystems) - these are coords from the nearest weather stations
nectar_coords <- data.frame(lat = rep(NA,9), lon = rep(NA,9))
rownames(nectar_coords) <- c("fitter", "harvard", "hubbard", "konzo", "luquillo", "mikesell", "niwot", "sevcore", "sevtrans")
nectar_coords["fitter",] <- c(51.77, 1.27) # This point is in the ocean
nectar_coords["fitter",] <- c(51.806760, 1.206607) # I found the nearest point on land manually
nectar_coords["harvard",] <- c(42.27, -71.87)
nectar_coords["hubbard",] <- c(44.03, -71.14)
nectar_coords["konza",] <- c(39.20, -96.58)
nectar_coords["luquillo",] <- c(18.43, -66.01)
nectar_coords["mikesell",] <- c(41.52, -84.15)
nectar_coords["niwot",] <- c(40.032136, -105.535874)
nectar_coords["sevcore",] <- c(32.25, -107.75)
nectar_coords["sevtrans",] <- c(32.25, -107.75)
nectar_coords_all <- mapply(match.coords, site_id = nectar_data$station, dataset = "nectar")
nectar_coords_all <- t(nectar_coords_all)
nectar_data$lat <- nectar_coords_all[,1]
nectar_data$lon <- nectar_coords_all[,2]
nectar_data$dataset <- "nectar"

# month_lag allows you to look at the var value in months before month of doy (if doy is provided)
# check_water allows you to check if the lat and lon are over an ocean = returns TRUE if over water, FALSE if not - overrides var output
get.terraclim <- function(lat, lon, year=2010, month=1, doy=NA, var="tmax", month_lag=0, check_water=FALSE, counter=FALSE, full_year=FALSE){
  if(counter) global.counter()
  if(sum(is.na(c(lat,lon,year,month,var)) > 0)){
    warning("NA in one or more of the arguments")
    return(NA)
  }
  if(year < 1958 | year > 2018){
    warning("Tried to get data before 1958 or after 2018")
    return(NA)
  }
  if(!is.na(doy)){
    # month_of_doy <- month(as.Date(104, paste(year,01,01,sep="-"))) # I had this earlier... this might have been a major oversight :O
    month_of_doy <- month(as.Date(doy, paste(year,01,01,sep="-")))
    month <- month_of_doy - month_lag # month overwritten by doy
  }
  file_loc <- paste0("raw_data/terraclimate/TerraClimate_",var,"_",year,".nc")
  nc <- nc_open(file_loc)
  lon_extent <- ncvar_get(nc, "lon")
  lat_extent <- ncvar_get(nc, "lat")
  #flat = which(!is.na(match(abs(lat_extent - lat) < 1/48, 1)))
  #flon = which(!is.na(match(abs(lon_extent - lon) < 1/48, 1)))
  flat <- which.min(abs(lat_extent - lat))
  flon <- which.min(abs(lon_extent - lon))
  #print(flat)
  #print(flon)
  start <- c(flon, flat, month)
  count <- c(1,1,1)
  if(full_year){
    start <- c(flon, flat, 1)
    count <- c(1,1,12)
  }
  output <- ncvar_get(nc, varid = "tmax", start = start, count = count, raw_datavals = F)
  nc_close(nc)
  #print(output)
  if(any(is.na(output))){
    if(check_water) return(TRUE)
    warning(paste("The point at (", lat, ", ", lon, ") may be over water. Consider using check_water=TRUE", sep=""))
  } 
  if(check_water) return(FALSE)
  return(output)
}


all_coords <- unique(rbind(nectar_data[,c("lat", "lon", "dataset")],
                           korea_data[,c("lat", "lon", "dataset")],
                           japan_data[,c("lat", "lon", "dataset")],
                           rmbl_data[,c("lat", "lon", "dataset")],
                           manomet_data[,c("lat", "lon", "dataset")],
                           #kivach_data[,c("lat", "lon", "dataset")],
                           rothamsted_data[,c("lat", "lon", "dataset")]))
all_coords <- as.data.frame(all_coords)
if(length(which(is.na(all_coords$lat) | is.na(all_coords$lon))) > 0) all_coords <- all_coords[-which(is.na(all_coords$lat) | is.na(all_coords$lon)),] #removes NAs if there are any
all_coords$over_water <- with(all_coords, mapply(get.terraclim, lat=lat, lon=lon, check_water=T))
all_coords[which(all_coords$over_water),] # none :)
all_coords[which(is.na(all_coords$lat) | is.na(all_coords$lon)),] # none :)

# global.counter(TRUE)
# korea_data$tmax <- with(korea_data, mapply(get.terraclim, lat=lat, lon=lon, year=year, doy=DOY, check_water=F, counter=T))
# japan_data$tmax <- with(japan_data, mapply(get.terraclim, lat=lat, lon=lon, year=year, doy=DOY, check_water=F, counter=T))
# rmbl_data$tmax <- with(rmbl_data, mapply(get.terraclim, lat=lat, lon=lon, year=year, doy=estimate, check_water=F, counter=T))
# nectar_data$tmax <- with(nectar_data, mapply(get.terraclim, lat=lat, lon=lon, year=year, doy=DOY, check_water=F, counter=T))
# manomet_data$tmax <- with(manomet_data, mapply(get.terraclim, lat=lat, lon=lon, year=year, doy=first_estimate, check_water=F, counter=T))
# #kivach_data$tmax <- with(kivach_data, mapply(get.terraclim, lat=lat, lon=lon, year=year(ymd(Date)), doy=doy, check_water=F, counter=T))
# rothamsted_data$tmax <- with(rothamsted_data, mapply(get.terraclim, lat=lat, lon=lon, year=year, doy=doy, check_water=F, counter=T))

# this should work faster than the above thanks to data.table
korea_data[, tmax := get.terraclim(lat=lat, lon=lon, year=year, doy=DOY, check_water=F),
           by = .(lat, lon, year, DOY)]
write.csv(korea_data, "clean_data/korea_with_temp.csv", row.names = FALSE)
rm(korea_data)

japan_data[, tmax := get.terraclim(lat=lat, lon=lon, year=year, doy=DOY, check_water=F),
           by = .(lat, lon, year, DOY)]
write.csv(japan_data, "clean_data/japan_with_temp.csv", row.names = FALSE)
rm(japan_data)

rmbl_data[, tmax := get.terraclim(lat=lat, lon=lon, year=year, doy=estimate, check_water=F),
           by = .(lat, lon, year, estimate)]
write.csv(rmbl_data, "clean_data/rmbl_with_temp.csv", row.names = FALSE)
rm(rmbl_data)

nectar_data[, tmax := get.terraclim(lat=lat, lon=lon, year=year, doy=DOY, check_water=F),
          by = .(lat, lon, year, DOY)]
write.csv(nectar_data, "clean_data/nectar_with_temp.csv", row.names = FALSE)
rm(nectar_data)

manomet_data[, tmax := get.terraclim(lat=lat, lon=lon, year=year, doy=first_estimate, check_water=F),
           by = .(lat, lon, year, first_estimate)]
write.csv(manomet_data, "clean_data/manomet_with_temp.csv", row.names = FALSE)
rm(manomet_data)

rothamsted_data[, tmax := get.terraclim(lat=lat, lon=lon, year=year, doy=doy, check_water=F),
           by = .(lat, lon, year, doy)]
write.csv(rothamsted_data, "clean_data/rothamsted_with_temp.csv", row.names = FALSE)
rm(rothamsted_data)

gc()

### The above gets temperature for the month in which a phenophase occured
# This might skew results
# So here I'm getting temperature for the typical (median) month of the phenophase across all years

korea_data <- fread("clean_data/korea_with_temp.csv")
japan_data <- fread("clean_data/japan_with_temp.csv")
manomet_data <- fread("clean_data/manomet_with_temp.csv")
rmbl_data <- fread("clean_data/rmbl_with_temp.csv")
nectar_data <- fread("clean_data/nectar_with_temp.csv")
#kivach_data <- fread("Edited Data/kivach/kivach_with_temp.csv")
rothamsted_data <- fread("clean_data/rothamsted_with_temp.csv")

korea_data$tmax_dynamic <- korea_data$tmax
japan_data$tmax_dynamic <- japan_data$tmax
manomet_data$tmax_dynamic <- manomet_data$tmax
rmbl_data$tmax_dynamic <- rmbl_data$tmax
nectar_data$tmax_dynamic <- nectar_data$tmax
#kivach_data$tmax_dynamic <- kivach_data$tmax
rothamsted_data$tmax_dynamic <- rothamsted_data$tmax

get.median.month <- function(dates){
  dates <- ymd(dates)
  months <- month(dates)
  return(floor(median(months, na.rm = T)))
}

#korea_data[1:100, get.median.month(Phenology), .(Station, Species, Stage)]
korea_data[, med_month := get.median.month(phenology), .(station, species, stage)]
japan_data[, med_month := get.median.month(phenology), .(station, species, stage)]
rmbl_data[, med_month := get.median.month(as.Date(round(estimate), paste(year,01,01,sep="-"))), .(station, species)]
nectar_data[, med_month := get.median.month(phenology), .(station, species, stage)]
manomet_data$first_estimate_date <- as.Date(round(manomet_data$first_estimate), paste(manomet_data$year,01,01,sep="-")) # I was having trouble getting 
manomet_data[, med_month := get.median.month(first_estimate_date), .(species)]
# kivach_data[, Date := mdy(Date)]
# kivach_data[, med_month := get.median.month(Date), .(Species, Event)]
rothamsted_data[, med_month := get.median.month(as.Date(round(doy), paste(year,01,01,sep="-"))), .(site, species)]

korea_data$tmax <- with(korea_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year, month=med_month, check_water=F))
japan_data$tmax <- with(japan_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year, month=med_month, check_water=F))
rmbl_data$tmax <- with(rmbl_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year, month=med_month, check_water=F))
nectar_data$tmax <- with(nectar_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year, month=med_month, check_water=F))
manomet_data$tmax <- with(manomet_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year, month=med_month, check_water=F))
#kivach_data$tmax <- with(kivach_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year(ymd(Date)), month=med_month, check_water=F))
rothamsted_data$tmax <- with(rothamsted_data, pbmapply(get.terraclim, lat=lat, lon=lon, year=year, month=med_month, check_water=F))

write.csv(korea_data, "clean_data/korea_with_temp.csv", row.names = FALSE)
write.csv(japan_data, "clean_data/japan_with_temp.csv", row.names = FALSE)
write.csv(manomet_data, "clean_data/manomet_with_temp.csv", row.names = FALSE)
write.csv(rmbl_data, "clean_data/rmbl_with_temp.csv", row.names = FALSE)
write.csv(nectar_data, "clean_data/nectar_with_temp.csv", row.names = FALSE)
write.csv(rothamsted_data, "clean_data/rothamsted_with_temp.csv", row.names = FALSE)

# cleanup
rm(japan_data, korea_data, japan_coords, nectar_data, korea_coords, rothamsted_data, rmbl_data, rothamsted_coords)
gc()
