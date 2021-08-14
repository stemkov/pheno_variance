### Cleaning USSR data

#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

data <- fread("raw_data/ussr/phenology.csv", stringsAsFactors = FALSE)
taxonomy <- fread("raw_data/ussr/taxonomy.csv", stringsAsFactors = FALSE)
sites <- fread("raw_data/ussr/studysites.csv", stringsAsFactors = FALSE)

# sort(table(data$eventtype))
# sort(table(data[taxonidentifier %in% taxonomy[class == "Aves", taxonidentifier],]$eventtype))

# just the taxa that I want to keep
taxonomy <- taxonomy[class %in% c("Insecta","Aves") | kindgom == "Plantae",]
taxonomy <- taxonomy[taxonomiclevel == "Species",]
weird_taxa <- taxonomy$taxon[which(grepl("s\\.l\\.|;|s\\. str\\.|Ã", taxonomy$taxon))] # I want to get rid of these

data <- data[taxonidentifier %in% taxonomy$taxonidentifier,]
data <- data[eventtype %in% c("onset of blooming", "1st occurrence", "onset of leaf unfolding", "start of vegetation", "")]
data <- data[quality == "OK",]
data <- data[taxon %!in% weird_taxa,]
data <- data[year >= 1958,]
setkey(sites, studysite)
setkey(data, studysite)
data <- sites[data]

sites$records <- sapply(sites$studysite, function(x) nrow(data[studysite == x,"dayofyear"]))
sites$location <- paste(sites$latitude, sites$longitude)
sites <- sites[sites[, .I[records == max(records)], by=location]$V1] # some sites are at the same coordinates, so I'm selecting the site that has the most records at each location - then I can really call them sites
paste("keeping", length(unique(sites$studysite)), "site of original", length(unique(data$studysite)))
data <- data[studysite %in% sites$studysite,] # gets rid of ~30,000 record, which isn't actually huge in the context of this dataset

get.median.month <- function(dates){
  dates <- ymd(dates)
  months <- month(dates)
  return(floor(median(months, na.rm = T)))
}

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

seasonality <- function(lat, lon, counter=F, tmax_mean=F, both=F){
  if(counter) global.counter()
  years <- 1958:2018
  seasonality_rec <- rep(NA,length(years))
  for(i in 1:length(years)){
    y <- years[i]
    temp <- get.terraclim(lat, lon, y, full_year = T)
    seasonality <- diff(c(min(temp), max(temp)))
    if(tmax_mean) seasonality <- mean(temp)
    seasonality_rec[i] <- seasonality
  }
  return(mean(seasonality_rec))
}

# I'm getting temp data and standardizing here
data[, med_month := get.median.month(as.Date(round(dayofyear), paste(year,01,01,sep="-"))), .(studysite, taxon)]
#data$tmax <- with(data, pbmapply(get.terraclim, lat=latitude, lon=longitude, year=year, month=med_month, check_water=F)) # you know you could make this faster if you do it by location/month combinations
n_combs <- uniqueN(data[,c("studysite", "med_month", "year")])
#data$tmax <- with(data, pbmapply(get.terraclim, lat=latitude, lon=longitude, year=year, month=med_month, check_water=F)) # you know you could make this faster if you do it by location/month combinations
data[, tmax := {cat(.GRP/n_combs*100,"%\n"); get.terraclim(latitude, longitude, year, med_month, check_water = F)}, .(studysite, med_month, year)]
n_combs <- uniqueN(data[,c("studysite")])
data[, seasonality := {cat(.GRP/n_combs*100,"% seasonality\n"); seasonality(unique(latitude), unique(longitude))}, by=studysite]
data[, mean_temp := {cat(.GRP/n_combs*100,"% mean_temp\n"); seasonality(unique(latitude), unique(longitude), tmax_mean=T)}, by=studysite]

#all_data[which(all_data$phenophase %in% first_leaf_phases), "phenophase"]
data[eventtype == "onset of blooming", "eventtype"] <- "first_flower"
data[eventtype %in% c("onset of leaf unfolding", "start of vegetation"), "eventtype"] <- "first_leaf"
data[eventtype == "1st occurence", "eventtype"] <- "first_occurence"

formatted_data <- data.table(site = data$studysite,
                             species = data$taxon,
                             phenophase = data$eventtype,
                             doy = data$dayofyear,
                             year = data$year,
                             tmax = data$tmax,
                             dataset = "ussr",
                             lat = data$latitude,
                             lon = data$longitude)

ussr_sites <- data[, .(records = .N, lat = unique(latitude), lon = unique(longitude), seasonality = mean(seasonality), mean_temp = mean(mean_temp)), by=studysite]
ussr_sites <- data.table(site = ussr_sites$studysite,
                        seasonality = ussr_sites$seasonality,
                        mean_temp = ussr_sites$mean_temp)

write.csv(formatted_data, "clean_data/ussr_clean.csv", row.names = FALSE)
write.csv(ussr_sites, "clean_data/ussr_sites.csv", row.names = FALSE) # doing this so that I don't have to recalculate seasonality and mean_temp later

rm(data, formatted_data)

