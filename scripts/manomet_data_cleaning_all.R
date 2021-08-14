###
# Cleaning up Manomet data for analysis
#
# This is fundementally different from the Korea/Japan data in that it has counts rather than "first sightings"
# So, it makes sense to use some kind of estimator, like phest
#
# Still needs doing: 
# Get sampling effort data and incorporate it to get weighted counts
#
# Outcomes:
# counts up the number of birds of a species caught at each location at a given date
# estimates first and last migrant of each species
# Formatted into columns: station, year, species, estimated onset, estimated end, CIs for the previous two, total counts, total unique observation dates
# 
#
###

#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

read.excel <- function(...) return(as.data.frame(read_excel(...)))

read.excel.allsheets <- function(filename,col_names,col_types=NULL) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) read.excel(filename, sheet = X, col_names=col_names, col_types=col_types))
  names(x) <- sheets
  return(x)
}

###loading in all data sheets

files <- list.files("raw_data/Manomet")
first_file <- read.excel.allsheets(paste("raw_data/Manomet/",files[1],sep=""),col_names = T,col_types = c("guess","numeric","numeric","numeric","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","date","skip","text","text","numeric","numeric","guess","guess","guess","guess","guess"))
manomet_all <- as.data.frame(matrix(NA,ncol=ncol(first_file[[1]]),nrow=0))
colnames(manoment_all) <- colnames(first_file[[1]])

# I excluded the time column because it wasn't loading in properly and I was just over it
for(f in files){
  file <- read.excel.allsheets(paste("raw_data/Manomet/",f,sep=""), col_names = T, col_types = c("guess","numeric","numeric","numeric","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","date","skip","text","text","numeric","numeric","guess","guess","guess","guess","guess"))
  sheets <- 1:length(file)
  for(s in sheets){
    this_sheet <- file[[s]]
    manomet_all <- rbind(manomet_all,this_sheet)
  }
}

get.counts <- function(x){
  #x <- manomet[which(manomet$Letter == "BGGN"),]
  if (length(unique(x$Letter)) != 1) stop("get.counts only accepts one species at a time")
  dates <- unique(x$Date)
  counts <- rep(NA,length(dates))
  for(i in 1:length(dates)){
    counts[i] <- nrow(x[which(x$Date == dates[i]),])
  }
  return(list(date = ymd(dates), count = counts, doy = yday(dates)))
}

years <- sort(unique(year(manomet_all$Date)))
species <- sort(unique(manomet_all$Letter))
stations <- sort(unique(manomet_all$Location))

# estimating phenologies using phest

full_length <- length(stations) * length(years) * length(species)

manomet_formatted <- data.frame(station = rep(NA,full_length),
                                year = rep(NA,full_length),
                                species = rep(NA,full_length),
                                first_estimate = rep(NA,full_length),
                                first_lower = rep(NA,full_length),
                                first_upper = rep(NA,full_length),
                                last_estimate = rep(NA,full_length),
                                last_lower = rep(NA,full_length),
                                last_upper = rep(NA,full_length),
                                count = rep(NA,full_length),
                                dates = rep(NA,full_length))
ticker <- 1
for (y in years){
  for(s in stations){
    for(sp in species){
      #y <- years[1]
      #s <- stations[1]
      #sp <- species[4]
      #sp <- "BADE"
      m <- manomet_all[which(year(manomet_all$Date) == y & manomet_all$Location == s & manomet_all$Letter == sp),]
      if(nrow(m) == 0){
        next
      }
      sp_counts <- get.counts(m)
      doy_counts <- sort(rep(sp_counts$doy,sp_counts$count))
      doy_counts_head <- head(doy_counts,50)
      doy_counts_tail <- tail(doy_counts,50)
      if(length(doy_counts) < 5) {
        first_estimate <- list(NA,NA,NA)
        last_estimate <- list(NA,NA,NA)
      } else {
        first_estimate <- weib.limit(doy_counts_head)
        last_estimate <- weib.limit(doy_counts_tail,upper=T)
      }
      
      manomet_formatted$year[ticker] <- y
      manomet_formatted$station[ticker] <- s
      manomet_formatted$species[ticker] <- sp
      manomet_formatted$first_estimate[ticker] <- first_estimate[[1]]
      manomet_formatted$first_lower[ticker] <- first_estimate[[2]]
      manomet_formatted$first_upper[ticker] <- first_estimate[[3]]
      manomet_formatted$last_estimate[ticker] <- last_estimate[[1]]
      manomet_formatted$last_lower[ticker] <- last_estimate[[3]]
      manomet_formatted$last_upper[ticker] <- last_estimate[[2]]
      manomet_formatted$count[ticker] <- sum(sp_counts$count)
      manomet_formatted$dates[ticker] <- length(sp_counts$count)
      
      ticker <- ticker + 1
      print(paste(y,"---",sp))
    }
  }
}

manomet_formatted <- manomet_formatted[1:ticker,]

# It's "almost" ready for analysis because I still need to get counts weighted by sampling effort
write.csv(manomet_formatted,"clean_data/manomet_all_almost_ready_for_analysis.csv")

manomet_formatted <- read.csv("clean_data/manomet_all_almost_ready_for_analysis.csv", stringsAsFactors = F)
head(manomet_formatted)

alpha_codes <- read.csv("raw_data/Bird_codes.csv", stringsAsFactors = F)

manomet_formatted$species <- mgsub(manomet_formatted$species, alpha_codes$SPEC, alpha_codes$SCINAME)
table(manomet_formatted$species) # problematic codes: BADE, BOBW, ETTI, FLIN, NSTS, SSTS, YWAR, WPWI
# I asked Sam Roberts about these codes and made the below adjustments
problem_codes <- c("BADE", "BOBW", "ETTI", "FLIN", "NSTS", "SSTS", "YWAR", "WPWI")
fixed_codes <- c("band destroyed", "Colinus virginianus", "Baeolophus bicolor", "Colaptes auratus", "Ammodramus nelsoni", "Ammodramus caudacutus", "Setophaga petechia", "Caprimulgus vociferus")
manomet_formatted$species <- mgsub(manomet_formatted$species, problem_codes, fixed_codes)

write.csv(manomet_formatted,"clean_data/manomet_ready_for_analysis.csv", row.names = F)

# cleanup
rm(manomet_all, first_file, file, this_sheet, alpha_codes)
