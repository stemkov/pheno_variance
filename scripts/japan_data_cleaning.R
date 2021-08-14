###
# Cleaning up Japan data for analysis
#
# Phenology codes:
# FFD: First flowering date
# FHD: First heard date (the date an insect or bird was first heard)
# FOD: First observation date (the date an insect or bird was first seen)
# LCD: Leaf coloring date
# LDD: Leaf drop date
# PFD: Peak flowering date
# FGD: *** uncertain. First gudding date. wtf is that
#
# Outcomes:
#
# Got rid of first column, which was just redundant sp./event code
# Formatted into columns: station, year, species, stage, phenology, DOY
# Assuming the phenology dates are represented as "MMDD". Need to check on this
# Also assuming that "0" phenology dates are actually "NA" ß
# I need to figure out what the "event" codes mean
#
# Some phenphases are in the winter, and dates are reported straddling then new year but without the year specifies.
# To deal with this, I'm just keeping the year as stated and converting pre-January dates into negative DOYs when there is evidence of split
#
###

#setwd("/Users/Michael/Documents/Grad School/Research Projects/Phenology variability")
#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")


# old version that doesn't take into account winter phenology problem
phen.format.japan <- function(data_row){
  #data_row <- japan_expanded[1,]
  phen <- data_row[5]
  if (is.na(phen)){
    data_row[5] <- NA
  } else if (as.numeric(phen) == 0){
    data_row[5] <- NA
  } else {
    year <- data_row[2]
    if (nchar(phen) == 4){
      month <- substr(phen,1,2)
      day <- substr(phen,3,4)
    } else if (nchar(phen) == 3) {
      month <- substr(phen,1,1)
      day <- substr(phen,2,3)
    } else {
      stop("Phenology is in an unknown format")
    }
    date <- ymd(paste(year,month,day,sep="-"))
    data_row[5] <- as.character(date)
  }

  return(data_row)
}

# new version that takes into account winter phenology

japan <- read.csv("raw_data/JMA_ALL_2011.csv")
unique(japan$event)
japan <- japan[,-c(1,5,65:ncol(japan))]
colnames(japan) <- gsub("X","",colnames(japan))

years <- colnames(japan)[4:ncol(japan)]

station_sp_event <- mapply(function(x,y,z) return(paste(x,y,z,sep="_")) , x = japan$Station, y = japan$sp.name, z = japan$event)
unique_combs <- unique(station_sp_event)

# station, year, species, stage, phenology, DOY
full_length <- length(years) * nrow(japan)
japan_expanded <- data.frame(station = rep(NA,full_length),
                             year = rep(NA,full_length),
                             species = rep(NA,full_length),
                             stage = rep(NA,full_length),
                             phenology = rep(NA,full_length),
                             DOY = rep(NA,full_length))

top <- 1

for (i in 1:length(years)){
  # i <- 1 
  # i <- 2
  # i <- 3
  year <- years[i]
  year_block <- japan[,c(1:3,(3+i))]
  colnames(year_block) <- c("species","stage","station","phenology")
  colnames(year_block)
  year_block$year <- year
  insert_range <- top:(top + nrow(year_block)-1)
  
  japan_expanded$station[insert_range] <- year_block$station
  japan_expanded$year[insert_range] <- year_block$year
  japan_expanded$species[insert_range] <- as.character(year_block$species)
  japan_expanded$phenology[insert_range] <- year_block$phenology
  japan_expanded$stage[insert_range] <- as.character(year_block$stage)
  
  top <- top + nrow(year_block)
}

head(japan_expanded,10)
dim(japan_expanded)
#hist(japan_expanded$phenology[which(japan_expanded$phenology != 0)], 365) # quite a few happening in Decemeber and January

#t(apply(japan_expanded[1:10,],1,phen.format.japan))
japan_formatted <- t(apply(japan_expanded,1,phen.format.japan))
#head(japan_formatted,10)
#dim(japan_formatted)

# to make it easier to work with
jdf <- as.data.frame(japan_formatted, stringsAsFactors=FALSE)
jdt <- data.table(jdf)

# identifying species/sites/stages that have both January and December phenology
is.split.phen <- function(x, month_range=1){
  x <- ymd(as.character(x))
  months <- lubridate::month(x)
  if(month_range == 1) if(12 %in% months & 1 %in% months) return(TRUE) #checks if there are January AND December phenologies
  if(month_range == 2) if((12 %in% months | 11 %in% months) & (1 %in% months | 2 %in% months)) return(TRUE) #checks if there are (January OR February) AND (December OR November) phenologies
  
  return(FALSE)
}

insert.split.phen <- function(x){
  split_sub <- split_phens_trues[which(split_phens_trues$station == x[1] & split_phens_trues$species == x[2] & split_phens_trues$stage == x[3]),]
  if(nrow(split_sub) == 0) return(FALSE)
  return(TRUE)
}

negative.doys <- function(x){
  if(is.na(x)) return(NA)
  date <- ymd(x)
  year <- year(x)
  month <- month(date)
  if(month <= 9) return(yday(date))
  end_of_year <- ymd(paste(year, 12, 31, sep="-"))
  return(as.numeric(date - end_of_year))
}

split_phens_1 <- jdt[, .(split_phen = is.split.phen(phenology, month_range=1)), by=.(station, species, stage)]
split_phens_2 <- jdt[, .(split_phen = is.split.phen(phenology, month_range=2)), by=.(station, species, stage)]
split_phens_1[which(split_phen),] 
split_phens_2[which(split_phen),] # thank god there are not problems beyond the January/December split

split_phens <- as.data.frame(apply(split_phens_1,2,as.character), stringsAsFactors=FALSE) # R just loves turning things into factors doesn't it
split_phens$split_phen <- gsub(" TRUE", "TRUE", split_phens$split_phen)
split_phens_trues <- split_phens[which(as.logical(split_phens$split_phen)),]

split_rows <- which(as.logical(apply(jdf[,c("station", "species", "stage")], 1, insert.split.phen))) # this is a dumb way to do it but I couldn't think of a better way
jdf$DOY <- yday(jdf$phenology)
jdf[split_rows,"DOY"] <- sapply(jdf[split_rows,"phenology"], negative.doys) # jdf now has phenology values corrected for split years

# I'm just keeping the year as stated and converting pre-January dates into negative DOYs
# A bit of justification for negative DOYs:
# 1. In areas with a winter snowpack, it makes most sense to look at "water years", which start at the first snow event
# 2. The calendar new year is after the winter solstice, so daylenghts are already increasing for the very early, late December phenology

japan_formatted <- jdf

j <- data.frame(station = as.factor(japan_formatted[,"station"]),
                    year = as.numeric(japan_formatted[,"year"]),
                    species = as.factor(japan_formatted[,"species"]),
                    stage = as.factor(japan_formatted[,"stage"]),
                    phenology = as.character(japan_formatted[,"phenology"]),
                    DOY = as.numeric(japan_formatted[,"DOY"]))


write.csv(j,"clean_data/japan_ready_for_analysis.csv")

#cleanup
rm(jdt, japan_formatted, jdf, j, japan_expanded, japan)
