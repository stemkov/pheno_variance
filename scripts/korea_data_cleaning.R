###
# Cleaning up data for analysis... this will take a while
#
# Outcomes:
#
# Combines all sheets into one
# stacks the many columns into the following: year, station, species, stage, phenology DOY
# substitute the korea characters with "missing"
# treats "missing" as same as NA, and removed these rows
# formats the dates into standard values
# * It's important to note that I'm making some assumptions about date formatting. In the raw excel files, there is "4. 1" and "4.1". I take these to mean "April 1" and " April 10", respectively
#
###

#setwd("/Users/Michael/Documents/Grad School/Research Projects/Phenology variability")
#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

# function borrowed from Jeromy Anglim on stack exchange 
read.excel.allsheets <- function(filename,col_names){
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) { readxl::read_excel(filename, sheet = X, col_names=col_names) })
  x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
}

# rounds numerics only
smart.round <- function(x,digits=2){
  if (is.numeric(x)){rounded <- round(x,digits)}
  else {rounded <- x}
  return(rounded)
}

# counts the number of occurances of a string in another string
g.occur <- function(pattern,x){
  sub_string <- gsub(pattern,"",x)
  difference <- nchar(x) - nchar(sub_string)
  return(difference)
}

# formats the phenology values
phen.format.korea <- function(data_row){
  phen <- data_row[5]
  if (is.na(phen)){
    data_row[5] <- NA
  } else if (g.occur("\\.",phen) == 1) { 
    year <- data_row[2]
    month <- gsub("\\..*$","",phen)
    day <- gsub("^*.\\.","",phen)
    # deals with the annoying thing where values are expanded
    if (nchar(day) > 2) {day <- substring(round(as.numeric(substr(day,1,3)),-1),1,2)} 
    date <- ymd(paste(year,month,day,sep="-"))
    data_row[5] <- as.character(date)
  } else if (g.occur("\\.",phen) == 2){
    date <- ymd(phen)
    data_row[5] <- as.character(date)
  } else {}
  
  return(data_row)
  
}

korea1 <- read.excel.allsheets("raw_data/Korea_part.1_species.xls",col_names=F)
korea2 <- read.excel.allsheets("raw_data/Korea_part.2_species.xls",col_names=F)

# bringing the two Koreas together in a peaceful reunification
one_korea <- c(korea1,korea2)
one_korea <- one_korea[which(names(one_korea) != "code")] # getting rid of "code" sheets
one_korea <- one_korea[order(names(one_korea))] # ordering the sheets. Is it too convoluted?
stations <- names(one_korea)


# getting the total number of rows and columns across all sheets
# and pulling out some lists
total_dims <- c(0,0)
max_dims <- c(0,0)
all_years <- c()
all_species <- c()
all_stages <- c()
for (i in 1:length(stations)){
  
  dims <- dim(one_korea[[i]])
  total_dims <- total_dims + dims
  
  if (dims[1] > max_dims[1]){max_dims[1] <- dims[1]}
  if (dims[2] > max_dims[2]){max_dims[2] <- dims[2]}
  
  st_years <- one_korea[[i]][3:dims[1],1]
  all_years <- as.numeric(append(all_years,st_years))
  
  st_species <- as.character(one_korea[[i]][1,2:dims[2]])
  # annoying thing in the data where the column names are supposed to span multiple columns but only exist in one.
  # this fixes that problem
  for (j in 1:length(st_species)){ if(is.na(one_korea[[i]][1,j+1])){ one_korea[[i]][1,j+1] <- one_korea[[i]][1,j] } }
  all_species <- append(all_species,st_species)
  
  st_stages <- one_korea[[i]][2,2:dims[2]]
  all_stages <- as.character(append(all_stages,st_stages))
  
}

years <- sort(unique(all_years))
species <- sort(unique(all_species))
stages <- sort(unique(all_stages))


# korea_mashed <- matrix(NA,nrow=max_dims[1],ncol=total_dims[2])
# korea_mash <- matrix(NA,nrow=max_dims[1],ncol=0)

korea_final <- data.frame(station=as.character(NA),year=as.numeric(NA),species=as.character(NA),stage=as.character(NA),phenology=as.numeric(NA))

korea_final <- matrix(NA,ncol=5,nrow=0)
colnames(korea_final) <- c("Station","Year","Species","Stage","Phenology")

# combining all of the sheets into one. Each sheet number is a station code ..... oh look at me writing this horrendous for loop before I knew anything :'D
for (i in stations){
   #i <- "100"
  this_stn <- one_korea[[i]]
  for (j in years){
     #j <- 1976
    for (k in species){
      #k <- "plume"
      for (l in stages){
        # l <- "First"
        # k <- "barn swallow"
        # i <- "100"
        # j <- "1974"
        this_yr_sp_sta <- this_stn[which(this_stn[,1] == j),which(this_stn[1,] == k & this_stn[2,] == l)]
        if (length(this_yr_sp_sta) > 0){
          if (is.na(this_yr_sp_sta)){
          } else if (this_yr_sp_sta == "결측") { this_yr_sp_sta <- "missing" }
          to_insert <- c(i,j,k,l,this_yr_sp_sta)
          korea_final <- rbind(korea_final,to_insert)
          } else {}
        }
    }
    print(c(i,j)) 
  }
}

write.csv(korea_final,"clean_data/korea_formatted.csv")

korea_formatted <- read.csv("clean_data/korea_formatted.csv")
korea_formatted <- korea_formatted[,-1]
korea_formatted[,5] <- as.character(korea_formatted[,5])

korea_formatted <- t(apply(korea_formatted,1,phen.format.korea))
k_f <- korea_formatted
#korea_formatted <- k_f

korea_formatted <- korea_formatted[-which(is.na(korea_formatted[,5]) | korea_formatted[,5] == "missing"),]
korea_formatted <- korea_formatted[-which(nchar(korea_formatted[,5]) != 10),]
korea_formatted[,"Stage"] <- gsub("First","first",korea_formatted[,"Stage"])
korea_formatted[,"Stage"] <- gsub("Last","last",korea_formatted[,"Stage"])
korea_formatted[,"Species"] <- gsub("Ginko leaves","ginko leaves",korea_formatted[,"Species"])
korea_formatted[,"Species"] <- gsub("Ginko","ginko leaves",korea_formatted[,"Species"])
korea_formatted[korea_formatted[,"Species"] == "maple","Species"] <- "maple leaves"
korea_formatted[,"Species"] <- gsub("Sakura","sakura",korea_formatted[,"Species"])

# head(korea_formatted)
# dim(korea_formatted)
# table(korea_formatted[,"Stage"])
# table(korea_formatted[,"Species"])

korea <- data.frame(station = as.factor(korea_formatted[,"Station"]), 
                    year = as.numeric(korea_formatted[,"Year"]), 
                    species = as.factor(korea_formatted[,"Species"]), 
                    stage = as.factor(korea_formatted[,"Stage"]), 
                    phenology = as.character(korea_formatted[,"Phenology"]),
                    DOY = yday(korea_formatted[,"Phenology"]))

write.csv(korea, "clean_data/korea_ready_for_analysis.csv")

# coming back to match up common names with latin names
korea <- read.csv("clean_data/korea_ready_for_analysis.csv", stringsAsFactors = F)
common_scientific <- read_excel("raw_data/List of speciesindex_Korea.xls")
common_scientific <- as.data.frame(common_scientific)[c(3:nrow(common_scientific)),c(3,4)]
colnames(common_scientific) <- c("common", "latin")
common_scientific # this is honestly so far from being matched up with the actual data that I'm just going to sub in the lat names manually
unique(korea$species)

common_latin <- data.frame(common = unique(korea$species), latin = unique(korea$species))
common_latin$latin <- c("Pieris rapae", "Hirundo rustica", "Oncotympana fuscata",
                        "Cuculus canorus", "Cosmos bipinnatus", "Rana nigromaculata",
                        "frost", "ice", "Pyrus communis",
                        "Alauda arvensis", "snow", "Orthetrum albistylum",
                        "Elaphe rufodorsata", "Rhododendron mucronulatum", "Robinia pseudoacacia",
                        "Prunus serrulata", "freezing river", "snow accumulation",
                        "Acer palmatum", "Forsythia koreana", "Acer palmatum",
                        "Prunus persica", "Prunus mume", "Ginkgo biloba",
                        "beach")
# I assumed that "plume" is "plum" and the latin name is "Prunus mume"
# I have no idea what "beach" is. So I'm just leaving it and assuming they didn't mean "beech"

korea$species <- mgsub(korea$species, as.character(common_latin$common), common_latin$latin)
unique(korea$species)
write.csv(korea, "clean_data/korea_ready_for_analysis.csv", row.names = FALSE)

#cleanup
rm(korea_formatted, k_f, one_korea, korea2, korea, korea1, korea_final)

