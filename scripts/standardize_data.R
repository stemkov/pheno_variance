# standardizing phenophase names and subsetting to just spring/onset events

korea <- fread("clean_data/korea_with_temp.csv", stringsAsFactors = FALSE)
japan <- fread("clean_data/japan_with_temp.csv", stringsAsFactors = FALSE)
manomet <- fread("clean_data/manomet_with_temp.csv", stringsAsFactors = FALSE)
rmbl <- fread("clean_data/rmbl_with_temp.csv", stringsAsFactors = FALSE)
nectar <- fread("clean_data/nectar_with_temp.csv", stringsAsFactors = FALSE)
#kivach <- read.csv("Edited Data/kivach/kivach_with_temp.csv", stringsAsFactors = FALSE)
rothamsted <- fread("clean_data/rothamsted_with_temp.csv", stringsAsFactors = FALSE)

all_data <- data.table(site = c(korea$station, japan$station, manomet$station, rmbl$station, nectar$station, rothamsted$site),
                       species = c(korea$species, japan$species, manomet$species, rmbl$species, nectar$species, rothamsted$species),
                       phenophase = c(korea$stage, japan$stage, rep("first_flight", nrow(manomet)), rep("first_flower", nrow(rmbl)), nectar$stage, rep("first_occurence", nrow(rothamsted))),
                       doy = c(korea$DOY, japan$DOY, manomet$first_estimate, rmbl$estimate, nectar$DOY, rothamsted$doy),
                       year = c(korea$year, japan$year, manomet$year, rmbl$year, nectar$year, rothamsted$year),
                       tmax = c(korea$tmax, japan$tmax, manomet$tmax, rmbl$tmax, nectar$tmax, rothamsted$tmax),
                       dataset = c(korea$dataset, japan$dataset, manomet$dataset, rmbl$dataset, nectar$dataset, rothamsted$dataset),
                       lat = c(korea$lat, japan$lat, manomet$lat, rmbl$lat, nectar$lat, rothamsted$lat),
                       lon = c(korea$lon, japan$lon, manomet$lon, rmbl$lon, nectar$lon, rothamsted$lon))


# excluding non-spring phenophases
# this removes only about 26% of the total data
# I'm assuming that "first gudding date" (FGD) from Japan is first leaf or bud
#sort(table(all_data$phenophase))

spring_phases <- c("onset of leaf unfolding", "opening", "firstleafbud", "firstleaf", "ffdbasket",
                   "onset of blooming", "1st occurrence", "fld", "first_flight", "first_flower",
                   "ffd", "FHD", "bud", "flower", "FOD", "first", "FFD", "FGD", "first_occurence")

all_data <- all_data[which(all_data$phenophase %in% spring_phases),]

# standardizing phase names into "first flower", "first leaf", and "first occurence" 
# the first two apply to plants, and the third is for all other taxa
# first occurence referes to the adult appearance

first_flower_phases <- c("ffdbasket", "opening", "onset of blooming", "first_flower", "ffd", "flower", "FFD")
first_leaf_phases <- c("onset of leaf unfolding", "firstleafbud", "firstleaf", "fld", "FGD", "bud")
first_occurence_phases <- c("1st occurrence", "first", "FOD", "FHD", "first_flight")

all_data[which(all_data$phenophase %in% first_flower_phases), "phenophase"] <- "first_flower"
all_data[which(all_data$phenophase %in% first_leaf_phases), "phenophase"] <- "first_leaf"
all_data[which(all_data$phenophase %in% first_occurence_phases), "phenophase"] <- "first_occurence"

# table(all_data$dataset)

# getting rid of abiotic "species" from Korea
abiotic_species <- c("frost", "ice", "snow", "snow accumulation", "freezing river")

`%!in%` <- Negate(`%in%`)

all_data <- all_data[which(all_data$species %!in% abiotic_species),]

write.csv(all_data, "clean_data/all_data_standardized.csv", row.names = FALSE)

# I added temp and standardized USSR data earlier, so I'm adding it here
# also Kivach is subsumed into USSR, so I'm getting rid of it
#all_data <- fread("Edited Data/all_data_standardized.csv", stringsAsFactors = FALSE)
#all_data <- all_data[dataset != "kivach",] # 218,901
ussr_data <- fread("clean_data/ussr_clean.csv", stringsAsFactors = FALSE) # 104,012
pep_data <- fread("clean_data/pep_clean.csv", stringsAsFactors = FALSE) # 136,312
all_data <- rbindlist(list(all_data, ussr_data, pep_data)) # 468,192
#all_data <- rbindlist(list(all_data, pep_data)) # 1,779,808

all_data[, site := paste(substr(dataset,1,1), site, sep="_"), by=.(site, dataset)]

write.csv(all_data, "clean_data/all_data_standardized.csv", row.names = FALSE)

