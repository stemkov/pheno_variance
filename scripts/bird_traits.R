# bird traits

library(plyr)

all_data <- fread("clean_data/all_data_standardized.csv")

# name reassignments to avoid problems
#problem_species <- read.csv("Edited Data/problem_species.csv", stringsAsFactors = FALSE)
problem_birds <- read.csv("raw_data/problem_birds.csv", stringsAsFactors = FALSE)
#all_data$species_synon <- mapvalues(all_data$species, problem_species$original_species, problem_species$fixed_species)
all_data$species_synon <- mapvalues(all_data$species, problem_birds$manomet_species, problem_birds$elton_species)

# limitting analysis to just 1958-2016, because that's the range of TerraClimate
all_data <- all_data[-which(is.na(all_data$tmax)),]

bird_species <- unique(all_data[which(all_data$dataset %in% c("manomet", "korea", "japan", "ussr")),c("species","species_synon")]) # species from datasets with potentially birds
bird_species$species_synon_clean <- sapply(bird_species$species_synon, name.cleanup)

# species <- fread("clean_data/species.csv", stringsAsFactors = F)
# bird_species <- species[tax_group == "bird",]
# bird_species$species_synon <- mapvalues(bird_species$species, problem_birds$manomet_species, problem_birds$elton_species)
# bird_species$species_synon_clean <- sapply(bird_species$species_synon, name.cleanup)

bird_trait_data <- fread("raw_data/traits/BirdFuncDat.txt")
write.csv(bird_trait_data, "raw_data/Eltontraits.csv", row.names = FALSE) #just to look at it more easily outside of R

bird_traits <- bird_trait_data[Scientific %in% bird_species$species | Scientific %in% bird_species$species_synon | Scientific %in% bird_species$species_synon_clean,]
bird_traits <- bird_traits[, c("Scientific","Diet-5Cat", "BodyMass-Value")] # I'm just pulling out mass & diet category

# paste0(round(nrow(bird_traits)/nrow(bird_species), 2)*100,"% of species have Elton trait data") #It's supposed to be 100% if you have only birds in the species list (Manomet)
# bird_species[which(bird_species$species_synon_clean %!in% bird_trait_data$Scientific)] # Ok, I got all species! (when looking at just Manomet - I threw in all species including non-birds into the species list for korea, japan, and kivach)

bird_species_traits <- data.frame(species = as.character(bird_species$species),
                                  species_synon = as.character(bird_species$species_synon),
                                  species_synon_clean = as.character(bird_species$species_synon_clean))
bird_species_traits$diet <- sapply(bird_species_traits$species_synon_clean, function(x) ifelse(x %!in% bird_traits$Scientific, NA , as.character(bird_traits[Scientific == x, "Diet-5Cat"])))
bird_species_traits$mass <- sapply(bird_species_traits$species_synon_clean, function(x) ifelse(x %!in% bird_traits$Scientific, NA , as.character(bird_traits[Scientific == x, "BodyMass-Value"])))
bird_species_traits[,c("species_synon", "diet", "mass")]
bird_species_traits <- bird_species_traits[which(complete.cases(bird_species_traits)),]
bird_species_traits <- bird_species_traits[-which(duplicated(bird_species_traits$species_synon_clean)),] # deleting a few duplicates that arose from the name cleaning and synonym finding

write.csv(bird_species_traits, "clean_data/bird_traits.csv", row.names = FALSE) # done!

