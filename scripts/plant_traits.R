# matching trait data to all species

library(plyr)

all_data <- fread("clean_data/all_data_standardized.csv")

# name reassignments to avoid problems
problem_species <- read.csv("raw_data/problem_species.csv", stringsAsFactors = FALSE) 
all_data$species <- mapvalues(all_data$species, problem_species$original_species, problem_species$fixed_species)


if(FALSE){
  name_candidates <- gnr_resolve(problem_species$original_species)
  gnr_resolve("Larix leptolepis")
  
  gnr_resolve("Diospyros kaki")
  
  gnr_resolve("Halictus conf")
  h_conf <- get_tsn("Halictus confusus", "scientific", accepted=TRUE)
  d_kaki <- get_tsn("Diospyros kaki", "scientific", accepted = T)
  
  tax_name("Platycodon grandiflorum", "family", accepted=TRUE)
  test <- get_tsn("Platycodon grandiflorum")
  test <- get_tsn("Halictus rubicundus")
  tax_name(test, "family")
  
  rotl_test_sp <- c("Platycodon grandiflorum", "Taraxacum officinale", "Diospyros kaki", "Tyto alba")
  resolved_names <- tnrs_match_names(rotl_test_sp)
  resolved_names
  my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
  plot(my_tree)
  
  bird_tree <- read.tree("Raw Data/phylogenies/bird_phylogeny/2012-03-04206D-MCC_trees/MCC_trees/Backbone_trees/BackboneEricsonMCC.trees")
  plant_tree <- read.tree("Raw Data/phylogenies/Vascular_Plants_rooted.dated.tre")
  plot(plant_tree)
  
  name_candidates
  ?get_tsn
  ?itis_acceptname
  ?classification
  summary(name_candidates$score)
  tax_name(problem_species$original_species)
}

# limitting analysis to just 1958-2016, because that's the range of TerraClimate
all_data <- all_data[-which(is.na(all_data$tmax)),]

### getting trait data
name.cleanup <- function(x, count=F){
  if(count) global.counter()
  x <- gsub("_", " ", x)
  x <- mgsub(x, c("\\(sp\\)", "\\(sp\\.\\)", "\\_spp", "spp\\.", " a\\. ", " h\\. "), c("sp\\.", "sp\\.", " sp\\.", "sp\\.", " ", " "))
  x <- gsub("sp$", "sp.", x)
  x <- gsub("spp$", "sp.", x)
  x <- gsub("NA$", "sp.", x)
  x <- gsub("/.*","",x) # removeds second specific epithet
  x <- gsub("sect\\..*","sp.",x) # simplifies sect. to sp.
  if(!grepl(" ", x)) x <- paste(x, "sp.")
  return(as.character(x))
}

get.kingdom <- function(species, dataset){
  if(dataset %in% c("rmbl", "nectar", "pep")) return("Plantae")
  if(dataset == "manomet") return("Animalia")
  tax_lookup <- tax_name(species, "kingdom")
  if(is.na(tax_lookup$kingdom)) return("unknown")
  return(tax_lookup$kingdom)
}

get.taxanom <- function(species, dataset){
  if(dataset %in% c("rmbl", "nectar", "pep")) return(c("Plantae", NA, NA, NA))
  if(dataset == "manomet") return(c("Animalia", "Chordata", "Aves", NA))
  tax_lookup <- tax_name(species, c("kingdom", "phylum", "class", "order"))
  if(is.na(tax_lookup$kingdom)) return(c("unknown", "unknown", "unknown", "unknown"))
  return(c(tax_lookup$kingdom, tax_lookup$phylum, tax_lookup$class, tax_lookup$order))
}


all_data$species <- sapply(all_data$species, name.cleanup)

all_species <- unique(all_data[which(all_data$dataset != "manomet") ,c("species", "dataset")]) # minus Manoment b/c I know those are bird
all_species_KJ <- all_species[which(all_species$dataset %in% c("korea", "japan", "ussr")),] # mixture of plant and others
all_species_RN <- all_species[which(all_species$dataset %in% c("rmbl", "nectar", "pep")),] # all plants

install.packages("curl", type= "source")
library(curl)

# looking up taxonomy, trying ITIS first, then NCBI
taxon_info_KJ1 <- tax_name(as.character(all_species_KJ$species)[1:400], c("kingdom", "phylum", "class", "order"), accepted=TRUE, ask=FALSE)
taxon_info_KJ2 <- tax_name(as.character(all_species_KJ$species)[401:800], c("kingdom", "phylum", "class", "order"), accepted=TRUE, ask=FALSE)
taxon_info_KJ3 <- tax_name(as.character(all_species_KJ$species)[801:length(all_species_KJ$species)], c("kingdom", "phylum", "class", "order"), accepted=TRUE, ask=FALSE)
taxon_info_KJ <- rbind(taxon_info_KJ1, taxon_info_KJ2, taxon_info_KJ3) # I was getting a http timeout error. Hense, this unsatisfying workaround

not_in_ITIS_sp <- taxon_info_KJ[which(is.na(taxon_info_KJ$kingdom)),"query"]
taxon_info_KJ[which(is.na(taxon_info_KJ$kingdom)),] <- tax_name(not_in_ITIS_sp, c("kingdom", "phylum", "class", "order"), db="ncbi", accepted=TRUE, ask=FALSE)
taxon_info_KJ[which(is.na(taxon_info_KJ$kingdom)),] # still not found

# datasets with known taxanomy
taxon_info_RN <- mapply(get.taxanom, species = all_species_RN$species, dataset = all_species_RN$dataset)
taxon_info_RN <- t(taxon_info_RN)

# combining
taxon_info <- data.table(species = c(taxon_info_KJ$query, rownames(taxon_info_RN)),
                         kingdom = c(taxon_info_KJ$kingdom, taxon_info_RN[,1]),
                         phylum = c(taxon_info_KJ$phylum, taxon_info_RN[,2]),
                         class = c(taxon_info_KJ$class, taxon_info_RN[,3]),
                         order = c(taxon_info_KJ$order, taxon_info_RN[,4]),
                         stringsAsFactors = FALSE)

sort(taxon_info$species)

# cleaning up differences between ITIS and NCBI taxonomy names
taxon_info$kingdom <- gsub("Metazoa", "Animalia", taxon_info$kingdom)
taxon_info$kingdom <- gsub("Viridiplantae", "Plantae", taxon_info$kingdom)


write.csv(taxon_info, "clean_data/species_list_taxonomy.csv", row.names = FALSE)

traits_of_interest <- c("flower pollination syndrome", "whole plant height", "whole plant growth form",
                        "whole plant dispersal syndrome", "whole plant sexual system", "whole plant vegetative phenology",
                        "seed mass", "leaf relative growth rate", "fruit type", "leaf area",
                        "leaf dry mass", "leaf life span", "maximum whole plant longevity",
                        "plant flowering duration", "diameter at breast height (1.3 m)",
                        "leaf area per leaf dry mass")


# averages together all of the available trait values and returns vector
# mean of numeric traits, mode of discrete traits
# if no species-level traits available, goes to genus... if no genus, then NAs returned
BIEN_trait_wrapper <- function(species){
  #species <- "Taraxa"
  #species <- "Taraxacum sp."
  #species <- "Taraxacum officinale"
  #species <- "Tilia vulgaris"
  #species <- "Pyrus communis"
  #global.counter()
  # checking if function was given genus
  if(grepl("sp\\.", species) | !grepl("\\s", species)){
    # genus passed to function
    query <- gsub("\\s.*$", "", species)
    traits_pull <- data.table(BIEN_trait_traitbygenus(query, traits_of_interest, all.taxonomy = FALSE))
    tax_level <- "genus"
  } else{
    traits_pull <- data.table(BIEN_trait_traitbyspecies(species, traits_of_interest, all.taxonomy = FALSE))
    tax_level <- "species"
  }
  
  # going to genus and then to family if no species data found
  if(nrow(traits_pull) == 0){
    genus <- gsub("\\s.*$", "", species)
    traits_pull <- data.table(BIEN_trait_traitbygenus(genus, traits_of_interest, all.taxonomy = FALSE))
    tax_level <- "genus"
    if(nrow(traits_pull) == 0){
      family <- BIEN_taxonomy_genus(genus)$scrubbed_family[1]
      if(is.null(family)) return(rep("not found", length(traits_of_interest)+1)) # no information found for species, genus, or family
      traits_pull <- data.table(BIEN_trait_traitbyfamily(family, traits_of_interest, all.taxonomy = FALSE))
      tax_level <- "family"
    }
  }
  
  if(nrow(traits_pull) == 0) return(rep("not found", length(traits_of_interest)+1)) # no information found for species, genus, or family
  trait_output <- rep(NA, length(traits_of_interest))
  names(trait_output) <- traits_of_interest
  trait_avgs <- traits_pull[, average.traits(trait_value), by = trait_name]
  if(FALSE){
    test_trait <- traits_pull[which(traits_pull$trait_name == "diameter at breast height (1.3 m)") ,]
    trait_avgs <- test_trait[, average.traits(trait_value)]
  }
  trait_output[trait_avgs$trait_name] <- trait_avgs$V1
  output <- c(trait_output, tax_level = tax_level)
  return(output)
  
  #BIEN_trait_mean(species, traits_of_interest)
  #BIEN_trait_mean(species, "flower pollination syndrome") # this doesn't work
}

# smart averaging of numeric and categorical traits - used by BIEN_trait_wrapper
average.traits <- function(trait_vals){
  #trait_vals <- c("Herb", "Herb", "Herb")
  #trait_vals <- test_trait$trait_value
  #trait_vals <- traits_pull$trait_value
  #trait_vals <- traits_pull$trait_value[2090:2095]
  #trait_vals <- c("17.272", "16.4592", "Tree", "15.24")
  #trait_vals <- c("17", "16.4592", "3", "15.24")
  #trait_vals <- c("17.272", "16.4592", "15.24")
  trait_vals <- unlist(strsplit(trait_vals, ", ")) #dealing with some trait values having two entries per element, separated by a comma
  if(sum(is.na(trait_vals)) >= length(trait_vals)) return(NA)
  
  if(all.is.numeric(trait_vals)) trait_vals <- as.numeric(trait_vals) # converting to numeric if it looks numeric
  if(is.numeric(trait_vals)) return(as.character(mean(trait_vals, na.rm=T))) # I'm sorry - data.table wants all outputs in the same data type, so I have to convert back to string
  if(is.character(trait_vals) | is.factor(trait_vals)){
    trait_vals <- trait_vals[which(!is.na(trait_vals))]
    trait_vals <- tolower(trait_vals)
    return(names(sort(table(trait_vals),decreasing=TRUE)[1]))
  } 
  return(NA)
}

just_plants <- taxon_info[which(taxon_info$kingdom == "Plantae"),]
#test <- sapply(just_plants$species[1:5], BIEN_trait_wrapper)
#test <- sapply(c("Taraxacum officinale", "Taraxacum officina", "Taraxa"), BIEN_trait_wrapper)
global.counter(T)
plant_traits <- sapply(just_plants$species, BIEN_trait_wrapper)
plant_traits <- as.data.frame(t(plant_traits))

write.csv(plant_traits, "clean_data/plant_traits.csv", row.names = TRUE)

# summary(plant_traits)
# most commonly available traits for my species:
# whole plant growth form (most)
# whole plant height (half)
# seed mass (half)
# leaf area (half)
# leaf area per leaf dry mass (half)

# skipping the hassle of earlier code
#plant_species <- unique(trends[tax_group == "plant", species])

plant_traits <- read.csv("clean_data/plant_traits.csv", stringsAsFactors = FALSE)
plant_traits <- plant_traits[,c("species", "tax_level", "whole.plant.growth.form", "whole.plant.height", "seed.mass", "leaf.area", "leaf.area.per.leaf.dry.mass")]
plant_traits$whole.plant.growth.form <- mgsub(plant_traits$whole.plant.growth.form,
                                              c("forb", "geophyte", "fern", "hemicryptophyte", "herb\\*"),
                                              c("herb", "herb", "herb", "herb", "herb")) # I checked that all of the hemicryptophytes are herbs
plant_traits$whole.plant.growth.form <- mgsub(plant_traits$whole.plant.growth.form,
                                              c("vine", "climbing legume", "non-woody epiphyte", "climber", "liana", "hemiepiphyte", "parasite"),
                                              c("dependant", "dependant", "dependant", "dependant", "dependant", "dependant", "dependant")) # This is all species that need others to grow properly
plant_traits$whole.plant.growth.form <- mgsub(plant_traits$whole.plant.growth.form,
                                              c("shrub\\*", "graminoid", "aquatic\\*", "3", "4", "cactus", "not found"),
                                              c("shrub", "grass", "NA", "NA", "NA", "NA", "NA")) # misc. groups. Getting rid of the 3 aqautic & cactus species
colnames(plant_traits) <- c("species", "tax_level", "growth_form", "height", "seed_mass", "leaf_area", "leaf_area_per_mass")
table(plant_traits$tax_level) # almost all down to species
plant_traits <- transform(plant_traits,
                          tax_level = as.factor(tax_level),
                          growth_form = as.factor(growth_form),
                          height = as.numeric(height),
                          seed_mass = as.numeric(seed_mass),
                          leaf_area = as.numeric(leaf_area),
                          leaf_area_per_mass = as.numeric(leaf_area_per_mass))
plant_traits$species <- gsub("\\.", " ", plant_traits$species)
plant_traits <- as.data.table(plant_traits)

plant_species <- read.csv("clean_data/species_list_taxonomy.csv", stringsAsFactors = FALSE)
missing_species <- plant_species[which(plant_species$species %!in% plant_traits$species), "species"]
new_plant_traits <- pbsapply(missing_species, BIEN_trait_wrapper)
new_plant_traits <- t(new_plant_traits)
new_plant_traits <- as.data.table(new_plant_traits)
new_plant_traits$species <- missing_species

colnames(new_plant_traits) <- gsub(" ", "\\.", colnames(new_plant_traits))
new_plant_traits <- new_plant_traits[,c("species", "tax_level", "whole.plant.growth.form", "whole.plant.height", "seed.mass", "leaf.area", "leaf.area.per.leaf.dry.mass")]
colnames(new_plant_traits) <- c("species", "tax_level", "growth_form", "height", "seed_mass", "leaf_area", "leaf_area_per_mass")

plant_traits <- rbind(plant_traits, new_plant_traits)
write.csv(plant_traits, "clean_data/plant_traits.csv", row.names = TRUE)
