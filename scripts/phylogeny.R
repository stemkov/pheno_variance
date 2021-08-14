# matching trait data to all species

bird_traits <- fread("clean_data/bird_traits.csv", stringsAsFactors = FALSE)
plant_traits <- fread("clean_data/plant_traits.csv", stringsAsFactors = FALSE) # this is just plants with traits

trends_focus <- fread("Edited Data/trends_focus.csv", stringsAsFactors = F)

plant_traits <- trends_focus[tax_group == "plant",] # all plant trends
#plant_traits <- plant_traits[pheno_group == "leaf",] # I forgot to do this earlier... I was lumping leaves and flowers together
plant_traits <- plant_traits[pheno_group == "flower",] # I forgot to do this earlier... I was lumping leaves and flowers together
plant_traits <- plant_traits[, .(mean_slope = mean(mean_slope, na.rm=T),
                                 mean_slope_temp = mean(mean_slope_temp, na.rm=T),
                                 sigma_slope = mean(sigma_slope, na.rm=T),
                                 sigma_slope_temp = mean(sigma_slope_temp, na.rm=T),
                                 mean_se = mean(mean_se, na.rm=T),
                                 mean_se_temp = mean(mean_se_temp, na.rm=T),
                                 sigma_se = mean(sigma_se, na.rm=T),
                                 sigma_se_temp = mean(sigma_se_temp, na.rm=T)),by=species]
# length(unique(plant_traits$species))
# sort(plant_traits$species)

#colnames(plant_traits)[1] <- "species"
# cleaning up plant species
plant_traits$old_species <- plant_traits$species
plant_traits$species <- gsub("\\b(\\w+)\\s+\\1\\b", "\\1", plant_traits$species) # removed duplicate genus names
#gsub("\\b(\\w+)\\s+\\1\\b", "\\1", "Psephellus Psephellus carbonatus")
#gsub("\\b(\\w+)\\s+\\1\\b", "\\1", "Psephellus psephellus")
plant_traits$species <- gsub("_"," ", plant_traits$species)
plant_traits$species <- gsub("([[:lower:]]) ([[:lower:]])","\\1_\\2", plant_traits$species)
plant_traits$species <- gsub("\\.1","", plant_traits$species)
plant_traits <- plant_traits[species %!in% c("na_na", "NA.sp."),]
plant_traits <- plant_traits[!grepl("_sp\\.", plant_traits$species),] # removing genera-level names

#plant_tree <- read.tree("Raw Data/phylogenies/Vascular_Plants_rooted.dated.tre") # old tree - 31749 tips - 33% sp missing
# plant_tree <- read.tree("Raw Data/phylogenies/GBMB.tre") # best tree, might be hard to use - 79874 tips - 34% sp missing
plant_tree <- read.tree("Raw Data/phylogenies/ALLOTB.tre") # super tree method - 353185 tips - 20% sp missing
#plant_tree_forced <- force.ultrametric(plant_tree) # ALLOTB is ultrametric, but there's a "rounding" error or something, so I need to force it
bird_tree <- read.tree("Raw Data/phylogenies/bird_tree.tre") # 9993 tips, full tree

### PLANTS

# which plant trait names aren't in the phylo tree?
not_in_phylo <- sort(plant_traits[species %!in% plant_tree$tip.label, species]) 
print(paste0(round((length(not_in_phylo)/nrow(plant_traits)),2)*100, "% of plant species not in phylogeny"))

genera_in_tree <- unique(sub("_.*","",plant_tree$tip.label))
genera_in_data <- unique(sub("_.*","",plant_traits$species))
missing_genera <- genera_in_data[which(genera_in_data %!in% genera_in_tree)] #24 genera in data but not in phylogeny - there may be synonyms
missing_genera_sp <- plant_traits$species[which(sub("_.*","",plant_traits$species) %in% missing_genera)] #species from missing genera to look for synonyms
missing_tsns <- get_tsn(missing_genera_sp, ask=FALSE)
missing_tsns_found <- missing_tsns[which(!is.na(missing_tsns))]
possible_synonym_genera <- sub(" .*","",itis_acceptname(missing_tsns_found)$acceptedname) # entered 1 at every prompt
genera_synonyms <- data.frame(old_name = sub("_.*","",missing_genera_sp[which(!is.na(missing_tsns))]),
                              new_name  = possible_synonym_genera,
                              stringsAsFactors = FALSE)
genera_synonyms <- genera_synonyms[which(complete.cases(genera_synonyms)),] # removes NAs for genera where synonym wasn't found
beep(4)

# replacing the synonym genera names, hoping these are in the phylogeny - assuming that the epithets are the same
plant_traits$genus <- sub("_.*","",plant_traits$species)
plant_traits$epithet <- sub(".*_","",plant_traits$species)
plant_traits$genus_synon <- mapvalues(plant_traits$genus, from = genera_synonyms$old_name, to = genera_synonyms$new_name) # don't worry, it returns warning because there's duplicates in the replacement
plant_traits$species_synon <- paste(plant_traits$genus_synon, plant_traits$epithet, sep="_")

# which plant names aren't in the phylo tree after looking for genera synonyms?
not_in_phylo <- sort(plant_traits[species_synon %!in% plant_tree$tip.label, species_synon]) 
#not_in_phylo <- sort(plant_traits$species_synon[which(plant_traits$species_synon %!in% plant_tree$tip.label)]) 
print(paste0(round((length(not_in_phylo)/nrow(plant_traits)),2)*100, "% of plant species not in phylogeny after looking for genera synon"))
genera_in_data <- unique(plant_traits$genus_synon)
missing_genera <- genera_in_data[which(genera_in_data %!in% genera_in_tree)]

plant_traits$species_synon

# I need to use force.ultrametric on the plant tree but it's too big, so I need to subset to just species in data and genera in which species are missing for congeneric.impute
plant_tree$genus.label <- sub("_.*","",plant_tree$tip.label)
not_in_phylo_sp <- sort(plant_traits$species_synon[which(plant_traits$species_synon %!in% plant_tree$tip.label)]) 
not_in_phylo_genus <- unique(sub("_.*","",not_in_phylo_sp))
genera_to_keep <- not_in_phylo_genus[which(not_in_phylo_genus %!in% missing_genera)] # I want to keep these genera in the pruned phylogeny for congeneric.impute - these are genera that have species that are present in data but missing from phylogeny 
plant_tree_pruned <- drop.tip(plant_tree, plant_tree$tip.label[which(plant_tree$genus.label %!in% genera_to_keep & plant_tree$tip.label %!in% plant_traits$species_synon)])
# plant_tree_pruned_forced <- force.ultrametric(plant_tree_pruned) # doesn't error immediately, but overloads memory and ends up erroring

# ^ this still results in too large of a tree - I'm going to try getting just one species from each genus for congeneric.impute
plant_tree$first_in_genus <- rep(FALSE, length(plant_tree$tip.label)) #write something that gives true or false for each tip
plant_tree$first_in_genus[match(genera_to_keep, plant_tree$genus.label)] <- TRUE
tips_to_keep <- which(plant_tree$first_in_genus | plant_tree$tip.label %in% plant_traits$species_synon)
plant_tree_pruned <- drop.tip(plant_tree, plant_tree$tip.label[-tips_to_keep])
plant_tree_pruned_forced <- force.ultrametric(plant_tree_pruned)

plant_traits_phylo_sub <- plant_traits[which(plant_traits$genus_synon %in% genera_in_tree),] # only loses 12 species of 1095
plant_tree_imputed <- congeneric.impute(plant_tree_pruned_forced, plant_traits_phylo_sub$species_synon)
plant_tree_imputed

not_in_phylo <- sort(plant_traits[species_synon %!in% plant_tree_imputed$tip.label, species_synon]) 
print(paste0(round((length(not_in_phylo)/nrow(plant_traits)),2)*100, "% of plant species not in phylogeny after imputing"))

saveRDS(plant_tree_imputed, "Edited Data/plant_tree_imputed.RDS")
write.csv(plant_traits, "Edited Data/plant_traits_from_phylo.csv", row.names = FALSE)

# svg("Figures/plant_phylogeny.svg", height=40, width=15)
# plot(plant_tree_imputed, "phylogram", show.tip.label=TRUE, cex = 0.2) # very minor changes to the overall structure
# dev.off()

# bootstrapping the above
# - picking random species to keep (instead of first)
# - just running congeneric.impute repeatedly

sample.match <- function(x,table) as.numeric(sapply(x, function(x) sample(which(table == x),1))) # random match value instead of first
# match(c("a","b","d"), rep(letters[1:4], each=5)) # first matching value
# sample.match(c("a","b","d"), rep(letters[1:4], each=5)) # random matching value

n_boot <- 100
for(i in 1:n_boot){
  plant_tree$first_in_genus <- rep(FALSE, length(plant_tree$tip.label))
  plant_tree$first_in_genus[sample.match(genera_to_keep, plant_tree$genus.label)] <- TRUE
  #match(genera_to_keep, plant_tree$genus.label)
  #sample.match(genera_to_keep, plant_tree$genus.label)
  ##### YOU're HERE. Replace the above code with your new shiny sample.match function, and proceed!
  tips_to_keep <- which(plant_tree$first_in_genus | plant_tree$tip.label %in% plant_traits$species_synon)
  plant_tree_pruned_boot <- drop.tip(plant_tree, plant_tree$tip.label[-tips_to_keep])
  plant_tree_pruned_boot_forced <- force.ultrametric(plant_tree_pruned)
  
  plant_traits_phylo_sub_boot <- plant_traits[which(plant_traits$genus_synon %in% genera_in_tree),] # only loses 12 species of 1095
  plant_tree_imputed_boot <- congeneric.impute(plant_tree_pruned_boot_forced, plant_traits_phylo_sub_boot$species_synon)
  plant_tree_imputed_boot
  
  #not_in_phylo <- sort(plant_traits_phylo_sub_boot[species_synon %!in% plant_tree_imputed_boot$tip.label, species_synon]) 
  #print(paste0(round((length(not_in_phylo)/nrow(plant_traits)),2)*100, "% of plant species not in phylogeny after imputing"))
  print(paste(i, "of", n_boot, "bootstraps"))
  saveRDS(plant_tree_imputed_boot, paste0("Edited Data/boot_trees/plant_tree_boot_",i,".RDS"))
}


### BIRDS

bird_tree <- read.tree("Raw Data/phylogenies/bird_tree.tre") # 9993 tips, full tree
bird_traits <- read.csv("Edited Data/bird_traits.csv", stringsAsFactors = FALSE)

bird_traits$species <- gsub(" ", "_", bird_traits$species)
bird_traits$species_synon <- gsub(" ", "_", bird_traits$species_synon)
bird_traits$species_synon_clean <- gsub(" ", "_", bird_traits$species_synon_clean)

# which bird trait names aren't in the phylo tree?
bird_traits$species[which(bird_traits$species_synon_clean %!in% bird_tree$tip.label)] # cool! I have them all
bird_traits <- bird_traits[which(complete.cases(bird_traits)),]
bird_tree_pruned <- keep.tip(bird_tree, bird_traits$species_synon_clean)
saveRDS(bird_tree_pruned, "Edited Data/bird_tree_pruned") # wow so much easier than the plant phylogeny

plot(bird_tree_pruned, cex=0.5)

