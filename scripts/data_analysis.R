### analyzing all of the phenology data for changes in mean & variance in relation to time and temperature

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

# tests for shift and change in variance
mean.sigma.test <- function(doy, explan, quant_reg=FALSE, plot=FALSE, counter=FALSE, ...){
  
  if(FALSE){
    test <- all_data[which(all_data$species == "Spinus tristis"),]
    doy <- test$doy
    explan <- test$year
    
    test <- all_data[species == "Pieris rapae" & site == 100 & phenophase == "first_occurence",]
    doy <- test$doy
    explan <- test$year
    
    test <- all_data[species == "Hirundo rustica" & site == 100 & phenophase == "first_occurence",]
    doy <- test$doy
    explan <- test$year
  }
  
  if(counter) global.counter()
  if(length(doy) != length(explan)) return(list(mean_slope = -8888, mean_se = -8888, mean_p = -8888, sigma_slope = -8888, sigma_se = -9999, sigma_p = -8888, bp_p = -8888, n=length(doy))) # doy and explan differ in lengths - this should never happen
  if(sum(!is.na(doy)) < 10 | sum(!is.na(explan)) < 10) return(list(mean_slope = -9999, mean_se = -9999, mean_p = -9999, sigma_slope = -9999, sigma_se = -9999, sigma_p = -9999, bp_p = -9999, n=length(doy))) # output for data.table all needs to be of the same type, so I can't output NA
  # if(sum(!is.na(doy)) != sum(!is.na(explan))) return(list(mean_slope = -7777, mean_p = -7777, sigma_slope = -7777, sigma_p = -7777, bp_p = -7777, n=length(doy))) # there in an NA in the doy or explan timeseries, but not both ... this was causing problems with manomet because some doys are NAs because phest couldn't estimate onset properly. I'm changing this to remove doys with NAs
  if(sum(!is.na(doy)) != sum(!is.na(explan))){
    warning("Unequal number of NAs in doy vs. explan") # not a problem if there are supposed to be some NAs in doy
    doy_explan <- data.frame(doy, explan)
    doy_explan <- doy_explan[which(complete.cases(doy_explan)),] # getting rid of years with either missing doy or explan data
    doy <- doy_explan$doy
    explan <- doy_explan$explan
  }

  mean_model <- lm(doy ~ explan)
  mean_summ <- summary(mean_model)

  if(quant_reg){
    res <- abs(residuals(mean_model))
    res_model <- rq(res ~ explan, tau=0.6827)
    res_summ <- summary(res_model, se="boot")
  } else{
    res <- abs(residuals(mean_model)) * sqrt(2)
    res_model <- lm(res ~ explan)
    res_summ <- summary(res_model)
  }
  
  bp_test <- bptest(mean_model)
  
  if(plot){
    par(mfrow=c(2,1), mar=c(2.5, 4.1, 1, 2.1))
    visreg(mean_model, ylab="Phenophase DOY", ...)
    par(mar=c(4.1, 4.1, 0, 2.1))
    visreg(res_model, ylab="Absolute residuals * sqrt(2)", ...)
    par(mfrow=c(1,1), mar=c(4.1, 4.1, 2.1, 2.1))
  }
  
  results <- list(mean_slope = as.numeric(mean_summ$coefficients[2]),
                  mean_p = as.numeric(mean_summ$coefficients[8]),
                  mean_se = as.numeric(mean_summ$coefficients[4]),
                  sigma_slope = as.numeric(res_summ$coefficients[2]),
                  sigma_se = as.numeric(res_summ$coefficients[4]),
                  sigma_p = as.numeric(res_summ$coefficients[8]),
                  bp_p = as.numeric(bp_test$p.value),
                  n = length(doy))
  
  return(results)
}

# scrubbing data
cleaning <- function(doy, explan, outlier_sd = NA, BIC_thresh = NA, outlier_range = 150, discont_thresh=0.25, min_final_year=NA, show=F){
  if(length(doy) != length(explan)) stop("doy and explan of unequal length") # this should never happen
  if(sum(!is.na(doy)) < 10 | sum(!is.na(explan)) < 10) return("too short")
  if(sum(!is.na(doy)) != sum(!is.na(explan))){
    warning("Unequal number of NAs in doy vs. explan") # not a problem if there are supposed to be some NAs in doy
    doy_explan <- data.frame(doy, explan)
    doy_explan <- doy_explan[which(complete.cases(doy_explan)),] # getting rid of years with either missing doy or explan data
    doy <- doy_explan$doy
    explan <- doy_explan$explan
  }
  
  if(!is.na(min_final_year)){
    if(max(explan, na.rm = T) < min_final_year) return("ends too soon")
  }
  
  # trying to root out erronious outliers
  if(nchar(round(explan[1])) == 4){
    if(length(unique(na.omit(explan))) < 10) return("not enough unique years")
    if(length(explan) < discont_thresh*length(min(explan):max(explan))) return("large discontinuity in years")
  }
  
  if(any(doy > 300) & any(doy < 50)) return("likely year split outlier")
  if(any(doy < 0) | any(doy > 365)) return("likely phest outlier")
  if(!is.na(outlier_sd)) if(sd(doy) > outlier_sd) return("likely data outlier based on large SD")
  if(!is.na(BIC_thresh)){
    mclust_model <- mclustBIC(doy, 1:3, verbose=FALSE)
    if(all(is.na(mclust_model[,]))) return("mclust failed")
    mclust_max <- max(mclust_model[1,], na.rm = TRUE)
    if(show) print(mclust_model)
    if(!all(is.na(mclust_model[2:3,]))){ # I need this because if v2 and v3 are all NA, then the data's probably ok
      if(any(mclust_model[2:3,] > (mclust_max + BIC_thresh), na.rm=T)) return("likely data outlier based on clustering")
      #if(any(mclust_model[2:3,] < mclust_max-BIC_thresh, na.rm=T)) return("likely data outlier based on clustering")
    }
  }
  if(!is.na(outlier_range)) if(diff(range(doy,na.rm=TRUE)) > outlier_range) return("likely data outlier based on large range")
  return("NA")
}

# auxillary plotting function - reporting number of outliers outside plotting range
report.outliers <- function(data, axis = "y", unified=FALSE, zero=FALSE, just_report=FALSE){
  if(any(!is.numeric(na.omit(data)))) stop("data must be numeric")
  if(axis %!in% c("y","x")) stop("axis must be 'y' or 'x'")
  if(axis == "y"){
    range <- par("usr")[3:4] 
    over_side <- 3; under_side <- 1
    if(zero) abline(h=0, lty=2)
  }
  if(axis == "x"){
    range <- par("usr")[1:2]
    over_side <- 4; under_side <- 2
    if(zero) abline(v=0, lty=2)
  }
  over <- sum(data > range[2])
  under <- sum(data < range[1])
  if(just_report){
    print(paste(over, "points over plotting range"))
    print(paste(under, "points under plotting range"))
  } else{
    if(unified){
      if(over+under > 0) mtext(paste(over+under, "points outside plotting range  "), 3, line=-1.2, adj=1, cex=0.8)
      
    } else{
      if(over > 0) mtext(paste(over, "points over plotting range  "), over_side, line=-1.2, adj=1, cex=0.8)
      if(under > 0) mtext(paste(under, "points under plotting range  "), under_side, line=-1.2, adj=1, cex=0.8)
    }
  }
}

all_data <- fread("clean_data/all_data_standardized.csv", stringsAsFactors = FALSE)

#all_data$problem <- NA
# preliminary flagging of problematic timeseries
# BIC_thresh=100... hardly any positives
# BIC_thresh=10... many false positives
# BIC_thresh=50... no false positives - probably missed some
# BIC_thresh=25... seems to be the sweet spot
n_combs <- uniqueN(all_data[,c("site", "species", "phenophase")])
all_data[, problem := cleaning(doy, year, BIC_thresh=25, outlier_range=90, min_final_year = 2000),
         .(site, species, phenophase)]


# for checking effectiveness of cluster cleaning
# cluster_problems <- all_data[problem == "likely data outlier based on clustering",]
# cluster_problems[, mean.sigma.test(doy, year, plot=T), by=.(species,site,phenophase)]
# discont_problems <- all_data[problem == "large discontinuity in years",]
# discont_problems[, mean.sigma.test(doy, year, plot=T), by=.(species,site,phenophase)]
# range_problems <- all_data[problem == "likely data outlier based on large range",]
# range_problems[, mean.sigma.test(doy, year, plot=T), by=.(species,site,phenophase)]
problem_table <- table(all_data$problem)
paste0(round(problem_table["likely data outlier based on clustering"]/nrow(all_data),4)*100, "% of data likely outliers based on clustering")
paste0(round(problem_table["likely phest outlier"]/nrow(all_data),4)*100, "% of data likely phest outlier")
paste0(round(problem_table["likely year split outlier"]/nrow(all_data),4)*100, "% of data likely year split outlier")
paste0(round(problem_table["large discontinuity in years"]/nrow(all_data),4)*100, "% of data large discontinuity in years")
paste0(round(problem_table["likely data outlier based on large range"]/nrow(all_data),4)*100, "% of data likely aseasonal based on large range")
paste0(round(problem_table["ends too soon"]/nrow(all_data),4)*100, "% of data ends too soon")

# why is so much NECTAR data being removed?
table(all_data[dataset=="nectar",problem]) # because it ends too soon or is to short

# getting the number of timeseries with problems, rather than number of data rows as above
problem_timeseries <- all_data[, unique(problem), .(site, species, phenophase)]
problem_timeseries_table <- table(problem_timeseries$V1)
paste0("Of the ", nrow(problem_timeseries), " total timeseries that were evaluated, the following were excluded:")
paste0(round(problem_timeseries_table["likely data outlier based on clustering"]/nrow(problem_timeseries),4)*100, "% of data likely outliers based on clustering (", problem_timeseries_table["likely data outlier based on clustering"],")")
paste0(round(problem_timeseries_table["likely phest outlier"]/nrow(problem_timeseries),4)*100, "% of data likely phest outlier (", problem_timeseries_table["likely phest outlier"],")")
paste0(round(problem_timeseries_table["likely year split outlier"]/nrow(problem_timeseries),4)*100, "% of data likely year split outlier (", problem_timeseries_table["likely year split outlier"],")")
paste0(round(problem_timeseries_table["large discontinuity in years"]/nrow(problem_timeseries),4)*100, "% of data large discontinuity in years (", problem_timeseries_table["large discontinuity in years"],")")
paste0(round(problem_timeseries_table["likely data outlier based on large range"]/nrow(problem_timeseries),4)*100, "% of data likely aseasonal based on large range (", problem_timeseries_table["likely data outlier based on large range"],")")
paste0(round(problem_timeseries_table["ends too soon"]/nrow(problem_timeseries),4)*100, "% of data ends too soon (", problem_timeseries_table["ends too soon"],")")



all_data <- all_data[problem == "NA",]

summary(all_data)
unique(all_data$dataset)
plot(table(all_data$dataset))

# limitting analysis to just 1958-2018, because that's the range of TerraClimate - removes 8% of data
all_data <- na.omit(all_data, "tmax")
hist(all_data$year,100)

# summary of timeseries data - before further pruning (trends_focus)
datasets <- table(all_data[,c("species", "phenophase", "site")])
nonzero_datasets <- datasets[which(datasets != 0)]
paste("There are", length(nonzero_datasets), "timeseries datasets")
paste("The median length of timeseries is", median(nonzero_datasets), "years")
paste("There are", length(unique(all_data$species)), "species across", length(unique(all_data$site)), "sites")



##############################
# Residual Method - step one #
##############################

### trend over YEARS (shift)

n_combs <- uniqueN(all_data[,c("site", "species", "phenophase")])
#trends <- all_data[, {cat(.GRP/n_combs*100,"% year trends\n"); mean.sigma.test(doy, year, counter=F)}, .(site, species, phenophase)]
#trends <- all_data[, {cat(.GRP/n_combs*100,"% year trends\n"); mean.sigma.test(log(doy+1), year, quant_reg=T)}, .(site, species, phenophase)] # checking for Lizzy Wolkovich
trends <- all_data[, {cat(.GRP/n_combs*100,"% year trends\n"); mean.sigma.test(doy, year, quant_reg=T)}, .(site, species, phenophase)]
all_data_unique <- all_data[, .(mean_doy = mean(doy),
                                min_year = min(year),
                                max_year = max(year),
                                mean_tmax = mean(tmax),
                                lat = mean(lat),
                                lon = mean(lon),
                                dataset = unique(dataset)), by = .(site, species, phenophase)]
all.equal(trends[,c("site", "species", "phenophase")], all_data_unique[,c("site", "species", "phenophase")]) # it's safe to merge them
trends <- cbind(trends, all_data_unique[,-c(1,2,3)])

length(which(trends$mean_slope == -9999)) # There are a few new "too short" datasets because the tmax trimming cut off some historic data
length(which(trends$mean_slope == -7777)) # There are a few new "too short" datasets because the tmax trimming cut off some historic data
trends <- trends[mean_p %!in% c(-9999,-7777),]
# trend over TEMP (sensitivity)
#temp_trends <- all_data[, {cat(.GRP/n_combs*100,"% temp trends\n"); mean.sigma.test(log(doy+1), tmax, quant_reg=T)}, .(site, species, phenophase)] # checking for Lizzy
temp_trends <- all_data[, {cat(.GRP/n_combs*100,"% temp trends\n"); mean.sigma.test(doy, tmax, quant_reg=T)}, .(site, species, phenophase)]
length(which(temp_trends$mean_slope == -9999)) # There are a few new "too short" datasets because the tmax trimming cut off some historic data
length(which(temp_trends$mean_slope == -7777))
temp_trends <- temp_trends[mean_slope %!in% c(-9999,-7777),]
# merging year and temp trends
all.equal(trends[,c("site", "species", "phenophase")], temp_trends[,c("site", "species", "phenophase")]) # it's safe to merge them
trends[,c("mean_slope_temp", "mean_se_temp", "mean_p_temp", "sigma_slope_temp", "sigma_se_temp", "sigma_p_temp", "bp_p_temp", "n_temp")] <- temp_trends[,c("mean_slope", "mean_se", "mean_p", "sigma_slope", "sigma_se", "sigma_p", "bp_p", "n")]


write.csv(trends, "clean_data/trends.csv", row.names = FALSE)
# trends <- fread("clean_data/trends.csv", stringsAsFactors = FALSE) ### to load in previously calculated trends

#temp over time - how is climate at the sites changing?
n_combs <- uniqueN(all_data[,c("site", "species", "phenophase")])
clim_trends <- all_data[, {cat(.GRP/n_combs*100,"% climate trends\n"); mean.sigma.test(tmax, year, counter=F)}, .(site, species, phenophase)]
all.equal(clim_trends[,c("site", "species", "phenophase")], all_data_unique[,c("site", "species", "phenophase")]) # it's safe to merge them
clim_trends <- cbind(clim_trends, all_data_unique[,-c(1,2,3)])
clim_trends <- clim_trends[mean_p %!in% c(-9999,-7777),]


#############################
# Explain Trends - step two #
#############################

scale.data <- function(x, log=FALSE){
  if(is.numeric(x)){
    if(log) x <- log(x)
    return(as.numeric(arm::rescale(x))) # rescale() divides by 2 SDs to make catagorical vars comparable to numeric vars
  }  
  else return(x)
}


### seasonality and mean temperature for every coords
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

# I already got seasonality and mean_temp for pep and ussr
trends[,`:=`(seasonality = as.numeric(NA), mean_temp = as.numeric(NA))]
pep_sites <- fread("clean_data/pep_sites.csv", stringsAsFactors = F); pep_sites[, site := paste("p", site, sep="_")]
ussr_sites <- fread("clean_data/ussr_sites.csv", stringsAsFactors = F); ussr_sites[, site := paste("u", site, sep="_")]
pep_ussr_sites <- rbindlist(list(pep_sites, ussr_sites))
trends[is.na(seasonality) & dataset %in% c("pep", "ussr"), "seasonality"] <- as.numeric(sapply(trends[is.na(seasonality) & dataset %in% c("pep", "ussr"),site], function(x) pep_ussr_sites[site == x, seasonality])) # this is inefficient, but older code was causing problems
trends[is.na(mean_temp) & dataset %in% c("pep", "ussr"), "mean_temp"] <- as.numeric(sapply(trends[is.na(mean_temp) & dataset %in% c("pep", "ussr"),site], function(x) pep_ussr_sites[site == x, mean_temp]))

# previously calculated seasonality and mean_temp values
if(file.exists("clean_data/sites.csv")){
  sites <- fread("clean_data/sites.csv", stringsAsFactors = F)
  trends[is.na(seasonality), "seasonality"] <- as.numeric(sapply(trends[is.na(seasonality),site], function(x) sites[site == x, seasonality])) # this is inefficient, but older commented code below was causing problems
  trends[is.na(mean_temp), "mean_temp"] <- as.numeric(sapply(trends[is.na(mean_temp),site], function(x) sites[site == x, mean_temp]))
} 

n_combs <- uniqueN(trends[is.na(seasonality), c("site")])
trends[is.na(seasonality), seasonality := {cat(.GRP/n_combs*100,"% seasonality\n"); seasonality(mean(lat), mean(lon), counter=F)}, by=site]
trends[is.na(mean_temp), mean_temp := {cat(.GRP/n_combs*100,"% mean_temp\n"); seasonality(mean(lat), mean(lon), tmax_mean=T)}, by=site]


if(!file.exists("clean_data_sites.csv")){ # saving a bit of time for reruns
  sites <- trends[, .(records = .N,
                      dataset = unique(dataset),
                      lat = unique(lat),
                      lon = unique(lon),
                      seasonality = unique(seasonality),
                      mean_temp = unique(mean_temp)), by=site]
  write.csv(sites, "clean_data/sites.csv", row.names = FALSE)
}


### world map
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- ne_coastline(scale = "medium", returnclass = "sf")
seasonality_palette <- colorRampPalette(c("orange","#85039c"), space="Lab")
sites$col <- seasonality_palette(50)[as.numeric(cut(sites$seasonality,breaks = 50))]
sites_sf <- st_as_sf(sites[order(dataset),], coords = c("lon", "lat"), crs = st_crs(world))

### getting average phenological position - this how many days away from the average doy (across all species) at that site
pheno.position <- function(species_raw, phenophase_raw, site_raw){
  all_sp <- all_data[phenophase == phenophase_raw & site == site_raw,]
  focal_sp <- all_data[species == species_raw & phenophase == phenophase_raw & site == site_raw,]
  
  all_mean <- mean(all_sp$doy, na.rm=T)
  focal_mean <- mean(focal_sp$doy, na.rm=T)
  position <- diff(c(all_mean, focal_mean))
  return(position)
}
n_combs <- uniqueN(trends[,c("species", "phenophase", "site")])
trends[, pheno_position := {cat(.GRP/n_combs*100,"% pheno_position\n"); pheno.position(species, phenophase, site)},
       by = .(species, phenophase, site)]


bird_traits <- read.csv("clean_data/bird_traits.csv", stringsAsFactors = FALSE)
# merging diet types into 3 basics: herbi, omni, and carni
bird_traits$diet <- mgsub(bird_traits$diet, c("PlantSeed", "FruiNect", "Invertebrate", "VertFishScav"), c("Herbivore", "Herbivore", "Carnivore", "Carnivore"))
bird_traits$diet <- as.factor(bird_traits$diet)

plant_traits <- plant_traits[tax_level != "not found",]

### broad taxonomic groups: plant, bird, insect, other
leftover_tax_groups <- read.csv("raw_data/species_tax_groups.csv", stringsAsFactors = F)
tax.group <- function(species, dataset){
  if(dataset == "manomet") return("bird")
  if(dataset == "rothamsted") return("insect")
  if(dataset %in% c("nectar","rmbl", "pep")) return("plant")
  if(dataset == "ussr"){
    sp <- species
    if(!exists("taxonomy")) taxonomy <- fread("raw_data/ussr/taxonomy.csv", stringsAsFactors = FALSE)
    tax_info <- taxonomy[taxon == sp,]
    if(tax_info$kindgom == "Plantae") return("plant")
    if(tax_info$class == "Aves") return("bird")
    if(tax_info$class == "Insecta") return("insect")
    return("NA")
  }
  
  if(species %in% leftover_tax_groups$species) return(leftover_tax_groups[which(leftover_tax_groups$species == species), "tax_group"])
  species <- gsub("_", " ", species)
  if(species %in% leftover_tax_groups$species) return(leftover_tax_groups[which(leftover_tax_groups$species == species), "tax_group"])
  if(species %in% plant_traits$species) return("plant")
  if(species %in% c(bird_traits$species,bird_traits$species_synon,bird_traits$species_synon_clean)) return("bird")
  
  tax_lookup <- tax_name(species, c("kingdom", "phylum", "class", "order"), accepted=TRUE, ask=FALSE)
  if(is.na(tax_lookup$kingdom)) return("NA")
  if(tax_lookup$class == "Insecta") return("insect")
  if(tax_lookup$class == "Aves") return("bird")
  if(tax_lookup$kingdom == "Plantae") return("plant")
  return("NA")
}

trends[, tax_group := tax.group(species, dataset), by=.(species,dataset)]
plot(table(trends$tax_group))

species <- trends[,.(timeseries = .N,
                     observations = sum(n),
                     datasets = paste(unique(dataset), sep=", "),
                     tax_group = unique(tax_group),
                     year_shift = mean(mean_slope),
                     temp_shift = mean(mean_slope_temp),
                     pheno_position = mean(pheno_position)), by=species]
write.csv(species, "clean_data/species.csv", row.names = FALSE)



# merging broad taxonomic group and phenophase - because I want plant flowers and plant leaves
trends$pheno_group <- paste(trends$tax_group, trends$phenophase)
trends$pheno_group <- mgsub(trends$pheno_group,
                            c("insect first_occurence", "insect 1st_occurrence", "insect 1st occurrence", "bird first_occurence", "bird 1st occurrence", "plant first_leaf", "plant first_flower"),
                            c("insect", "insect", "insect", "bird", "bird", "leaf", "flower")) # making the pheno groups more readable
trends <- trends[pheno_group != "plant first_occurence",] # first occurence of Acer palmatum is first red - didn't catch it in data cleaning
sort(table(trends$pheno_group))

### subsetting to just groups I want
trends_focus <- trends[tax_group %in% c("plant", "insect", "bird"),]
trends_focus$pheno_group <- factor(trends_focus$pheno_group, levels = c("leaf", "flower", "insect", "bird"))
write.csv(trends_focus, "clean_data/trends_focus_quantreg.csv", row.names = FALSE)

# trends_focus <- fread("clean_data/trends_focus_quantreg.csv") # load this to just to do second-level analysis

### overall distribution statistics - data not scaled, so effect sizes are in original units
mean_shift_model <- lm(mean_slope ~ 1 ,data = trends_focus)
summary(mean_shift_model)
mean_sens_model <- lm(mean_slope_temp ~ 1 ,data = trends_focus)
summary(mean_sens_model)
sigma_shift_model <- lm(sigma_slope ~ 1 ,data = trends_focus)
summary(sigma_shift_model)
sigma_sens_model <- lm(sigma_slope_temp ~ 1 ,data = trends_focus)
summary(sigma_sens_model)

### summary by pheno_group
trends_summary <- trends_focus[, .(mean_min = mean(min_year),
                                   mean_max = mean(max_year),
                                   sd_min = sd(min_year),
                                   sd_max = sd(max_year),
                                   min25 = quantile(min_year, 0.25),
                                   min75 = quantile(min_year, 0.75),
                                   max25 = quantile(max_year, 0.25),
                                   max75 = quantile(max_year, 0.75),
                                   n = .N,
                                   mean_lat = mean(lat),
                                   n_species = length(unique(species))),
                               by=pheno_group]
pheno_group_pallete <- brewer.pal(length(levels(as.factor(trends_focus$pheno_group))), "Dark2")
trends_summary[, col := pheno_group_pallete[c(3,4,2,1)]]
trends_summary[, group_name := c("Insects", "Flowers", "Leaves", "Birds")]

# statistics for the contours plot relationship
model <- lm(mean_slope ~ mean_slope_temp, data=trends_focus)
summary(model)

### summary statistics
trends_focus
paste("There are", nrow(trends_focus), "timeseries datasets")
paste("The median length of timeseries is", median(trends_focus$n), "years")
paste("Total number of observations is", sum(trends_focus$n))
paste("There are", length(unique(trends_focus$species)), "species across", length(unique(trends_focus$site)), "sites")
paste(nrow(trends_focus[pheno_group == "leaf",]), "leaves",
      nrow(trends_focus[pheno_group == "flower",]), "flowers",
      nrow(trends_focus[pheno_group == "bird",]), "birds",
      nrow(trends_focus[pheno_group == "insect",]), "insects time-series")
paste(length(unique(trends_focus[tax_group == "plant", species])), "plant",
      length(unique(trends_focus[pheno_group == "bird", species])), "bird",
      length(unique(trends_focus[pheno_group == "insect", species])), "insect species")
range(trends$mean_temp)
range(trends$seasonality)

paste0(round(nrow(trends_focus[pheno_group == "leaf" & mean_slope_temp < 0,])/ nrow(trends_focus[pheno_group == "leaf",]),4), "% leaves ",
      round(nrow(trends_focus[pheno_group == "flower" & mean_slope_temp < 0,])/ nrow(trends_focus[pheno_group == "flower",]),4), "% flowers ",
      round(nrow(trends_focus[pheno_group == "bird" & mean_slope_temp < 0,])/ nrow(trends_focus[pheno_group == "bird",]),4), "% birds ",
      round(nrow(trends_focus[pheno_group == "insect" & mean_slope_temp < 0,])/ nrow(trends_focus[pheno_group == "insect",]),4), "% insects with negative mean sensitivity")

paste(round(nrow(trends_focus[sigma_slope < 0,])/ nrow(trends_focus),4), " overall with negative sigma slope")


### example species

# variance sensitivity
clear_trends_variance_temp <- trends_focus[, .(sigma_slope = median(sigma_slope_temp), N = .N, perc_negative = sum(sigma_slope_temp < 0)/.N), by = .(species, phenophase)]
clear_trends_variance_temp[N > 10 & sigma_slope %between% c(-2,-0.5),][order(-N)]

# hist(trends_focus[species == "Acer platanoides", sigma_slope_temp])
# mapply(check.raw.data,
#        trends_focus[species == "Acer platanoides",site],
#        MoreArgs = list(species_raw = "Acer platanoides", phenophase_raw = "first_flower", plot=T, explan="temp"))
# 
# hist(trends_focus[species == "Hydrangea_macrophylla", sigma_slope_temp])
# mapply(check.raw.data,
#        trends_focus[species == "Hydrangea_macrophylla",site],
#        MoreArgs = list(species_raw = "Hydrangea_macrophylla", phenophase_raw = "first_flower", plot=T, explan="temp"))

# variance shift
clear_trends_variance <- trends_focus[, .(sigma_slope = median(sigma_slope), N = .N, perc_negative = sum(sigma_slope < 0)/.N), by = .(species, phenophase)]
clear_trends_variance[N > 10 & sigma_slope %between% c(-1,-0.1),][order(-N)]

# hist(trends_focus[species == "Vaccinium myrtillus", sigma_slope], 10)
# mapply(check.raw.data,
#        trends_focus[species == "Vaccinium myrtillus",site],
#        MoreArgs = list(species_raw = "Vaccinium myrtillus", phenophase_raw = "first_flower", plot=T, explan="year"))


### datasets summary table
datasets_summary <- trends_focus[, .(`n time-series` = .N,
                                     sites = length(unique(site)),
                                     species = length(unique(species)),
                                     `n obs.` = .N*mean(n),
                                     `mean time-series length` = round(mean(n)),
                                     `mean lat.` = round(mean(lat),1),
                                     `mean lon.` = round(mean(lon),1),
                                     `mean annual max. temp.` = round(mean(mean_temp), 1),
                                     `mean seaslty` = round(mean(seasonality), 1),
                                     `data range` = paste(min(min_year),max(max_year),sep="-")),
                                 by=dataset]

datasets_summary <- datasets_summary[order(`n time-series`, decreasing=T),]
datasets_summary$dataset <- c("PEP", "CNC", "Japan", "Korea", "RMBL", "Rothamsted", "NECTAR", "Manomet")
datasets_summary$`pheno. groups` <- c("leaves, flowers", "all", "all", "all", "flowers", "insects", "leaves, flowers", "birds")

datasets_table <- flextable(datasets_summary, cwidth = 0.5)
datasets_table <- fontsize(datasets_table, size=9)
datasets_table <- fontsize(datasets_table, size=9, part="header")
datasets_table <- theme_vanilla(datasets_table)
datasets_table <- width(datasets_table, j=2, 0.3)
save_as_html(datasets_table, path="manuscript/datasets_table.html")

### how many countries does the data come from?
countries <- ne_countries(scale = "medium", returnclass = "sf")
sites_sf <- st_as_sf(sites[order(dataset),], coords = c("lon", "lat"), crs = st_crs(world))
site_countries <- st_intersects(countries, sites_sf)
represented_countries <- countries$name[which(apply(site_countries, 1, any))]
paste("These countries are represented:", paste(represented_countries, collapse = ", "))
paste(length(represented_countries), "countries total")


### Main models
hist(1/trends$mean_se, 100)
hist(1/trends$mean_se_temp, 100)
hist(1/trends$sigma_se, 100)
hist(1/trends$sigma_se_temp, 100)

trends_scaled <- trends_focus
cols_to_scale <- c("seasonality", "mean_temp", "pheno_position")
trends_scaled[, (cols_to_scale) := lapply(.SD, scale.data), .SDcols=cols_to_scale]

year_mean_model <- lmer(mean_slope ~ 0 + pheno_group + seasonality + mean_temp + pheno_position + (1|site:dataset) + (1|species),
                                  data = trends_scaled,
                                  weights=1/mean_se)
summary(year_mean_model)
anova(year_mean_model)
r.squaredGLMM(year_mean_model)

temp_mean_model <- lmer(mean_slope_temp ~ 0 + pheno_group + seasonality + mean_temp + pheno_position + (1|dataset:site) + (1|species),
                                  data = trends_scaled,
                                  weights = 1/mean_se_temp)
summary(temp_mean_model)
anova(temp_mean_model)
r.squaredGLMM(temp_mean_model)

year_sigma_model <- lmer(sigma_slope ~ 0 + pheno_group + seasonality + mean_temp + pheno_position + (1|dataset:site) + (1|species),
                        data = trends_scaled,
                        weights=1/sigma_se)
summary(year_sigma_model)
anova(year_sigma_model)
r.squaredGLMM(year_sigma_model)

temp_sigma_model <- lmer(sigma_slope_temp ~ 0 + pheno_group + seasonality + mean_temp + pheno_position + (1|dataset:site) + (1|species),
                        data = trends_scaled,
                        weights = 1/sigma_se_temp)
summary(temp_sigma_model)
anova(temp_sigma_model)
r.squaredGLMM(temp_sigma_model)


### Trait analysis
# plants
plant_traits <- read.csv("clean_data/plant_traits.csv", stringsAsFactors = FALSE)
plant_traits <- as.data.table(plant_traits)
#plant_traits <- plant_traits[,c("X", "tax_level", "whole.plant.growth.form", "whole.plant.height", "seed.mass", "leaf.area", "leaf.area.per.leaf.dry.mass")]
plant_traits$growth_form <- mgsub(plant_traits$growth_form,
                                              c("forb", "geophyte", "fern", "hemicryptophyte", "herb\\*"),
                                              c("herb", "herb", "herb", "herb", "herb")) # I checked that all of the hemicryptophytes are herbs
plant_traits$growth_form <- mgsub(plant_traits$growth_form,
                                              c("vine", "climbing legume", "non-woody epiphyte", "climber", "liana", "hemiepiphyte", "parasite"),
                                              c("dependant", "dependant", "dependant", "dependant", "dependant", "dependant", "dependant")) # This is all species that need others to grow properly
plant_traits$growth_form <- mgsub(plant_traits$growth_form,
                                              c("shrub\\*", "graminoid", "aquatic\\*", "3", "4", "cactus", "not found"),
                                              c("shrub", "grass", "NA", "NA", "NA", "NA", "NA")) # misc. groups. Getting rid of the 3 aqautic & cactus species
#colnames(plant_traits) <- c("species", "tax_level", "growth_form", "height", "seed_mass", "leaf_area", "leaf_area_per_mass")
table(plant_traits$tax_level) # almost all down to species
table(plant_traits$growth_form)
plant_traits <- transform(plant_traits,
                          tax_level = as.factor(tax_level),
                          growth_form = as.factor(growth_form),
                          height = as.numeric(height),
                          seed_mass = as.numeric(seed_mass),
                          leaf_area = as.numeric(leaf_area),
                          leaf_area_per_mass = as.numeric(leaf_area_per_mass))
plant_traits$species <- gsub("\\.", " ", plant_traits$species)
# I've pulled out just traits that have more than 50% coverage
# - growth form (collapsed into herb, shrub, tree, grass, dependant)
# - height
# - seed mass
# - leaf area/mass
plant_traits <- plant_traits[tax_level != "not found",]
summary(plant_traits)

# matching up trait data to pheno records
plant_data <- trends_focus[species %in% plant_traits$species,]
plant_data <- plant_data[tax_group == "plant",] # just making sure there aren't any hangers-on
plant_data$growth_form <- plant_traits[match(plant_data$species, plant_traits$species),"growth_form"]
plant_data$height <- plant_traits[match(plant_data$species, plant_traits$species),"height"]
plant_data$seed_mass <- plant_traits[match(plant_data$species, plant_traits$species),"seed_mass"]
plant_data$leaf_area <- plant_traits[match(plant_data$species, plant_traits$species),"leaf_area"]
plant_data$leaf_area_per_mass <- plant_traits[match(plant_data$species, plant_traits$species),"leaf_area_per_mass"]

# stupid fix to an "NA" vs NA problem
plant_data[which(plant_data$growth_form == "NA"), "growth_form"] <- NA
plant_data <- droplevels(plant_data)

cor(plant_traits[,c("height", "seed_mass", "leaf_area", "leaf_area_per_mass")], use="complete.obs") # no major correlation between traits - biggest is seed mass:height
table(plant_data[,.(phenophase, growth_form)]) # there are mostly just trees with leaf data

flower_data_scaled <- plant_data[phenophase == "first_flower",] # to look at traits
tree_data_scaled <- plant_data[growth_form %in% c("tree", "shrub", "herb"),] # to look at phenophases
flower_data_scaled[,c("height", "seed_mass", "leaf_area", "leaf_area_per_mass")] <- lapply(flower_data_scaled[,c("height", "seed_mass", "leaf_area", "leaf_area_per_mass")], scale.data, log=T)
tree_data_scaled[,c("height", "seed_mass", "leaf_area", "leaf_area_per_mass")] <- lapply(tree_data_scaled[,c("height", "seed_mass", "leaf_area", "leaf_area_per_mass")], scale.data, log=T)


### flower phenophase only - traits model

# mean shift model over time
flower_traits_year_mean <- lmer(mean_slope ~ growth_form + height + seed_mass + leaf_area_per_mass + (1|dataset:site) + (1|species),
                                data = flower_data_scaled, weights=(1/mean_se))
summary(flower_traits_year_mean)
anova(flower_traits_year_mean)
vif(flower_traits_year_mean)
r.squaredGLMM(flower_traits_year_mean)

# mean over temp
flower_traits_temp_mean <- lmer(mean_slope_temp ~ growth_form + height + seed_mass + leaf_area_per_mass + (1|dataset:site) + (1|species),
                                data = flower_data_scaled, weights=(1/mean_se_temp))
summary(flower_traits_temp_mean)
anova(flower_traits_temp_mean)
vif(flower_traits_temp_mean)
r.squaredGLMM(flower_traits_temp_mean)

# sigma over year
flower_traits_year_sigma <- lmer(sigma_slope ~ growth_form + height + seed_mass + leaf_area_per_mass + (1|dataset:site) + (1|species),
                               data = flower_data_scaled, weights=(1/sigma_se))
summary(flower_traits_year_sigma)
anova(flower_traits_year_sigma)
vif(flower_traits_year_sigma)
r.squaredGLMM(flower_traits_year_sigma)

# sigma over year
flower_traits_temp_sigma <- lmer(sigma_slope_temp ~ growth_form + height + seed_mass + leaf_area_per_mass + (1|dataset:site) + (1|species),
                                data = flower_data_scaled, weights=(1/sigma_se_temp))
summary(flower_traits_temp_sigma)
anova(flower_traits_temp_sigma)
vif(flower_traits_temp_sigma)
r.squaredGLMM(flower_traits_temp_sigma)


### leaves vs flowers - just trees, shrubs, and herbs

# mean shift model over time
plant_phase_year_mean <- lmer(mean_slope ~ growth_form*phenophase + (1|site) + (1|species),
                               data = tree_data_scaled, weights=(1/mean_se))
summary(plant_phase_year_mean)
anova(plant_phase_year_mean)
#visreg(plant_phase_year_mean, "phenophase", by="growth_form")
plant_phase_temp_mean <- lmer(mean_slope_temp ~ growth_form*phenophase + (1|site) + (1|species),
                              data = tree_data_scaled, weights=(1/mean_se_temp))
summary(plant_phase_temp_mean)
anova(plant_phase_temp_mean)
#visreg(plant_phase_temp_mean, "phenophase", by="growth_form", points.par=list(col="#f0a21d"))
plant_phase_year_sigma <- lmer(sigma_slope ~ growth_form*phenophase + (1|site) + (1|species),
                              data = tree_data_scaled, weights=(1/sigma_se))
summary(plant_phase_year_sigma)
anova(plant_phase_year_sigma)
#visreg(plant_phase_year_sigma, "phenophase", by="growth_form")
plant_phase_temp_sigma <- lmer(sigma_slope_temp ~ growth_form*phenophase + (1|site) + (1|species),
                              data = tree_data_scaled, weights=(1/sigma_se_temp))
summary(plant_phase_temp_sigma)
anova(plant_phase_temp_sigma)
#visreg(plant_phase_temp_sigma, "phenophase", by="growth_form")


### Bird traits
bird_traits <- read.csv("clean_data/bird_traits.csv", stringsAsFactors = FALSE)
# merging diet types into 3 basics: herbi, omni, and carni
bird_traits$diet <- mgsub(bird_traits$diet, c("PlantSeed", "FruiNect", "Invertebrate", "VertFishScav"), c("Herbivore", "Herbivore", "Carnivore", "Carnivore"))
bird_traits$diet <- as.factor(bird_traits$diet)

#bird_data <- trends[which(trends$species %in% bird_traits$species), ]
bird_data <- trends[which(trends$species %in% bird_traits$species | trends$species %in% bird_traits$species_synon | trends$species %in% bird_traits$species_synon_clean), ]
bird_data$diet <- sapply(bird_data$species, function(x) bird_traits[which(bird_traits$species == x | bird_traits$species_synon == x | bird_traits$species_synon_clean == x),"diet"])
bird_data$mass <- sapply(bird_data$species, function(x) bird_traits[which(bird_traits$species == x | bird_traits$species_synon == x | bird_traits$species_synon_clean == x),"mass"])
bird_data$species_synon_clean <- sapply(bird_data$species, function(x) bird_traits[which(bird_traits$species == x | bird_traits$species_synon == x | bird_traits$species_synon_clean == x),"species_synon_clean"])
bird_data$log_mass <- log(bird_data$mass)

#bird_data <- bird_data[which(bird_data$mean_slope_temp < 50 & bird_data$mean_slope_temp > -50),]
bird_data_scaled <- bird_data
bird_data_scaled$mass <- scale.data(bird_data$mass)
bird_data_scaled$log_mass <- scale.data(bird_data$log_mass)

bird_traits_model_mean_time <- lmer(mean_slope ~ diet + log_mass + (1|dataset:site) + (1|species),
                                    data = bird_data_scaled, weights=1/mean_se)
summary(bird_traits_model_mean_time)
anova(bird_traits_model_mean_time)
vif(bird_traits_model_mean_time)

bird_traits_model_mean_temp <- lmer(mean_slope_temp ~ diet + log_mass + (1|dataset:site) + (1|species),
                                    data = bird_data_scaled, weights=1/mean_se_temp)
summary(bird_traits_model_mean_temp)
anova(bird_traits_model_mean_temp)
vif(bird_traits_model_mean_temp)

bird_traits_model_sigma_time <- lmer(sigma_slope ~ diet + log_mass + (1|dataset:site) + (1|species),
                                     data = bird_data_scaled, weights=1/sigma_se)
summary(bird_traits_model_sigma_time)
anova(bird_traits_model_sigma_time)
vif(bird_traits_model_sigma_time)

bird_traits_model_sigma_temp <- lmer(sigma_slope_temp ~ diet + log_mass + (1|dataset:site) + (1|species),
                                    data = bird_data_scaled, weights=1/sigma_se_temp)
summary(bird_traits_model_sigma_temp)
anova(bird_traits_model_sigma_temp)
vif(bird_traits_model_sigma_temp)


beep(4)

### model coefficient tables

# main model
mean_shift<- data.table(summary(year_mean_model)$coefficients, keep.rownames = "Coefficient"); mean_shift[, Metric := "μ shift"]
mean_sens<- data.table(summary(temp_mean_model)$coefficients, keep.rownames = "Coefficient"); mean_sens[, Metric := "μ sensitivity"]
sigma_shift<- data.table(summary(year_sigma_model)$coefficients, keep.rownames = "Coefficient"); sigma_shift[, Metric := "σ shift"]
sigma_sens<- data.table(summary(temp_sigma_model)$coefficients, keep.rownames = "Coefficient"); sigma_sens[, Metric := "σ sensitivity"]

coefs_export <- rbindlist(list(mean_shift, mean_sens, sigma_shift, sigma_sens))
setcolorder(coefs_export, c(7,1:6)); setnames(coefs_export, c("t value", "Pr(>|t|)"), c("t-value", "p-value"))
coefs_export[, c("Estimate", "Std. Error", "df", "t-value", "p-value") := lapply(.SD, round, 3), .SDcols = c("Estimate", "Std. Error", "df", "t-value", "p-value")]
coefs_export[, `p-value` := as.character(`p-value`)]
coefs_export[, `p-value` := ifelse(`p-value` == "0", "< 0.001", `p-value`)]
coefs_export[, Coefficient := gsub("pheno_group", "", Coefficient)]
coefs_export <- flextable(coefs_export, cwidth = 1); coefs_export <- theme_vanilla(coefs_export)
coefs_export <- merge_at(coefs_export, i = c(1:7), j = 1); coefs_export <- merge_at(coefs_export, i = c(8:14), j = 1); coefs_export <- merge_at(coefs_export, i = c(15:21), j = 1); coefs_export <- merge_at(coefs_export, i = c(22:28), j = 1) # there must be a more efficient way to do this
coefs_export
save_as_html(coefs_export, path="manuscript/main_model_table.html")

# plant traits model (just flowers)
mean_shift<- data.table(summary(flower_traits_year_mean)$coefficients, keep.rownames = "Coefficient"); mean_shift[, Metric := "μ shift"]
mean_sens<- data.table(summary(flower_traits_temp_mean)$coefficients, keep.rownames = "Coefficient"); mean_sens[, Metric := "μ sensitivity"]
sigma_shift<- data.table(summary(flower_traits_year_sigma)$coefficients, keep.rownames = "Coefficient"); sigma_shift[, Metric := "σ shift"]
sigma_sens<- data.table(summary(flower_traits_temp_sigma)$coefficients, keep.rownames = "Coefficient"); sigma_sens[, Metric := "σ sensitivity"]

coefs_export <- rbindlist(list(mean_shift, mean_sens, sigma_shift, sigma_sens))
setcolorder(coefs_export, c(7,1:6)); setnames(coefs_export, c("t value", "Pr(>|t|)"), c("t-value", "p-value"))
coefs_export[, c("Estimate", "Std. Error", "df", "t-value", "p-value") := lapply(.SD, round, 3), .SDcols = c("Estimate", "Std. Error", "df", "t-value", "p-value")]
coefs_export[, `p-value` := as.character(`p-value`)]
coefs_export[, `p-value` := ifelse(`p-value` == "0", "< 0.001", `p-value`)]
coefs_export[, Coefficient := gsub("growth_form", "", Coefficient)]
coefs_export[, Coefficient := gsub("leaf_area_per_mass", "SLA", Coefficient)]
coefs_export <- flextable(coefs_export, cwidth = 1); coefs_export <- theme_vanilla(coefs_export)
coefs_export <- merge_at(coefs_export, i = c(1:8), j = 1); coefs_export <- merge_at(coefs_export, i = c(9:16), j = 1); coefs_export <- merge_at(coefs_export, i = c(17:24), j = 1); coefs_export <- merge_at(coefs_export, i = c(25:32), j = 1) # there must be a more efficient way to do this
coefs_export
save_as_html(coefs_export, path="manuscript/flower_traits_table.html")

# plant phenophase model
mean_shift<- data.table(summary(plant_phase_year_mean)$coefficients, keep.rownames = "Coefficient"); mean_shift[, Metric := "μ shift"]
mean_sens<- data.table(summary(plant_phase_temp_mean)$coefficients, keep.rownames = "Coefficient"); mean_sens[, Metric := "μ sensitivity"]
sigma_shift<- data.table(summary(plant_phase_year_sigma)$coefficients, keep.rownames = "Coefficient"); sigma_shift[, Metric := "σ shift"]
sigma_sens<- data.table(summary(plant_phase_temp_sigma)$coefficients, keep.rownames = "Coefficient"); sigma_sens[, Metric := "σ sensitivity"]

coefs_export <- rbindlist(list(mean_shift, mean_sens, sigma_shift, sigma_sens))
setcolorder(coefs_export, c(7,1:6)); setnames(coefs_export, c("t value", "Pr(>|t|)"), c("t-value", "p-value"))
coefs_export[, c("Estimate", "Std. Error", "df", "t-value", "p-value") := lapply(.SD, round, 3), .SDcols = c("Estimate", "Std. Error", "df", "t-value", "p-value")]
coefs_export[, `p-value` := as.character(`p-value`)]
coefs_export[, `p-value` := ifelse(`p-value` == "0", "< 0.001", `p-value`)]
coefs_export[, Coefficient := gsub("growth_form", "", Coefficient)]
coefs_export[, Coefficient := gsub("phenophase", "", Coefficient)]
coefs_export <- flextable(coefs_export, cwidth = 1); coefs_export <- theme_vanilla(coefs_export)
coefs_export <- merge_at(coefs_export, i = c(1:4), j = 1); coefs_export <- merge_at(coefs_export, i = c(5:8), j = 1); coefs_export <- merge_at(coefs_export, i = c(9:12), j = 1); coefs_export <- merge_at(coefs_export, i = c(13:16), j = 1) # there must be a more efficient way to do this
coefs_export
save_as_html(coefs_export, path="manuscript/phenophase_table.html")

# bird traits model
mean_shift<- data.table(summary(bird_traits_model_mean_time)$coefficients, keep.rownames = "Coefficient"); mean_shift[, Metric := "μ shift"]
mean_sens<- data.table(summary(bird_traits_model_mean_temp)$coefficients, keep.rownames = "Coefficient"); mean_sens[, Metric := "μ sensitivity"]
sigma_shift<- data.table(summary(bird_traits_model_sigma_time)$coefficients, keep.rownames = "Coefficient"); sigma_shift[, Metric := "σ shift"]
sigma_sens<- data.table(summary(bird_traits_model_sigma_temp)$coefficients, keep.rownames = "Coefficient"); sigma_sens[, Metric := "σ sensitivity"]

coefs_export <- rbindlist(list(mean_shift, mean_sens, sigma_shift, sigma_sens))
setcolorder(coefs_export, c(7,1:6)); setnames(coefs_export, c("t value", "Pr(>|t|)"), c("t-value", "p-value"))
coefs_export[, c("Estimate", "Std. Error", "df", "t-value", "p-value") := lapply(.SD, round, 3), .SDcols = c("Estimate", "Std. Error", "df", "t-value", "p-value")]
coefs_export[, `p-value` := as.character(`p-value`)]
coefs_export[, `p-value` := ifelse(`p-value` == "0", "< 0.001", `p-value`)]
coefs_export[, Coefficient := gsub("diet", "", Coefficient)]
coefs_export[, Coefficient := gsub("log_mass", "log(mass)", Coefficient)]
coefs_export <- flextable(coefs_export, cwidth = 1); coefs_export <- theme_vanilla(coefs_export)
coefs_export <- merge_at(coefs_export, i = c(1:4), j = 1); coefs_export <- merge_at(coefs_export, i = c(5:8), j = 1); coefs_export <- merge_at(coefs_export, i = c(9:12), j = 1); coefs_export <- merge_at(coefs_export, i = c(13:16), j = 1) # there must be a more efficient way to do this
coefs_export
save_as_html(coefs_export, path="manuscript/bird_traits_table.html")





