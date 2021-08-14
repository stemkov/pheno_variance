### Cleaning Rothamsted data

#setwd("/home/michael/Documents/Grad School/Research Projects/Phenology variability")

data <- fread("raw_data/rothamsted/AphidData.csv", stringsAsFactors = FALSE)
data <- data[!is.na(`5% flight`),]

# cleaning up site names to match the coordinates spreadsheet
data$Trap <- mgsub(data$Trap,
                   c("Brooms' Barn S", "Elgin One", "Hereford \\(Rosemaund\\)", "Rothamsted", "Silwood", "Writtle S"),
                   c("Brooms Barn", "Elgin", "Rosemaund, Hereford", "Rothamsted Tower", "Silwood Park", "Writtle"))

# difference between 1st and 5%
data$diff <- data$`5% flight` - data$`1st flight`
data_scrubbed <- data[diff < 90]
#(nrow(data)-nrow(data_scrubbed))/nrow(data) # this scrubbing removes 1.5% of the data

#plot(`5% flight` ~ `1st flight`, data=data, xlab="1st flight DOY", ylab="5% flight DOY")
#plot(`5% flight` ~ `1st flight`, data=data[diff < 90], xlab="1st flight DOY", ylab="5% flight DOY")

formatted_data <- data.table(species = data_scrubbed$Binomial,
                             site = data_scrubbed$Trap,
                             year = data_scrubbed$Season,
                             doy = data_scrubbed$`1st flight`, 
                             dataset = "rothamsted")

write.csv(formatted_data, "clean_data/rothamsted_clean.csv", row.names = FALSE)

#cleanup
rm(data, data_scrubbed, formatted_data)

