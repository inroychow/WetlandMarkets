#-------------------------------------------#
#- Saving geometry for outlier large banks  #
#-------------------------------------------#

library(sf)

bank_acres = read.csv("Footprints/Bank Acres.csv")
bank_acres$Name <- gsub(" ", "_", bank_acres$Name)

banks = readRDS("Footprints/footprints_and_buffers.rds")
banks <- merge(banks, bank_acres, by.x = "Name", by.y = "Name", all.x = TRUE)
banks <- banks[!is.na(banks$Bank.Acres), ]
banks = banks[c("Name", "Bank.Acres")]

# Select banks that are erroneously large
banks$geom_area <- expanse(banks, unit = "m")
quantile(banks$geom_area, c(0.95, 0.99))
huge_banks <- banks[banks$geom_area > quantile(banks$geom_area, 0.995), ]
plot(huge_banks)

selected_banks <- banks[banks$Name %in% huge_banks$Name, ]
centroids <- centroids(selected_banks)

acre_to_m2 <- 4046.86
radii <- sqrt((huge_banks$Bank.Acres * acre_to_m2) / pi)

# Create buffers from centroids using those radii
buffers <- buffer(centroids, width = radii)

buffers$Name <- huge_banks$Name
buffers$Bank.Acres <- huge_banks$Bank.Acres

saveRDS(buffers, "Bank Footprints/outlier_banks.rds")

test = huge_banks[huge_banks$Name == "Upper_Coldwater_Mitigation_Bank"]
plot(test)
