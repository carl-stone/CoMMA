## code to prepare `genome_sites_bed` dataset goes here
# Reorder columns
genome_sites_bed <- genome_sites[, c("Site", "Left", "Right", "Type", "Strand")]

# Add a dummy score column as required by the BED format (set all scores to 0)
genome_sites_bed$Score <- 0

# Save as BED file
write.table(genome_sites_bed, file = "genome_sites.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

usethis::use_data(genome_sites_bed, overwrite = TRUE)
