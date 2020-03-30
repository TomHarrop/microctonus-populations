#!/usr/bin/env Rscript

library(data.table)
library(seqinr)

RevComp <- function(x){
    toupper(c2s(rev(comp(s2c(x)))))
}


sample_sheet_raw <- fread("setup/SampleSheet_MiSeq_E7600.csv", skip = 21)
all_indivs <- fread("setup/Microctonus_indiv_species_location.csv")
final_pool <- fread("setup/C818-microctonus-samples.csv")

# mung the barcodes in sample_sheet
munged_barcodes <- rbind(
    sample_sheet_raw[, .(index_id = I7_Index_ID,
                         index = index)],
    unique(sample_sheet_raw[, .(index_id = I5_Index_ID,
                                index = RevComp(index2)),
                            by = index2][, .(index_id, index)]))

# set up barcode sequences
samples_only <- final_pool[, .(
    sample = sample_id,
    i7_primer,
    i5_primer)]


samples_i7 <- merge(samples_only,
                    munged_barcodes[, .(index_id, i7_index = index)],
                    by.x = "i7_primer", by.y = "index_id")
samples_both <- merge(samples_i7,
                      munged_barcodes[, .(index_id, i5_index = index)],
                      by.x = "i5_primer", by.y = "index_id")
samples_both[, barcode_sequence := paste(i7_index, i5_index, sep = "+")]

# set up metadata
merged_data <- merge(samples_both,
                     all_indivs[, .(
                         sample = indiv,
                         para_name = species,
                         location = location
                     )],
                     by = "sample",
                     all.x = TRUE,
                     all.y = FALSE)

sample_table <- merged_data[, .(
    sample = paste(
        sample,
        gsub("_", "", para_name),
        location,
        sep = "_"),
    barcode = barcode_sequence,
    r1_path = paste0("data/reads/", sample, "/", sample, "_R1.fq.gz"),
    r2_path = paste0("data/reads/", sample, "/", sample, "_R2.fq.gz"),
    metadata = paste(
        gsub("_", "", para_name),
        location,
        sep = "_")
)]

fwrite(sample_table, "data/samples.csv")
