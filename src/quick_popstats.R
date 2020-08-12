
library(adegenet)
library(vcfR)
library(data.table)
library(ggplot2)


vcf_file <- "output/020_filtered/hyp/pruned.vcf.gz"


vcf <- read.vcfR(vcf_file)

vcf_gl <- vcfR2genlight(vcf)

# drop aeth?
snp_data <- vcf_gl[grepl("mhyp", vcf_gl$ind.names),]

# main scaffolds only
vcf_scaf <- snp_data[, startsWith(snp_data$loc.names, "PGA_scaffold")]

# impute
na_means <- tab(vcf_scaf, NA.method = "mean")
snps_imputed <- new("genlight",
                    na_means,
                    ploidy = 2)
pop(snps_imputed) <- sapply(snps_imputed$ind.names, function(x)
    paste(unlist(strsplit(x, split = "_"))[2:3], collapse = "_"))


glPlot(snps_imputed)


pops <- sapply(colnames(vcf@gt), function(x)
    paste(unlist(strsplit(x, split = "_"))[2:3], collapse = "_"))

my_diff <- genetic_diff(vcf, as.factor(pops[-1]))
diff_dt <- data.table(my_diff)

#melt
pd <- melt(diff_dt,
     id.vars = c("CHROM", "POS"),
     measure.vars = c("Hs_maethm_lincoln", "Hs_mhyp_lincoln", "Hs_mhyp_ruakura", "Ht"))

ggplot(pd, aes(x = variable, y = value)) +
    geom_boxplot()


## No id variables; using all as measure variables

p <- ggplot(diff_dt, aes(x=variable, y=Ht)) + geom_violin(fill="#2ca25f", adjust = 1.2)

# [startsWith(CHROM, "PGA_scaffold")]

# set up coords for dirty plot
# chromosome coordinates
chr_coords <- diff_dt[, .(chr_length = as.numeric(max(POS))),
                      by = CHROM]
setorder(chr_coords, CHROM)
chr_coords[, chr_end := cumsum(chr_length)]
chr_coords[, chr_start := chr_end - chr_length + 1]
chr_coords[, lab_pos := chr_start + round(mean(c(chr_start, chr_end)), 0), by = Chr]
pos_to_label = chr_coords[, seq(1, max(chr_end), length.out = n_to_label)]

label_positions <- sapply(pos_to_label, function(x)
    chr_coords[, .I[which.min(abs(x - lab_pos))]])

chr_coords[label_positions, x_lab := Chr]
chr_coords[is.na(x_lab), x_lab := ""]

g_coords <- merge(diff_dt,
                   chr_coords,
                   by = "CHROM",
                   all.x = TRUE,
                   all.y = FALSE)

g_coords[, POS := as.numeric(POS)]
setorder(g_coords, -chr_length, CHROM, POS)
g_coords[, bp_coord := POS + chr_start - 1]





ggplot(g_coords, aes(x = bp_coord, y = Gprimest)) +
    geom_point()

