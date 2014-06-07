# compare_gene_rankings.R
# calculates a z score for each gene based on various factors (does same thing as rank_genes_zscores.pl)

# ranks the genes based on several features, ranks by total number of z scores away from mean
# ranks based on cds length, exon length, intron length, codon usage
#
# DON'T USE THIS

setwd("~/Desktop/Lab/DataFiles/")

# reads in files
cds_length <- read.table("cds_lengths.txt", sep = ",")
exon_length <- read.table("exon_lengths.txt", sep = ",")
intron_length <- read.table("intron_lengths.txt", sep = ",")
codon_usage <- read.table("gene_codon_usage_relative_entropies.txt", sep = ",")
donor <- read.table("splice_site_donor_usage.txt", sep = ",")
acceptor <- read.table("splice_site_acceptor_usage.txt", sep = ",")
intron_usage <- read.table("intron_usage.txt", sep = ",")

# log transformation
cds_length[[4]] <- log(cds_length[[1]])
exon_length[[4]] <- log(exon_length[[1]])
intron_length[[4]] <- log(intron_length[[1]])
codon_usage[[4]] <- log(codon_usage[[1]])
donor[[4]] <- log(donor[[1]])
acceptor[[4]] <- log(acceptor[[1]])
intron_usage[[4]] <- log(intron_usage[[1]])

# calculates z score
cds_length[[5]] <- (cds_length[[4]] - mean(cds_length[[4]])) / sd(cds_length[[4]])
exon_length[[5]] <- (exon_length[[4]] - mean(exon_length[[4]])) / sd(exon_length[[4]])
intron_length[[5]] <- (intron_length[[4]] - mean(intron_length[[4]])) / sd(intron_length[[4]])
codon_usage[[5]] <- (codon_usage[[4]] - mean(codon_usage[[4]])) / sd(codon_usage[[4]])
donor[[5]] <- (donor[[4]] - mean(donor[[4]])) / sd(donor[[4]])
acceptor[[5]] <- (acceptor[[4]] - mean(acceptor[[4]])) / sd(acceptor[[4]])
intron_usage[[5]] <- (intron_usage[[4]] - mean(intron_usage[[4]])) / sd(intron_usage[[4]])

cds_length$V2 <- as.character(cds_length$V2)
  
# exports file
write.table(cds_length, "log_cds_lengths.csv", row.names = FALSE, col.names = TRUE, sep = ",")

# merges the files into one

parameters <- c(cds_length, exon_length, intron_length, codon_usage, donor, acceptor, intron_usage)
