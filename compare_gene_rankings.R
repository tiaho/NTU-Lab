# compare_gene_rankings.R
# ranks the genes based on several features, ranks by total number of z scores away from mean
# ranks based on cds length, exon length, intron length, codon usage
#
# need to avg the z scores for exons and introns

setwd("~/Desktop/Lab/DataFiles/")

# reads in files
cds_length <- read.table("cds_lengths.txt", sep = ",")
exon_length <- read.table("exon_lengths.txt", sep = ",")
intron_length <- read.table("intron_lengths.txt", sep = ",")
codon_usage <- read.table("gene_codon_usage_relative_entropies.txt", sep = ",")

# puts the data frames of the parameters into an array
parameters <- c(cds_length, exon_length, intron_length, codon_usage)

# renames the columns
colnames(cds_length) <- c('length', 'chromosome', 'gene')
colnames(exon_length) <- c('length', 'chromosome', 'gene')
colnames(intron_length) <- c('length', 'chromosome', 'gene')
colnames(codon_usage) <- c('length', 'chromosome', 'gene')

# calculates z score
cds_length$zscore <- (cds_length$length - mean(cds_length$length)) / sd(cds_length$length)
cds_length[[4]] <- (cds_length[[1]] - mean(cds_length[[1]])) / sd(cds_length[[1]])



# summaries
summary(cds_length)
summary(exon_length)
summary(intron_length)
summary(codon_usage)

