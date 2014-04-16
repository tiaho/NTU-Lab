setwd('~/Desktop/Lab/DataFiles/')

# imports exons
exon_lengths <- read.table('exon_lengths.txt', sep = ",")
colnames(exon_lengths) <- c('length', 'chromosome', 'cdsID', 'type')
summary(exon_lengths$length)

# imports introns
intron_lengths <- read.table('intron_lengths.txt', sep = ",")
colnames(intron_lengths) <- c('length', 'chromosome', 'cdsID', 'type')
summary(intron_lengths$length)

# combines the two data frames
exon_intron_lengths <- rbind(exon_lengths, intron_lengths)

library(ggplot2)
qplot(length, data = exon_intron_lengths, geom = 'density', color = type, xlim = c(0, 2000),
      main = 'Exon and Intron Length Distribution', xlab = 'Length in base pairs', ylab = 'Density')
