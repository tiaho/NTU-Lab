# cds_codon_composition_relative_entropies.R
# plots the distribution of the KL distances for the codon composition in all genes

setwd('~/Desktop/Lab/DataFiles/')

codon_usage <- read.table('gene_codon_usage_relative_entropies.txt', sep = ",")
colnames(codon_usage) <- c('entropy', 'id')
codon_usage <- codon_usage[order(codon_usage$entropy, decreasing = FALSE), ] # sorts descending numerically

summary(codon_usage$entropy)

library(lattice)
densityplot(codon_usage$entropy,
            main= list('Density of Relative Entropy of Gene Codon Usage', cex = 2),
            xlab = list('Relative Entropy', cex = 1.5),
            ylab = list(cex = 1.5),
            scales = list(cex = 1.2))
densityplot(codon_usage$entropy,
            main= list('Density of Relative Entropy of Gene Codon Usage', cex = 2),
            xlab = list('Relative Entropy', cex = 1.5),
            ylab = list(cex = 1.5),
            scales = list(cex = 1.2),
            xlim = c(0, 1.5))

library(ggplot2)
qplot(entropy, data = codon_usage, geom = 'density')
