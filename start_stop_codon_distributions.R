# start_stop_codon_distributions.R
# plots the distribution of the start and stop codon usage

setwd('~/Desktop/Lab/DataFiles/')

start_codon <- read.csv('start_codon_distribution.csv', header = FALSE)
colnames(start_codon) <- c('codon', 'freq')
num_of_codons <- sum(start_codon$freq)
# start_codon$freq <- start_codon$freq/num_of_codons

stop_codon <- read.csv('stop_codon_distribution.csv', header = FALSE) 
colnames(stop_codon) <- c('codon', 'freq')
# stop_codon$freq <- stop_codon$freq/num_of_codons

library(lattice)
xyplot(start_codon$freq ~ start_codon$codon,
       main = 'Scatterplot of the Start Codon Distribution',
       xlab = 'Codon',
       ylab = 'Frequency',
       col = ifelse(start_codon$codon == 'atg', 'red', 'blue'))



xyplot(stop_codon$freq ~ stop_codon$codon,
       main = 'Scatterplot of the Stop Codon Distribution',
       xlab = 'Codon',
       ylab = 'Frequency')

