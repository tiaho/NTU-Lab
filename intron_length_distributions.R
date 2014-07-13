# intron_length_distribution.R
# plots the distribution of intron lengths

setwd('~/Desktop/Lab/DataFiles/')

intron_lengths <- read.table('intron_lengths.txt', sep = ",")
colnames(intron_lengths) <- 'lengths'

summary(intron_lengths$lengths)
IQR <- quantile(intron_lengths$lengths, .75) - quantile(intron_lengths$lengths, .25)
upper_bound <- quantile(intron_lengths$lengths, .75) + 1.5 * IQR
lower_bound <- quantile(intron_lengths$lengths, .25) - 1.5 * IQR

library(lattice)
densityplot(intron_lengths$lengths,
            main = list('Density Plot of Intron Lengths', cex = 2),
            xlab= list('Intron Lengths', cex = 1.5),
            ylab= list(cex = 1.5),
            scales = list(cex = 1.2),
            key = list(text = list(c('Mean', 'Upper Bound'), cex = 1.5),
                       lines = list(type = 'l', lty = 1, col = c('green', 'red')),
                       border = TRUE,
                       columns = 2),
            panel = function(...){
              panel.densityplot(...)
              panel.abline(v = upper_bound, col = 'red')
              panel.abline(v = mean(intron_lengths$lengths), col = 'green')
            }
)
