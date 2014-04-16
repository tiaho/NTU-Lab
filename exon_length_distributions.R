setwd('~/Desktop/Lab/DataFiles/')

exon_lengths <- read.table('exon_lengths.txt', sep = ",")
colnames(exon_lengths) <- 'lengths'

summary(exon_lengths$lengths)
IQR <- quantile(exon_lengths$lengths, .75) - quantile(exon_lengths$lengths, .25)
upper_bound <- quantile(exon_lengths$lengths, .75) + 1.5 * IQR
lower_bound <- quantile(exon_lengths$lengths, .25) - 1.5 * IQR

library(lattice)
densityplot(exon_lengths$lengths,
            main = list('Density Plot of Exon Lengths', cex = 2),
            xlab= list('Exon Lengths', cex = 1.5),
            ylab= list(cex = 1.5),
            scales = list(cex = 1.2),
            key = list(text = list(c('Mean', 'Upper Bound'), cex = 1.5),
                       lines = list(type = 'l', lty = 1, col = c('green', 'red')),
                       border = TRUE,
                       columns = 2),
            panel = function(...){
              panel.densityplot(...)
              panel.abline(v = upper_bound, col = 'red')
              panel.abline(v = mean(exon_lengths$lengths), col = 'green')
            }
)