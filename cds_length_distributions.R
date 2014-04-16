setwd('~/Desktop/Lab/DataFiles/')

cds_lengths <- read.table('cds_lengths.txt', sep = ",")
colnames(cds_lengths) <- 'lengths'

summary(cds_lengths$lengths)
IQR <- quantile(cds_lengths$lengths, .75) - quantile(cds_lengths$lengths, .25)
upper_bound <- quantile(cds_lengths$lengths, .75) + 1.5 * IQR
lower_bound <- quantile(cds_lengths$lengths, .25) - 1.5 * IQR

library(lattice)
densityplot(cds_lengths$lengths,
            main = list('Density Plot of CDS Lengths', cex = 2),
            xlab= list('CDS Lengths', cex = 1.5),
            ylab= list(cex = 1.5),
            scales = list(cex = 1.2),
            key = list(text = list(c('Median', 'Mean', 'Upper Bound'), cex = 1.5),
                       lines = list(type = 'l', lty = 1, col = c('orange', 'green', 'red')),
                       border = TRUE,
                       columns = 3),
            panel = function(...){
              panel.densityplot(...)
              panel.abline(v = upper_bound, col = 'red')
              panel.abline(v = mean(cds_lengths$lengths), col = 'green')
              panel.abline(v = quantile(cds_lengths$lengths, 0.5), col = 'orange')
            }
)
