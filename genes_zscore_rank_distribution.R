# genes_zscore_rank_distribution.R
# plots the distribution of the ranking of the genes (z-score method)

setwd("~/Desktop/Lab/DataFiles/")

zscores = read.table("log_zscore_rankings.txt", sep = ",")
colnames(zscores) = c("zscore", "chromosome", "gene")

summary(zscores$zscore)
IQR <- quantile(zscores$zscore, .75) - quantile(zscores$zscore, .25)
upper_bound <- quantile(zscores$zscore, .75) + 1.5 * IQR
lower_bound <- quantile(zscores$zscore, .25) - 1.5 * IQR

library(lattice)
densityplot(zscores$zscore,
            main = list('Density Plot of Zscores of Gene Rankings', cex = 2),
            xlab= list('Zscores', cex = 1.5),
            ylab= list(cex = 1.5),
            scales = list(cex = 1.2),
            key = list(text = list(c('Mean', 'Upper Bound'), cex = 1.5),
                       lines = list(type = 'l', lty = 1, col = c('green', 'red')),
                       border = TRUE,
                       columns = 2),
            panel = function(...){
              panel.densityplot(...)
              panel.abline(v = upper_bound, col = 'red')
              panel.abline(v = mean(zscores$zscore), col = 'green')
            }
)

# zoomed in
densityplot(zscores$zscore,
            main = list('Density Plot of Zscores of Gene Rankings', cex = 2),
            xlab= list('Zscores', cex = 1.5),
            ylab= list(cex = 1.5),
            xlim = c(0, 2),
            scales = list(cex = 1.2),
            key = list(text = list(c('Mean', 'Upper Bound'), cex = 1.5),
                       lines = list(type = 'l', lty = 1, col = c('green', 'red')),
                       border = TRUE,
                       columns = 2),
            panel = function(...){
              panel.densityplot(...)
              panel.abline(v = upper_bound, col = 'red')
              panel.abline(v = mean(zscores$zscore), col = 'green')
            }
)