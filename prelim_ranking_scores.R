setwd("~/Desktop/Lab/DataFiles/")

scores <- read.table("prelim_ranking_scores.txt", sep = " ")
colnames(scores) <- c("gene", "abs", "rel")

library(lattice)
plot(as.numeric(as.character(scores$abs)), as.numeric(as.character(scores$rel)))

smoothScatter(as.numeric(as.character(scores$abs)), as.numeric(as.character(scores$rel)))
