setwd('~/Desktop/Lab/DataFiles/')

# reads in the codon usage counts for all coding sequences
codon_usage <- read.table('codon_usage_counts.txt', sep = ",")
colnames(codon_usage) <- c('codon', 'count', 'type')
codon_usage <- codon_usage[order(codon_usage$codon, decreasing = FALSE), ] # sorts alphabetically

# converts the counts to frequency
total_num_codons <- sum(codon_usage$count)
codon_usage$count <- codon_usage$count / total_num_codons


# reads in the codon usage counts for a particular gene and adds it to the data frame
setwd('~/Desktop/Lab/DataFiles/CDSCodonCounts/')
gene_list <- c('T24B8.3c', 'F31B12.1a', 'ZK381.62')
for (gene in gene_list){
  gene <- read.table(gene, sep = ",")
  colnames(gene) <- c('codon', 'count', 'type')
  gene_total <- sum(gene$count)
  gene$count <- gene$count / gene_total
  codon_usage <- rbind(codon_usage, gene)
}

# plots
library(ggplot2)
ggplot(codon_usage, aes(codon, count)) + geom_point(aes(color = factor(type), size = 2)) + ggtitle("Codon Usage Frequency")
