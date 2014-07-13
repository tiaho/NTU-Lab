# modeling_significant_reads.R
# Tiffany Ho
# 07/09/2014

# models the NGS reads to determine what amount of reads is "significant"

# reads in data
ngs <- read.table("~/Desktop/NTU Lab/NGS data.csv", sep = ",", header = T)

# need to remove kingdom, genus, species reads
data <- subset(ngs, Level == "Phylum" | Level == "Class" | Level == "Order" | Level == "Family")

# makes a vector with the sample names
sample.names <- unique(data$SampleName)
# test.names <- sample.names[1:2]

# makes a vector with the level names
data.levels <- unique(data$Level)
# test.levels <- data.levels[1:2]

# makes a vector with the % reads
percent.of.reads <- seq(0.01, 0.50, by = 0.005)
# test.percent <- seq(0.49, 0.50, by = 0.005)

# summary vectors
summary.sample <- vector()
summary.level <- vector()
summary.percent.reads <- vector()
summary.stdev <- vector()

# 4 graphs on a page
par(mfrow = c(2,2))


# runs the simulation for each sample
# for (i in "ZT2-C2"){
for (i in sample.names){
  sample <- subset(data, SampleName == i)  # gets the sample
#   print(sum(sample$Reads))

  # makes 2 vectors to keep track of % reads and % of level represented
  reads.data <- vector()
  level.data <- vector()

  # runs the simulation for each of the levels (phylum, class, order, family)
#   for (j in "Phylum"){
  for (j in data.levels){
    sample.level <- subset(sample, Level == j) # only wants the data from individual levels
#     print(sum(sample.level$Reads))
#     print(sample.level.sorted <- sample.level[order(sample.level$Group),])

    # make a new data frame to sample from
    sample.level.data <- vector()

    # fills in the data frame to sample from
    for (k in 1:nrow(sample.level)){
      tmp <- rep(as.character(sample.level$Group[k]), sample.level$Reads[k]) # puts the reps for each group into a temp vector
      sample.level.data <- c(sample.level.data, tmp) # combines the existing vector with the temp vector
#       print(sample.level.data)
    }

    # runs the simulation for each of the % read values
#     for (m in c(0.50, 0.10)){
    for (m in percent.of.reads){
      if ((as.integer(m * length(sample.level.data))) < 1) {next} # skips this % of reads sample if it results in no samples
      samples <- sample(sample.level.data, as.integer(m * length(sample.level.data)), replace = FALSE) # samples from the data frame
#       print(paste(i, j, m * 100, sep = " - "))
      tmptable <- table(samples)/length(samples) * 100 # creates a table of % of each group represented
#       print(tmptable)

      # fills in 2 vectors with the % read values and % of each group represented
      for (n in 1:length(tmptable)){
        reads.data <- c(reads.data, m * length(sample.level.data))
        level.data <- c(level.data, tmptable[n])
      }

      # data frame with: sample, level, % read, standard deviation
      summary.sample <- c(summary.sample, i)
      summary.level <- c(summary.level, j)
      summary.percent.reads <- c(summary.percent.reads, m)
      summary.stdev <- c(summary.stdev, sd(tmptable))
    }

    # makes a data frame for the data to be plotted
    graph.data <- data.frame(reads.data, level.data)
    plot(graph.data, main = paste(i, j), xlab = "Number of Reads", ylab = paste("% of", j, "represented"))
  }
}

# makes a data frame of the summary
summary.df <- data.frame(summary.sample, summary.level, summary.percent.reads, summary.stdev)

# writes out the table to a csv file
# write.table(summary.df, file = "~/Desktop/NTU Lab/stdev.csv", sep = ",", row.names = F)
