# modeling_significant_reads.R
# Tiffany Ho
# 07/09/2014

# models the NGS reads to determine what amount of reads is "significant"
  
# changes the working directory
setwd("~/Desktop/NTU Lab/")

# reads in data
ngs <- read.table("NGS data.csv", sep = ",", header = T)

# need to remove kingdom, genus, species reads
data <- subset(ngs, Level == "Phylum" | Level == "Class" | Level == "Order" | Level == "Family")

# makes a vector with the sample names
sample.names <- unique(data$SampleName)

# makes a vector with the level names
data.levels <- unique(data$Level)

# makes a vector with the % reads
percent.of.reads <- seq(0.10, 0.90, by = 0.10)

# summary vectors
summary.sample <- vector()
summary.level <- vector()
summary.percent.reads <- vector()
summary.number.of.reads <- vector()
summary.group <- vector()
summary.stdev <- vector()
summary.mean <- vector()

# summary data frame
summary.df <- data.frame()

# 4 graphs on a page
par(mfrow = c(2,2))


# runs the simulation for each sample
for (i in sample.names){
  sample <- subset(data, SampleName == i)  # gets the sample

  # makes 2 vectors to keep track of % reads and % of level represented
    level.data <- vector()
    percentage.data <- vector()
  
  # runs the simulation for each of the levels (phylum, class, order, family)
  for (j in data.levels){
    sample.level <- subset(sample, Level == j) # only wants the data from individual levels

    # make a new data frame to sample from
    sample.level.data <- vector()

    # fills in the data frame to sample from
    for (k in 1:nrow(sample.level)){
      tmp <- rep(as.character(sample.level$Group[k]), sample.level$Reads[k]) # puts the reps for each group into a temp vector
      sample.level.data <- c(sample.level.data, tmp) # combines the existing vector with the temp vector
    }

    # runs the simulation for each of the % read values
    for (m in percent.of.reads){
      if ((as.integer(m * length(sample.level.data))) < 1) {next} # skips this % of reads sample if it results in no samples
      
      tmptable <- vector()
    
      # runs the simulation 100 times
      for (n in 1:100){
        samples <- sample(sample.level.data, as.integer(m * length(sample.level.data)), replace = TRUE) # samples from the data frame
        tmptable <- table(samples)/length(samples) * 100

        # fills in 2 vectors with the level and % of each group represented
        level.data <- c(level.data, names(tmptable))
        percentage.data <- c(percentage.data, as.vector(tmptable))

        # prints the progress of the script
        print(paste(i, j, m, n))
      }
            
      # makes a data frame out of the level and percentage data
      level.percentage.data <- data.frame(level.data, percentage.data)

      # vector with the names of the different grousp for this level
      level.names <- unique(sample.level$Group)

      # calculates the mean and standard deviations for each group of the level at each % of total reads sampling point
      for (v in 1:length(level.names)){
        tmpdataframe <- subset(level.percentage.data, level.data == as.character(level.names[v])
                               )
      # vectors with: sample, level, % read, standard deviation
        summary.sample <- c(summary.sample, i)
        summary.level <- c(summary.level, j)
        summary.percent.reads <- c(summary.percent.reads, m * 100)
        summary.number.of.reads <- c(summary.number.of.reads, as.integer(m * length(sample.level.data)))
        summary.group <- c(summary.group, as.character(level.names[v]))
        summary.stdev <- c(summary.stdev, sd(tmpdataframe[,2]))
        summary.mean <- c(summary.mean, mean(tmpdataframe[,2]))
      }
    }
  }
}

# makes a data frame of the summary
summary.df <- data.frame(summary.sample, summary.level, summary.percent.reads, summary.number.of.reads, summary.group, summary.stdev, summary.mean)

# plots
for (a in 1:length(unique(summary.df$summary.sample))){
  tmpdfsample <- subset(summary.df, summary.sample == unique(summary.df$summary.sample)[a])
  

  for (b in 1:length(unique(summary.level))){
    tmpdflevel <- subset(tmpdfsample, summary.level == unique(tmpdfsample$summary.level)[b])
    # plots mean
    plot(x = tmpdflevel$summary.number.of.reads, y = tmpdflevel$summary.mean, col = as.factor(tmpdflevel$summary.group),
         main = paste(unique(summary.df$summary.sample)[a], unique(summary.level)[b]), ylim = c(0, max(tmpdfsample$summary.mean, na.rm = T)),
         xlab = "Number of Reads", ylab = paste("Mean of % of", unique(summary.level)[b], "represented"))
    # plots standard deviation
    plot(x = tmpdflevel$summary.number.of.reads, y = tmpdflevel$summary.stdev, col = as.factor(tmpdflevel$summary.group),
         main = paste(unique(summary.df$summary.sample)[a], unique(summary.level)[b]), ylim = c(0, max(tmpdfsample$summary.stdev, na.rm = T)),
         xlab = "Number of Reads" , ylab = paste("Standard Deviation of % of", unique(summary.level)[b], "represented"), cex.lab = 0.75)
  }
}

# changes the colummn names so they are more readable
colnames(summary.df) <- c("Sample", "Level", "Percent of Total Reads", "Number of Reads", "Group", "Standard Deviation", "Mean")

# writes out the table to a csv file
write.table(summary.df, file = "summary.csv", sep = ",", row.names = F)