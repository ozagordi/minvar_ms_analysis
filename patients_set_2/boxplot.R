library(ggplot2)
reads_count <- read.csv("reads_count.csv", header=TRUE,
                        colClasses = c('factor', 'numeric'))
p = ggplot(reads_count, aes(subtype, percent_mapped))
p = p + geom_boxplot(outlier.colour = "green", outlier.size = 3)
p = p + geom_jitter()
ggsave(filename='boxplot.pdf')
