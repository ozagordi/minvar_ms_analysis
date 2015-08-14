library(dplyr)
library(ggplot2)
mix1=10
mix2=18

m1 = paste0('mix_', mix1, '.csv') %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
    rename(freq1 = freq)

m2 = paste0('mix_', mix2, '.csv') %>%
    read.csv(stringsAsFactors=FALSE) %>%
    filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
    rename(freq2 = freq)

all = full_join(m1, m2, by=c('gene', 'pos', 'mut')) %>%
    mutate(freq1=ifelse(is.na(freq1), 0, freq1)) %>%
    mutate(freq2=ifelse(is.na(freq2), 0, freq2))

p = ggplot(data=all, aes(x=freq1, y=freq2, color=gene)) + geom_point()
paste0('scatterplot_', mix1, '_', mix2, '.pdf') %>%
    ggsave()
