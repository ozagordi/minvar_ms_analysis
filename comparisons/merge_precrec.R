library(xtable)
library(dplyr)
library(tidyr)
library(ggplot2)
# The palette with grey:
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
mixes = c(9, 10, 11, 12, 17, 18, 19, 20)

# here precision and recall
mv_res = data.frame()
for(mix in mixes){
    x = paste0('MinVar_analysis/summary/precrec_', mix, '.csv') %>%
        read.csv(sep='\t')
    mv_res = bind_rows(mv_res, x)
}
mv_res$tool = rep('MinVar', length(mixes))
mv_res$platform = rep('Miseq', length(mixes))

vv_res = data.frame()
for(mix in mixes){
    x = paste0('VirVarSeq_analysis/summary/precrec_', mix, '.csv') %>%
        read.csv(sep='\t')
    vv_res = bind_rows(vv_res, x)
}
vv_res$tool = rep('VirVarSeq', length(mixes))
vv_res$platform = rep('Miseq', length(mixes))

all = bind_rows(mv_res, vv_res) %>%
    gather(measure, value, c(precision, recall))
all$mix = factor(all$mix)

p = ggplot(data=all, aes(x=mix, y=value, fill=tool)) +
    facet_grid(. ~ measure) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=cbPalette)
ggsave('precision_recall.pdf')


res_454 = data.frame()
for(mix in mixes){
    x = paste0('MinVar_analysis_454/summary/precrec_', mix, '.csv') %>%
        read.csv(sep='\t')
    res_454 = bind_rows(res_454, x)
}
res_454$tool = rep('MinVar', length(mixes))
res_454$mix = factor(res_454$mix)

res_454 %>%
    gather(measure, value, c(precision, recall)) %>%
    ggplot(aes(x=mix, y=value, fill=tool)) +
    facet_grid(. ~ measure) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=cbPalette)
ggsave('precision_recall_454.pdf')

# tb = xtable(vv_res)
# digits(tb) = 3
# print(tb, include.rownames=FALSE)
