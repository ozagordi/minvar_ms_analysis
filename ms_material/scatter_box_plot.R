library(dplyr)
library(ggplot2)

theme_agile <- function(base_size = 12, base_family = "Arial", plot.type = "formal", lines.lwd = 0.25, ticks.type = "outer", plot.grid = TRUE, axis.font = base_family, title.size = base_size*1.2, legend.size = base_size,
                        bg.col = ifelse(plot.type == "formal", "white", "#F0F0F0"), title.font = base_family , base.col = "black", axis.lines = TRUE,
                        minor.grid = ifelse(plot.grid, TRUE, FALSE), vert.grid = ifelse(plot.grid, TRUE, FALSE), ticks.length = ifelse(ticks.type == "outer", 0.2, -0.2), horz.grid = ifelse(plot.grid, TRUE, FALSE), alpha.leg = 0.1, bord.size = 0,
                        legend.bg = ifelse(plot.type == "formal", "white", "#F0F0F0"), strip.bg = ifelse(plot.type == "formal", "white", "grey80")){
  theme_bw()+
    ggplot2::theme(
      # Plot margins and finally line annotations
      plot.margin = grid::unit(c(1, 1, .5, .7), "cm"),

      text = ggplot2::element_text(family = base_family, size = base_size),
      axis.line =  element_line(size = ifelse(axis.lines, grid::unit(lines.lwd, "mm"),0), color = "black"),
      axis.ticks.length = grid::unit(ticks.length, "cm"),
      axis.ticks.margin = grid::unit(ifelse(ticks.length > 0,0.25, -ticks.length + 0.25) , "cm"),
      axis.text.x = ggplot2::element_text(size = base_size, colour = base.col, family = axis.font),
      axis.text.y = ggplot2::element_text(size = base_size, colour = base.col, family = axis.font),
      axis.title.y = ggplot2::element_text(size =  base_size, colour = base.col, vjust = 1.5, family = axis.font),
      axis.title.x = ggplot2::element_text(size = base_size,colour = base.col,vjust = -.5, family = axis.font),
      panel.background = ggplot2::element_rect(fill = bg.col),
      plot.background = ggplot2::element_rect(fill = bg.col),
      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size = bord.size),
      panel.grid.major.x = ggplot2::element_line(colour = ifelse(vert.grid, "grey60",bg.col), size = ifelse(vert.grid,0.45, 0)),
      panel.grid.minor.x = ggplot2::element_line(colour = ifelse(vert.grid, ifelse(minor.grid, "grey80",bg.col),bg.col), size = ifelse(vert.grid,0.35, 0)),
      panel.grid.major.y = ggplot2::element_line(colour = ifelse(horz.grid, "grey60",bg.col), size = ifelse(horz.grid,0.45, 0)),
      panel.grid.minor.y = ggplot2::element_line(colour = ifelse(horz.grid, ifelse(minor.grid, "grey80",bg.col),bg.col), size = ifelse(horz.grid,0.35, 0)),
      panel.grid.major = ggplot2::element_line(colour = "grey40", size=0.45),
      plot.title = ggplot2::element_text(face="bold",hjust = ifelse(plot.type == "formal", 0.5, 0) ,vjust = 2, colour = base.col, size = title.size, family = title.font),
      legend.background = ggplot2::element_rect(fill = scales::alpha(legend.bg, alpha.leg)), legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend.size),
      legend.title = element_blank(),
      strip.background = ggplot2::element_rect(fill = strip.bg),
      strip.text.x = ggplot2::element_text(size = base_size + 1),
      strip.text.y = ggplot2::element_text(size = base_size + 1)
    )
}

args <- commandArgs(trailingOnly = TRUE)
platform=as.character(args[1])

# mix is given as argument, this retrieves the relevant pair of mix numbers
mix_in = as.character(args[2])
mixes = c(9, 10, 11, 12, 13)
names(mixes) = c('GEO', '20-20', '60-10', '80-5', '90-2.5')
mix_numbers = c(mixes[mix_in], mixes[mix_in] + 8)

all_m = data.frame()
all_t = data.frame()
for(mix in mix_numbers){
    muts = paste0('../mixes_', platform, '_MinVar/minvar_results_', mix, '/annotated_mutations.csv') %>%
        read.csv() %>%
        filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
        rename(measured = freq)
    m = nrow(muts)
    muts$mix = rep(mix, m)

    truth = paste0('../true_mixes/mix_', mix, '.csv') %>%
        read.csv() %>%
        rename(expected = freq)
    m = nrow(truth)
    truth$mix = rep(mix, m)

    all_m = rbind(all_m, muts)
    all_t = rbind(all_t, truth)
}
# add the viral load  to the titers
mix_labels = paste0(mix_in, c('_10E5', '_10E4'))
df_results = full_join(all_t, all_m, by=c("mix", "gene", "pos", "mut")) %>%
  mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
  mutate(expected=ifelse(is.na(expected), 0, expected)) %>%
  mutate(mix=ifelse(mix == mix_numbers[1], mix_labels[1], mix_labels[2]))

write.csv(df_results, paste0('merged_truth_table_', platform, '_', mix_in, '.csv'), row.names=FALSE)

p = ggplot(df_results, aes(x=expected, y=measured, color=gene)) +
  geom_point(aes(shape=factor(gene)),
             position=position_jitter(width=0.02, height=0.0),
             size=1.8, solid=FALSE) +
  facet_grid(. ~ mix) +
  scale_x_continuous(limits=c(-0.03, 1.03), breaks=seq(0, 1, 0.1)) +
  scale_y_continuous(limits=c(-0.03, 1.03), breaks=seq(0, 1, 0.1)) +
  scale_shape(solid=FALSE) +
  #guides(color=FALSE) +
  guides(color=FALSE, shape=FALSE) +
  theme_bw()

gname = paste0("scatter_plot_mix_", platform, "_", mix_in, ".pdf")
ggsave(file=gname, width=200, height=120, unit='mm')

# Now for mutations identified as unique

unique_muts = read.csv('../true_mixes/unique_mutations.csv') %>%
  left_join(df_results, by=c("gene", "pos", "mut"))

write.csv(unique_muts, paste0('merged_truth_table_unique_', platform, '_', mix_in, '.csv'), row.names=FALSE)

p = ggplot(unique_muts, aes(x=clone, y=measured)) +
  geom_boxplot() + geom_jitter() +
  facet_grid(. ~ mix) +
  scale_shape(solid=FALSE) +
  guides(shape=FALSE) +
  theme_bw()

gname = paste0("boxplot_mix_unique_", platform, "_", mix_in, ".pdf")
ggsave(file=gname, width=200, height=120, unit='mm')
