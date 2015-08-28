library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if(is.null(args[1])) {
    gene_in = 'RT'
    pos_in = 67
    mut_in = '-'
} else {
    gene_in = args[1]
    pos_in = as.integer(args[2])
    mut_in = args[3]
}

mixes = c(9, 10, 11, 12, 13, 17, 18, 19, 20, 21)
# exp = rep(c(3.2, 20, 10, 5, 2.5), 2)
# names(exp) = as.character(mixes)
paste('mix', 'Miseq', '454', 'virus titre', sep='\t') %>%
    write(file='')
for (mix in mixes){
mut_miseq = paste0('../mixes_Miseq_MinVar/minvar_results_', mix, '/annotated_mutations.csv') %>%
    read.csv() %>%
    filter(gene==gene_in, pos==pos_in, mut==mut_in)
mut_454 = paste0('../mixes_454_MinVar/minvar_results_', mix, '/annotated_mutations.csv') %>%
    read.csv() %>%
    filter(gene==gene_in, pos==pos_in, mut==mut_in)

freq_miseq = round(100 * mut_miseq$freq, 1)
freq_454 = round(100 * mut_454$freq, 1)
vt = ifelse(mix < 15, 100000L, 10000L)
paste(mix, freq_miseq, freq_454, vt, sep='\t') %>%
    write(file='')
}
