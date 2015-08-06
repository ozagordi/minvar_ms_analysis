library(dplyr)
library(ggplot2)
library(assertthat)
#library(testthat)

args <- commandArgs(trailingOnly = TRUE)
mix_number=as.numeric(args[1])

minvarfile = Sys.glob(file.path(paste0('minvar_results_', mix_number, '*','/annotated_DRM.csv')))
assert_that(length(minvarfile) == 1)

minvar = minvarfile %>%
    read.csv() %>%
    filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
    rename(measured = freq)

truth =  paste0('../true_mixes/mix_', mix_number, '.csv') %>%
    read.csv() %>%
    rename(expected = freq)

df_results = full_join(truth, minvar, by=c("gene", "pos", "mut")) %>%
  mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
  mutate(expected=ifelse(is.na(expected), 0, expected))

p = ggplot(df_results, aes(x=expected, y=measured, color=gene)) +
    geom_point(position=position_jitter(width=0.02, height=0.00), size=2.5) +
    scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1))

gname = paste0("images/scatter_plot_mix_", mix_number, ".pdf")
ggsave(file=gname)

# single_prec_rec = function(sv){
#   # true positives, false positives, false negatives
#   tp = sum(sv$expected & sv$measured)
#   fp = sum(!sv$expected & sv$measured)
#   fn = sum(sv$expected & !sv$measured)
#   # false discovery rate (1 - precision), true positive rate (aka sensitivity or recall)
#   # false negative rate (1 - sensitivity), positive predictive value (aka precision)
#   if(fp > 0){
#     fdr = fp / (fp + tp)
#   }
#   else{
#     fdr = 0.0
#   }
#   if(tp > 0){
#     ppv = tp / (fp + tp)
#     tpr = tp / (tp + fn)
#   }
#   else{
#     ppv = 0.0
#     tpr = 0.0
#   }
#   if(fn > 0){
#     fnr = fn / (tp + fn)
#   }
#   else{
#     fnr = 0.0
#   }
#
#   return(c(ppv, tpr))
# }
#
# main_summary = function(mix_number){
# mir = parse_mix(mix_number) %>%
#   single_scatter %>%
#   print
# pr = parse_mix(mix_number) %>%
#   filter(platform=='454') %>%
#   single_prec_rec
# sprintf("454 - precision: %s, recall: %s", pr[1], pr[2]) %>%
#   print
# pr = parse_mix(mix_number) %>%
#   filter(platform=='miseq') %>%
#   single_prec_rec
#   sprintf("miseq - precision: %s, recall: %s", pr[1], pr[2]) %>%
#     print
# }
#
# sanger_results = function(mix_number){
#   mix_name = paste0("Mix_", mix_number, ".cons")
#   sanger_muts = read.csv("../Mix_1-24_CS_Aligned.csv") %>%
#     filter(mix == mix_name)
#
#   truth =  paste0('../mix_', mix_number, '.csv') %>%
#     read.csv() %>%
#     rename(freq_real = freq)
#
#   tp_muts = semi_join(sanger_muts, truth, by=c("gene", "pos", "mut"))
#   fp_muts = anti_join(sanger_muts, truth, by=c("gene", "pos", "mut"))
#   fn_muts = anti_join(truth, sanger_muts, by=c("gene", "pos", "mut"))
#   assert_that(nrow(tp_muts) + nrow(fp_muts) == nrow(sanger_muts))
#
#   # false discovery rate (1 - precision), true positive rate (aka sensitivity or recall)
#   # false negative rate (1 - sensitivity), positive predictive value (aka precision)
#   fdr = nrow(fp_muts) / (nrow(fp_muts) + nrow(tp_muts))
#   ppv = nrow(tp_muts) / (nrow(fp_muts) + nrow(tp_muts))
#   tpr = nrow(tp_muts) / (nrow(tp_muts) + nrow(fn_muts))
#   #fnr = fn / (tp + fn)
#   sprintf("Sanger - precision: %s, recall: %s", ppv, tpr) %>%
#     print
# #  return(c(ppv, tpr))
# }
#
#
#
#
#
# # simple case, totally wrong
# df = data.frame(expected=c(1, 0), measured=c(0,1))
# expect_equal(single_prec_rec(df), c(0.0, 0.0))
# # totally right
# df = data.frame(expected=c(0.1), measured=c(0.1))
# expect_equal(single_prec_rec(df), c(1.0, 1.0))
# # 1 data, false positive
# df = data.frame(expected=c(0.0), measured=c(0.1))
# expect_equal(single_prec_rec(df), c(0.0, 0.0))
# # 1 TP, 1 FP
# df = data.frame(expected=c(0, 0, 0, 1), measured=c(0, 0, 1, 1))
# expect_equal(single_prec_rec(df), c(0.5, 1.0))
# # 1 TP, 1 FN
# df = data.frame(expected=c(0, 0, 1, 1), measured=c(0, 0, 0, 1))
# expect_equal(single_prec_rec(df), c(1.0, 0.5))
