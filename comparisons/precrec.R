library(tidyr)
library(dplyr)
library(ggplot2)
library(assertthat)
#library(testthat)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
single_prec_rec = function(sv){
  # true positives, false positives, false negatives
  tp = sum(sv$expected & sv$measured)
  fp = sum(!sv$expected & sv$measured)
  fn = sum(sv$expected & !sv$measured)
  # false discovery rate (1 - precision), true positive rate (aka sensitivity or recall)
  # false negative rate (1 - sensitivity), positive predictive value (aka precision)
  if(fp > 0){
    fdr = fp / (fp + tp)
  }
  else{
    fdr = 0.0
  }
  if(tp > 0){
    ppv = tp / (fp + tp)
    tpr = tp / (tp + fn)
  }
  else{
    ppv = 0.0
    tpr = 0.0
  }
  if(fn > 0){
    fnr = fn / (tp + fn)
  }
  else{
    fnr = 0.0
  }

  return(c(ppv, tpr))
}

df = data.frame(mix=numeric(), platform=character(), tool=character(),
    precision=numeric(), recall=numeric())
for(mix_number in c(9, 10, 11, 12, 17, 18, 19, 20)){
    print("")
    cat(paste("Doing now", mix_number, "  ---  "))
    # parse files

    minvar = paste0('minvar_expanded_', mix_number, '.csv') %>%
        read.csv(header=TRUE) %>%
        filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
        rename(measured = freq)
    cat("MinVar done ")

    minvar_454 = paste0('minvar_454_expanded_', mix_number, '.csv') %>%
        read.csv() %>%
        filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
        rename(measured = freq)
    cat("MinVar 454 done ")

    virvar = paste0('virvar_expanded_', mix_number, '.csv') %>%
        read.csv() %>%
        filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
        rename(measured = freq)
    cat("VirVar done")

    # Sanger mix 10 did not work
    if(mix_number != 10){
        sanger = paste0('../mixes_Sanger/Mix_', mix_number, '_mutations.csv') %>%
            read.csv() %>%
            filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
            rename(measured = freq)
        cat("Sanger done")
    }

    truth =  paste0('truth_expanded_', mix_number, '.csv') %>%
        read.csv() %>%
        rename(expected = freq)
    cat("Truth done")

    # join results with expected
    mv_results = full_join(truth, minvar, by=c("gene", "pos", "mut")) %>%
      mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
      mutate(expected=ifelse(is.na(expected), 0, expected))

    mv_454_results = full_join(truth, minvar_454, by=c("gene", "pos", "mut")) %>%
      mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
      mutate(expected=ifelse(is.na(expected), 0, expected))

    vv_results = full_join(truth, virvar, by=c("gene", "pos", "mut")) %>%
    mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
    mutate(expected=ifelse(is.na(expected), 0, expected))
    if(mix_number != 10){
      sn_results = full_join(truth, sanger, by=c("gene", "pos", "mut")) %>%
        mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
        mutate(expected=ifelse(is.na(expected), 0, expected))
    }

    # compute precision and recall
    pr = single_prec_rec(mv_results)
    print(pr)
    row = data.frame(mix=mix_number, platform='Miseq', tool='MinVar', precision=pr[1], recall=pr[2])
    df = rbind(df, row)

    pr = single_prec_rec(mv_454_results)
    print(pr)
    row = data.frame(mix=mix_number, platform='454', tool='MinVar', precision=pr[1], recall=pr[2])
    df = rbind(df, row)

    pr = single_prec_rec(vv_results)
    print(pr)
    row = data.frame(mix=mix_number, platform='Miseq', tool='VirVarSeq', precision=pr[1], recall=pr[2])
    df = rbind(df, row)

    if(mix_number != 10){
        pr = single_prec_rec(sn_results)
        print(pr)
        row = data.frame(mix=mix_number, platform='Sanger', tool='manual', precision=pr[1], recall=pr[2])
        df = rbind(df, row)
    }
}

df2 = df %>%
    mutate(descr=paste(platform, tool, sep=' / ')) %>%
    mutate(F1=sqrt(precision * recall))

df2 %>%
    select(mix, descr, precision, recall, F1) %>%
    write.table(file='summary.csv', sep="\t", row.names=FALSE)
df2$mix = factor(df2$mix)

df2 %>%
    gather(measure, value, c(precision, recall, F1)) %>%
    filter(measure == "F1") %>%
    ggplot(aes(x=mix, y=value, fill=descr)) +
    geom_bar(stat="identity", position="dodge") +
    ylim(-0.01, 1.01) +
    theme_bw() +
    scale_fill_manual(values=cbPalette) +
    ylab("F1") +
    guides(fill = guide_legend(title = "platform / tool")) +
ggsave("F1_all.pdf")

df2 %>%
    gather(measure, value, c(precision, recall, F1)) %>%
    filter(measure != "F1") %>%
    ggplot(aes(x=mix, y=value, fill=descr)) +
    facet_grid(. ~ measure) +
    geom_bar(stat="identity", position="dodge") +
    ylim(-0.01, 1.01) +
    theme_bw() +
    scale_fill_manual(values=cbPalette) +
    guides(fill = guide_legend(title = "platform / tool")) +

ggsave("precision_recall_all.pdf", width=200, height=120, unit='mm')
warnings()
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
