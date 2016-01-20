
library(tidyr)
library(dplyr)
library(ggplot2)
library(assertthat)
#library(testthat)

cbPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
minfreq = rep(c(0.20, 0.10, 0.05, 0.025), 2)
names(minfreq) = c(10, 11, 12, 13, 18, 19, 20, 21)
#mix_names = c('GEO_1E5', '20-20_1E5', '60-10_1E5', '80-5_1E5', '90-2.5_1E5',
#              'GEO_1E4', '20-20_1E4', '60-10_1E4', '80-5_1E4', '90-2.5_1E4')
#names(mix_names) = c(9, 10, 11, 12, 13, 17, 18, 19, 20, 21)

f_mix_names = function(mix_n){
  mix_names = c('GEO_1E5', '20-20_1E5', '60-10_1E5', '80-5_1E5', '90-2.5_1E5',
                'GEO_1E4', '20-20_1E4', '60-10_1E4', '80-5_1E4', '90-2.5_1E4')
  names(mix_names) = c(9, 10, 11, 12, 13, 17, 18, 19, 20, 21)
  mix_names[as.character(mix_n)]
}

graph_results = function(all_res, suffix_name, meas_to_plot){
    all_res %>%
        gather(measure, value, c(precision, recall, F1)) %>%
        filter(measure == meas_to_plot) %>%
        mutate(mix_name = f_mix_names(mix)) %>%
        separate(mix_name, into=c('mix_id', 'viral_load'), sep='_') %>%
        ggplot(aes(x=mix_id, y=value)) +
        facet_grid(descr ~ viral_load) +
        geom_bar(stat="identity", position="dodge", fill="lightgrey", colour="black") +
        ylim(-0.01, 1.01) +
        scale_fill_manual(values=cbPalette) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1., size=12)) +
        xlab('mix type') + ylab(meas_to_plot) +
        guides(fill=guide_legend(title="platform / tool"))
    ggsave(paste0(meas_to_plot, "_", suffix_name, ".eps"), width=120, height=200, unit='mm')

    # all_res %>%
    #     gather(measure, value, c(precision, recall, F1)) %>%
    #     filter(measure != "F1") %>%
    #     ggplot(aes(x=mix, y=value)) +
    #     facet_grid(descr ~ measure) +
    #     geom_bar(stat="identity", position="dodge", fill="lightgrey", colour="black") +
    #     ylim(-0.01, 1.01) +
    #     scale_fill_manual(values=cbPalette) +
    #     #scale_x_discrete(labels=mix_names) +
    #     theme_bw() +
    #     theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1., size=12)) +
    #     guides(fill = guide_legend(title = "platform / tool"))
    #
    # ggsave(paste0("precision_recall_", suffix_name, ".pdf"), width=200, height=180, unit='mm')
    warnings()
}

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

# data frame for keeping all precision recall results
df = data.frame(mix=numeric(), platform=character(), tool=character(),
    precision=numeric(), recall=numeric())
# and the one for keeping only those derived from low frequency
lf_df = data.frame(mix=numeric(), platform=character(), tool=character(),
    precision=numeric(), recall=numeric())
for(mix_number in c(9, 10, 11, 12, 13, 17, 18, 19, 20, 21)){
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

    sanger = paste0('../mixes_Sanger/Mix_', mix_number, '_mutations.csv') %>%
        read.csv() %>%
        filter(gene != "GagPolTF", pos < 336, gene != 'integrase') %>%
        rename(measured = freq)
    cat("Sanger done")

    truth =  paste0('truth_expanded_', mix_number, '.csv') %>%
        read.csv() %>%
        rename(expected = freq)
    cat("Truth done")

    # join results with expected
    # MinVar Miseq
    mv_results = full_join(truth, minvar, by=c("gene", "pos", "mut")) %>%
      mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
      mutate(expected=ifelse(is.na(expected), 0, expected))
    # MinVar 454
    mv_454_results = full_join(truth, minvar_454, by=c("gene", "pos", "mut")) %>%
      mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
      mutate(expected=ifelse(is.na(expected), 0, expected))
    # VirVar
    vv_results = full_join(truth, virvar, by=c("gene", "pos", "mut")) %>%
    mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
    mutate(expected=ifelse(is.na(expected), 0, expected))

    # Sanger
    sn_results = full_join(truth, sanger, by=c("gene", "pos", "mut")) %>%
    mutate(measured=ifelse(is.na(measured), 0, measured)) %>%
    mutate(expected=ifelse(is.na(expected), 0, expected))

    # compute precision and recall on all variants
    pr = single_prec_rec(mv_results)
    #print(pr)
    row = data.frame(mix=mix_number, platform='Miseq', tool='MinVar', precision=pr[1], recall=pr[2])
    df = rbind(df, row)

    pr = single_prec_rec(mv_454_results)
    #print(pr)
    row = data.frame(mix=mix_number, platform='454', tool='MinVar', precision=pr[1], recall=pr[2])
    df = rbind(df, row)

    pr = single_prec_rec(vv_results)
    #print(pr)
    row = data.frame(mix=mix_number, platform='Miseq', tool='VirVarSeq', precision=pr[1], recall=pr[2])
    df = rbind(df, row)

    pr = single_prec_rec(sn_results)
    #print(pr)
    row = data.frame(mix=mix_number, platform='Sanger', tool='manual', precision=pr[1], recall=pr[2])
    df = rbind(df, row)


    # compute precision and recall on low freq variants only
    if(as.character(mix_number) %in% names(minfreq)){
        freq_here = minfreq[as.character(mix_number)]
        #print("frequency here is")
        #print(freq_here)

        # MinVar on Miseq
        lf_mv_results = mv_results %>% filter(expected <= freq_here)
        pr = single_prec_rec(lf_mv_results)
        row = data.frame(mix=mix_number, platform='Miseq', tool='MinVar',
                         precision=pr[1], recall=pr[2])
        lf_df = rbind(lf_df, row)

        # MinVar on 454
        lf_mv_454_results = mv_454_results %>% filter(expected <= freq_here)
        pr = single_prec_rec(lf_mv_454_results)
        row = data.frame(mix=mix_number, platform='454', tool='MinVar',
                         precision=pr[1], recall=pr[2])
        lf_df = rbind(lf_df, row)

        # VirVar on Miseq
        lf_vv_results = vv_results %>%
          filter(expected <= freq_here)
        pr = single_prec_rec(lf_vv_results)
        row = data.frame(mix=mix_number, platform='Miseq', tool='VirVarSeq',
                         precision=pr[1], recall=pr[2])
        lf_df = rbind(lf_df, row)

        lf_sn_results = sn_results %>%
          filter(expected <= freq_here)
        pr = single_prec_rec(lf_sn_results)
        row = data.frame(mix=mix_number, platform='Sanger', tool='manual',
                         precision=pr[1], recall=pr[2])
        lf_df = rbind(lf_df, row)

    }

}
# add F1 and save to file
df2 = df %>%
    mutate(descr=paste(platform, tool, sep=' / ')) %>%
    mutate(F1=2 * precision * recall / (precision + recall))
lf_df2 = lf_df %>%
    mutate(descr=paste(platform, tool, sep=' / ')) %>%
    mutate(F1=2 * precision * recall / (precision + recall))

df2 %>%
    select(mix, descr, precision, recall, F1) %>%
    write.table(file='summary_all.csv', sep="\t", row.names=FALSE)
lf_df2 %>%
    select(mix, descr, precision, recall, F1) %>%
    write.table(file='summary_low_freq.csv', sep="\t", row.names=FALSE)

df2$mix = factor(df2$mix)
lf_df2$mix = factor(lf_df2$mix)

print("Plotting all")
graph_results(df2, "all", "F1")
graph_results(df2, "all", "precision")
graph_results(df2, "all", "recall")
print("Plotting low freq")
graph_results(lf_df2, "low_freq", "F1")
graph_results(lf_df2, "low_freq", "precision")
graph_results(lf_df2, "low_freq", "recall")

warnings()
