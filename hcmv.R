library(sleuth) # make the sleuth package available to your current R session #

sample.tab <- read.table('hcmv.metadata', header = TRUE) # read the metadata table already constructed #

sleuth.obj <- sleuth_prep(sample.tab) # prepare the data table for DGE analysis (defaults to the 'full' model) #

sleuth.obj <- sleuth_fit(sleuth.obj, ~condition, 'condition.fit') # fit the data to the error model using the condition as the source of comparison #

sleuth.obj <- sleuth_fit(sleuth.obj, ~1, 'condition.re') # fit the data to a reduced fit model for a source of comparison between the two conditions #

sleuth.obj <- sleuth_lrt(sleuth.obj, 'condition.re', 'condition.fit') # perform the likelihood ratio test to determine differential gene expression #

library(dplyr) # make the dplyr package available to your current R session #
sleuth.res <- sleuth_results(sleuth.obj, 'condition.re:condition.fit', 'lrt', show_all = FALSE) # extract the significant results before FDR of the model # 

sleuth.filt <- dplyr::filter(sleuth.res, qval <= 0.05) |> dplyr::arrange(pval) # filter out q-values after FDR that are less than 0.05 and order the remaining pvals in ascending order #
sleuth.fin <- dplyr::select(sleuth.filt, target_id, test_stat, pval, qval) # take only the target_id, test_stat, pval, and qval columns of the data frame #

write.table(sleuth.fin, './dge/hcmv_sigs.tsv', sep = '\t', row.names = FALSE) # write the data frame as a tsv file #

