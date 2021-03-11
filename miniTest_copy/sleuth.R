library(sleuth)
library(dplyr)


#read sample_info.txt created after kallisto
stab <- read.table("sample_info.txt",header=TRUE,stringsAsFactors=FALSE)

#initialize sleuth object
so <- sleuth_prep(stab)


#fit a model comparing the two conditions 
so <- sleuth_fit(so, ~condition, 'full')

#fit the reduced model to compare in the likelihood ratio test 
so <- sleuth_fit(so, ~1, 'reduced')

#perform the likelihood ratio test for differential expression between conditions 
so <- sleuth_lrt(so, 'reduced', 'full')

#extract the test results from the sleuth object 
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant <- dplyr::filter(sleuth_table, qval < 0.05) %>% dplyr::arrange(pval) 

#select to output to log file
sleuth_log <- (dplyr::select(sleuth_significant, target_id, test_stat, pval, qval))

#write to log file, if Rscript is run via python by creating temp sleuth_log file
write.table(sleuth_log, file="sleuth_log.txt",quote = FALSE,row.names = FALSE,sep="\t")

#write to log file, this line will only work if Rscript is run separately from python
#write.table(sleuth_log, file="logfile.txt",quote = FALSE,row.names = FALSE,sep="\t",append=TRUE)







#print top 10 transcripts
#head(sleuth_significant, n=10)

#write FDR < 0.05 transcripts to file
#write.table(sleuth_significant, file="fdr_results.txt",quote = FALSE,row.names = FALSE)

#just show transcript, pval, qval (select by column header names) 
#head(dplyr::select(sleuth_significant, target_id, test_stat, pval, qval), n=10)

