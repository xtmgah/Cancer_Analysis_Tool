library('doParallel')
library("BSgenome.Hsapiens.UCSC.hg19")

motifs_list <- c()
for(mut in c('CA','CG','CT','TA','TC','TG')) {
  for(left1 in c('A','C','G','T')) {
    for(right1 in c('A','C','G','T')) {
      for(left2 in c('A','C','G','T')) {
        for(right2 in c('A','C','G','T')) {
          ref_c <- paste0(left2,left1,substr(mut,1,1),right1,right2)
          alt <-substr(mut,2,2)
          motif <- paste0(ref_c,alt)
          motif_eq <- paste0(reverseComplement(DNAString(ref_c)),complement(DNAString(alt)))
          motifs_list <- c(motifs_list,motif,motif_eq)
        }
      }
    }
  }
}
motif_id <- sort(rep(c(1:(6*16*16)),2))
names(motif_id) <- motifs_list

count_mutations <- function(tt)
{
  
  CHROM <- as.character(tt$chromosome)
  POS <- as.character(tt$position)
  REF <- as.character(tt$reference)
  ALT <- as.character(tt$mutated_to)
  
  MT_indxs <- which(CHROM=='MT')
  if(length(MT_indxs) != 0){CHROM[MT_indxs] <- 'M'}
  
  L1.L2.R.R1.R2.A <- paste0(BSgenome::getSeq(Hsapiens,paste0('chr',CHROM),as.numeric(POS)-2,as.numeric(POS)+2),ALT) # 5-mer context
  
  motifs_ids <- as.numeric(motif_id[L1.L2.R.R1.R2.A])
  
  tm <- table(motifs_ids)
  counts <- rep(0,(6*16*16))
  counts[as.numeric(names(tm))] <- as.numeric(tm)
  
  return(counts)
}


setProgress(0.2, detail = 'Opening input table...')

catalog <- read.table(paste0('result/',input_table_file_name,'.csv'),sep = ',',header = T)


setProgress(0.3, detail = 'Extracting tables of each sample...')
cat('\nExtracting tables of each sample...')
sample_names <- unique(catalog$sample_id)
# ---------------------
write.table(sample_names,file='result/Samples_IDs.txt',row.names=F,col.names=F)
# ---------------------
sample_indxs <- c(1:length(sample_names))  
separator <- function(i) {catalog[which(catalog$sample_id == sample_names[i]),]}

#---------------------------------------------------------------------------------------------------------------
if(number_of_cpu_cores > 1)
{
  cat('prepareing for parallelization...')
  st_prepare <- Sys.time()
  cl <- makeCluster(number_of_cpu_cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  se_prepare <- Sys.time()
  cat(paste0('Done (after ',format(se_prepare-st_prepare),')...'))
  
  cat('parallelizing...')
  sample_tables <- foreach(i=1:length(sample_indxs)) %dopar% {separator(i)}
  stopCluster(cl)
} else {
  sample_tables <- lapply(sample_indxs,separator)
}
#---------------------------------------------------------------------------------------------------------------
cat('All tables extracted.\n')



# setProgress(0.6, detail = 'Counting mutations...\nIt might take several minutes to finish this step.')

cat('\nCounting mutations...')

#---------------------------------------------------------------------------------------------------------------
if(number_of_cpu_cores > 1)
{
  cat('prepareing for parallelization...')
  st_prepare <- Sys.time()
  cl <- makeCluster(number_of_cpu_cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  se_prepare <- Sys.time()
  cat(paste0('Done (after ',format(se_prepare-st_prepare),')...'))
  
  cat('parallelizing...')
  counts <- foreach(i=1:length(sample_tables),.packages = c("BSgenome.Hsapiens.UCSC.hg19")) %dopar% {count_mutations(sample_tables[[i]])}
  stopCluster(cl)
} else {
  counts <- lapply(sample_tables,count_mutations)
}
#---------------------------------------------------------------------------------------------------------------

cat('All mutations counted.\n')


setProgress(0.9, detail = 'Saving the preprocessing results...')




# For 5-mer:
M_5mer <- simplify2array(counts)
#M_5mer <- t(t(M_5mer)/colSums(M_5mer)) # Normalize catalog
if(save_M5mer == T){write.table(data.frame(M_5mer),file=paste0('result/',M5mer_file_name,'.csv'),sep=",",col.names=FALSE,row.names=FALSE)}

# For 3-mer:
M_3mer <- t(simplify2array(lapply(c(0:95),function(i) colSums(M_5mer[c((i*16+1):(i*16+16)),]))))
#M_3mer <- t(t(M_3mer)/colSums(M_3mer)) # Normalize catalog
if(save_M3mer == T){write.table(data.frame(M_3mer),file=paste0('result/',M3mer_file_name,'.csv'),sep=",",col.names=FALSE,row.names=FALSE)}


setProgress(1, detail = 'Preprocessing finished!')




