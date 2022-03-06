setwd("~/BLASTResults/")

library(dplyr)
sampleNames <- c("B73", "B107", "F888", "FC1890", "Lo1056", "Lo1270", "Lo1290", "PHG83")

header <- c("subject_id","query_sequence_id" ,"subject_length", "s_start", "s_end", "evalue", "bit_score", "score", "alignment_length", "identity", "query_coverage")

dat.all <- list()

#mRNA1 transcript analysis

best.genes.Zm00001d049443_mRNA1 <- list()
for ( j in sampleNames){
f <- grep(dir(), pattern = paste0( j, ".Zm00001d049443_mRNA1.id_gene"), value = T)

for ( i in 1:(length(f)/2)){
  file <- paste( j, ".Zm00001d049443_mRNA1.id_gene", i, ".fa.blastn.txt", sep = "")
  dat <- read.table(file, sep = "\t")
  colnames(dat) <- header
  
  file <- paste( j, ".Zm00001d049443_mRNA1.id_gene", i, ".fa.blastn.reverse.txt", sep = "")
  dat.rev <- read.table(file, sep = "\t")
  colnames(dat.rev) <- paste(header, "_rev", sep = "")
  
  dat.all[[i]] <- cbind(dat, dat.rev)
  
}

dat.best <- lapply(dat.all, function(x){ x[order(x$alignment_length, decreasing = T)[1],] }) %>%
  do.call( what = rbind)

dat.select <- subset(dat.best, query_coverage > 90 & query_coverage_rev > 90)
row.names(dat.select) <- j



     best.genes.Zm00001d049443_mRNA1[[j]] <-dat.select
}

#mRNA2 transcript analysis

sampleNames <- c("B73", "B107", "F888", "FC1890", "Lo1056", "Lo1270", "Lo1290", "PHG83")


best.genes.Zm00001d049443_mRNA2 <- list()
for ( j in sampleNames){
  f <- grep(dir(), pattern = paste0( j, ".Zm00001d049443_mRNA2.id_gene"), value = T)
  
  for ( i in 1:(length(f)/2)){
    file <- paste( j, ".Zm00001d049443_mRNA2.id_gene", i, ".fa.blastn.txt", sep = "")
    dat <- read.table(file, sep = "\t")
    colnames(dat) <- header
    
    file <- paste( j, ".Zm00001d049443_mRNA2.id_gene", i, ".fa.blastn.reverse.txt", sep = "")
    dat.rev <- read.table(file, sep = "\t")
    colnames(dat.rev) <- paste(header, "_rev", sep = "")
    
    dat.all[[i]] <- cbind(dat, dat.rev)
    
  }
  
  dat.best <- lapply(dat.all, function(x){ x[order(x$alignment_length, decreasing = T)[1],] }) %>%
    do.call( what = rbind)
  
  dat.select <- subset(dat.best, query_coverage > 95 & query_coverage_rev > 90)
  row.names(dat.select) <- j
  
  
  
  best.genes.Zm00001d049443_mRNA2[[j]] <-dat.select
}
