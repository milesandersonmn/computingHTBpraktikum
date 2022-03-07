setwd("~/BLASTResults/")

library(dplyr)
sampleNames <- c("B73", "B106", "B107", "F888", "FC1890", "Lo1056", "Lo1270", "Lo1290", "PHG83")

header <- c("subject_id","query_sequence_id" ,"subject_length", "s_start", "s_end", 
            "evalue", "bit_score", "score", "alignment_length", "identity", 
            "query_coverage")

dat.fwd.all <- list()
dat.rev.all <- list()

######################## mRNA1 transcript analysis ###########################

best.genes.Zm00001d049443_mRNA1 <- list()
for ( j in sampleNames){

gene_model <- "Zm00001d049443_mRNA1"  
 
  
  f <- grep(dir(), pattern = paste0( j, ".Zm00001d049443_mRNA1.id_gene"), value = T)

for ( i in 1:(length(f)/2)){
  file <- paste( j, ".Zm00001d049443_mRNA1.id_gene", i, ".fa.blastn.txt", sep = "")
  dat.fwd <- read.table(file, sep = "\t")
  colnames(dat.fwd) <- header
  #rownames(dat) <- paste0( j, ".Zm00001d049443_mRNA1")
  
  file <- paste( j, ".Zm00001d049443_mRNA1.id_gene", i, ".fa.blastn.reverse.txt", sep = "")
  dat.rev <- read.table(file, sep = "\t")
  colnames(dat.rev) <- paste(header, "_rev", sep = "")
  
  dat.fwd.all[[i]] <- dat.fwd
  dat.rev.all[[i]] <- dat.rev
}

for ( s in 1:length(dat.fwd.all)){
  if ( nrow(dat.fwd.all[[s]]) > 1 ){
    longest_alignment <- order(dat.fwd.all[[s]]$alignment_length, decreasing = T)[1]
    dat.fwd.all[[s]] <- dat.fwd.all[[s]][longest_alignment, ]
  } 
}  
for ( s in 1:length(dat.rev.all)){
  if ( nrow(dat.rev.all[[s]]) > 1 ){
    longest_alignment <- order(dat.rev.all[[s]]$alignment_length_rev, decreasing = T)[1]
    dat.rev.all[[s]] <- dat.rev.all[[s]][longest_alignment, ]
  } 
}
  
dat.fwd.best <- lapply(dat.fwd.all, function(x){ x[order(x$subject_id, decreasing = T)[1],] }) %>%
  do.call( what = rbind)

dat.fwd.select <- subset(dat.fwd.best, query_coverage > 85)

dat.rev.best <- lapply(dat.rev.all, function(x){ x[order(x$alignment_length, decreasing = T)[1],] }) %>%
  do.call( what = rbind)

dat.rev.select <- subset(dat.rev.best, query_coverage_rev > 85)

dat.rev.select <- dat.rev.select[match(dat.fwd.select$subject_id, dat.rev.select$query_sequence_id_rev),]

dat.select <- cbind(dat.fwd.select, dat.rev.select)
dat.select['sample']=j
dat.select['gene_model']=gene_model

     best.genes.Zm00001d049443_mRNA1[[j]] <-dat.select
}

best.genes <- bind_rows(best.genes.Zm00001d049443_mRNA1)


############### mRNA2 transcript analysis ####################

best.genes.Zm00001d049443_mRNA2 <- list()
for ( j in sampleNames){
  
  gene_model <- "Zm00001d049443_mRNA2"  
  
  
  f <- grep(dir(), pattern = paste0( j, ".Zm00001d049443_mRNA2.id_gene"), value = T)
  
  for ( i in 1:(length(f)/2)){
    file <- paste( j, ".Zm00001d049443_mRNA2.id_gene", i, ".fa.blastn.txt", sep = "")
    dat.fwd <- read.table(file, sep = "\t")
    colnames(dat.fwd) <- header
    #rownames(dat) <- paste0( j, ".Zm00001d049443_mRNA2")
    
    file <- paste( j, ".Zm00001d049443_mRNA2.id_gene", i, ".fa.blastn.reverse.txt", sep = "")
    dat.rev <- read.table(file, sep = "\t")
    colnames(dat.rev) <- paste(header, "_rev", sep = "")
    
    dat.fwd.all[[i]] <- dat.fwd
    dat.rev.all[[i]] <- dat.rev
  }
  
  for ( s in 1:length(dat.fwd.all)){
    if ( nrow(dat.fwd.all[[s]]) > 1 ){
      longest_alignment <- order(dat.fwd.all[[s]]$alignment_length, decreasing = T)[1]
      dat.fwd.all[[s]] <- dat.fwd.all[[s]][longest_alignment, ]
    } 
  }  
  for ( s in 1:length(dat.rev.all)){
    if ( nrow(dat.rev.all[[s]]) > 1 ){
      longest_alignment <- order(dat.rev.all[[s]]$alignment_length_rev, decreasing = T)[1]
      dat.rev.all[[s]] <- dat.rev.all[[s]][longest_alignment, ]
    } 
  }
  
  dat.fwd.best <- lapply(dat.fwd.all, function(x){ x[order(x$subject_id, decreasing = T)[1],] }) %>%
    do.call( what = rbind)
  
  dat.fwd.select <- subset(dat.fwd.best, query_coverage > 85)
  
  dat.rev.best <- lapply(dat.rev.all, function(x){ x[order(x$alignment_length, decreasing = T)[1],] }) %>%
    do.call( what = rbind)
  
  dat.rev.select <- subset(dat.rev.best, query_coverage_rev > 85)
  
  dat.rev.select <- dat.rev.select[match(dat.fwd.select$subject_id, dat.rev.select$query_sequence_id_rev),]
  
  dat.select <- cbind(dat.fwd.select, dat.rev.select)
  dat.select['sample']=j
  dat.select['gene_model']=gene_model
  
  best.genes.Zm00001d049443_mRNA2[[j]] <-dat.select
  
}

x <- bind_rows(best.genes.Zm00001d049443_mRNA2)
best.genes <- rbind(best.genes, x)
