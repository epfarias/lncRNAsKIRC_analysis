gdcRNAMerge2 <- function(metadata, path, data.type, organized=FALSE) {
  
  #if (endsWith(path, '/')) {
  #  path = substr(path, 1, nchar(path)-1)
  #}
  
  if (organized==TRUE) {
    filenames <- file.path(path, metadata$file_name, 
                           fsep = .Platform$file.sep)
  } else {
    filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                           fsep = .Platform$file.sep)
  }
  
  if (data.type=='RNAseq') {
    message ('############### Merging RNAseq data ################\n',
             '### This step may take a few minutes ###\n')
    
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.table(gzfile(fl), fill = TRUE)$V2))
    rownames(rnaMatrix) <- read.table(gzfile(filenames[1]), fill =  TRUE)$V1
    rownames(rnaMatrix) <- unlist(lapply(strsplit(rownames(rnaMatrix), 
                                                  '.', fixed=TRUE), function(gene) gene[1]))
    #colnames(rnaMatrix) <- metadata$sample
    colnames(rnaMatrix) <- metadata$submitter_id
    
    #rnaMatrix <- rnaMatrix[biotype$ensemblID,]
    
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    
    return (rnaMatrix)
  } else if (data.type=='pre-miRNAs') {
    message ('############### Merging pre-miRNAs data ################\n',
             '### This step may take a few minutes ###\n')
    
    rnaMatrix <- do.call("cbind", lapply(filenames, function(fl) 
      read.delim(fl)$read_count))
    rownames(rnaMatrix) <- read.delim(filenames[1])$miRNA_ID
    
    colnames(rnaMatrix) <- metadata$sample
    
    nSamples = ncol(rnaMatrix)
    nGenes = nrow(rnaMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of genes: ', nGenes, '\n', sep=''))
    
    return (rnaMatrix)
    
    
  } else if (data.type=='miRNAs') {
    message ('############### Merging miRNAs data ###############\n')
    
    mirMatrix <- lapply(filenames, function(fl) cleanMirFun(fl))
    #mirs <- sort(unique(names(unlist(mirMatrix))))
    mirs <- rownames(mirbase)
    mirMatrix <- do.call('cbind', lapply(mirMatrix, 
                                         function(expr) expr[mirs]))
    
    rownames(mirMatrix) <- mirbase$v21[match(mirs,rownames(mirbase))]
    colnames(mirMatrix) <- metadata$sample
    
    mirMatrix[is.na(mirMatrix)] <- 0
    
    nSamples = ncol(mirMatrix)
    nGenes = nrow(mirMatrix)
    
    message (paste('Number of samples: ', nSamples, '\n', sep=''),
             paste('Number of miRNAs: ', nGenes, '\n', sep=''))
    
    return (mirMatrix)
  } else {
    return ('error !!!')
  }
}
