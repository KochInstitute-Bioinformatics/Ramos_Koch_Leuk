## pca
full_pca_deseq <- function(dds, plot = FALSE, ntop = 500){
  vsd <- vst(dds, blind=FALSE)
  
  # plotPCA(vsd, intgroup="together")
  
  intgroup = "condition"; returnData = FALSE # ntop is usually 500
  
  # pca <- prcomp(vsd)
  # loadings <- as.data.frame(pca$rotation)
  
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  
  # Elbow Plot
  var <- pca$sdev^2
  variance <- var/sum(var)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  pca$percentVar <- percentVar
  # 
  # plot(1:length(variance), variance, type='l')
  # points(1:length(variance), variance)
  # 
  # head(pca$rotation[order(abs(pca$rotation[,1]), decreasing = TRUE),1])
  # head(pca$rotation[order(abs(pca$rotation[,2]), decreasing = TRUE),2])
  # 
  # PLOT
  # if (plot){
  #   PC <- data.frame(pca$x, var=interaction(meta$Treatment, meta$Gender, meta$TCPBOP))
  #   ggplot(PC,aes(x=PC1,y=PC2, col=as.factor(var)))+
  #     geom_point(size=3)+
  #     theme_classic()
  # }

  return(pca)
}

full_pca <- function(data, count_cols = '.l2tpm$', var_col = 'l2tpmVar', gene_col = 'MGI.symbol', ntop = 500){
  # vsd <- vst(dds, blind=FALSE)
  rownames(data) <- data[,gene_col]
  cols <- grep(count_cols, colnames(data))
  # intgroup = "condition"; returnData = FALSE # ntop is usually 500
  
  rv <- data[,var_col]
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(data[select, cols]))
  
  # Elbow Plot
  var <- pca$sdev^2
  variance <- var/sum(var)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  pca$percentVar <- percentVar
  
  return(pca)
}


## identify markers
identify_markers_from_compar_data <- function(compar_data, lfc_threshold = 1, adp_threshold = 0.05){
  
  markers <- c()
  col_idxs <- grep(pattern = 'logFC|adp', x = colnames(compar_data))
  compar_data <- compar_data[,col_idxs]
  
  comparisons <- gsub(pattern = '_logFC|_adp', replacement = '', x = colnames(compar_data))
  comparisons <- unique(comparisons)
  
  for (i in comparisons){
    col_idxs <- grep(pattern = i, x = colnames(compar_data))
    compar_i_data <- compar_data[,col_idxs]
    
    lfc_idx <- grep('logFC', colnames(compar_i_data))
    adp_idx <- grep('adp', colnames(compar_i_data))
    
    fwd_markers <- which((compar_i_data[,lfc_idx] >= lfc_threshold) & (compar_i_data[,adp_idx] <= adp_threshold))
    rev_markers <- which((compar_i_data[,lfc_idx] <= -1*lfc_threshold) & (compar_i_data[,adp_idx] <= adp_threshold))
    
    fwd_markers <- rownames(compar_i_data)[fwd_markers]
    rev_markers <- rownames(compar_i_data)[rev_markers]
    
    compar_groups <- strsplit(i, split = '_v_')[[1]]
    markers[compar_groups[1]] <- list(fwd_markers)
    markers[compar_groups[2]] <- list(rev_markers)
  }

  return(markers)
}

selective_pca <- function(gene_counts, genes = c(), metadata, design, metadata_cols = c(),
                         intgroup = "condition", ntop = 500, returnData = FALSE){

  ## basic setup
  
  if (length(genes) > 0){
    genes <- which(genes %in% rownames(gene_counts))
    if (length(genes) < 1){
      print('no valid genes!')
      return()
    }
    gene_counts <- gene_counts[genes,]
  }
  
  # if (ntop > nrow(gene_counts)){
  #   ntop <- nrow(gene_counts)
  # }
  
  if (length(metadata_cols) == 0){
    metadata_cols <- colnames(metadata)
  }
  
  ## Set up dds object
  dds <- DESeqDataSetFromMatrix(gene_counts, colData=metadata, design)
  dds_complex <- DESeq(dds)
  complex_res <- results(dds_complex)
  
  ## run PCA
  vsd <- varianceStabilizingTransformation(dds) #vst(dds, blind=FALSE)
  # intgroup = "condition"; ntop = 500; returnData = FALSE
  
  # pca <- prcomp(vsd)
  # loadings <- as.data.frame(pca$rotation)
  
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  
  # Elbow Plot
  var <- pca$sdev^2
  variance <- var/sum(var)
  
  groups <- c()
  for (group in metadata_cols){
    groups <- paste(groups, metadata[,group], sep = '.')
  }
  
  # PLOT
  PC <- data.frame(pca$x, var=groups)
  # ggplot(PC,aes(x=PC2,y=PC4, col=as.factor(var)))+
  #   geom_point(size=3)+
  #   theme_classic()
  # 
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  percentVar <- round(100 * percentVar)
  pca$percentVar <- percentVar
  
  return(list(PC=PC, variance=variance, data=pca))
}

## add new data points to existing pca transformation
add_samples_to_pca <- function(pca, new_samples){
  # extract factors used in pca
  factors <- rownames(pca$rotation)
  
  # subset data in new samples
  factors <- factors[which(factors %in% rownames(new_samples))]
  new_samples <- new_samples[factors,]
  rotation <- pca$rotation[factors,]
  
  # prep output
  output <- as.data.frame(t(pca$x))
  
  for (sample in colnames(new_samples)){
    loading <- c()
    for (PC in colnames(rotation)){
      new_loading <- sum(unlist(Map('*',new_samples[,sample],rotation[,PC])))
      loading <- c(loading, new_loading)
    }
    output[,sample] <- loading
  }
  output <- t(output)
  return(output)
}

arrange_matrix_for_ggplot <- function(mtx, dim = 1){
  values <- c()
  names <- c()
  
  if (dim == 1){
    for (i in 1:nrow(mtx)){
      new_values <- mtx[i,]
      new_names <- rep(rownames(mtx)[i], ncol(mtx))
      
      values <- c(values, new_values)
      names <- c(names, new_names)
    }
  }
  
  if (dim == 2){
    for (i in 1:ncol(mtx)){
      new_values <- mtx[,i]
      new_names <- rep(colnames(mtx)[i], nrow(mtx))
      
      values <- c(values, new_values)
      names <- c(names, new_names)
    }
  }
  
  df <- data.frame(values = values, names = names)
  return(df)
}

arrange_matrix_for_ggplot <- function(mtx, name_col = 1, value_cols = c()){
  values <- c()
  value_names <- c()
  names <- c()
  
  for (i in value_cols){
    new_names <- as.character(mtx[,name_col])
    new_values <- mtx[,i]
    new_value_names <- rep(colnames(mtx)[i], nrow(mtx))
    
    values <- c(values, new_values)
    value_names <- c(value_names, new_value_names)
    names <- c(names, new_names)
  }
  
  df <- data.frame(names = names, value_names = value_names, values = values)
  return(df)
}


remove_suffix <- function(string){
  idxs <- unlist(gregexpr(pattern ='_',string))#which(strsplit(levels[1], "")[[1]]=="_")
  idx <- idxs[length(idxs)]
  sub <- substr(string, 1, idx-1)
  return(sub)
}
# identify_misread_genes <- function()

calculate_zscore <- function(data){
  for (i in 1:nrow(data)){
    data[i,] <- scale(data[i,])
  }
  
  return(data)
}

##

regroup_metadata <- function(samples, column, targets = c(), queries = c(), 
                             target_name, query_name, background = '0'){
  
  target_idxs <- which(samples[,column] %in% targets)
  query_idxs <- which(samples[,column] %in% queries)
  
  samples[,column] <- as.character(samples[,column])  # background
  samples[target_idxs,column] <- as.character(target_name)
  samples[query_idxs,column] <- as.character(query_name)
  
  samples[,column] <- factor(samples[,column])
  samples[,column] <- relevel(samples[,column], target_name)
  
  return(samples)
}

targeted_dds <- function(samples, column, targets = c(), queries = c(), 
                          target_name, query_name, background = '0', print_samples = FALSE){
  
  rownames <- row.names(samples)
  samples <- regroup_metadata(samples, column, targets, queries, target_name, query_name, background)
  samples <- as.data.frame(samples[,column], row.names = rownames)
  colnames(samples) <- column
  
  if (print_samples){print(as.data.frame(samples))}
  
  design <- as.formula(paste('~', column, sep = ''))
  
  dds <- DESeqDataSetFromMatrix(data, colData=samples, design)
  dds <- DESeq(dds)
  
  return(dds)
  
}

## DE_analysis

DE_analysis <- function(dds, coeff_idx, output_directory = ''){
  
  i <- coeff_idx
  print(resultsNames(dds)[i])
  res_norm <- lfcShrink(dds, coef=i, type="normal")
  res_ape <- lfcShrink(dds, coef=i, type="apeglm")
  
  normData <- as.data.frame(res_norm)
  apeData <- as.data.frame(res_ape)
  
  colnames(normData) <- paste('norm.', colnames(normData), sep = '')
  colnames(apeData) <- paste('ape.', colnames(apeData), sep = '')
  
  # normData <- normData %>% rename(
  #   norm.baseMean = baseMean,
  #   norm.logFC = log2FoldChange,
  #   norm.lfcSE = lfcSE,
  #   norm.stat = stat,
  #   norm.p = pvalue,
  #   norm.adp = padj
  # )
  # 
  # apeData <- apeData %>% rename(
  #   ape.baseMean = baseMean,
  #   ape.logFC = log2FoldChange,
  #   ape.lfcSE = lfcSE,
  #   ape.p = pvalue,
  #   ape.adp = padj
  # )
  # 
  compar_name <- resultsNames(dds)[i]
  
  results <- merge(normData, apeData, by=0, all=TRUE)
  write.csv(as.data.frame(res_norm), file=paste(output_directory, compar_name, '_norm.csv', sep = ''))
  write.csv(as.data.frame(res_ape), file=paste(output_directory, compar_name, '_ape.csv', sep = ''))
  write.table(as.data.frame(results), file=paste(output_directory, compar_name, '_results.tsv', sep = ''), quote=FALSE, sep="\t", row.names = FALSE)
}

base_comparison <- function(x, base){
  idx <- which(x == base)
  x[c(1,idx)] <- x[c(idx,1)]
  return(x)
}
