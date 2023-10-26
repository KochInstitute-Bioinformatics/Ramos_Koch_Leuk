library(readxl)
# library(preprocessCore)

# tximport tools

salmon_process <- function(dir = getwd(), samples_file = "samples.xlsx",  tx_to_gene_file, DGE = FALSE, output = TRUE,
                           output_data_file = "summarized.txt", output_session_file = "sessionInfo.txt",
                           salmon_data_directory_key = "quant.sf", txi_type = "salmon"){
  
  log_count_type <- ".l2tpm"
  if (DGE){
    log_count_type <- ".l2cpm"
  }
  
  # list.files(dir,"d")
  samples <- read.xlsx(samples_file)
  
  if (!(c('Sample', 'Folder', 'Condition') %in% colnames(samples) %>% all)){
    print('samples file missing necessary fields')
    return()
  }
  
  files <- file.path(dir,samples$Folder, salmon_data_directory_key)
  names(files) <- paste0(samples$Condition,".",samples$Sample)
  if (!all(file.exists(files))){
    print('missing quant files')
    return()
  }
  
  tx2gene <- read_xlsx(tx_to_gene_file)
  txi <- tximport(files, type = txi_type, tx2gene = tx2gene)
  
  intCt <- round((txi$counts),0)
  colnames(intCt)<-paste0(colnames(intCt),".intCt")
  
  logCt <- log((txi$abundance+1),2)
  colnames(logCt)<-paste0(colnames(logCt),log_count_type)
  dumpDat <- merge(intCt, logCt, by=0, all=TRUE)
  
  # abundance <- txi$abundance
  # colnames(abundance)<-paste0(colnames(abundance),".abnd")
  # dumpDat <- merge(intCt, abundance, by=0, all=TRUE)
  # 
  colnames(dumpDat)[1] <- 'geneid'
  
  if (output){
    write.table(dumpDat, sep='\t',col.names=TRUE, quote=FALSE, row.names=FALSE, file=output_data_file)
    writeLines(capture.output(sessionInfo()), output_session_file)
  }
  
  return(list(data=dumpDat,txi=txi))
}


label_and_reduce <- function(file, key=NULL, drop_extra = TRUE, label = TRUE){
  
  keep = c("gene_id")
  extra_info <- c("transcript_id.s.","length","effective_length")
  
  # drop <- c("transcript_id.s.","length","effective_length")
  # drop <- drop[which(!(drop %in% keep))]
  sample_name <- strsplit(file, split = '\\.')[[1]][1]
  sample_name <- strsplit(sample_name, split = '/')[[1]]
  sample_name <- sample_name[length(sample_name)]
  
  # print(sample_name)
  
  data <- read.table(file, header = TRUE)
  if (drop_extra){
    drop_idxs <- which(colnames(data) %in% extra_info)
    data <- data[,-drop_idxs]
  }
  
  # process counts
  # print(colnames(data))
  data$l2TPM <- log2(data$pme_TPM+1) # used to be pme_TPM
  # data$l2TPM <- normalize.quantiles(log2(round(data$pme_TPM)+1))
  # data$IntCt <- round(data$expected_count)
  data$IntCt <- round(data$posterior_mean_count)
  #
  print(sample_name)
  
  sample_key <- sample_name
  if (!is.null(key)){
    sample_key <- as.character(key[sample_name,])
  }
  
  sample_cols <- which(!(colnames(data) %in% c(keep, extra_info)))
  colnames(data)[sample_cols] <- paste(sample_key, colnames(data)[sample_cols], sep = '_')
  
  return(as.data.frame(data))
}

reduce_and_merge <- function(files, key = NULL, drop_extra = TRUE, label = TRUE){
  
  data <- label_and_reduce(files[1], key, drop_extra = drop_extra, label = label)
  
  if (length(files) == 1){
    return(data)
  }
  
  for (i in 2:length(files)){
    file <- files[i]
    new_data <- label_and_reduce(file, key, label = label)
    data <- merge(data, new_data, by="gene_id")
  }
  
  return(data)
}

# process_and_merge <- function(files){
#   keep <- c("Gene")
#   extra_info <- c("Transcript","Length","Effective_Length")
#   
#   for (file in files){
#     data <- read.table(file, header = TRUE)
#     
#     drop_idxs <- which(colnames(data) %in% extra_info)
#     data <- data[,-drop_idxs]
#     
#     counts_idxs <- grep('counts', colnames(data))
#     
#     for (i in counts_idxs){
#       
#     }
#     
#     # log2(round(data$pme_TPM)+1)
#   }
# }

countsMatrix <- function(data, col_ids = 'l2TPM', row_id = 'gene_id', 
                         normalize_magnitude=TRUE, normalize_quantiles=TRUE){
  col_idxs <- grep(paste(col_ids, collapse = '|'), colnames(data))
  counts <- data[,col_idxs]
  col_names <- colnames(counts)
  row_names <- data[,row_id]
  rownames(counts) <- row_names
  
  if (normalize_quantiles){
    counts <- normalize.quantiles(as.matrix(counts))
    counts <- as.data.frame(counts)
    rownames(counts) <- row_names
    colnames(counts) <- col_names
  }
  
  if (normalize_magnitude){
    counts <- t(apply(counts, 1, function(x){x/norm(x, type = '2')}))
  }
  
  return(counts)
}

intCountsMatrix <- function(data, col_ids = 'IntCt', row_id = 'gene_id', 
                            normalize_magnitude=TRUE, normalize_quantiles=TRUE){
  col_idxs <- grep(paste(col_ids, collapse = '|'), colnames(data))
  counts <- data[,col_idxs]
  col_names <- colnames(counts)
  row_names <- data[,row_id]
  rownames(counts) <- row_names
  
  return(counts)
}


add_annotation <- function(data, annotation_file, merge_col = 'geneid'){
  ann <- as.data.frame(read_excel(annotation_file))
  ann <- ann[-which(duplicated(ann$`Gene stable ID`)),] # remove duplicate gene ids?
  rownames(ann) <- ann$`Gene stable ID`
  # ann <- ann[,-c(1)]
  colnames(ann)[1] <- merge_col 
  # print(dim(ann))
  # gene_id <- data$gene_id #rownames(data)
  
  # ann <- ann[as.character(data$gene_id),]
  # ann$gene_id <- rownames(ann)
  output_data <- merge(data, ann, by = merge_col)
  
  # data$gene_type <- ann$`Gene type`
  # data$MGI <- ann$`MGI symbol` #`HGNC symbol`#`MGI symbol`
  
  return(output_data)
}

calculate_MPG <- function(count_data, gene_col, count_cols = 'l2tpm', mpg_col = 'MPG'){
  col_idxs <- grep(count_cols, colnames(count_data))
  
  # avg <- rowMeans(count_data[,col_idxs])
  avg <- rowMeans(apply(count_data[,col_idxs], 2, as.numeric))
  # var <- rowVars(as.matrix(count_data[,col_idxs]))
  var <- rowVars(apply(count_data[,col_idxs], 2, as.numeric))
  
  mpg <- rep('No', nrow(count_data))
  for (i in unique(count_data[,gene_col])){
    if (is.na(i)){next}
    idxs <- which(count_data[,gene_col] == i)
    max_idx <- which(avg[idxs] == max(avg[idxs]))[1]
    
    mpg[idxs[max_idx]] <- 'Yes'
  }
  
  count_data$Avg <- avg
  count_data$Var <- var
  count_data$MPG <- mpg
  
  colnames(count_data)[which(colnames(count_data) =='MPG')] <- mpg_col
  
  return(count_data)
}

add_orthology <- function(data, orthology_file, match_column = 'MGI symbol', orthology_name = 'HuSym',
                          ref_col = 2, orth_col = 1){ # which column in the orthology file is which
  orth <- as.data.frame(read_xlsx(orthology_file))
  # orth <- orth[-which(duplicated(orth[,ref_col])),] # remove duplicate gene ids?
  # rownames(orth) <- orth[,orth_idx]
  
  colnames(orth)[ref_col] <- match_column
  colnames(orth)[orth_col] <- orthology_name
  
  # orth <- orth[data[,data_idx],]
  # orth <- orth[,-c(orth_idx)]
  # 
  # data$orthology <- orth#[data[,data_idx]]
  
  data <- merge(data, orth, by = match_column)
  
  return(data)
}


# add_orthology <- function(data, data_idx, orthology_file, orth_idx){
#   orth <- as.data.frame(read_xlsx(orthology_file))
#   orth <- orth[-which(duplicated(orth[,orth_idx])),] # remove duplicate gene ids?
#   rownames(orth) <- orth[,orth_idx]
#   orth <- orth[data[,data_idx],]
#   orth <- orth[,-c(orth_idx)]
#   
#   data$orthology <- orth#[data[,data_idx]]
#   return(data)
# }

add_MPG <- function(data, counts_identifier = 'IntCt', gene_column = ''){
  genes <- as.character(data[,gene_column])
  
  MPG <- matrix(data = 'No', nrow = length(genes), ncol = 1)
  
  count_data <- data[,grep(pattern = counts_identifier, x = colnames(data))]
  
  for (gene in unique(genes)){
    idxs <- which(genes == gene)
    if (length(idxs) == 1){
      MPG[idxs] <- 'Yes'
      next
    }
    if (length(idxs) > 1){
      # print(gene)
      # print(idxs)
      count_avgs <- unlist(lapply(idxs, function(x){mean(as.numeric(count_data[x,]))}))
      # print(count_avgs)
      max_idx <- idxs[which(count_avgs == max(count_avgs))]
      # print(which(count_avgs == max(count_avgs)))
      # if (length(max_idx) == 1){}
      MPG[max_idx[1]] <- 'Yes'
    }
  }
  
  data$MPG <- MPG
  
  return(data)
}

read_QC_table <- function(files, rows=c(1:8)){
  qc <- data.frame()
  
  for (file in files){
    key_file <- as.data.frame(readxl::read_xlsx(file))
    rmv <- which(is.na(key_file$ShortSamp))
    key_file <- key_file[-rmv,rows]
    # print(colnames(key_file))
    # key <- as.data.frame(key_file)
    qc <- rbind(qc, key_file)
    key_file <- NULL
  }
  return(qc)
}

QC_key <- function(files){
  
  qc <- data.frame()
  
  for (file in files){
    key_file <- as.data.frame(readxl::read_xlsx(file))
    rmv <- which(is.na(key_file$ShortSamp))
    if (length(rmv) > 0){
      key_file <- key_file[-rmv,]
    }
    key <- data.frame(key = key_file$Key, name = key_file$ShortSamp)
    qc <- rbind(qc, key)
  }
  
  row_names <- qc$name
  qc <- data.frame(key = qc$key)
  rownames(qc) <- row_names
  
  return(qc)
}


format_deseq2_comparison_data <- function(filename, columns = c('logFC', 'p', 'adp', 'stat', 'result')){
  data <- read.table(filename, header = TRUE)
  
  index <- interaction(data[,'Comparison'], data[,'geneid'])
  rm <- which(duplicated(index))
  if (length(rm) >= 1){
    data <- data[-rm,]
  }
  
  files <- unique(data$Comparison)
  
  # print(files)
  
  gene_idx <- which(colnames(data) %in% c('geneid'))
  col_idxs <- which(colnames(data) %in% columns)
  
  output <- data[which(data$Comparison == files[1]), c(gene_idx, col_idxs)]
  # col_idxs <- grep(paste(columns, collapse = '|'), colnames(data))
  colnames(output)[2:6] <- paste(files[1], colnames(output)[2:6], sep = '_')
  # output <- 
  
  if (length(files) == 1){return(output)}
  
  for (i in files[2:length(files)]){
    new_data <- data[which(data$Comparison == i), c(gene_idx, col_idxs)]
    colnames(new_data)[2:6] <- paste(i, colnames(new_data)[2:6], sep = '_')
    # print(head(new_data))
    output <- merge(output, new_data, by = 'geneid', all=TRUE)
  }
  
  return(output)
}

# need to make column names a parameter
prepare_GSEA_rnk_files <- function(compar_data, count_data, return_data = TRUE, output = TRUE, sym = 'HuSym', mpg = 'HuMPG', split_stat = '_stat'){
  # annotation <- as.data.frame(read_excel(annotation_file))
  
  # stat <- read.table(stat_file, header = TRUE, sep = ',')
  # stat_cols <- grep('stat', colnames(stat))
  
  # colnames(count_data)[1] <- 'geneid'
  cols <- grep(pattern = paste(mpg, 'geneid', 'Gene type', sym, sep = '|'), colnames(count_data))
  count_data <- count_data[, cols]
  
  stat <- merge(count_data, compar_data, by = 'geneid', all=TRUE)
  keep_rows <- which(stat[,'Gene type'] == 'protein_coding' & stat[,mpg] == 'Yes')
  stat <- stat[keep_rows,]
  # print(head(stat))
  stat_cols <- grep('stat', colnames(stat))
  
  for (i in stat_cols){
    name <- strsplit(colnames(stat)[i], split = split_stat)[[1]][1]
    # print(name)
    
    # stat <- read.table(i, header = TRUE, sep = ',')
    # colnames(stat)[1] <- 'Gene stable ID'
    # stat <- merge(stat, annotation, by = 'Gene stable ID')
    
    rnk_data <- data.frame(sym = stat[,sym], rank = stat[,i])
    # print(head(rnk_data))
    colnames(rnk_data) <- c(paste('#', sym, sep = ''), name)
    
    rm <- union(which(is.na(rnk_data[,1])), which(is.na(rnk_data[,2])))
    rnk_data <- rnk_data[-rm,]
    
    if (output){
      filename <- paste(name, '.rnk', sep = '')
      print(filename)
      write.table(rnk_data, file = filename,  sep = '\t', quote = FALSE, row.names = FALSE)
    }
  }
  
  if (return_data){return(stat)}
}

prepare_GSEA_rnk_files_direct <- function(compar_data, output = TRUE, stat_identifer = 'stat$', 
                                          gene_sym = '#', gene_id_col = 1){
  stat_cols <- grep(pattern = stat_identifer, colnames(compar_data))
  
  for (i in stat_cols){
    name <- colnames(compar_data)[i] #strsplit(colnames(compar_data)[i], split = '_')[[1]][1]
    rnk_data <- data.frame(sym = compar_data[,gene_id_col], rank = compar_data[,i])
    colnames(rnk_data) <- c(gene_sym, name)
    
    # QC
    rm <- union(which(is.na(rnk_data[,1])), which(is.na(rnk_data[,2])))
    # print(rm)
    if (length(rm) > 0){
      rnk_data <- rnk_data[-rm,]
    }
    
    if (output){
      filename <- paste(name, '.rnk', sep = '')
      write.table(rnk_data, file = filename,  sep = '\t', quote = FALSE, row.names = FALSE)
    }
  }
  
  
}

format_IPA_data <- function(data, lfc_col = 'logFC', gene_col = 'HGNC symbol', p_col = 'adp'){
  cols <- c(gene_col, lfc_col, p_col)
  col_idxs <- grep(pattern = paste(cols, collapse = '|'), colnames(data))
  
  ipa_data <- data[,col_idxs]
  
  p_idxs <- grep(p_col, colnames(ipa_data))
  for (i in p_idxs){
    idxs <- which(ipa_data[,i] > 1 | is.na(ipa_data[,i]))
    ipa_data[idxs,i] <- 1
  }
  
  fc_idxs <- grep(lfc_col, colnames(ipa_data))
  for (i in fc_idxs){
    idxs <- which(is.na(ipa_data[,i]))
    ipa_data[idxs,i] <- 0
  }
  
  return(ipa_data)
}

# write_rnk_files <- function(stat_data, stat_col = '.stat', gene_symbol_idx = 1){
#   stat_col_idxs <- grep(stat_col, colnames(stat_data))
#   
#   
#   for (i in stat_col_idxs){
#     output <- stat_data[,c(gene_symbol_idx, i)]
#     file_name <- strsplit(colnames(stat_data)[i], split = stat_col)[[1]][1]
#     file_name <- paste(file_name, '.rnk', sep = '')
#     write.table(output, file = file_name, sep = '\t', quote = FALSE, row.names = FALSE)
#   }
# }

summarize_GSEA_output <- function(directory, comparisons, gene_collections = c('h')){
  files <- list.files(directory, recursive = TRUE)
  all_files <- files[intersect(grep(pattern = 'gsea_report_', files), grep(pattern = '\\.xls|\\.tsv', files))]
  all_files <- paste(directory, all_files, sep = '')
  df <- data.frame(NAME = c('NA'), Collection = c('NA'))
  
  for (compar in comparisons){
    # print(paste(compar, '.', sep = ''))
    compar_files <- all_files[grep(compar, all_files)]
    # print(head(compar_files))
    
    compar_df <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
    colnames(compar_df) <- c('NAME', paste(compar, c('_pos_NES', '_pos_FDR.q.val'), sep = ''),
                             paste(compar, c('_neg_NES', '_neg_FDR.q.val'), sep = ''),
                             'Collection')
    
    for (collection in gene_collections){
      files <- compar_files[grep(paste('.', collection, '.GseaPreranked', sep = ''), compar_files)]
      if (length(files) < 2){
        print(paste(collection,'missing files!'))
        next
      }
      print(files)
      
      pos <- read.table(files[1],
                        header = TRUE, sep = '\t')
      # print(head(pos,2))
      pos <- pos[,c(1,6,8)]
      colnames(pos)[2:3] <- paste(compar, colnames(pos)[2:3], sep = '_')
      neg <- read.table(files[2],
                        header = TRUE, sep = '\t')
      neg <- neg[,c(1,6,8)]
      colnames(neg)[2:3] <- paste(compar, colnames(neg)[2:3], sep = '_')
      
      data <- rbind(pos,neg)
      
      # data <- merge(pos, neg, by = 'NAME', all = TRUE)
      # print(head(data))
      data$Collection <- collection
      
      # print(head(data,2))
      compar_df <- rbind(compar_df, data)
    }
    
    df<- merge(df, compar_df, by = c('NAME', 'Collection'), all = TRUE)
    
    # print(head(df, 2))
  }
  
  # reorganize columns
  
  nes_idxs <- grep('NES', colnames(df))
  FDR_idxs <- grep('FDR', colnames(df))
  
  df <- df[,c(1,2,nes_idxs, FDR_idxs)]
  df <- df[-which(df$NAME == 'NA'),]
  return(df)
}

append_GSEA_data <- function(base_data, addn_data){
  
  base_data_names <- paste(gsea_summary$NAME, gsea_summary$Collection, sep = '-')
  addn_data_names <- paste(addn_gsea$NAME, addn_gsea$Collection, sep = '-')
  
  rownames(base_data) <- base_data_names
  rownames(addn_data) <- addn_data_names
  
  sets <- union(base_data_names, addn_data_names) #union()
  
  full_data <- matrix(nrow = length(sets), ncol = ncol(base_data))
  colnames(full_data) <- colnames(base_data)
  rownames(full_data) <- sets
  
  # full_data[,1] <- unlist(strsplit(sets, split = '-')[[1]])
  # full_data[,2] <- unlist(strsplit(sets, split = '-')[[1]])
  
  for (i in 1:nrow(base_data)){
    row <- rownames(base_data)[i]
    for (j in 1:2){
      data <- base_data[row, j]
      if (is.na(data)){next}
      full_data[row,j] <- as.character(data)
    }
    for (j in 3:ncol(base_data)){
      data <- base_data[row, j]
      if (is.na(data)){next}
      full_data[row,j] <- data
    }
  }
  
  for (i in 1:nrow(addn_data)){
    row <- rownames(addn_data)[i]
    for (j in 1:2){
      data <- addn_data[row, j]
      if (is.na(data)){next}
      full_data[row,j] <- as.character(data) #data
    }
    for (j in 3:ncol(addn_data)){
      data <- addn_data[row, j]
      if (is.na(data)){next}
      full_data[row,j] <- data
    }
  }
  
  return(full_data)
}

clean_gene_sets <- function(sets){
  for (i in 1:length(sets)){
    sets[[i]] <- unique(sets[[i]])
  }
  
  return(sets)
}

# sets is a list of lists
write_gene_sets_to_gmx <- function(sets, filename = 'gene_sets.gmx'){
  sets <- clean_gene_sets(sets)
  
  len <- max(unlist(lapply(sets, length)))
  
  mtx <- matrix(data = NA, nrow = len+1, ncol = length(sets))
  colnames(mtx) <- names(sets)
  
  for (i in 1:length(sets)){
    # print(head(sets[[i]]))
    mtx[2:(length(sets[[i]])+1),i] <- as.character(sets[[i]])
  }
  
  mtx[1,] <- rep('na', ncol(mtx))
  
  write.table(mtx, file = filename, quote = FALSE, sep = '\t', na = '', row.names = FALSE)
}

read_gmx <- function(mtx, check_formatting = TRUE, remove_notes = TRUE){
  
  # if (remove_notes){
  #   mtx <- mtx[-c(1),]
  # }
  
  markers <- c()
  
  for (i in 1:ncol(mtx)){
    name <- colnames(mtx)[i] %>% gsub(pattern = '/', replacement = '.', x = .)
    idxs <- which(!is.na(mtx[,i]) & (mtx[,i] != ''))
    elems <- mtx[idxs,i]
    if (remove_notes){
      elems <- elems[-c(1)]
    }
    elems <- unique(elems) # check for redundancy
    
    if (check_formatting){ # check for misread genes (genes read as dates by excel)
      misread <- which(!is.na(as.numeric(elems)))
      if (length(misread) > 0){
        print(paste('misread genes likely!', name, ':', elems[misread], '(', misread, ')'))
      }
    }
    
    markers[name] <- list(elems)
  }
  
  return(markers)
}

identify_misread_genes <- function(sets){
  for (i in 1:length(sets)){
    set <- sets[[i]]
    misread <- which(!is.na(as.numeric(set)))
    if (length(misread) > 0){
      print(paste('misread genes!', names(sets)[i], ':', set[misread], '(', misread, ')'))
    }

  }
}

read_gmx_csv <- function(filename, delimiter = '\t'){
  mtx <- read.csv(filename, sep = delimiter, header = TRUE)
  # mtx <- mtx[-c(1),]
  
  # markers <- c()
  # 
  # for (i in 1:ncol(mtx)){
  #   name <- colnames(mtx)[i]
  #   idxs <- which(!is.na(mtx[,i]) & (mtx[,i] != ''))
  #   elems <- mtx[idxs,i]
  #   markers[name] <- list(elems)
  # }
  # print(head(mtx))
  markers <- read_gmx(mtx)
  return(markers)
}

read_gmx_file_sheets <- function(file_path, sheets = c(), 
                                 check_gene_formatting = TRUE, remove_notes = TRUE){
  
  if (length(sheets) == 0){
    all_sheets <- excel_sheets(file_path)
    sheets <- all_sheets
  }
  
  sets <- c()
  
  for (sheet in sheets){
    print(paste('sheet:', sheet))
    mtx <- as.matrix(read_excel(file_path, sheet = sheet))
    new_sets <- read_gmx(mtx, check_formatting = TRUE, remove_notes = TRUE)
    sets <- c(sets, new_sets)
  }
  
  return(sets)
}

## concatenate matrix into dataframe for ggplot

concat_matrix <- function(mtx, valuename = 'Expression'){
  levs <- colnames(mtx)
  
  data <- c()
  labels <- c()
  
  n <- nrow(mtx)
  
  for (i in 1:ncol(mtx)){
    data <- c(data, mtx[,i])
    labels <- c(labels, rep(levs[i], ))
  }
}

alter_duplicate_names <- function(names){
  dupl <- which(duplicated(names))
  dupl <- unique(names[dupl])
  
  for (i in dupl){
    idxs <- which(names == i)
    for (j in 1:length(idxs)){
      idx <- idxs[j]
      
      names[idx] <- paste(i, j, sep = '_')
    }
  }
  
  return(names)
}

organize_named_csv <- function(file = 'named_csv', sigFC = 1, sigP = 0.05){
  data <- read.table(file, header = FALSE, sep = ',')
  
  start_idxs <- which(data$V2 == 'baseMean')
  
  names <- as.character(data$V1[start_idxs])
  names <- unlist(lapply(names, 
                      function(x){strsplit(x, split = '\\.')[[1]][1]}))
  
  # for (i in 2:6){
  #   data[,i] <- as.numeric(data[,i])
  # }
  
  genes <- unique(unlist(lapply(as.character(data$V1), function(x){strsplit(x, split = '\t')[[1]][2]})))
  genes <- genes[which(!is.na(genes))]
  
  output <- data.frame(geneid = genes)
  
  for (i in 1:length(start_idxs)){
    start <- start_idxs[i]+1
    end <- start+1
    if (i < length(start_idxs)){
      end <- start_idxs[i+1]-1
    }
    else {
      end <- dim(data)[1]
    }
    sub_mtx <- data[start:end,]
    sub_mtx$V1 <- unlist(lapply(as.character(sub_mtx$V1), 
                                function(x){strsplit(x, split = '\t')[[1]][2]}))
    colnames(sub_mtx) <- c('geneid', 
                           paste(names[i], c('baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj'),
                                 sep = '_'))
    results <- rep('n', length(genes))
    up_idxs <- which(as.numeric(sub_mtx$log2FoldChange) >= sigFC)  
    if (length(up_idxs) > 0){results[up_idxs] <- 'up'}
    dn_idxs <- which(as.numeric(sub_mtx$log2FoldChange) <= -sigFC)
    if (length(dn_idxs) > 0){results[dn_idxs] <- 'dn'}
    sub_mtx$results <- results
    
    output <- merge(output, sub_mtx, by = 'genes')
  }
  return(output)
  # unlist(lapply(head(as.character(data$V1)), function(x){strsplit(x, split = '\t')[[1]][2]}))
}


## replace charlie's named csv function with named_tsv

named_tsv_process <- function(filepath = 'named_tsv', output = 'named_tsv.txt', lfc_algorithm = 'ape', p_algorithm = 'ape', adp_algorithm = 'ape',
                              lfc_name = 'log2FoldChange', adp_name = 'padj', p_name = 'pvalue'){
  
  # filepath = 'named_tsv', output = 'named_tsv.txt', lfc_algorithm = 'ape', p_algorithm = 'ape', adp_algorithm = 'ape',
  # lfc_name = 'logFC', adp_name = 'adp', p_name = 'p'
  
  sigp = 0.05
  sigFC = 1
  
  out_lines <- c("Comparison\tgeneid\tlogFC\tp\tadp\tstat\tresult")
  # writeLines("Comparison\tgeneid\tlogFC\tp\tadp\tstat\tresult\n", out)
  
  
  con = file(filepath, "r")
  columns = readLines(con, n = 1)
  columns <- strsplit(columns, '\t|,')[[1]]
  lfc_col <- grep(paste(lfc_algorithm, lfc_name, sep = '.'), columns)
  p_col <-grep(paste(p_algorithm, p_name, sep = '.'), columns)
  adp_col <- grep(paste(adp_algorithm, adp_name, sep = '.'), columns)
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    # print(line)
    
    line <- gsub('\"', '', line)
    cols <- strsplit(line, '\t|,')[[1]]
    # print(cols)

    if (!any(grepl('baseMean', cols))){
      if (cols[4] != 'NA' & cols[10] != 'NA'){
        comparison <- gsub('.tsv', '', cols[1])
        gene_id <- cols[2]
        lfc <- cols[lfc_col] #10
        p <- gsub('NA', '1.01', cols[p_col]) #12
        adp <- gsub('NA', '1.01', cols[adp_col]) #13
        stat <- cols[6]

        result <- 'n'

        if (as.numeric(lfc) >= sigFC & as.numeric(adp) <= sigp){
          result <- 'up'
        }
        if (as.numeric(lfc) <= -1*sigFC & as.numeric(adp) <= sigp){
          result <- 'down'
        }

        write_line <- paste(comparison,gene_id,lfc,p,adp,stat,result, collapse = '\t')
        # print(write_line)
        
        out_lines <- c(out_lines, write_line)
        
        # out <- file(output, 'w')
        # write(write_line,out, append = TRUE)
        # close(out)
      }
    }

  }
  
  close(con)
  
  out <- file(output, 'w')
  writeLines(out_lines,out)
  close(out)
}

## UPDATED VERSION

named_tsv_process <- function(filepath = 'named_tsv', output = 'named_tsv.txt', lfc_algorithm = 'ape', p_algorithm = 'ape', adp_algorithm = 'ape',
                              lfc_name = 'log2FoldChange', adp_name = 'padj', p_name = 'pvalue'){
  
  # filepath = 'named_tsv', output = 'named_tsv.txt', lfc_algorithm = 'ape', p_algorithm = 'ape', adp_algorithm = 'ape',
  # lfc_name = 'logFC', adp_name = 'adp', p_name = 'p'
  
  sigp = 0.05
  sigFC = 1
  
  out_lines <- c("Comparison\tgeneid\tlogFC\tp\tadp\tstat\tresult")
  # writeLines("Comparison\tgeneid\tlogFC\tp\tadp\tstat\tresult\n", out)
  
  
  con = file(filepath, "r")
  columns = readLines(con, n = 1)
  columns <- strsplit(columns, '\t|,')[[1]]
  lfc_col <- grep(lfc_name, columns)
  p_col <-grep(p_name, columns)
  adp_col <- grep(adp_name, columns)
  
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    # print(line)
    
    line <- gsub('\"', '', line)
    cols <- strsplit(line, '\t|,')[[1]]
    # print(cols)
    
    if (!any(grepl('baseMean', cols))){
      if (cols[4] != 'NA'){
        comparison <- gsub('.tsv', '', cols[1])
        gene_id <- cols[2]
        lfc <- cols[lfc_col] #10
        p <- gsub('NA', '1.01', cols[p_col]) #12
        adp <- gsub('NA', '1.01', cols[adp_col]) #13
        stat <- cols[6]
        
        result <- 'n'
        
        if (as.numeric(lfc) >= sigFC & as.numeric(adp) <= sigp){
          result <- 'up'
        }
        if (as.numeric(lfc) <= -1*sigFC & as.numeric(adp) <= sigp){
          result <- 'down'
        }
        
        write_line <- paste(comparison,gene_id,lfc,p,adp,stat,result, collapse = '\t')
        # print(write_line)
        
        out_lines <- c(out_lines, write_line)
        
        # out <- file(output, 'w')
        # write(write_line,out, append = TRUE)
        # close(out)
      }
    }
    
  }
  
  close(con)
  
  out <- file(output, 'w')
  writeLines(out_lines,out)
  close(out)
}

aggregate_deseq_results <- function(result_files, count_data, sigp = 0.05, sigFC = 1){
  for (file in result_files){
    result_name <- file %>% gsub(pattern = '_results.tsv', replacement = '', x = .)
    res <- read.csv(file, sep = '\t')
    res <- res[,-c(2)]
    
    res$result <- 'n'
    res$result[res$log2FoldChange >= sigFC & res$padj <= sigp] <- 'up' 
    res$result[res$log2FoldChange <= -1*sigFC & res$padj <= sigp] <- 'dn' 
    
    colnames(res)[2:ncol(res)] <- paste(result_name, colnames(res)[2:ncol(res)], sep = '_')
    count_data <- merge(count_data, res, by = 'geneid')
  }
  return(count_data)
}

arrange_deseq_results <- function(data, columns = c('_log2FoldChange$', '_pvalue$','_padj$','_stat$','_result$'), 
                                  merge_col = 'geneid', base_cols){
  base <- as.data.frame(data[,merge_col])
  colnames(base) <- merge_col
  
  # col_idxs <- grep(merge_col, colnames(data))
  col_idxs <- base_cols
  
  for (i in columns){
    # print(i)
    idxs <- grep(i, colnames(data))
    # print(idxs)
    col_idxs <- c(col_idxs, idxs)
    # add_data <- data[,c(merge_col,colnames(data)[col_idxs])]
    # base <- merge(base, add_data, by = merge_col)
  }
  return(data[,col_idxs])
  
  # return(base)
}

## read broken xls

read_broken_xls <- function(file){
  lines <- readLines(file)
  # print(length(lines))
  columns <- strsplit(lines[1], split = '\t')[[1]]
  
  mtx <- matrix(nrow = (length(lines)-1), ncol = length(columns))
  colnames(mtx) <- columns
  for (i in 2:length(lines)){
    line <- strsplit(lines[i], split = '\t')[[1]]
    mtx[i-1,] <- line
  }
  data <- as.data.frame(mtx)
  return(data)
}

calculate_MPG <- function(count_data, gene_col, avg_col = '', var_col = '', count_cols = '.l2tpm', mpg_col = 'MPG'){
  
  if (avg_col == ''){
    col_idxs <- grep(count_cols, colnames(count_data))
    avg <- rowMeans(apply(count_data[,col_idxs], 2, as.numeric))
    count_data$Avg <- avg
    avg_col <- 'Avg'
  }
  if (var_col == ''){
    col_idxs <- grep(count_cols, colnames(count_data))
    var <- rowVars(apply(count_data[,col_idxs], 2, as.numeric))
    count_data$Var <- var
    var_col <- 'Var'
  }
  
  avg <- count_data[,avg_col]
  var <- count_data[,var_col]
  mpg <- rep('No', nrow(count_data))
  for (i in unique(count_data[,gene_col])){
    if (is.na(i)){next}
    idxs <- which(count_data[,gene_col] == i)
    max_idx <- which(avg[idxs] == max(avg[idxs]))[1]
    mpg[idxs[max_idx]] <- 'Yes'
    # if (i == "7420426K07Rik"){print(idxs[max_idx])}
  }
  # print(mpg[which(count_data[,gene_col] == "7420426K07Rik")])
  count_data[,mpg_col] <- mpg
  # print(count_data[which(count_data[,gene_col] == "7420426K07Rik"),mpg_col])
  # colnames(count_data)[which(colnames(count_data) =='new_mpg')] <- mpg_col
  # print(count_data[which(count_data[,gene_col] == "7420426K07Rik"),mpg_col])
  return(count_data)
}

add_orthology <- function(data, orthology_file, match_column = 'MGI.symbol', orthology_name = 'HuSym',
                          ref_col = 2, orth_col = 1){ # which column in the orthology file is which
  
  orth <- orthology_file[-which(duplicated(orthology_file[,ref_col])),]
  rownames(orth) <- orth[,ref_col]
  HuSym_orthologs <- orth[data[,match_column],orth_col]
  data[,orthology_name] <- HuSym_orthologs
  
  return(data)
}


## metadata tools

name_cols <- function(table){
  table <- data.frame(matrix(unlist(table), ncol=length(table$V1), byrow=T) %>% t,stringsAsFactors=FALSE)
  colnames(table) <- table[1,]
  table <- table[-1,]
  return(table)
}

is_outlier <- function(x) {
  out <- x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
  return(out)
}

write_metadata_record <- function(design, metadata_samples, annotation_file = '', infosite_file = '', orthology_file = '', 
                                  sample_file_name = 'deseq_samples.csv', file_name = 'run_record.txt'){
  metadata_record <- c()
  metadata_record$time <- Sys.time()
  metadata_record$design <- design
  metadata_record$annotation_file <- annotation_file
  metadata_record$orthology_file <- orthology_file
  metadata_record$QC_file <- infosite_file
  
  # date <- gsub(pattern = '-', replacement = '_', x = Sys.Date())
  # dir.create(file.path('.', date))
  # setwd(date)
  
  write.csv(metadata_samples, file = paste0(sample_file_name), row.names = TRUE)
  
  for (i in names(metadata_record)){
    string <- paste(as.character(metadata_record[[i]]), collapse = ' ')
    # print(string)
    write(string, file = file_name, append = TRUE)
  }
  
}

