library(gplots)

enrichment_plot <- function(rankings, components, bin=0.01, 
                            title = '',xtitle = '', ytitle = ''){
  rank_order <- abs(rankings[order(abs(rankings), decreasing = TRUE)])
  thresh <- seq(from = 0, to = max(rank_order), by = bin)
  
  y <- c()
  
  for (i in thresh){
    set <- names(rank_order[which((rank_order >= i) & (rank_order < (i+bin)))])
    enrich <- length(intersect(components,set))
    total <- length(set)
    
    frac <- enrich/total
    # print(frac)
    if (total == 0){frac <- 1.1}
    y <- c(y, frac)
  }
  
  plot(thresh, y, type='l', main=title, xlab=xtitle,ylab=ytitle)
}

## Volcano Plot

volcano_plot <- function(data, comparison, adp_identifier = 'adp', lfc_identifier = 'logFC', gene_column = "HGNC symbol", 
                         genes = c(), title = '', label_y_axis_space = 5,
                         lfc_threshold = 1.0, pval_threshold = 0.05,
                         num_upregulated = 0, num_downregulated = 0, num_low_padj = 0, output = FALSE,
                         label_limit = 25,
                         output_file = ''){
  if (title == ''){title <- comparison}
  exp <- as.data.frame(data)
  exp <- exp[,grep(comparison, colnames(exp))]
  
  adp_idx <- grep(adp_identifier, colnames(exp))
  lfc_idx <- grep(lfc_identifier, colnames(exp))
  
  exp$log_padj <- -log10(exp[,adp_idx])
  exp$padj <- exp[,adp_idx]
  exp$log2FoldChange <- exp[,lfc_idx]
  rownames(exp) <- data[, gene_column]
  
  idx <- which(rownames(exp) %in% genes)
  exp$set <- 'other'
  exp$set[idx] <- genes
  exp$opacity <- 0.25
  exp$opacity[idx] <- 5
  exp$size <- 0.1
  exp$size[idx] <- 0.2
  exp$gene_names <- as.character(rownames(exp))
  
  # significant fold change with low p-val
  hi_vals <- which((exp$log2FoldChange >= lfc_threshold) & (exp$padj <= pval_threshold))
  lo_vals <- which((exp$log2FoldChange <= -1*lfc_threshold) & (exp$padj <= pval_threshold))
  
  # top fold change
  top_hi_vals <- head(sort(exp$log2FoldChange, decreasing = TRUE), num_upregulated)
  top_idxs <- which(exp$log2FoldChange %in% top_hi_vals)
  # bottom fold change
  top_lo_vals <- head(sort(exp$log2FoldChange, decreasing = FALSE), num_downregulated)
  bot_idxs <- which(exp$log2FoldChange %in% top_lo_vals)
  
  lo_padj <- head(sort(exp$padj, decreasing = FALSE), num_low_padj) # low p-adj values
  lo_idxs <- which(exp$padj %in% lo_padj)
  # print(union(which(lo_idxs %in% top_idxs), which(lo_idxs %in% bot_idxs)))
  # rm_dup <- union(which(lo_idxs %in% top_idxs), which(lo_idxs %in% bot_idxs))
  # if (length(rm_dup) > 0){
  #   lo_idxs <- lo_idxs[-union(which(lo_idxs %in% top_idxs), which(lo_idxs %in% bot_idxs))]
  # }
  
  if (output){
    print('upregulated genes:')
    print(exp[top_idxs, c('padj', 'log2FoldChange')])
    print('downregulated genes:')
    print(exp[bot_idxs, c('padj', 'log2FoldChange')])
    print('low p-values:')
    print(exp[lo_idxs, c('padj', 'log2FoldChange')])
  }
  
  hi_label <- hi_vals
  if (length(hi_label) > label_limit){
    sort_vals <- exp$log2FoldChange[hi_label]*exp$log_padj[hi_label]
    hi_label_vals <- head(sort(sort_vals, decreasing = TRUE), label_limit)
    hi_label <- hi_label[which(sort_vals %in% hi_label_vals)]
  }
  
  lo_label <- lo_vals
  if (length(lo_label) > label_limit){
    sort_vals <- exp$log2FoldChange[lo_label]*exp$log_padj[lo_label]*-1
    lo_label_vals <- head(sort(sort_vals, decreasing = TRUE), label_limit)
    lo_label <- lo_label[which(sort_vals %in% lo_label_vals)]
  }
  
  label_idxs <- unique(c(lo_idxs, top_idxs, bot_idxs, idx, hi_label, lo_label))
  
  gp <- ggplot() + 
    geom_point(data=exp, 
               mapping = aes(y=log_padj, x=log2FoldChange), color='gray') +
    ggtitle(title) +
    ylab('-log10(adp)') +
    scale_x_continuous(breaks = round(seq(-6, 6, by = 1),1)) +
    # geom_text(aes(label = genes, y = 0.98, x = -0.25)) +
    geom_point(data=exp[hi_vals,], 
               mapping = aes(y=log_padj, x=log2FoldChange), color='firebrick1') +
    geom_point(data=exp[lo_vals,], 
               mapping = aes(y=log_padj, x=log2FoldChange), color='dodgerblue2') +
    geom_point(data=exp[idx,], 
               mapping = aes(y=log_padj, x=log2FoldChange), color='purple') + 
    geom_label_repel(aes(label = exp$gene_names[label_idxs], y = (exp$log_padj[label_idxs]),
                         x = exp$log2FoldChange[label_idxs])) #+
  # geom_label_repel(aes(label = exp$gene_names[top_idxs], y = (exp$log_padj[top_idxs]), 
  # x = exp$log2FoldChange[top_idxs])) +
  # geom_label_repel(aes(label = exp$gene_names[bot_idxs], y = (exp$log_padj[bot_idxs]), 
  # x = exp$log2FoldChange[bot_idxs])) +
  # geom_label_repel(aes(label = exp$gene_names[idx], y = (exp$log_padj[idx]), 
  # x = exp$log2FoldChange[idx]))
  if (output_file != ''){
    ggsave(filename = output_file, width = 10, height = 10)
  }
  
  return(gp)
}

## lfc compar

lfc_compar_plot <- function(data, comparisons, adp_identifier = 'adp', lfc_identifier = 'logFC', gene_column = "HGNC symbol", 
                         genes = c(), title = '', label_y_axis_space = 5,
                         lfc_threshold = 1.0, pval_threshold = 0.05,
                         num_upregulated = 10, num_downregulated = 10, num_low_padj = 0, output = FALSE,
                         label_limit = 25){
  
  if (length(comparisons) > 2){
    print('Only 2 Axes')
    return()
  }
  
  exp <- as.data.frame(data)
  col_1 <- grep(comparisons[1], colnames(exp))
  col_2 <- grep(comparisons[2], colnames(exp))
  exp <- exp[,c(col_1, col_2)]
  
  adp_idxs <- grep(adp_identifier, colnames(exp))
  lfc_idxs <- grep(lfc_identifier, colnames(exp))
  
  colnames(exp)[lfc_idxs[1]] <- 'x_lfc'
  colnames(exp)[lfc_idxs[2]] <- 'y_lfc'
  
  colnames(exp)[adp_idxs[1]] <- 'x_adp'
  colnames(exp)[adp_idxs[2]] <- 'y_adp'
  
  exp[,'x_lfc'] <- as.numeric(exp[,'x_lfc'])
  exp[,'y_lfc'] <- as.numeric(exp[,'y_lfc'])
  exp[,'x_adp'] <- as.numeric(exp[,'x_adp'])
  exp[,'y_adp'] <- as.numeric(exp[,'y_adp'])
  
  # print(head(exp))
  
  significant_both  <- which(exp$x_adp <= pval_threshold & exp$y_adp <= pval_threshold & abs(exp$x_lfc) >= lfc_threshold & abs(exp$y_lfc) >= lfc_threshold)
  significant_one <- which(exp$x_adp <= pval_threshold & !(exp$y_adp <= pval_threshold & abs(exp$y_lfc) >= lfc_threshold) & abs(exp$x_lfc) >= lfc_threshold)
  significant_two <- which(exp$y_adp <= pval_threshold & !(exp$x_adp <= pval_threshold & abs(exp$x_lfc) >= lfc_threshold) & abs(exp$y_lfc) >= lfc_threshold)
  
  rownames(exp) <- data[, gene_column]
  
  idx <- which(rownames(exp) %in% genes)
  exp$set <- 'other'
  exp$set[idx] <- genes
  exp$opacity <- 0.25
  exp$opacity[idx] <- 5
  exp$size <- 0.1
  exp$size[idx] <- 0.2
  exp$gene_names <- as.character(rownames(exp))
  
  exp$significant <- 'not significant'
  exp$significant[significant_both] <- 'significant in both'
  exp$significant[significant_one] <- paste('significant in', comparisons[1])
  exp$significant[significant_two] <- paste('significant in', comparisons[2])

  # label idxs 
  sig_one_up <- head(order(exp$x_lfc[significant_one], decreasing = TRUE), num_upregulated)
  sig_one_up_idxs <- significant_one[sig_one_up]
  sig_two_up <- head(order(exp$y_lfc[significant_two], decreasing = TRUE), num_upregulated)
  sig_two_up_idxs <- significant_two[sig_two_up]
  
  sig_one_down <- head(order(exp$x_lfc[significant_one], decreasing = FALSE), num_upregulated)
  sig_one_down_idxs <- significant_one[sig_one_down]
  sig_two_down <- head(order(exp$y_lfc[significant_two], decreasing = FALSE), num_upregulated)
  sig_two_down_idxs <- significant_two[sig_two_down]
  
  xy_sig <- exp$y_lfc*exp$x_lfc
  minus_idxs <- which(exp$x_lfc < 0 & exp$y_lfc < 0)
  plus_idxs <- which(exp$x_lfc > 0 & exp$y_lfc > 0)
  xy_sig[minus_idxs] <- -1*xy_sig[minus_idxs]
  xy_sig[plus_idxs] <- abs(xy_sig[plus_idxs])
  
  sig_both_up <- head(order(xy_sig[significant_both], decreasing = TRUE), num_upregulated)
  sig_both_up_idxs <- significant_both[sig_both_up]
  sig_both_down <- head(order(xy_sig[significant_both], decreasing = FALSE), num_upregulated)
  sig_both_down_idxs <- significant_both[sig_both_down]
  
  label_idxs <- unique(c(sig_one_up_idxs, sig_one_down_idxs, sig_two_up_idxs, sig_two_down_idxs, sig_both_up_idxs, sig_both_down_idxs))
  
  # exp$sig <- 'none'
  # exp$sig[c(sig_one_up_idxs, sig_one_down_idxs)] <- ''
  
  ggplot() + 
    geom_point(data=exp, 
               mapping = aes(y=y_lfc, x=x_lfc), color='gray') +
    ggtitle(title) +
    xlab(paste(comparisons[1], 'log2FoldChange')) +
    ylab(paste(comparisons[2], 'log2FoldChange')) +
    scale_x_continuous(breaks = round(seq(-6, 6, by = 1),1)) +
    # geom_point(data=exp[idx,], 
    #            mapping = aes(y=y_lfc, x=x_lfc), color='purple') + 
    # # geom_text(aes(label = genes, y = 0.98, x = -0.25)) +
    geom_point(data=exp[significant_both,],
               mapping = aes(y=y_lfc, x=x_lfc), color='orchid3') +
    geom_point(data=exp[significant_one,],
               mapping = aes(y=y_lfc, x=x_lfc), color='firebrick1') +
    geom_point(data=exp[significant_two,],
               mapping = aes(y=y_lfc, x=x_lfc), color='dodgerblue2') +
    geom_label_repel(aes(label = exp$gene_names[label_idxs], y = (exp$y_lfc[label_idxs]), 
                         x = exp$x_lfc[label_idxs]))
  # geom_label_repel(aes(label = exp$gene_names[top_idxs], y = (exp$log_padj[top_idxs]), 
  # x = exp$log2FoldChange[top_idxs])) +
  # geom_label_repel(aes(label = exp$gene_names[bot_idxs], y = (exp$log_padj[bot_idxs]), 
  # x = exp$log2FoldChange[bot_idxs])) +
  # geom_label_repel(aes(label = exp$gene_names[idx], y = (exp$log_padj[idx]), 
  # x = exp$log2FoldChange[idx]))
}


##

plot_heatmap <- function(mtx, color_ramp = c('red','white', 'blue'), margins = c(5,5), rownames = c()){
  col_func <- colorRampPalette(color_ramp)
  
  if (length(rownames) == 0){
    rownames <- rownames(mtx)
  }
  
  heatmap.2(mtx, col = col_func(100), margins = margins, trace = 'none', cexRow = 0.6, cexCol = 0.6, labRow = rownames)
  
}

## GSEA enrichment plot

GSEA_enrichment_plot <- function(directory, comparison, gene_sets, rank_thresh = NULL,
                                 title = ''){
  
  all_files <- list.files(directory, recursive = TRUE)
  all_files <- grep(comparison, all_files, value = TRUE)
  
  enrichment_files <- paste(gene_sets, '.xls', sep = '')
  
  files <- c()
  for (i in enrichment_files){
    new_files <- grep(i, all_files, value = TRUE)
    files <- c(files, new_files)
  }
  
  files <- paste(directory, files, sep = '')
  
  df <- data.frame(ES = c(NA), rank = c(NA), gene_set = c(NA), color = c(NA))
  
  for (i in 1:length(files)){
    file <- files[i]
    enrich_data <- read_broken_xls(file)
    new_data <- data.frame(ES = enrich_data$`RUNNING ES`, rank = enrich_data$`RANK METRIC SCORE`, gene_set = gene_sets[i], color = cols[i])
    df <- rbind(df, new_data)
  }
  
  df <- df[-1,]
  
  df$ES <- as.numeric(df$ES)
  df$rank <- as.numeric(df$rank)
  
  if (!is.null(rank_thresh)){
    df <- df[df$rank >= rank_thresh,]
  }
  
  ggplot(df, aes(x = rank, y = ES, color = gene_set)) + 
    geom_point() + ggtitle(title)
    
}


GSEA_double_enrichment_plot <- function(directory, comparisons, gene_sets, rank_thresh = NULL,
                                 title = ''){
  
  if (length(comparisons) > 2 | length(gene_sets) > 2){
    print('wrong number of inputs')
    return()
  }
  
  all_files <- list.files(directory, recursive = TRUE)
  all_files <- grep(comparison, all_files, value = TRUE)
  
  enrichment_files <- paste(gene_sets, '.xls', sep = '')
  
  files <- c()
  for (i in enrichment_files){
    new_files <- grep(i, all_files, value = TRUE)
    files <- c(files, new_files)
  }
  
  files <- paste(directory, files, sep = '')
  
  df <- data.frame(ES = c(NA), rank = c(NA), gene_set = c(NA), color = c(NA))
  
  for (i in 1:length(files)){
    file <- files[i]
    enrich_data <- read_broken_xls(file)
    new_data <- data.frame(ES = enrich_data$`RUNNING ES`, rank = enrich_data$`RANK METRIC SCORE`, gene_set = gene_sets[i], color = cols[i])
    df <- rbind(df, new_data)
  }
  
  df <- df[-1,]
  
  df$ES <- as.numeric(df$ES)
  df$rank <- as.numeric(df$rank)
  
  if (!is.null(rank_thresh)){
    df <- df[df$rank >= rank_thresh,]
  }
  
  ggplot(df, aes(x = rank, y = ES, color = gene_set)) + 
    geom_point() + ggtitle(title)
  
}

replotGSEA_multple_enrichment <- function(path, gene.sets, metric.range,
                                          enrichment.score.range, title = "Enrichment") {
  
  colors <- c('red', 'dodgerblue1', 'forestgreen', 'gold1', 'darkorchid1', 'cyan2', 'chartreuse2', 'deeppink2')
  
  # if (missing(path)) {
  #   stop("Path argument is required")
  # }
  # if (!file.exists(path)) {
  #   stop("The path folder could not be found. Please change the path")
  # }
  # if (missing(gene.set)) {
  #   stop("Gene set argument is required")
  # }
  # 
  
  hit_indices <- c()
  es_profile <- c()
  rnk <- c()
  
  for (i in gene.sets){
    path.rnk <- list.files(path = file.path(path, "edb"),
                           pattern = ".rnk$", full.names = TRUE)
    gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
    # print(path.rnk)
    colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
    # if (missing(metric.range)) {
    #   metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
    # }  
    # metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
    # rnk[i] <- list(gsea.rnk$metric) #list(gsea.metric)
    ## Load .edb data
    path.edb <- list.files(path = file.path(path, "edb"),
                           pattern = ".edb$", full.names = TRUE)
    gsea.edb <- read.delim(file = path.edb,
                           header = FALSE, stringsAsFactors = FALSE)
    gsea.edb <- unlist(gsea.edb)
    gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
    gsea.metric <- unlist(strsplit(gsea.metric, " "))
    gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
    gsea.metric <- gsub("METRIC=", "", gsea.metric)
    gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
    # Select the right gene set
    if (length(gsea.edb) == 0) {
      stop(paste("The gene set name was not found, please provide",
                 "a correct name"))
    }
    if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", i), " "), gsea.edb)) > 1) {
      warning(paste("More than 1 gene set matched the gene.set",
                    "argument; the first match is plotted"))
    }
    gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", i), " "), gsea.edb)[1]]
    if (is.na(gsea.edb[[1]])){next}
    # Get template name
    gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
    gsea.edb <- unlist(strsplit(gsea.edb, " "))
    gsea.template <- gsea.edb[1]
    
    # Get gene set name
    gsea.gene.set <- gsea.edb[2]
    gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
    
    # Get enrichment score
    gsea.enrichment.score <- gsea.edb[3]
    gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
    
    # Get gene set name
    gsea.normalized.enrichment.score <- gsea.edb[4]
    gsea.normalized.enrichment.score <- gsub("NES=", "",
                                             gsea.normalized.enrichment.score)
    
    # Get nominal p-value
    gsea.p.value <- gsea.edb[5]
    gsea.p.value <- gsub("NP=", "", gsea.p.value)
    gsea.p.value <- as.numeric(gsea.p.value)
    
    # Get FDR
    gsea.fdr <- gsea.edb[6]
    gsea.fdr <- gsub("FDR=", "", gsea.fdr)
    gsea.fdr <- as.numeric(gsea.fdr)
    
    # Get hit indices
    gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
    gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
    gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
    gsea.hit.indices <- as.integer(gsea.hit.indices)
    
    # Get ES profile
    gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
    gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
    gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
    gsea.es.profile <- as.numeric(gsea.es.profile)
    
    hit_indices[i] <- list(gsea.hit.indices)
    es_profile[i] <- list(gsea.es.profile)
    rnk[i] <- list(gsea.rnk$metric)
  }
  # enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
  es_max <- max(unlist(es_profile))
  es_min <- min(unlist(es_profile))
  enrichment.score.range <- c(es_min, es_max)
  
  # metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
  rnk_min <- min(unlist(rnk))
  rnk_max <- max(unlist(rnk))
  metric.range <- c(rnk_min, rnk_max)
  
  ## Create GSEA plot
  # Save default for resetting
  # def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  # dev.new(width = 10, height = 20)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(seq(1, (2*length(hit_indices)+2))), heights = c(1.7, rep(c(0.2, 0.5), length(hit_indices)), 0.5))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, hit_indices[[1]], length(gsea.rnk$metric)),
       c(0, es_profile[[1]], 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", ylab = "Enrichment score (ES)",
       ylim = enrichment.score.range,
       main = list(title, font = 1, cex = 2),
       xlab = "Rank in ordered gene list", xlim = c(0, length(rnk[[1]])),
       # panel.first = abline(h = seq(metric.range[1] / 2,
       #                              metric.range[2] - metric.range[1] / 4,
       #                              metric.range[2] / 2), col = "gray95", lty = 2),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  for (i in 2:length(hit_indices)){
    print(i)
    lines(c(1, hit_indices[[i]], length(rnk[[i]])),
          c(0, es_profile[[i]], 0), col = colors[i])
  }
  plot.coordinates <- par("usr")
  for (i in 1:(length(hit_indices)-1)){
    par(mar = c(0, 5, 0, 2))
    # plot()
    plot(0,type='n',axes=FALSE, ylab = '')
    plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = list(gene.sets[i], font = 1, cex = 0.3), yaxt = "n",
         ylab = '', xlim = c(1, length(gsea.rnk$metric)))
    abline(v = hit_indices[[i]], lwd = 0.75, col = colors[i])
    mtext(gene_sets[i], side = 3, cex = 0.7)
  }
  for (i in length(hit_indices)){
    par(mar = c(0, 5, 0, 2))
    plot(0,type='n',axes=FALSE, ylab = '')
    plot(0, type = "n", xaxs = "i", xlab = "Rank in ordered gene list", yaxt = "n",
         ylab = '', xlim = c(1, length(gsea.rnk$metric)))
    abline(v = hit_indices[[i]], lwd = 0.75, col = colors[i])
    mtext(gene_sets[i], side = 3, cex = 0.7)
  }
  mtext('Rank in ordered gene list', side = 1, padj = 2.5)
  # for (i in length(hit_indices)){
  #   plot(gsea.rnk$metric, type = "n", xaxs = "i",
  #        xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
  #        ylim = metric.range, yaxs = "i",
  #        )#,
  #        # panel.first = abline(h = seq(metric.range[1] / 2,
  #        #                              metric.range[2] - metric.range[1] / 4,
  #        #                              metric.range[2] / 2), col = "gray95", lty = 2))
  #   abline(v = hit_indices[[i]], lwd = 0.75, col = colors[i])
  # }
  
  
}

## GSEA dot plot

GSEA_dot_plot <- function(directory, gene.sets, comparisons, title = ''){
  
  # hit_indices <- c()
  compar <- c()
  sets <- c()
  NES <- c()
  padj <- c()
  # rnk <- c()
  
  for (i in comparisons){
    path <- paste0(directory, i)
    paths <- list.files(path = directory, pattern = i)
    
    for (j in paths){
      neg_path <- list.files(path = paste0(directory, '/', j, '/'), pattern = 'gsea_report_for_na_neg')
      pos_path <- list.files(path = paste0(directory, '/', j, '/'), pattern = 'gsea_report_for_na_pos')
      
      neg_path <- grep(pattern = '.xls$', neg_path, value = TRUE)
      pos_path <- grep(pattern = '.xls$', pos_path, value = TRUE)
      
      # print(paste(neg_path, pos_path))
      
      neg_report <- read.table(paste0(directory, '/', j, '/', neg_path), sep = '\t', header = TRUE)
      pos_report <- read.table(paste0(directory, '/', j, '/', pos_path), sep = '\t', header = TRUE)
      # pos_report <- as.data.frame(read_excel(paste0(directory, '/', j, '/', pos_path)))
      
      full_report <- rbind(neg_report, pos_report)
      
      idxs <- which(full_report$NAME %in% gene.sets)
      
      if (length(idxs) < 1){next}
      
      report <- full_report[idxs,]
      
      # print(report$NES)
      # print(report$FDR.q.val)
      # 
      compar <- c(compar, rep(i, length(idxs)))
      sets <- c(sets, as.character(report$NAME))
      NES <- c(NES, report$NES)
      padj <- c(padj, report$FDR.q.val)
      # print(head(report))
    }
    # print(paths)
  }
  
  gene_set_df <- data.frame(comparison = compar, gene_set = sets, NES = NES, pval=padj)
  # gene_set_df
  
  color_ramp = c('dodgerblue','firebrick1','firebrick2','firebrick3','firebrick4') #c('blue','red')
  col_func <- colorRampPalette(color_ramp)#colorRampPalette(color_ramp)
  
  col_map <- col_func(100)
  # scaled_pval <- log10(gene_set_df$pval+1)
  scaled_pval <- gene_set_df$pval*100 #scales::rescale(c(0,0.2, gene_set_df$pval), to = c(0,100))  #scales::rescale(gene_set_df$pval, to = c(0,1))
  
  # scaled_pval[which(scaled_pval == 0)]
  
  pval_colors <- col_map[as.integer(scaled_pval)+1] #unlist(lapply(scales::rescale(gene_set_df$pval, to = c(0,1)), function(x){col_func(as.numeric(x))}))
  
  ggplot(data=gene_set_df, aes(y=gene_set, x=comparison, size=pval, color=NES)) + 
    geom_point() + scale_color_gradient2(midpoint=0, low="blue", mid="white",
                                         high="red", space ="Lab" )+
    scale_size_continuous(name="p-value",range=c(8,1))+
    ggtitle(title) + theme(legend.position="right") + ylab('Gene Set') + xlab('Comparison')
  
}

##

GSEA_multiple_enrichment_data <- function(path, gene.sets, metric.range,
                                          enrichment.score.range) {
  
  # colors <- c('red', 'dodgerblue1', 'forestgreen', 'gold1', 'darkorchid1', 'cyan2', 'chartreuse2', 'deeppink2')
  
  # if (missing(path)) {
  #   stop("Path argument is required")
  # }
  # if (!file.exists(path)) {
  #   stop("The path folder could not be found. Please change the path")
  # }
  # if (missing(gene.set)) {
  #   stop("Gene set argument is required")
  # }
  # 
  
  hit_indices <- c()
  es_profile <- c()
  rnk <- c()
  es_score <- c()
  pval <- c()
  
  for (i in gene.sets){
    path.rnk <- list.files(path = file.path(path, "edb"),
                           pattern = ".rnk$", full.names = TRUE)
    gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
    # print(path.rnk)
    colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
    # if (missing(metric.range)) {
    #   metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
    # }  
    # metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
    # rnk[i] <- list(gsea.rnk$metric) #list(gsea.metric)
    ## Load .edb data
    path.edb <- list.files(path = file.path(path, "edb"),
                           pattern = ".edb$", full.names = TRUE)
    gsea.edb <- read.delim(file = path.edb,
                           header = FALSE, stringsAsFactors = FALSE)
    gsea.edb <- unlist(gsea.edb)
    gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
    gsea.metric <- unlist(strsplit(gsea.metric, " "))
    gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
    gsea.metric <- gsub("METRIC=", "", gsea.metric)
    gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
    # Select the right gene set
    if (length(gsea.edb) == 0) {
      stop(paste("The gene set name was not found, please provide",
                 "a correct name"))
    }
    if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", i), " "), gsea.edb)) > 1) {
      warning(paste("More than 1 gene set matched the gene.set",
                    "argument; the first match is plotted"))
    }
    gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", i), " "), gsea.edb)[1]]
    if (is.na(gsea.edb[[1]])){next}
    # Get template name
    gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
    gsea.edb <- unlist(strsplit(gsea.edb, " "))
    gsea.template <- gsea.edb[1]
    
    # Get gene set name
    gsea.gene.set <- gsea.edb[2]
    gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
    
    # Get enrichment score
    gsea.enrichment.score <- gsea.edb[3]
    gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
    
    # Get gene set name
    gsea.normalized.enrichment.score <- gsea.edb[4]
    gsea.normalized.enrichment.score <- gsub("NES=", "",
                                             gsea.normalized.enrichment.score)
    
    # Get nominal p-value
    gsea.p.value <- gsea.edb[5]
    gsea.p.value <- gsub("NP=", "", gsea.p.value)
    gsea.p.value <- as.numeric(gsea.p.value)
    
    # Get FDR
    gsea.fdr <- gsea.edb[6]
    gsea.fdr <- gsub("FDR=", "", gsea.fdr)
    gsea.fdr <- as.numeric(gsea.fdr)
    
    # Get hit indices
    gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
    gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
    gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
    gsea.hit.indices <- as.integer(gsea.hit.indices)
    
    # Get ES profile
    gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
    gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
    gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
    gsea.es.profile <- as.numeric(gsea.es.profile)
    
    hit_indices[i] <- list(gsea.hit.indices)
    es_profile[i] <- list(gsea.es.profile)
    rnk[i] <- list(gsea.rnk$metric)
    es_score[i] <- gsea.normalized.enrichment.score
    pval[i] <- gsea.fdr
  }
  
  data <- c()
  data$hit_indices <- hit_indices
  data$es_profile <- es_profile
  data$rnk <- rnk
  data$es_score <- es_score
  data$pval <- pval
  
  return(data)
}

## plot

GSEA_multiple_enrichment_plot <- function(data, metric.range,
                                          enrichment.score.range, title = "Enrichment") {
  
  colors <- c('red', 'dodgerblue1', 'forestgreen', 'gold1', 'darkorchid1', 'cyan2', 'chartreuse2', 'deeppink2', 'darkorange','darkgreen')
  
  
  hit_indices <- data$hit_indices
  es_profile <- data$es_profile
  rnk <- data$rnk
  es_score <- data$es_score
  pval <- data$pval
  
  gsea.rnk <- rnk[[1]]
  gene.sets <- names(hit_indices)
  
  # enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
  es_max <- max(unlist(es_profile))
  es_min <- min(unlist(es_profile))
  enrichment.score.range <- c(es_min, es_max)
  
  # metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
  rnk_min <- min(unlist(rnk))
  rnk_max <- max(unlist(rnk))
  metric.range <- c(rnk_min, rnk_max)
  
  ## Create GSEA plot
  # Save default for resetting
  # def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  # dev.new(width = 10, height = 20)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(seq(1, (2*length(hit_indices)+2))), heights = c(5, rep(c(0.2, 0.5), length(hit_indices)), 0.5))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, hit_indices[[1]], length(rnk[[1]])),
       c(0, es_profile[[1]], 0), type = "p", pch = 16, col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", ylab = "Enrichment score (ES)",
       ylim = enrichment.score.range,
       main = list(title, font = 1, cex = 2),
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk)),
       # panel.first = abline(h = seq(metric.range[1] / 2,
       #                              metric.range[2] - metric.range[1] / 4,
       #                              metric.range[2] / 2), col = "gray95", lty = 2),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  for (i in 2:length(hit_indices)){
    print(i)
    points(c(1, hit_indices[[i]], length(rnk[[i]])),
          c(0, es_profile[[i]], 0), col = colors[i], pch = 16)
  }
  plot.coordinates <- par("usr")
  for (i in 1:(length(hit_indices)-1)){
    subtitle <- paste(gene.sets[i], ' [ES SCore: ', es_score[i],', FDR q-val: ', pval[i],']', sep = '')
    par(mar = c(0, 5, 0, 2))
    # plot()
    plot(0,type='n',axes=FALSE, ylab = '')
    plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = list(gene.sets[i], font = 1, cex = 0.3), yaxt = "n",
         ylab = '', xlim = c(1, length(gsea.rnk)))
    abline(v = hit_indices[[i]], lwd = 0.75, col = colors[i])
    mtext(subtitle, side = 3, cex = 0.7)
  }
  for (i in length(hit_indices)){
    subtitle <- paste(gene.sets[i], ' [ES SCore: ', es_score[i],', FDR q-val: ', pval[i],']', sep = '')
    par(mar = c(0, 5, 0, 2))
    plot(0,type='n',axes=FALSE, ylab = '')
    plot(0, type = "n", xaxs = "i", xlab = "Rank in ordered gene list", yaxt = "n",
         ylab = '', xlim = c(1, length(gsea.rnk)))
    abline(v = hit_indices[[i]], lwd = 0.75, col = colors[i])
    mtext(subtitle, side = 3, cex = 0.7)
  }
  mtext('Rank in ordered gene list', side = 1, padj = 2.5)
  # for (i in length(hit_indices)){
  #   plot(gsea.rnk$metric, type = "n", xaxs = "i",
  #        xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
  #        ylim = metric.range, yaxs = "i",
  #        )#,
  #        # panel.first = abline(h = seq(metric.range[1] / 2,
  #        #                              metric.range[2] - metric.range[1] / 4,
  #        #                              metric.range[2] / 2), col = "gray95", lty = 2))
  #   abline(v = hit_indices[[i]], lwd = 0.75, col = colors[i])
  # }
  
  
}

