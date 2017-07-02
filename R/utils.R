# prepare annotation file---------------------------------------
prepare_ann <- function(gene_ann_path, line_skip, sep, gene_id) {
  # read gtf file into data frame for processing
  ann <- read.table(file = gene_ann_path, header = FALSE, sep = sep, skip = line_skip,
                    stringsAsFactors = FALSE)
  ann$V4 <- as.integer(ann$V4)
  ann$V5 <- as.integer(ann$V5)
  # select exons and eliminate useless columns (2, 3, 6, 8)
  ann <- ann[ann$V3 == 'exon', c(1, 4, 5, 7, 9)]
  ann$V9 <- stringr::str_extract(ann$V9, paste(gene_id, ' [^;]*;', sep = ''))
  ann$V9 <- stringr::str_replace_all(ann$V9, gene_id, '')
  ann$V9 <- stringr::str_replace_all(ann$V9, ';', '')
  ann$V9 <- stringr::str_trim(ann$V9, 'both')
  ann
}

# summarize single end read for one sam file-------------------------------
summarize_read_single_end_each <- function(input_sam_file_path, block_ann, gene_range, minOverlap) {
  read_type_df <- create_read_type_cpp(input_sam_file_path, block_ann, gene_range, minOverlap)
  read_type_dt <- data.table::data.table(read_type_df)
  temp <- read_type_dt[, length(read_type), by = "read_id"]
  read_type_dt <- merge(read_type_dt, temp, by = "read_id", all.x = TRUE)
  data.table::setkey(read_type_dt, V1)
  read_type_dt <- read_type_dt[list(1), list(read_id, read_gene, read_type)]
  read_type_dt <- read_type_dt[, length(read_id), keyby = c("read_gene", "read_type")]
  read_type_dt
}

# summarize paired end read for one sam file-------------------------------
summarize_read_paired_end_each <- function(input_sam_file_path, block_ann, gene_range, minOverlap) {
  read_type_df <- create_read_type_cpp(input_sam_file_path, block_ann, gene_range, minOverlap)
  read_type_dt <- data.table::data.table(read_type_df)
  temp <- read_type_dt[, length(read_type), by = "read_id"]
  read_type_dt <- merge(read_type_dt, temp, by = "read_id", all.x = TRUE)
  data.table::setkey(read_type_dt, V1)
  read_type_dt <- read_type_dt[list(2), list(read_id, read_gene, read_type)]
  temp <- read_type_dt[, length(read_type), by = c("read_id", "read_gene")]
  read_type_dt <- merge(read_type_dt, temp, by = c("read_id", "read_gene"), all.x = TRUE)
  data.table::setkey(read_type_dt, V1)
  read_type_dt <- read_type_dt[list(2), list(read_id, read_gene, read_type)]
  read_type_dt <- read_type_dt[, length(read_id), keyby = c("read_gene", "read_type")]
  read_type_dt
}

# TMM normalization
TMM <- function(count_df, lib_size)
{
  trim_M <- 0.3
  trim_A <- 0.1
  res <- rep(NA, length(lib_size))
  res_len <- dim(count_df)[2]
  ref_idx <- which(lib_size == median(lib_size))
  for(i in 1 : res_len)
  {
    count_df_curr <- count_df[c(ref_idx, i)]
    count_df_curr <- count_df_curr[count_df_curr[[1]] > 0 & count_df_curr[[2]] > 0, ]
    M <- log((count_df_curr[[2]] / lib_size[i]) / (count_df_curr[[1]] / lib_size[ref_idx]))
    A <- abs((1 / 2) * log((count_df_curr[[2]] / lib_size[i]) * (count_df_curr[[1]] / lib_size[ref_idx])))
    count_df_curr <- cbind(count_df_curr, M, A)
    count_df_curr <- count_df_curr[order(count_df_curr$M, count_df_curr$A), ]
    l <- dim(count_df_curr)[1]
    trim_size <- as.integer(l * trim_M)
    count_df_curr <- count_df_curr[(trim_size + 1) : (l - trim_size), ]
    count_df_curr <- count_df_curr[order(count_df_curr$A), ]
    l <- dim(count_df_curr)[1]
    trim_size <- as.integer(l * trim_A)
    count_df_curr <- count_df_curr[(trim_size + 1) : (l - trim_size), ]
    weight <- (lib_size[i] - count_df_curr[[2]]) / (lib_size[i] * count_df_curr[[2]]) +
              (lib_size[ref_idx] - count_df_curr[[1]]) / (lib_size[ref_idx] * count_df_curr[[1]])
    weight <- 1 / weight
    logTMM <- (weight %*% count_df_curr[[3]]) / sum(weight)
    TMM_factor <- exp(logTMM)
    res[i] <- TMM_factor
  }
  res
}


