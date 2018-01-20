#' Modify gene annotation file
#'
#' This function creates 3 files in .RData format, which are block annotation,
#' gene range lookup table, and gene size lookup table, from the gene annotation
#' file in gtf format downloaded from Ensembl for later analysis.
#'
#' The block annotation is a list has three levels. First is a list of chromosomes. Second is a list of
#' genes. Third is a data frame of block annotations.
#'
#' The gene range lookup table is a list has two levels. First is a list of chromosomes.
#' Second is a data frame of gene range information.
#'
#' The gene size lookup table is a list has one level, which is a list of gene size information.
#'
#'@param gene_ann_path A physical path to gene annotation file in gtf format downloaded from Ensembl.
#'@param line_skip The number of lines that need to be skipped when reading the gene annotation file. Default is 5.
#'@param sep The separation used in the gene annotation file. Default is tab.
#'@param gene_id A fixed start string of gene id in the gene annotation file column 9.
#'Default is 'gene_id'.
#'@param block_ann_path A physical path to save block annotation in .RData format.
#'@param gene_range_path A physical path to save gene range lookup table in .RData format.
#'@param gene_size_path A physical path to save gene size lookup table in .RData format.
#'
#'@examples
#'modify_ann("./Homo.gtf", 5, "./block_ann.RData", "./gene_range.RData", "./gene_size.RData")
#'
#'@return none.
#'
#'@export
# modify gene annotation file--------------------------------------------
modify_ann <- function(gene_ann_path, line_skip = 5, sep = '\t', gene_id = 'gene_id',
                       block_ann_path, gene_range_path, gene_size_path){

  # prepare data structure
  ann <- prepare_ann(gene_ann_path, line_skip, sep, gene_id)
  ann <- ann[order(ann$V1, ann$V9), ]

  # create block annotation
  res <- create_block_cpp(ann)
  block_ann <- res[['block']]
  gene_range <- res[['range']]
  gene_size <- res[['size']]

  # save result
  save(block_ann, file = block_ann_path)
  save(gene_range, file = gene_range_path)
  save(gene_size, file = gene_size_path)
}



#' Summarize RNA-seq single end read
#'
#' This function summarize the RNA-seq single end read into read type.
#'
#' @param input_sam_folder_path_0 A physical path to a folder of mapped RNA-seq samples in .sam format
#' under condition 0.
#' @param input_sam_folder_path_1 A physical path to a folder of mapped RNA-seq samples in .sam format
#' under condition 1.
#' @param block_ann_ls The block annotation created by function \code{\link{modify_ann}}.
#' If set to NULL, the function will use the block annotation saved in the package for human. Default is NULL.
#' @param gene_range_ls The gene range lookup table created by function \code{\link{modify_ann}}.
#' If set to NULL, the function will use the gene range lookup table saved in the package for human. Default is NULL.
#' @param run.parallel A logical variable. If TRUE, run in parallel mode. Default is TRUE.
#' @param core.num CPU core numbers will be used if run.parallel is TRUE. Default is the number of all the cores.
#' @param min_overlap minimum overlap needed for summarizing read to read type, default is 1.
#' @examples
#' summarize_read_single_end("./group0", "./group1")
#' summarize_read_single_end("./group0", "./group1", block_ann, gene_range)
#'
#' @return A data.frame saving the summarized read.
#'
#' @export
# summarize RNA-seq single end read-------------------------------------------------
summarize_read_single_end <- function(input_sam_folder_path_0, input_sam_folder_path_1,
                                      block_ann_ls = NULL, gene_range_ls = NULL, run.parallel = TRUE,
                                      core.num = parallel::detectCores(),
                                      min_overlap = 1) {

  # get sam file path for group 0 and 1
  input_sam_file_path_0 <- list.files(path = input_sam_folder_path_0, pattern = '\\.sam', full.names = TRUE)
  input_sam_file_path_1 <- list.files(path = input_sam_folder_path_1, pattern = '\\.sam', full.names = TRUE)

  # summarize read
  read_summarized_list <- list()
  if (is.null(block_ann_ls) && !is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_single_end_each(x, block_ann, gene_range_ls, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_single_end_each(input_sam_file_path_0[i], block_ann, gene_range_ls, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_single_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann, gene_range_ls, min_overlap)
      }
    }
  }
  else if (!is.null(block_ann_ls) && is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_single_end_each(x, block_ann_ls, gene_range, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_single_end_each(input_sam_file_path_0[i], block_ann_ls, gene_range, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_single_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann_ls, gene_range, min_overlap)
      }
    }
  }
  else if (is.null(block_ann_ls) && is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_single_end_each(x, block_ann, gene_range, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_single_end_each(input_sam_file_path_0[i], block_ann, gene_range, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_single_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann, gene_range, min_overlap)
      }
    }
  }
  else if (!is.null(block_ann_ls) && !is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_single_end_each(x, block_ann_ls, gene_range_ls, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_single_end_each(input_sam_file_path_0[i], block_ann_ls, gene_range_ls, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_single_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann_ls, gene_range_ls, min_overlap)
      }
    }
  }

  # merge summarized read
  temp <- read_summarized_list[[1]]
  names(temp)[3] <- "group0_sample1"
  read_summarized_dt <- temp
  for (i in 2 : length(input_sam_file_path_0)) {
    temp <- read_summarized_list[[i]]
    names(temp)[3] <- paste("group0_sample", i, sep = '')
    read_summarized_dt <- merge(read_summarized_dt, temp, by = c("read_gene", "read_type"), all = TRUE)
  }
  for (i in 1 : length(input_sam_file_path_1)) {
    temp <- read_summarized_list[[i + length(input_sam_file_path_0)]]
    names(temp)[3] <- paste("group1_sample", i, sep = '')
    read_summarized_dt <- merge(read_summarized_dt, temp, by = c("read_gene", "read_type"), all = TRUE)
  }

  # output summarized data
  read_summarized_df <- as.data.frame(read_summarized_dt)
  temp <- read_summarized_df[-c(1, 2)]
  temp[is.na(temp)] <- 0
  read_summarized_df[-c(1, 2)] <- temp
  read_summarized_df
}



#' Summarize RNA-seq paired end read
#'
#' This function summarize the RNA-seq paired end read into read type.
#'
#' @param input_sam_folder_path_0 A physical path to a folder of mapped RNA-seq sample in .sam format
#' under condition 0.
#' @param input_sam_folder_path_1 A physical path to a folder of mapped RNA-seq sample in .sam format
#' under condition 1.
#' @param block_ann_ls The block annotation created by function \code{\link{modify_ann}}.
#' If set to NULL, the function will use the block annotation saved in the package for human. Default is NULL.
#' @param gene_range_ls The gene range lookup table created by function \code{\link{modify_ann}}.
#' If set to NULL, the function will use the gene range lookup table saved in the package for human. Default is NULL.
#' @param run.parallel A logical variable. If TRUE, run in parallel mode. Default is TRUE.
#' @param core.num CPU core numbers will be used if run.parallel is TRUE. Default is the number of all the cores.
#' @param min_overlap minimum overlap needed for summarizing read to read type, default is 1.
#' @examples
#' summarize_read_paired_end("./group0", "./group1")
#' summarize_read_paired_end("./group0", "./group1", block_ann, gene_range)
#'
#' @return A data.frame saving the summarized read.
#'
#' @export
# summarize RNA-seq paired end read
summarize_read_paired_end <- function(input_sam_folder_path_0, input_sam_folder_path_1,
                                      block_ann_ls = NULL, gene_range_ls = NULL, run.parallel = TRUE,
                                      core.num = parallel::detectCores(),
                                      min_overlap = 1) {

  # get sam file path for group 0 and 1
  input_sam_file_path_0 <- list.files(path = input_sam_folder_path_0, pattern = '\\.sam', full.names = TRUE)
  input_sam_file_path_1 <- list.files(path = input_sam_folder_path_1, pattern = '\\.sam', full.names = TRUE)

  # summarize read
  read_summarized_list <- list()
  if (is.null(block_ann_ls) && !is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_paired_end_each(x, block_ann, gene_range_ls, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_paired_end_each(input_sam_file_path_0[i], block_ann, gene_range_ls, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_paired_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann, gene_range_ls, min_overlap)
      }
    }
  }
  else if (!is.null(block_ann_ls) && is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_paired_end_each(x, block_ann_ls, gene_range, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_paired_end_each(input_sam_file_path_0[i], block_ann_ls, gene_range, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_paired_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann_ls, gene_range, min_overlap)
      }
    }
  }
  else if (is.null(block_ann_ls) && is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_paired_end_each(x, block_ann, gene_range, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_paired_end_each(input_sam_file_path_0[i], block_ann, gene_range, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_paired_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann, gene_range, min_overlap)
      }
    }
  }
  else if (!is.null(block_ann_ls) && !is.null(gene_range_ls)) {
    if (run.parallel == TRUE) {
      read_summarized_list <- parallel::mclapply(c(input_sam_file_path_0, input_sam_file_path_1),
                                                 function(x) summarize_read_paired_end_each(x, block_ann_ls, gene_range_ls, min_overlap),
                                                 mc.cores = core.num)
    }
    else {
      for (i in 1 : length(input_sam_file_path_0)) {
        print(i)
        read_summarized_list[[i]] <- summarize_read_paired_end_each(input_sam_file_path_0[i], block_ann_ls, gene_range_ls, min_overlap)
      }
      for (i in 1 : length(input_sam_file_path_1)) {
        print(i)
        read_summarized_list[[i + length(input_sam_file_path_0)]] <- summarize_read_paired_end_each(input_sam_file_path_1[i],
                                                                                                    block_ann_ls, gene_range_ls, min_overlap)
      }
    }
  }

  # merge summarized read
  temp <- read_summarized_list[[1]]
  names(temp)[3] <- "group0_sample1"
  read_summarized_dt <- temp
  for (i in 2 : length(input_sam_file_path_0)) {
    temp <- read_summarized_list[[i]]
    names(temp)[3] <- paste("group0_sample", i, sep = '')
    read_summarized_dt <- merge(read_summarized_dt, temp, by = c("read_gene", "read_type"), all = TRUE)
  }
  for (i in 1 : length(input_sam_file_path_1)) {
    temp <- read_summarized_list[[i + length(input_sam_file_path_0)]]
    names(temp)[3] <- paste("group1_sample", i, sep = '')
    read_summarized_dt <- merge(read_summarized_dt, temp, by = c("read_gene", "read_type"), all = TRUE)
  }

  # output summarized data
  read_summarized_df <- as.data.frame(read_summarized_dt)
  temp <- read_summarized_df[-c(1, 2)]
  temp[is.na(temp)] <- 0
  read_summarized_df[-c(1, 2)] <- temp
  read_summarized_df
}



#' Fit SGNB model
#'
#' This function estiamte the mle and returns p-value and estimated gene expression
#' level in a data.frame.
#'
#' @param read_summarized_df The returned data.frame from function \code{\link{summarize_read_single_end}} or
#'                           \code{\link{summarize_read_paired_end}}.
#' @param gene_size_ls The gene size file. If set it to NULL, the function will use the gene_size file saved in the
#'                     package data folder, which is for human. Default is NULL.
#' @param min_reduce A number between 0 and 1 indicating the minimum percentage of the number of parameters
#'                   needed to be reduced. If the parameter reduction procedure fail to reduce this amount,
#'                   then group all read types into one read type. Default is 0.3.
#' @param tol A numerical value to control estimation accuracy. Default is 0.001.
#' @param times A numerical value to control estimation loop times. Default is 100.
#'
#' @examples
#' fit_SGNB(read_summarized_df)
#' fit_SGNB(read_summarized_df, gene_size, min_reduce, tol = 0.01, times = 500)
#'
#' @return A data.frame saving gene_id, p-value, likelihood test statistics, degree of freedom and gene expression
#'         for condition 0 and 1.
#'
# fit SGNB model-----------------------------------------------------------
fit_SGNB <- function(read_summarized_df, gene_size_ls = NULL, min_reduce = 0.3, tol = 0.001, times = 100) {
  read_summarized_df <- read_summarized_df[order(read_summarized_df$read_gene, read_summarized_df$read_type), ]

  # detect sample number under two conditions
  sample_num_0 <- sum(stringr::str_detect(names(read_summarized_df), pattern = 'group0_'))
  sample_num_1 <- sum(stringr::str_detect(names(read_summarized_df), pattern = 'group1_'))
  group_sample_num <- c(sample_num_0, sample_num_1)

  # calculate TMM normalization factor
  TMM_dt <- data.table::data.table(read_summarized_df)
  TMM_dt[, read_type := NULL]
  TMM_dt <- TMM_dt[, lapply(.SD, sum), by = read_gene]
  TMM_dt[, read_gene := NULL]
  lib_size <- colSums(TMM_dt)
  tmm_norm_factor <- TMM(as.data.frame(TMM_dt), lib_size)
  lib_size_norm <- lib_size * tmm_norm_factor

  # group read types
  read_type_group_df <- create_read_type_group_cpp(unique(read_summarized_df$read_gene), read_summarized_df$read_gene,
                                                   read_summarized_df$read_type, min_reduce)
  read_count_dt <- data.table::data.table(read_summarized_df)
  read_type_group_dt <- data.table::data.table(read_type_group_df)
  read_count_dt <- merge(read_count_dt, read_type_group_dt, by = c("read_gene", "read_type"))
  read_count_dt[, read_type := NULL]
  read_count_dt <- read_count_dt[, lapply(.SD, sum), by = c("read_gene", "read_type_group")]
  read_count_df <- as.data.frame(read_count_dt)

  # fit model
  if (is.null(gene_size_ls)) {
    res <- fit_SGNB_cpp(unique(read_count_df$read_gene), read_count_df$read_gene,
                        as.matrix(read_count_df[-c(1, 2)]), lib_size_norm, group_sample_num,
                        gene_size, tol, times)
  }
  else if (!is.null(gene_size_ls)) {
    res <- fit_SGNB_cpp(unique(read_count_df$read_gene), read_count_df$read_gene,
                        as.matrix(read_count_df[-c(1, 2)]), lib_size_norm, group_sample_num,
                        gene_size_ls, tol, times)
  }
  res
}


#' Fit exact SGNB model
#'
#'
#' @param read_summarized_df The returned data.frame from function \code{\link{summarize_read_single_end}} or
#'                           \code{\link{summarize_read_paired_end}}.
#' @param simplify Whether apply model simplification. Default is TRUE.
#' @param combine Method for calculating global p-value. Could be "Fisher" or "Bonferroni". Default is "Fisher".
#' @param common_disp If set to TRUE, then estimate a common dispersion for all the isoforms of a gene. Defaultis FALSE.
#' @param min_reduce A number between 0 and 1 indicating the minimum percentage of the number of parameters
#'                   needed to be reduced. If the parameter reduction procedure fail to reduce this amount,
#'                   then group all read types into one read type. Default is 0.
#' @param tol A numerical value to control estimation accuracy. Default is 0.001.
#' @param times A numerical value to control estimation loop times. Default is 100.
#'
#' @examples
#' fit_SGNB_exact(read_summarized_df)
#'
#' @return A data.frame having gene_id, p-value and relative gene expression for condition 0 and 1.
#'
#' @export
# fit exact SGNB model-----------------------------------------------------------
fit_SGNB_exact <- function(read_summarized_df, simplify = TRUE, combine = "Bonferroni", common_disp = FALSE,
                           min_reduce = 0, gene_size_ls = NULL, tol = 0.001, times = 100) {
  read_summarized_df <- read_summarized_df[order(read_summarized_df$read_gene, read_summarized_df$read_type), ]
  # detect sample number under two conditions
  sample_num_0 <- sum(stringr::str_detect(names(read_summarized_df), pattern = 'group0_'))
  sample_num_1 <- sum(stringr::str_detect(names(read_summarized_df), pattern = 'group1_'))
  group_sample_num <- c(sample_num_0, sample_num_1)
  # calculate TMM normalization factor
  TMM_dt <- data.table::data.table(read_summarized_df)
  TMM_dt[, read_type := NULL]
  TMM_dt <- TMM_dt[, lapply(.SD, sum), by = read_gene]
  TMM_dt[, read_gene := NULL]
  lib_size <- colSums(TMM_dt)
  tmm_norm_factor <- TMM(as.data.frame(TMM_dt), lib_size)
  lib_size_norm <- lib_size * tmm_norm_factor
  # group read types
  if (simplify == TRUE) {
    read_type_group_df <- create_read_type_group_cpp(unique(read_summarized_df$read_gene), read_summarized_df$read_gene,
                                                     read_summarized_df$read_type, min_reduce)
    read_count_dt <- data.table::data.table(read_summarized_df)
    read_type_group_dt <- data.table::data.table(read_type_group_df)
    read_count_dt <- merge(read_count_dt, read_type_group_dt, by = c("read_gene", "read_type"))
    read_count_dt[, read_type := NULL]
    read_count_dt <- read_count_dt[, lapply(.SD, sum), by = c("read_gene", "read_type_group")]
    read_count_df <- as.data.frame(read_count_dt)
  } else if (simplify == FALSE) {
    read_count_df <- read_summarized_df
    read_count_df$read_type <- seq(1, dim(read_count_df)[1])
    names(read_count_df)[2] <- "read_type_group"
  }
  # fit model
  if (common_disp == FALSE) {
    temp <- fit_SGNB_exact_cpp(read_count_df$read_gene, read_count_df$read_type_group, as.matrix(read_count_df[-c(1, 2)]),
                               lib_size_norm, group_sample_num, tol, times)
  } else if (common_disp == TRUE) {
    temp <- fit_SGNB_exact_common_cpp(unique(read_count_df$read_gene), read_count_df$read_gene,
                                      read_count_df$read_type_group, as.matrix(read_count_df[-c(1, 2)]),
                                      lib_size_norm, group_sample_num, tol, times)
  }
  # summarize results
  temp_result <- temp$results
  #temp_pseudo_data <- temp$`pseudo data`
  if (combine == "Fisher") {
    res_pvalue <- aggregate(p_value ~ gene_id, temp_result,
                            function(x) pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE))
  } else if (combine == "Bonferroni") {
    res_pvalue <- aggregate(p_value ~ gene_id, temp_result, function(x) min(1, min(x) * length(x)))
  }
  res_theta0 <- aggregate(theta0 ~ gene_id, temp_result, sum)
  res_theta1 <- aggregate(theta1 ~ gene_id, temp_result, sum)
  res <- merge(res_pvalue, res_theta0, by.all = "gene_id", all = TRUE)
  res <- merge(res, res_theta1, by.all = "gene_id", all = TRUE)
  if (is.null(gene_size_ls)) {
    gene_size_df <- data.frame(gene_id = names(gene_size), size = unlist(gene_size))
    res <- merge(res, gene_size_df, by = 'gene_id', all.x = TRUE)
    res$theta0 <- res$theta0 / res$size
    res$theta1 <- res$theta1 / res$size
  }
  else if (!is.null(gene_size_ls)) {
    gene_size_df <- data.frame(gene_id = names(gene_size_ls), size = unlist(gene_size_ls))
    res <- merge(res, gene_size_df, by = 'gene_id', all.x = TRUE)
    res$theta0 <- res$theta0 / res$size
    res$theta1 <- res$theta1 / res$size
  }
  res <- res[, colnames(res) != "size"]
  res
}

















