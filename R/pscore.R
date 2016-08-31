
# Package implementation of laplacian score for variable ranking ----------


# check exists of packages needed
check_pkgs = function(check_pkg_list) {
  exist_pkg_list = installed.packages()[, 1]
  for (pkg in check_pkg_list) {
    cat("Check package:", pkg, "... ")
    if (is.element(pkg, exist_pkg_list)) {
      cat("Exists\n")
    } else {
      cat("Not exists/Install\n")
      install.packages(pkg)
    }
  }
}


# calculate mean absolute deviation
meanabsdev = function(y) {
  return( mean(abs(y - mean(y, na.rm = TRUE)), na.rm = TRUE) )
}


# calculate distance matrix
calc_dist_mat = function(x, method = "euclidean", stand = FALSE) {
  if (stand) {
    x = scale(x = x, scale = apply(X = x, MARGIN = 2, FUN = meanabsdev))
  }
  return( as.matrix(x = dist(x = x, method = method)) )
}


# calculate affinity matrix
calc_affin_mat = function(x, stand = FALSE, affin_method = "euclidean", rbf_gamma) {
  return( exp(-calc_dist_mat(x = x, method = affin_method, stand = stand) / rbf_gamma) )
}


# calculate graph laplacian matrix
calc_lapl_mat = function(affin_mat) {
  return( diag(as.vector(affin_mat %*% rep(1, ncol(affin_mat)))) - affin_mat )
}


#' Laplacian scoring for variables
#' 
#' @param x a data frame or matrix, row-wised input varible data for scoring
#' @param stand a logical value indicating if the varibles to be scaled
#' @param affin_method distance measure for affinity matrix calculation
#' @param max_core_num integer, number of cores/threads used for parallelized clustering
#' @return a data frame with significance scores for the variables
#' @export
laplacian_score = function(x, stand = FALSE, affin_method = "euclidean", 
                           max_core_num = getOption("max_core_num", 2L)) {
  cat("-------- Parallelized laplacian scoring --------\n")
  
  # validate input parameters
  cat("Validate input parameters\n")
  max_core_num = as.integer(max_core_num)
  if (nrow(x) < 1 || ncol(x) < 1) {
    stop("Error on the parameter 'x': Input feature matrix is NULL")
  }
  if (max_core_num < 1) {
    stop("Error on the parameter 'max_core_num': Number of cores used must be greater than 0")
  }
  
  # check exists of necessary packages
  check_pkg_list = c("parallel", "Matrix", "data.table")
  check_pkgs(check_pkg_list = check_pkg_list)
  # load packages
  cat("Load packages\n")
  lapply(check_pkg_list, require, character.only = TRUE)
  cat("\n")
  
  # prepare parameter
  x = as.matrix(x)
  one_in_mat = matrix(data = rep(1, nrow(x)), nrow = nrow(x), ncol = 1) # column vector
  
  cat("Calculate affinity matrix\n")
  ptm_lapl_step = proc.time()
  affin_mat = calc_affin_mat(x = x, stand = stand, rbf_gamma = 100)
  print(proc.time() - ptm_lapl_step)
  
  cat("Calculate graph laplacian\n")
  ptm_lapl_step = proc.time()
  ar_sum_mat = diag(as.vector(affin_mat %*% one_in_mat))
  lapl_mat = ar_sum_mat - affin_mat 
  print(proc.time() - ptm_lapl_step)
  
  cat("Calculate laplacian scores\n")
  lapl_score_fea = function(fea_idx) {
    fea_in_mat = matrix(data = x[, fea_idx], nrow = length(x[, fea_idx]), ncol = 1) # column vector
    norm_fea = ((t(fea_in_mat) %*% ar_sum_mat %*% one_in_mat) / (t(one_in_mat) %*% ar_sum_mat %*% one_in_mat))
    norm_fea = fea_in_mat - as.numeric(norm_fea) * one_in_mat
    lapl_score = (t(norm_fea) %*% lapl_mat %*% norm_fea) / (t(norm_fea) %*% ar_sum_mat %*% norm_fea)
    lapl_score = 1 - lapl_score
    
    lapl_score[is.na(lapl_score)] = 0
    lapl_score_df = data.frame("laplacian.score" = lapl_score, 
                               "adj.laplacian.score" = ifelse(lapl_score > 0, lapl_score, 0))
    return(lapl_score_df)
  }
  cat("In scoring ... ")
  ptm_lapl_step = proc.time()
  laplacian_score_list = mclapply(X = c(1:ncol(x)), FUN = lapl_score_fea, mc.cores = max_core_num)
  cat("Done in", (proc.time() - ptm_lapl_step)[3], "seconds", "\n")
  
  cat("Reshape scoring results\n")
  ptm_lapl_step = proc.time()
  laplacian_scores = as.data.frame(rbindlist(laplacian_score_list))
  remove("laplacian_score_list")
  print(proc.time() - ptm_lapl_step)
  
  # return laplacian scores
  cat("-------- Parallelized laplacian scoring --------\n")
  return(laplacian_scores)
}


#' Bagging laplacian scoring for variables
#' 
#' @param x a data frame or matrix, row-wised input varible data for scoring
#' @param stand a logical value indicating if the varibles to be scaled
#' @param affin_method distance measure for affinity matrix calculation
#' @param bag_vol integer, sample volume in each bag iteration
#' @param max_core_num integer, number of cores/threads used for parallelized clustering
#' @return a data frame with significance scores for the variables
#' @export
bag_lapl_score = function(x, stand = FALSE, affin_method = "euclidean", bag_vol, 
                          max_core_num = getOption("max_core_num", 2L)) {
  cat("-------- Bagging parallelized laplacian scoring --------\n")
  
  # validate input parameters
  cat("Validate input parameters\n")
  bag_vol = as.integer(bag_vol)
  max_core_num = as.integer(max_core_num)
  if (nrow(x) < 1 || ncol(x) < 1) {
    stop("Error on the parameter 'x': Input feature matrix is NULL")
  }
  if ( (bag_vol < 1) || (bag_vol > nrow(x)) ) {
    stop("Error on the parameter 'bag_vol': Number of bag volume must be greater than 0")
  }
  if (max_core_num < 1) {
    stop("Error on the parameter 'max_core_num': Number of cores used must be greater than 0")
  }
  
  # calculate bag number/last bag vol/sample index for bagging
  bag_num = nrow(x) %/% bag_vol
  last_bag_vol = nrow(x) %% bag_vol
  perm_sample_idx = sample(x = 1:nrow(x), replace = FALSE)
  
  # process on bags
  bag_num = ifelse(last_bag_vol>0, bag_num+1, bag_num)
  laplacian_score_df = as.data.frame(matrix(0, ncol = bag_num, nrow = ncol(x)))
  for (i in 1:bag_num) {
    cat("Scoring on bag", i, "/", bag_num, ": ")
    bag_start_idx = (i - 1) * bag_vol + 1
    bag_end_idx = ifelse(i < bag_num, i * bag_vol, length(perm_sample_idx))
    bag_sample_idx = perm_sample_idx[bag_start_idx:bag_end_idx]
    
    cat("with", length(bag_sample_idx), "samples\n")
    laplacian_score_df[, i] = laplacian_score(x = x[bag_sample_idx, ], stand = stand, affin_method = affin_method, 
                                              max_core_num = max_core_num)$adj.laplacian.score
    cat("\n")
  }
  
  # bagging measurement
  laplacian_score_df = cbind(laplacian_score_df, 
                             "adj.laplacian.score" = apply(X = laplacian_score_df, MARGIN = 1, FUN = median))
  
  # return laplacian scores
  cat("-------- Bagging parallelized laplacian scoring --------\n")
  return(laplacian_score_df)
}


#' Fisher scoring for variables
#' 
#' @param x a data frame or matrix, row-wised input varible data for scoring
#' @param y an integer vector, response vector with one label for each row/compoonent of x
#' @param stand a logical value indicating if the varibles to be scaled
#' @param max_core_num integer, number of cores/threads used for parallelized clustering
#' @return a vector with significance scores for the variables
#' @export
fisher_score = function(x, y, stand = FALSE, max_core_num = getOption("max_core_num", 2L)) {
  cat("-------- Parallelized fisher scoring --------\n")
  
  # validate input parameters
  cat("Validate input parameters\n")
  max_core_num = as.integer(max_core_num)
  if (nrow(x) < 1 || ncol(x) < 1) {
    stop("Error on the parameter 'x': Input feature matrix is NULL")
  }
  if (!is.vector(y) || length(y) < 1 || nrow(x) != length(y)) {
    stop("Error on the parameter 'y': Input lable vector is inconsistent to feature matrix x")
  }
  if (max_core_num < 1) {
    stop("Error on the parameter 'max_core_num': Number of cores used must be greater than 0")
  }
  
  # check exists of necessary packages
  check_pkg_list = c("parallel", "data.table")
  check_pkgs(check_pkg_list = check_pkg_list)
  # load packages
  cat("Load packages\n")
  lapply(check_pkg_list, require, character.only = TRUE)
  cat("\n")
  
  # prepare parameter
  class_label_vec = unique(y)
  if (length(class_label_vec) <= 1) {
    stop("Error on the parameter 'y': Number of classes must be greater than 1")
  }
  class_num_vec = simplify2array(mclapply(X = class_label_vec, 
                                          FUN = function(one_label) { return(length(which(y == one_label))) }, 
                                          mc.cores = max_core_num))
  class_num_vec = class_num_vec / length(y)
  class_info_df = data.frame("class_label" = class_label_vec, "class_num" = class_num_vec)
  remove("class_label_vec", "class_num_vec")
  
  # finisher scoring in parallel
  cat("Calculate fisher scores\n")
  fisher_score_fea = function(fea_idx) {
    all_mean_vec = rep(x = mean(x = x[, fea_idx], na.rm = TRUE), times = length(class_info_df$class_label))
    class_mean_var_list = lapply(X = class_info_df$class_label, 
                                 FUN = function(one_class_label) {
                                   return(data.frame("mean" = mean(x = x[y == one_class_label, fea_idx], na.rm = TRUE),
                                                     "var" = var(x = x[y == one_class_label, fea_idx], na.rm = TRUE)))})
    class_mean_var_df = as.data.frame(rbindlist(class_mean_var_list))
    remove("class_mean_var_list")
    fisher_score = sum((class_mean_var_df$mean - all_mean_vec) ^ 2 * class_info_df$class_num) / 
      sum(class_mean_var_df$var * class_info_df$class_num)
    return(fisher_score)
  }
  cat("In scoring ... ")
  ptm_lapl_step = proc.time()
  fisher_scores_vec = simplify2array(mclapply(X = c(1:ncol(x)), FUN = fisher_score_fea, mc.cores = max_core_num))
  cat("Done in", (proc.time() - ptm_lapl_step)[3], "seconds", "\n")
  
  # return laplacian scores
  cat("-------- Parallelized fisher scoring --------\n")
  return(fisher_scores_vec)
}

