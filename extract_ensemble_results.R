# This is the script to extract the results from ensemble learning
library(SummarizedExperiment)
library(dplyr)
command_args <- commandArgs(trailingOnly=TRUE)
outfile_name <- command_args[1]

files <- list.files(pattern = ".RDS")
get_AUC_from_scores <- function(testWithScore_single_set, combineSignatureName, 
                                include.last = FALSE) {
  cv_num <- length(testWithScore_single_set)
  re <- lapply(1:cv_num, function(i) {
    x <- testWithScore_single_set[[i]]
    n <- length(x)
    if (!include.last) {
      n <- length(x) - 1    
    }
    final <- lapply(1:n, function(j) {
      dt <- x[[j]]
      re_roc <- pROC::roc(dt$TBStatus, dt[, combineSignatureName])
      both_cutoff <- pROC::coords(re_roc, "best", best.method = "youden", 
                                  transpose = TRUE)
      if (!is.null(dim(both_cutoff))) {
        both_cutoff <- both_cutoff[, 1] # Maximize sensitivity
      }
      auc <- max(re_roc$auc, 1 - re_roc$auc)
      sens <- both_cutoff["sensitivity"]
      spec <- both_cutoff["specificity"]
      data.frame(AUC = auc, Sensitivity = round(sens, 2), 
                 Specificity = round(spec, 2),
                 Study = names(x)[j], sample_size = length(dt$TBStatus))
    }) |> 
      dplyr::bind_rows()
    data.frame(final, cv_num = paste0("CV_", i))
  }) |> 
    dplyr::bind_rows()
  return(re)
}

out_all <- lapply(files, function(file) {
  out <- readRDS(file) |> 
    get_AUC_from_scores("combine", FALSE)
  
  file_name <- gsub(".RDS", "", file)
  out$Set <- file_name
  out
}) |> 
  dplyr::bind_rows()

write.table(out_all,file = paste0(outfile_name, ".txt"), 
            quote = FALSE, sep = '\t')
# Example for import table
# out_all <- read.delim("out_all.txt")



