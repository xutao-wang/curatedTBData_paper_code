# This is the script to generate figures for curatedTBData paper
library(curatedTBData)
library(ggplot2)
boxplotTBSig <- function(object_list, annotationColName, signatureColNames) {
    .check_input(object_list)
    sig_list <- .signature_filter(object_list, annotationColName, signatureColNames)
    ## Combine list of data frame
    rbindx <- function(dfs) {
        ns <- lapply(dfs, colnames) |>
            unlist() |>
            unique()
        do.call(rbind, lapply(dfs, function(x) {
            for (n in ns[!ns %in% colnames(x)]) {
                x[[n]] <- NA
            }
            x
        }))
    }
    sig_data <- rbindx(sig_list)
    ## Check whether annotationColName or signatureColNames are all NA's
    if (all(is.na(sig_data[, signatureColNames]))) {
        stop(sprintf("Gene signature: %s is not found from the entire input.",
                     signatureColNames))
    } else if (all(is.na(sig_data[, annotationColName]))) {
        stop(sprintf("Annotation name: %s is not found from the entire input.",
                     annotationColName))
    }
    study_name <- unique(sig_data$Study)
    p_boxplot <- lapply(study_name, function(x, signatureColNames) {
        sig_data_gse <- sig_data |>
            dplyr::filter(.data$Study == x)
        sig_data_gse$anno_names <- sig_data_gse[, annotationColName] |>
            factor()
        is_signatureColNames_na <- sig_data_gse |>
            dplyr::select(signatureColNames) |>
            is.na() |>
            all()
        if (is_signatureColNames_na) {
            sprintf("Gene signature: %s not found from the input list.",
                    signatureColNames) |>
                paste("NULL is returnted") |>
                message()
            return(NULL)
        }
        sig_data1 <-
            SummarizedExperiment::SummarizedExperiment(colData = sig_data_gse)
        ## Create a custom color scale to deal with different factors
        anno_levels <- levels(sig_data1$anno_names)
        n <- length(anno_levels)
        if (n > 9L) {
            sprintf("The number of levels under %s",
                    annotationColName) |>
                paste("is greater than 9.")  |>
                paste("Only first 9 levels is included.") |>
                message()
            n <- 9
            ## Remove observations for the last level
            sig_data1 <- sig_data1[, sig_data1$anno_names %in%
                                       anno_levels[seq_len(n)]]
            sig_data1$anno_names <- factor(sig_data1$anno_names)
        } else if (n <= 3L) {
            n <- 3
        }
        myColors <- RColorBrewer::brewer.pal(n, "Set1")
        names(myColors) <- levels(sig_data1$anno_names)
        p <-
            TBSignatureProfiler::signatureBoxplot(inputData = sig_data1,
                                                  name = x,
                                                  signatureColNames = signatureColNames,
                                                  annotationColName = "anno_names",
                                                  rotateLabels = FALSE,
                                                  fill_colors = myColors)
        p1 <- p +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 12,
                                                              face = "bold"),
                           legend.position = "none",
                           axis.title.x = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_text(colour = "Black",
                                                               size = 12,
                                                               hjust = 0.5,
                                                               face = "bold"),
                           axis.text.y = ggplot2::element_text(size = 12,
                                                               angle = 0,
                                                               hjust = 0.5))
        p1
    }, signatureColNames)
    ## Remove empty element from list
    p_boxplot <- p_boxplot[!vapply(p_boxplot, is.null, TRUE)]
    n_sqrt <- sqrt(length(p_boxplot))
    p_combine <- do.call(gridExtra::"grid.arrange",
                         c(p_boxplot, ncol = floor(n_sqrt)))
    return(p_combine)
}

.signature_filter <- function(object_list, annotationColName,
                              signatureColNames) {
    obj_name <- names(object_list)
    object_list_seq <- seq_len(length(object_list))
    sig_list1 <- lapply(object_list_seq, function(i, signatureColNames) {
        x <- object_list[[i]]
        col_info <- SummarizedExperiment::colData(x)
        GSE <- rep(names(object_list[i]), nrow(col_info))
        index <- match(signatureColNames, colnames(col_info)) |>
            stats::na.omit()
        anno <- col_info[, annotationColName]
        if (length(index) == 0L) {
            sprintf("Gene signature: %s not found in study: %s,",
                    signatureColNames, obj_name[i]) |>
                paste("NA is returned.") |>
                message()
            result <- data.frame(anno, NA, GSE)
        } else {
            result <- data.frame(anno, col_info[, index], GSE)
        }
        colnames(result) <- c(annotationColName, signatureColNames, "Study")
        result
    }, signatureColNames)
    names(sig_list1) <- obj_name
    return(sig_list1)
}

combine_auc <- function(SE_scored_list, annotationColName, signatureColNames,
                        num.boot = NULL, percent = 0.95, AUC.abs = FALSE,
                        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
    .check_input(SE_scored_list)
    bpparam <- BPPARAM
    if (is.null(num.boot)) {
        paste("\"num.boot\" is NULL",
              "Bootstrap Confidence Interval is not computed.") |>
            message()
    }
    aucs_result <- BiocParallel::bplapply(SE_scored_list, function(x) {
        .get_auc_stats(x, annotationColName, signatureColNames, num.boot,
                       percent, AUC.abs)
    }, BPPARAM = bpparam)
    aucs_result_dat <- do.call(rbind, aucs_result)
    if (nrow(aucs_result_dat) == 0) {
        msg <- sprintf(" \"signatureColNames\": %s is/are not found in the list",
                       paste(signatureColNames, collapse = ", "))
        paste(msg, "in the study. Check \"signatureColNames\".") |>
            stop(call. = FALSE)
    }
    ## Re-order data based on their median AUC (from largest to smallest)
    ## Remove NA value
    aucs_result_dat_median <- aucs_result_dat |>
        dplyr::filter(!is.na(.data$AUC)) |>
        dplyr::group_by(.data$Signature) |>
        dplyr::summarise_all(stats::median) |>
        dplyr::arrange(dplyr::desc(.data$AUC))
    ## Order signatures based on median AUC values
    Signature_order <- as.character(aucs_result_dat_median$Signature)
    Sig_NA <- aucs_result_dat |>
        dplyr::filter(is.na(.data$AUC)) |>
        dplyr::select(.data$Signature) |>
        unlist(use.names = FALSE) |>
        unique()
    ## Re-order gene signature
    ## re-level this step is to let ridge plot ordered based on median value
    sig_levels <- unique(c(Signature_order, Sig_NA))
    aucs_result_dat$Signature <- factor(aucs_result_dat$Signature,
                                        levels = sig_levels)
    ## Label name of each data under column 'Study'
    aucs_result_dat$Study <- gsub("\\..*", "", row.names(aucs_result_dat))
    row.names(aucs_result_dat) <- NULL
    return(aucs_result_dat)
}

.get_auc_stats <- function(SE_scored, annotationColName, signatureColNames,
                           num.boot, percent, AUC.abs) {
    ## Check signatureColNames
    col_info <- SummarizedExperiment::colData(SE_scored)
    index <- match(signatureColNames, colnames(col_info)) |>
        stats::na.omit()
    if (length(index) == 0) {
        msg <- sprintf(" \"signatureColNames\": %s is/are not found",
                       paste(signatureColNames, collapse = ", "))
        paste(msg, "in the study. NULL is returned.\n") |>
            message()
    }
    signatureColNames <- colnames(col_info)[index]
    ## Check annotationColName
    index_anno <- match(annotationColName, colnames(col_info))
    if (is.na(index_anno)) {
        sprintf("\"annotationColName\": %s is not found from the study.\n",
                annotationColName) |>
            stop(call. = FALSE)
    }
    annotationData <- col_info[annotationColName][, 1] |>
        as.character() |>
        as.factor()
    ## Check levels of annotationData
    anno_level_len <- length(unique(annotationData))
    if (anno_level_len != 2L) {
        paste("Annotation data should have exactly two levels.",
              "The number of input levels is:",
              anno_level_len, ".\n") |>
            stop(call. = FALSE)
    }
    ## Get AUC value for each signature along with corresponding datasets
    if (is.null(num.boot)) {
        sig_result <- lapply(signatureColNames,
                             function(i, SE_scored, annotationData) {
                                 score <- col_info[i][, 1] |>
                                     as.vector()
                                 ## Deal with scores that have constant value (e.g. Sloot_HIV_2)
                                 if (length(unique(score)) == 1L) {
                                     sprintf(paste("Constant score found for siganture: %s,",
                                                   "results will be NA.\n"), i) |>
                                         message()
                                     dat <- data.frame(Signature = i, P.value = NA,
                                                       neg10xP.value = NA, AUC = NA)
                                     return(dat)
                                 }
                                 pvals <- stats::t.test(score ~ annotationData)$p.value
                                 neg10log <- -1 * log(pvals + 1e-4)
                                 out <- .get_stats(annotationData, score)
                                 # pred <- ROCit::rocit(score, annotationData)
                                 # if (AUC.abs) {
                                 #     aucs <- pred$AUC
                                 # } else {
                                 #     aucs <- max(pred$AUC, 1 - pred$AUC)
                                 # }
                                 ## Create data frame for With signature, P.value, AUC
                                 # data.frame(Signature = i, P.value = round(pvals, 4),
                                 #            neg10xP.value = round(neg10log, 4),
                                 #            AUC = round(aucs, 4))
                                 cbind(Signature = i, P.value = round(pvals, 4),
                                       neg10xP.value = round(neg10log, 4),
                                       out)
                             }, SE_scored, annotationData)
        result <- do.call(rbind, sig_result) |>
            as.data.frame()
        row.names(result) <- NULL
        return(result)
    } else {
        sig_result <- lapply(signatureColNames,
                             function(i, SE_scored, annotationData, percent) {
                                 score <- col_info[i][, 1]
                                 ## Get lower and upper quantile
                                 lower <- (1 - percent) / 2
                                 upper <- 1 - lower
                                 ## Deal with PLAGE that have constant score (e.g. Sloot_HIV_2)
                                 if (length(unique(score)) == 1L) {
                                     sprintf(paste("Constant score found for siganture: %s,",
                                                   "results will be NA.\n"), i) |>
                                         message()
                                     dat <- data.frame(i, NA, NA, NA, NA, NA)
                                     colnames(dat) <- c("Signature", "P.value", "neg10xP.value",
                                                        "AUC",
                                                        paste0("CI lower.", lower * 100, "%"),
                                                        paste0("CI upper.", upper * 100, "%"))
                                     return(dat)
                                 }
                                 pvals <- stats::t.test(score ~ annotationData)$p.value
                                 neg10log <- -1 * log(pvals + 1e-4)
                                 out <- .get_stats(annotationData, score)
                                 
                                 ## Calculate bootstrapped AUC confidence interval.
                                 ## Repeated sampling scores and annotationData
                                 ## compute the AUC for the sampled pairs
                                 bootCI <- lapply(seq_len(num.boot),
                                                  function(j, score, annotationData) {
                                                      index <- sample(seq_len(length(score)), replace = TRUE)
                                                      tmp_score <- score[index]
                                                      tmp_annotationData <- annotationData[index]
                                                      ## Consider when re-sampling only has 1 cases, remove it
                                                      if (length(unique(tmp_annotationData)) == 2L) {
                                                          .get_stats(tmp_annotationData, tmp_score)
                                                      } else {
                                                          data.frame(AUC = NA, 
                                                                     Sensitivity = NA,
                                                                     Specificity = NA)
                                                      }
                                                  }, score, annotationData)
                                 bootCI <- bootCI |>
                                     dplyr::bind_rows()
                                 bootCI <- bootCI[complete.cases(bootCI), ]
                                      
                                 LowerAUC <- stats::quantile(bootCI$AUC, prob = lower, na.rm = TRUE)
                                 UpperAUC <- stats::quantile(bootCI$AUC, prob = upper, na.rm = TRUE)
                                 
                                 LowerSens <- stats::quantile(bootCI$Sensitivity, prob = lower, na.rm = TRUE)
                                 UpperSens <- stats::quantile(bootCI$Sensitivity, prob = upper, na.rm = TRUE)
                                 
                                 LowerSpec <- stats::quantile(bootCI$Specificity, prob = lower, na.rm = TRUE)
                                 UpperSpec <- stats::quantile(bootCI$Specificity, prob = upper, na.rm = TRUE)
                                 
                                 dat <- data.frame(i, round(pvals, 4),
                                                   round(neg10log, 4), 
                                                   round(out$AUC, 4),
                                                   round(LowerAUC, 4),
                                                   round(UpperAUC, 4),
                                                   round(out$Sensitivity, 4),
                                                   round(LowerSens, 4),
                                                   round(UpperSens, 4),
                                                   round(out$Specificity, 4),
                                                   round(LowerSpec, 4),
                                                   round(UpperSpec, 4))
                                 colnames(dat) <- c("Signature", "P.value",
                                                    "neg10xP.value", 
                                                    "AUC",
                                                    paste0("AUC CI lower.", lower * 100, "%"),
                                                    paste0("AUC CI upper.", upper * 100, "%"),
                                                    "Sensitivity",
                                                    paste0("Sens CI lower.", lower * 100, "%"),
                                                    paste0("Sens CI upper.", upper * 100, "%"),
                                                    "Specificity",
                                                    paste0("Spec CI lower.", lower * 100, "%"),
                                                    paste0("Spec CI upper.", upper * 100, "%"))
                                 dat
                             }, SE_scored, annotationData, percent)
        result <- do.call(rbind, sig_result)
        row.names(result) <- NULL
        return(result)
    }
}

.get_stats <- function(annotationData, score) {
    myroc <- suppressMessages(pROC::roc(annotationData, score))
    both_cutoff <- pROC::coords(myroc, "best", 
                                best.method = "youden", 
                                transpose = TRUE)
    if (!is.null(dim(both_cutoff))) {
        both_cutoff <- both_cutoff[, 1] # Maximize sensitivity
    }
    sens <- both_cutoff["sensitivity"]
    spec <- both_cutoff["specificity"]
    
    data.frame(AUC = round(myroc$auc, 4), 
               Sensitivity = round(sens, 4), 
               Specificity = round(spec, 4))
}
.check_input <- function(object_list) {
    ## Check if the input is a list
    if (!is.list(object_list)) {
        paste("Function only supports a list of",
              "SummarizrdExperiment objetcs.") |>
            stop(call. = FALSE)
    }
    check_element_class <- vapply(object_list, function(x)
        methods::is(x, "SummarizedExperiment"), TRUE)
    if (!all(check_element_class)) {
        msg <- sprintf("Elements(s) %s ",
                       paste0(which(!check_element_class), collapse = ", "))
        paste("Function only supports class: SummarizedExperiment.", msg,
              "in the list is/are not SummarizedExperiment object.") |>
            stop(call. = FALSE)
    }
    ## Check valid list names
    list_name <- names(object_list)
    if (is.null(list_name)) {
        paste("names of the input list should not be NULL.",
              "Add unique name for each element from the list.") |>
            stop(call. = FALSE)
    } else if (!is.na(match("", list_name))) {
        paste("Names of the input contains \"\".",
              "Replace  \"\" with unique character.") |>
            stop(call. = FALSE)
    }
}

get_mean_auc <- function(df, column_name_variable, column_name_value,
                         method = c("percentile", "empirical"),
                         num.boot = 100, percent = 0.95) {
    ## Select signatures and associated AUC
    ## split them into list based on signature
    df_list <- df |>
        dplyr::select(c(column_name_variable, column_name_value)) |>
        dplyr::group_split(.data[[column_name_variable]])
    sigInfo <- df_list |>
        vapply(function(x) as.character(unlist(x[column_name_variable])[1]),
               character(1))
    ## Get summarized table and
    ## bootstrap 95% Confidence Interval for the mean AUC
    sigInfo <- data.frame(Signature = sigInfo)
    if (missing(method)) {
        method <- "empirical"
        paste("Missing method argument.",
              "Using the default method: empirical") |>
            message()
    } else {
        method <- match.arg(method)
    }
    sprintf("Use %s method to compute Bootstrap Confidence Interval",
            method) |>
        message()
    meanAUC_list <- lapply(df_list, function(x)
        .bootstrap_mean_CI(x, column_name_value, method, num.boot, percent))
    meanAUC <- do.call(rbind, meanAUC_list)
    return(cbind(sigInfo, meanAUC))
}

.bootstrap_mean_CI <- function(df, column_name_value, method,
                               num.boot, percent = 0.95) {
    lower <- (1 - percent) / 2
    upper <- 1 - lower
    x <- unlist(df[, column_name_value], use.names = FALSE)
    ## Remove NA's (e.g. in PLAGE method)
    x <- stats::na.omit(x)
    n <- length(x)
    if (n == 1L) {
        xbar <- x
        ci <- data.frame(round(xbar, 4), NA, NA)
        colnames(ci) <- c("MeanAUC", paste0("CI lower.", lower * 100, "%"),
                          paste0("CI upper.", upper * 100, "%"))
        row.names(ci) <- NULL
        return(ci)
    }
    ## Sample mean
    xbar <- mean(x)
    ## Random re-samples from x
    bootstrapsample <- lapply(seq_len(num.boot), function(i)
        sample(x, n, replace = TRUE))
    bootstrapsample <- do.call(cbind, bootstrapsample)
    ## Compute the means x∗
    bsmeans <- colMeans(bootstrapsample)
    if (method == "empirical") {
        ## Compute deltastar for each bootstrap sample
        deltastar <- bsmeans - xbar
        ## Find the 0.025 and 0.975 quantile for deltastar
        d <- stats::quantile(deltastar, c(lower, upper), na.rm = TRUE)
        ## Calculate the confidence interval for the mean.
        ci <- xbar - c(d[2], d[1])
    } else if (method == "percentile") {
        ci <- stats::quantile(bsmeans, c(lower, upper), na.rm = TRUE)
    }
    ## Set upper and lower bound for the confidence interval
    lower_ci <- round(ci[1], 4)
    lower_ci <- ifelse(lower_ci <= 0.5, 0.5, lower_ci)
    upper_ci <- round(ci[2], 4)
    upper_ci <- ifelse(upper_ci >= 1, 1, upper_ci)
    ci <- data.frame(round(xbar, 4), lower_ci, upper_ci)
    colnames(ci) <- c("MeanAUC", paste0("CI lower.", lower * 100, "%"),
                      paste0("CI upper.", upper * 100, "%"))
    row.names(ci) <- NULL
    return(ci)
}

heatmap_auc <- function(combine_dat, GSE_sig = NULL,
                        facet = TRUE, clustering = TRUE, show_avg = FALSE) {
    ## Subset input data.frame with the desired column names
    ## check whether column contains 'Signature', 'Study', 'AUC'
    expect_name <- c("Signature", "Study", "AUC")
    index_name <- match(expect_name, colnames(combine_dat))
    if (any(is.na(index_name))) {
        stop(sprintf("Column with name(s): %s is/are missing.",
                     paste0(expect_name[is.na(index_name)], collapse = ", ")))
    }
    dat <- combine_dat[, index_name]
    if (!is.factor(dat$Study)) {
        dat$Study <- factor(dat$Study)
    }
    data_wide <- reshape2::dcast(dat, stats::formula("Study ~ Signature"),
                                 value.var = "AUC")
    row.names(data_wide) <- data_wide$Study
    ## remove study name column
    dat_input <- as.matrix(data_wide[, -1])
    dat_input[is.na(dat_input)] <- NA
    if (length(unique(dat$Study)) > 1L && clustering == TRUE) {
        ## Clustering AUC values if the number of studies is greater than 1
        ## Transform form long to wide data:
        ## First column is the study names and column names is signatures
        ## This step is necessary for clustering
        dd <- stats::dist(dat_input)
        hc <- stats::hclust(dd)
        dat_input <- dat_input[hc$order, ]
    }
    if (show_avg) {
        ## Get mean AUC for each study across multiple gene signatures
        Avg <- rowMeans(dat_input, na.rm = TRUE)
        ## Transform back into long format
        datta <- cbind(dat_input, Avg = Avg) |>
            reshape2::melt()
    } else {
        datta <- dat_input |> 
            reshape2::melt()
    }
    if (!is.null(GSE_sig)) {
        GSE_sig <- .expand_study(GSE_sig)
        index <- lapply(seq_len(nrow(GSE_sig)), function(i) {
            kk <- datta[grep(GSE_sig$TBSignature[i], datta$Var2), ]
            kk$indx <- row.names(kk)
            indx <- kk[which(as.character(kk$Var1) %in%
                                 GSE_sig$Study[i]), "indx"]
            indx
        }) |>
            unlist(use.names = FALSE)
    } else {
        paste("GSE_sig not provided.",
              "Training data information is not available for the output") |>
            message()
        index <- NULL
    }
    ## Label signature type based on the input signatures
    signatureColNames <- as.character(unique(dat$Signature))
    sig_type_temp <- vapply(strsplit(signatureColNames, "_"),
                            function(x) x[2], character(1))
    sig_type_index <- sig_type_temp |>
        as.numeric() |>
        is.na() |>
        suppressWarnings()
    ## Get signature type
    sig_type <- sig_type_temp[sig_type_index] |>
        unique()
    ## Assign category:
    ## Disease for those do not have signature type: e.g. Anderson_42
    ## Get training data position index
    datta$trian <- FALSE
    datta$sig_typek <- "Disease"
    for (i in sig_type) {
        datta$sig_typek[grep(i, datta$Var2)] <- i
    }
    if (show_avg) {
        datta$sig_typek[grep("Avg", datta$Var2)] <- "Avg"
        datta$sig_typek <- factor(datta$sig_typek,
                                  levels = c("Avg", sig_type, "Disease"))
    } else {
        datta$sig_typek <- factor(datta$sig_typek,
                                  levels = c(sig_type, "Disease"))
    }
    
    datta[as.numeric(index), "trian"] <- TRUE
    ## Subset datta with training study and its associated signature(s)
    frames <- datta[datta$trian, c("Var1", "Var2", "sig_typek")]
    p <- ggplot2::ggplot(data = datta,
                         ggplot2::aes(x = .data$Var1, y = .data$Var2,
                                      fill = .data$value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_distiller("AUC values", palette = "RdPu", trans = "reverse") +
        ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2)),
                           cex = 3.5)
    if (facet) {
        p <- p + ggplot2::facet_grid(.data$sig_typek ~ ., switch = "y",
                                     scales = "free", space = "free")
        frame_facet <- .facet_rect_position(datta, frames)
        if (!nrow(frame_facet) == 0L) {
            p <- p + ggplot2::geom_rect(data = frame_facet,
                                        ggplot2::aes(xmin = .data$Var1 - 0.5,
                                                     xmax = .data$Var1 + 0.5,
                                                     ymin = .data$Var2 -  0.5,
                                                     ymax = .data$Var2 + 0.5),
                                        size = 1, fill = NA, colour = "black",
                                        inherit.aes = FALSE)
        }
    } else {
        frames$Var1 <- as.integer(frames$Var1)
        frames$Var2 <- as.integer(frames$Var2)
        p <- p +
            ggplot2::geom_rect(data = frames,
                               ggplot2::aes(xmin = .data$Var1 - 0.5,
                                            xmax = .data$Var1 + 0.5,
                                            ymin = .data$Var2 - 0.5,
                                            ymax = .data$Var2 + 0.5),
                               size = 1, fill = NA, colour = "black")
    }
    p <- p +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45,
                                                           vjust = 1,
                                                           size = 12,
                                                           hjust = 1),
                       axis.text.y = ggplot2::element_text(size = 12))
    return(p)
}

.facet_rect_position <- function(datta, frames) {
    # Split data frame into list based on different signature type
    frames_list <- frames |>
        dplyr::group_split(.data$sig_typek)
    names(frames_list) <- vapply(frames_list,
                                 function(x) as.character(x$sig_typek[1]),
                                 character(1))
    datta_list <- datta |>
        dplyr::group_split(.data$sig_typek)
    names(datta_list) <- vapply(datta_list,
                                function(x) as.character(x$sig_typek[1]),
                                character(1))
    ## Get the correct index in for training study change
    ## sig_type levels from sub list based on characters in the full list
    level_all <- levels(datta$Var2)
    frame_facet1 <- lapply(names(frames_list), function(i) {
        num_Var1 <- frames_list[[i]]$Var1 |>
            as.integer()
        frame_sig <- frames_list[[i]]$Var2
        datta_sig <- unique(datta_list[[i]]$Var2)
        num_Var2 <- factor(frame_sig, 
                           levels = level_all[which(level_all %in% datta_sig)]) |> 
            as.integer()
        data.frame(Var1 = num_Var1, Var2 = num_Var2, sig_typek = i)
    }) |> 
        dplyr::bind_rows() |> 
        as.data.frame()
    return(frame_facet1)
}

.expand_study <- function(GSE_sig) {
    n <- nrow(GSE_sig)
    col_name <- colnames(GSE_sig)
    ## check for column names:TBSignature and Study
    expect_name <- c("TBSignature", "Study")
    index_name <- match(expect_name, colnames(GSE_sig))
    if (any(is.na(index_name))) {
        stop(sprintf("Column with name(s): %s is/are missing.",
                     paste0(expect_name[is.na(index_name)], collapse = ", ")))
    }
    GSE_sig <- GSE_sig[, index_name]
    data_list <- lapply(seq_len(n), function(i) {
        study_vector <- strsplit(GSE_sig[, col_name[2]][i], split = "&")
        df <- data.frame(GSE_sig[, col_name[1]][i], study_vector)
        colnames(df) <- col_name
        df
    })
    re <- do.call(rbind, data_list) |>
        as.data.frame()
    return(re)
}

update_genenames <- function(siglist) {
    newgenes <- suppressMessages(suppressWarnings(
        HGNChelper::checkGeneSymbols(siglist,
                                     unmapped.as.na = FALSE)))$Suggested.Symbol
    ind <- grep("//", newgenes)
    if (length(ind) != 0) newgenes[ind] <- strsplit(newgenes[ind],
                                                    " /// ")[[1]][1]
    # if (any(newgenes != siglist))
    #    message("One or more gene names were altered.")
    return(newgenes)
}

extract_CI <- function(dat_combine) {
    signature_names <- unique(dat_combine$Signature)
    out <- lapply(signature_names, function(x) {
        dat_sub <- dat_combine |> 
            dplyr::filter(Signature == x)
        
        weighted_AUC <- weighted.mean(dat_sub$AUC, dat_sub$Size)
        AUC_CI_lower <- weighted.mean(dat_sub$`AUC CI lower.2.5%`, dat_sub$Size)
        AUC_CI_upper <- weighted.mean(dat_sub$`AUC CI upper.97.5%`, dat_sub$Size)
        AUC_CI <- sprintf("%.2f (%.2f - %.2f)", 
                          weighted_AUC, AUC_CI_lower, AUC_CI_upper)
        
        weighted_sens <- weighted.mean(dat_sub$Sensitivity, dat_sub$Size)
        Sens_CI_lower <- weighted.mean(dat_sub$`Sens CI lower.2.5%`, dat_sub$Size)
        Sens_CI_upper <- weighted.mean(dat_sub$`Sens CI upper.97.5%`, dat_sub$Size)
        Sens_CI <- sprintf("%.2f (%.2f - %.2f)", 
                           weighted_sens, Sens_CI_lower, Sens_CI_upper)
        
        weighted_spec <- weighted.mean(dat_sub$Specificity, dat_sub$Size)
        Spec_CI_lower <- weighted.mean(dat_sub$`Spec CI lower.2.5%`, dat_sub$Size)
        Spec_CI_upper <- weighted.mean(dat_sub$`Spec CI upper.97.5%`, dat_sub$Size)
        Spec_CI <- sprintf("%.2f (%.2f - %.2f)", 
                           weighted_spec, Spec_CI_lower, Spec_CI_upper)
        
        data.frame(Signature = x, AUC = AUC_CI, 
                   Sensitivity = Sens_CI, Specificity = Spec_CI)
    }) |> 
        dplyr::bind_rows()
    return(out)
}

PTB_Control_edit <- function(gsea_list, filter_size, out_file_path = NULL, 
                             save_file = TRUE) {
    gsea_list <- gsea_list[!names(gsea_list) %in% c("GSE74092", "GSE6112")]
    study_size <- lapply(gsea_list, ncol) |>
        unlist()
    gsea_list <-  gsea_list[study_size >= filter_size]
    if (save_file) {
        saveRDS(gsea_list, out_file_path) 
    }
    return(gsea_list)
}

PTB_LTBI_edit <- function(gsea_list, filter_size, out_file_path = NULL,
                          save_file = TRUE) {
    gsea_list <- gsea_list[names(gsea_list) != "GSE6112"]
    study_size <- lapply(gsea_list, ncol) |>
        unlist()
    gsea_list <-  gsea_list[study_size >= filter_size]
    if (save_file) {
        saveRDS(gsea_list, out_file_path)
    }
    return(gsea_list)
}

get_AUC_for_each_signature <- function(testWithScore_single_set, combineSignatureName,
                                       percent = 0.95, num.boot = NULL) {
    test_score_list_flat <- do.call(c, testWithScore_single_set)
    # Remove name with model and get unique studies
    unique_names <- names(test_score_list_flat) |> 
        unique()
    unique_names <- unique_names[unique_names != "Model"]
    test_score_list_input <- test_score_list_flat[unique_names]
    candidate_sig <- colnames(test_score_list_input[[1]])
    input_sig <- candidate_sig[-which(candidate_sig %in% c(combineSignatureName, "TBStatus"))]
    final_df <- lapply(1:length(test_score_list_input), function(i) {
        x <- test_score_list_input[[i]]
        single_study_auc <- lapply(input_sig, function(sig) {
            re_roc <- pROC::roc(x[, "TBStatus"], x[, sig])
            both_cutoff <- pROC::coords(re_roc, "best", best.method = "youden", 
                                        transpose = TRUE)
            if (!is.null(dim(both_cutoff))) {
                both_cutoff <- both_cutoff[, 1] # Maximize sensitivity
            }
            auc <- max(re_roc$auc, 1 - re_roc$auc)
            sens <- both_cutoff["sensitivity"]
            spec <- both_cutoff["specificity"]
            data.frame(Signature = sig, AUC = auc,
                       Sensitivity = round(sens, 2), 
                       Specificity = round(spec, 2))
        })
        cbind(Study = names(test_score_list_input)[i], 
              do.call(rbind, single_study_auc), sample_size = nrow(x))
    }) |> 
        dplyr::bind_rows()
    return(final_df)
}

extract_signature_score <- function(gsea_list, signature_names) {
    lapply(gsea_list, function(x) {
        colData(x) |> 
            as.data.frame() |> 
            dplyr::select(signature_names, "TBStatus") |> 
            list()
    })
}

make_ensembleSignaturesList <- function(signatureNames, n = 30, size = 10) {
    ensembleSignaturesList <- lapply(1:n, function(x) {
        set.seed(x)
        index <- sample(length(signatureNames), size)
        signatureNames[index]
    })
    names(ensembleSignaturesList) <- paste0("Set", 
                                            seq_len(length(ensembleSignaturesList)))
    return(ensembleSignaturesList)
}


find_sig_with_max_AUC <- function(gsea_AUC, gsea_ensl, ensemble_sig_list, 
                                  gsea_method) {
    gsea_ensl_list <- split(gsea_ensl, gsea_ensl$Set)
    test_study <- gsea_ensl_list[[1]] |> 
        dplyr::select(Study)
    # Get mean AUC for each signature using gene set scores
    gsea_AUC_mean <- gsea_AUC |>
        dplyr::group_by(Signature) |>
        dplyr::summarise(mean_AUC = weighted.mean(AUC, sample_size)) |>
        dplyr::mutate(ranks = dense_rank(-mean_AUC))
    # gsea_AUC_mean <- gsea_AUC |>
    #     dplyr::group_by(Signature) |>
    #     dplyr::summarise(mean_AUC = mean(AUC)) |>
    #     dplyr::mutate(ranks = dense_rank(-mean_AUC))
    n <- length(ensemble_sig_list)
    # Get signature with max AUC based on their mean
    signature_with_max_AUC <- lapply(1:n, function(i) {
        x <- ensemble_sig_list[[i]]
        sig_max <- gsea_AUC_mean |> 
            dplyr::filter(Signature %in% x) |> 
            dplyr::filter(ranks == min(ranks)) |> 
            dplyr::pull(Signature)

        data.frame(Set = names(ensemble_sig_list)[i], Signature = sig_max)
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::mutate(Signature_type = "Single") |> 
        dplyr::mutate(Method = gsea_method)
  
    # Get max siganture's AUC value
    max_signature_AUCs <- lapply(1:n, function(i) {
        sig_max <- signature_with_max_AUC |> 
            dplyr::filter(Set == names(ensemble_sig_list)[i]) |> 
            dplyr::pull(Signature)
        gsea_AUC_sig_max <- gsea_AUC |> 
            dplyr::filter(Signature == sig_max) |> 
            dplyr::right_join(test_study, by = "Study")
        gsea_AUC_sig_max <- gsea_AUC_sig_max[complete.cases(gsea_AUC_sig_max), ]
        weighted_mean_AUC <- vapply(1:1000, function(j) {
            set.seed(j^2)
            index <- sample(seq_len(nrow(gsea_AUC_sig_max)), replace = TRUE)
            weighted.mean(gsea_AUC_sig_max[index, "AUC"],
                          gsea_AUC_sig_max[index, "sample_size"])
            # mean(gsea_AUC_sig_max[index, "AUC"])
        }, numeric(1))
        
        data.frame(AUC = weighted_mean_AUC, Signature = sig_max) |> 
            dplyr::mutate(Set = names(ensemble_sig_list)[i]) |> 
            dplyr::mutate(Signature_type = "Single")
    }) |> 
        dplyr::bind_rows()
    
    gsea_ensl_sub <- lapply(1:n, function(i) {
        gsea_ensl_set <- gsea_ensl |> 
            dplyr::filter(Set == names(ensemble_sig_list)[i])
        mean_AUC <- vapply(1:1000, function(j) {
            set.seed(j^2)
            index <- sample(seq_len(nrow(gsea_ensl_set)), replace = TRUE)
            weighted.mean(gsea_ensl_set[index, "AUC"],
                          gsea_ensl_set[index, "sample_size"])
            # mean(gsea_ensl_set[index, "AUC"])
        }, numeric(1))
        data.frame(AUC = mean_AUC, Signature = names(ensemble_sig_list)[i]) |> 
            dplyr::mutate(Set = Signature) |> 
            dplyr::mutate(Signature_type = "Ensemble")
    }) |> 
        dplyr::bind_rows()
    
    sig_list_names_edit <- gsub("Set", "Set ", names(ensemble_sig_list))
    
    AUC_final <- rbind(gsea_ensl_sub, max_signature_AUCs) |> 
        dplyr::mutate(Method = gsea_method) |> 
        dplyr::mutate(Set = gsub("Set", "Set ", Set))
    AUC_final$Set <- factor(AUC_final$Set, levels = sig_list_names_edit)
    
    
    signature_with_max_AUC$Set <- gsub("Set", "Set ", signature_with_max_AUC$Set) |> 
        factor(levels = sig_list_names_edit)
        
    final_list <- list(AUC_final = AUC_final, 
                       signature_with_max_AUC = signature_with_max_AUC)
    return(final_list)
}
find_sig_with_max_AUC1 <- function(gsea_AUC, gsea_ensl, ensemble_sig_list, 
                                  gsea_method) {
    gsea_ensl_list <- split(gsea_ensl, gsea_ensl$Set)
    test_study <- gsea_ensl_list[[1]] |> 
        dplyr::select(Study)
    # Get mean AUC for each signature using gene set scores
    gsea_AUC_mean <- gsea_AUC |>
        dplyr::group_by(Signature) |>
        dplyr::summarise(mean_AUC = weighted.mean(AUC, sample_size)) |>
        dplyr::mutate(ranks = dense_rank(-mean_AUC))
    # gsea_AUC_mean <- gsea_AUC |>
    #     dplyr::group_by(Signature) |>
    #     dplyr::summarise(mean_AUC = mean(AUC)) |>
    #     dplyr::mutate(ranks = dense_rank(-mean_AUC))
    n <- length(ensemble_sig_list)
    # Get signature with max AUC based on their mean
    signature_with_max_AUC <- lapply(1:n, function(i) {
        x <- ensemble_sig_list[[i]]
        sig_max <- gsea_AUC_mean |> 
            dplyr::filter(Signature %in% x) |> 
            dplyr::filter(ranks == min(ranks)) |> 
            dplyr::pull(Signature)

        data.frame(Set = names(ensemble_sig_list)[i], Signature = sig_max)
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::mutate(Signature_type = "Single") |> 
        dplyr::mutate(Method = gsea_method)
  
    # Get max siganture's AUC value
    max_signature_AUCs <- lapply(1:n, function(i) {
        sig_max <- signature_with_max_AUC |> 
            dplyr::filter(Set == names(ensemble_sig_list)[i]) |> 
            dplyr::pull(Signature)
        gsea_AUC_sig_max <- gsea_AUC |> 
            dplyr::filter(Signature == sig_max) |> 
            dplyr::right_join(test_study, by = "Study")
        gsea_AUC_sig_max <- gsea_AUC_sig_max[complete.cases(gsea_AUC_sig_max), ]
        gsea_AUC_sig_max <- gsea_AUC_sig_max[!duplicated(gsea_AUC_sig_max), ]
        weighted_mean_AUC <- vapply(1:1000, function(j) {
            set.seed(j^2)
            index <- sample(seq_len(nrow(gsea_AUC_sig_max)), replace = TRUE)
            weighted.mean(gsea_AUC_sig_max[index, "AUC"],
                          gsea_AUC_sig_max[index, "sample_size"])
            # mean(gsea_AUC_sig_max[index, "AUC"])
        }, numeric(1))
        
        data.frame(AUC = weighted_mean_AUC, Signature = sig_max) |> 
            dplyr::mutate(Set = names(ensemble_sig_list)[i]) |> 
            dplyr::mutate(Signature_type = "Single")
    }) |> 
        dplyr::bind_rows()
    
    gsea_ensl_sub <- lapply(1:n, function(i) {
        gsea_ensl_set <- gsea_ensl |> 
            dplyr::filter(Set == names(ensemble_sig_list)[i]) |> 
            dplyr::group_by(Study) |> 
            dplyr::summarise(AUC = max(AUC), 
                             sample_size = unique(sample_size))
        mean_AUC <- vapply(1:1000, function(j) {
            set.seed(j^2)
            index <- sample(seq_len(nrow(gsea_ensl_set)), replace = TRUE)
            weighted.mean(gsea_ensl_set[index, "AUC"],
                          gsea_ensl_set[index, "sample_size"])
            # mean(gsea_ensl_set[index, "AUC"])
        }, numeric(1))
        data.frame(AUC = mean_AUC, Signature = names(ensemble_sig_list)[i]) |> 
            dplyr::mutate(Set = Signature) |> 
            dplyr::mutate(Signature_type = "Ensemble")
    }) |> 
        dplyr::bind_rows()
    
    sig_list_names_edit <- gsub("Set", "Set ", names(ensemble_sig_list))
    
    AUC_final <- rbind(gsea_ensl_sub, max_signature_AUCs) |> 
        dplyr::mutate(Method = gsea_method) |> 
        dplyr::mutate(Set = gsub("Set", "Set ", Set))
    AUC_final$Set <- factor(AUC_final$Set, levels = sig_list_names_edit)
    
    signature_with_max_AUC$Set <- gsub("Set", "Set ", signature_with_max_AUC$Set) |> 
        factor(levels = sig_list_names_edit)
        
    final_list <- list(AUC_final = AUC_final, 
                       signature_with_max_AUC = signature_with_max_AUC)
    return(final_list)
}

GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

get_results_for_ridge <- function(df_ensl, df_single, set_name, threshold = 0.8,
                                  ensemble_sig_info_list = ensemble_sig_list) {
    df_ensl_mean <- df_ensl |>
        dplyr::group_by(Study, Set) |>
        dplyr::summarise(AUC = max(AUC),
                         sample_size = mean(sample_size)) |>
        dplyr::mutate(Signature = Set) |>
        dplyr::select(-Set)
    # df_ensl_mean <- df_ensl |>
    #     dplyr::group_by(Study, Set) |>
    #     dplyr::summarise(AUC = mean(AUC), 
    #                      sample_size = mean(sample_size)) |>
    #     dplyr::mutate(Signature = Set) |>
    #     dplyr::select(-Set)

    set_ensl <- df_ensl_mean |> 
        dplyr::filter(Signature == set_name) |> 
        dplyr::select(Study, Signature, AUC, sample_size)
    
    studies <- set_ensl |> 
        dplyr::filter(AUC >= threshold) |> 
        dplyr::pull(Study)
    
    df_single_sig <- df_single |> 
        dplyr::filter(Signature %in% ensemble_sig_info_list[[set_name]])
    sig_order <- df_single_sig |> 
        dplyr::group_by(Signature) |> 
        dplyr::summarise(mean_AUC = weighted.mean(AUC, sample_size)) |> 
        dplyr::arrange(-1 * mean_AUC) |> 
        dplyr::pull(Signature) |> 
        as.character()
    df_combine <- df_single_sig |> 
        dplyr::select(Study, Signature, AUC, sample_size) |> 
        rbind(set_ensl) |> 
        dplyr::filter(Study %in% studies)
    
    df_combine$Signature <- factor(df_combine$Signature,
                                   levels = c(set_name,sig_order))
    return(df_combine)
}

# mean_aa <- ssgsea_PTB_LTBI_ensl |> 
#     dplyr::group_by(Set, Study) |> 
#     dplyr::summarise(AUC_mean = mean(AUC))
# median_aa <- ssgsea_PTB_LTBI_ensl |> 
#     dplyr::group_by(Set, Study) |> 
#     dplyr::summarise(AUC_median = median(AUC))
# wtd_mean_aa <- ssgsea_PTB_LTBI_ensl |> 
#     dplyr::group_by(Set, Study) |> 
#     dplyr::summarise(AUC_wtd_mean = weighted.mean(AUC, sample_size))
# 
# max_aa <- ssgsea_PTB_LTBI_ensl |> 
#     dplyr::group_by(Set, Study) |> 
#     dplyr::summarise(AUC_max = max(AUC))
# compare_re_tmp <- mean_aa |> 
#     dplyr::inner_join(max_aa, by = c("Set", "Study")) |> 
#     dplyr::inner_join(median_aa, by = c("Set", "Study")) |> 
#     dplyr::inner_join(wtd_mean_aa, by = c("Set", "Study")) |> 
#     dplyr::mutate(Diff_mean = AUC_max - AUC_mean) |> 
#     dplyr::mutate(Diff_median = AUC_max - AUC_median) |> 
#     dplyr::mutate(Diff_wtd_mean = AUC_max - AUC_wtd_mean) |> 
#     mutate(Diff_mean_median = AUC_mean - AUC_median) |> 
#     mutate(Diff_mean_wtd_mean = AUC_mean - AUC_wtd_mean)
