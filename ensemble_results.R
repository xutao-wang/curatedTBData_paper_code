library(SummarizedExperiment)
library(dplyr)
# This is the script to generate ensemble learning results for curatedTBData paper
# All the results and code are stored in "/restricted/projectnb/tuberculosis/data/TB_reprocess_2019/Ensemble_learning_curatedTBData_2023"
#### Functions ####
make_ensembleSignaturesList <- function(signatureNames, n, size) {
    ensembleSignaturesList <- lapply(1:n, function(x) {
        set.seed(x)
        index <- sample(length(signatureNames), size)
        signatureNames[index]
    })
    names(ensembleSignaturesList) <- paste0("Set", 
                                            seq_len(length(ensembleSignaturesList)))
    return(ensembleSignaturesList)
}

transform_input_list <- function(input_list, signatureNames, 
                                 scale_input, norm_input) {
    out_list <- lapply(input_list, function(x) {
        col_info <- colData(x)
        col_info_sub <- col_info[, signatureNames]
        if (norm_input) {
            # normalization code from common sense [0, 1]
            col_info_sub <- apply(col_info_sub, 2, function(x)
                (x - min(x)) / (max(x) - min(x))) %>%
                DataFrame()
            
            # normalization code from gsva code
            # col_info_sub <- apply(col_info_sub, 2, 
            #                       function(x) x / (range(x)[2] - range(x)[1])) |> 
            #     DataFrame()
        } else if (scale_input) {
            col_info_sub <- apply(col_info_sub, 2, scale) %>%
                DataFrame()
        }
        # if (scale_input) {
        #     col_info_sub <- apply(col_info_sub, 2, scale) %>%
        #         DataFrame()
        # } else if (norm_input) {
        #     col_info_sub <- apply(col_info_sub, 2, function(x) 
        #         (x - min(x)) / (max(x) - min(x))) %>%
        #         DataFrame()
        # }
        row.names(col_info_sub) <- row.names(col_info)
        col_info_sub$TBStatus <- col_info$TBStatus
        data.frame(col_info_sub)
    })
    return(out_list)
}

ensemble_train_LASSO <- function(trainData, GSEATestResultScaled, intercept, 
                                 batchCorrection, combineSignatureName, 
                                 weights, runModel) {
    inputMatrix <- as.matrix(trainData[, which(names(trainData) != "TBStatus")])
    outcome1 <- ifelse(trainData$TBStatus == "PTB", "PTB", "Others")
    outcome <- as.integer(factor(outcome1), levels = c("Others", "PTB"))
    # By default alpha = 1, and Lasso regression is called
    signatureModel <- glmnet::cv.glmnet(x = inputMatrix, y = outcome, weights = weights,
                                        intercept = intercept, family = "binomial")
    results <- lapply(GSEATestResultScaled, function (x) {
        dataTest <- as.matrix(x[, match(colnames(inputMatrix), colnames(x))])
        if (batchCorrection) {
            dataTest <- make_batch_corrected_dt(inputMatrix, dataTest)
        }
        # lambda.min : λ of minimum mean cross-validated error
        # lambda.1se : largest value of λ such that error is within 1 standard error 
        # of the cross-validated errors for lambda.min.
        predScore <- stats::predict(signatureModel$glmnet.fit, 
                                    s = signatureModel$lambda.min,
                                    newx = dataTest, na.action = na.omit, 
                                    type = "response")
        x[, combineSignatureName] <- as.vector(predScore)
        x
    })
    if (runModel) {
        results$Model <- signatureModel
    }
    return(results)
}

multiStudyEnsemble_single <- function(GSEATrainResult, GSEATestResult, 
                                      signatureNames, combineSignatureName, 
                                      training_method, intercept, 
                                      scale_input, norm_input, batchCorrection,
                                      runModel) {
    
    GSEATrainResultScaled <- transform_input_list(GSEATrainResult, 
                                                  signatureNames,
                                                  scale_input, norm_input)
    GSEATestResultScaled <- transform_input_list(GSEATestResult, 
                                                 signatureNames,
                                                 scale_input, norm_input)
    # Stack data together
    trainData <- do.call(rbind, GSEATrainResultScaled)
    # inputMatrix <- as.matrix(trainData[, which(names(trainData) != "TBStatus")])
    if (training_method == "LASSO") {
        
        message("Calling LASSO")
        results <- ensemble_train_LASSO(trainData, GSEATestResultScaled, 
                                        intercept = intercept, 
                                        batchCorrection, combineSignatureName,
                                        weights = NULL, runModel)
    } else if (training_method == "glm_with_weight") {
        message("Calling glm_with_weight")
        w1 <- vapply(GSEATrainResultScaled, nrow, numeric(1))
        w <- rep(1/w1, w1)
        results <- ensemble_train_glm_with_weight(trainData, GSEATestResultScaled,
                                                  batchCorrection, combineSignatureName,
                                                  weights = w, runModel)
    } else if (training_method == "gbm") {
        message("calling gbm")
        results <- ensemble_train_gbm(trainData, GSEATestResultScaled, 
                                      batchCorrection, combineSignatureName,
                                      runModel)
    } else if (training_method == "nnet") {
        results <- ensemble_train_neural_network(trainData, GSEATestResultScaled, 
                                                 batchCorrection, 
                                                 combineSignatureName,
                                                 runModel)
    } else if (training_method == "svm") {
        message("calling svm")
        results <- ensemble_train_svm(trainData, GSEATestResultScaled, 
                                      batchCorrection, combineSignatureName,
                                      runModel)
    }
    return(results)
}

get_ensemble_result_cv <- function(theObjectList, 
                                   signatureNames, combineSignatureName, 
                                   training_method, 
                                   cross_study = FALSE, times = 1, size = 0, 
                                   leave_one_cv = FALSE,
                                   intercept = TRUE, 
                                   scale_input = FALSE, 
                                   norm_input = FALSE,
                                   batchCorrection = FALSE,
                                   runModel = FALSE) {
    if (!leave_one_cv) {
        message("Not leave one test study, use multiple datasets for testing")
        index <- lapply(1:times, function(x) {
            set.seed(x^2)
            sample(length(theObjectList), size = size)
        })
    } else {
        message("Perform leave one test")
        index <- mapply(list, 1:length(theObjectList))
    }
    
    predScoreList <- lapply(1:length(index), function(i) {
        GSEATrainResult <- theObjectList[-index[[i]]] # List of studies
        GSEATestResult <- theObjectList[index[[i]]] # Single studies
        if (cross_study) {
            message("cross_study is TRUE")
            stop("Do not have cross study feature")
        } else {
            results <- multiStudyEnsemble_single(
                GSEATrainResult, GSEATestResult, 
                signatureNames = signatureNames, 
                combineSignatureName = combineSignatureName,
                training_method = training_method,
                intercept = intercept,
                batchCorrection = batchCorrection,
                scale_input = scale_input,
                norm_input = norm_input,
                runModel = runModel)
        }
        results
    })
    return(predScoreList)
}

#### Analysis ####

# Remove Hoang_OD_3, GSE6112, and study with sample size < filter_size for PTB vs. LTBI 

sigatures_names <- c("Maertzdorf_4", "Maertzdorf_15", "LauxdaCosta_OD_3",
                     "Verhagen_10", "Jacobsen_3", "Sambarey_HIV_10",
                     "Leong_24", "Berry_OD_86", "Berry_393", "Bloom_OD_144",
                     "Suliman_RISK_4", "Zak_RISK_16", "Leong_RISK_29",
                     "Anderson_42", "Anderson_OD_51", "Kaforou_27",
                     "Kaforou_OD_44", "Kaforou_OD_53", "Sweeney_OD_3",
                     "Blankley_5", "Singhania_OD_20", "Thompson_9", "Esmail_82",
                     "Thompson_FAIL_13", "Tornheim_71", "Darboe_RISK_11", 
                     "Hoang_OD_20", "Roe_3", "Gong_OD_4",
                     "Gjoen_10", "Duffy_23", "Jenum_8", "Qian_OD_17", "Estevez_133")
ensemble_sig_list <- make_ensembleSignaturesList(sigatures_names, n = 30, size = 5)

# import score list
command_args <- commandArgs(trailingOnly=TRUE)
infile_path <- command_args[1]
index <- command_args[2] %>% 
  as.numeric()
outfile_path <- command_args[3]

gsea_list <- readRDS(infile_path)


merge_result <- get_ensemble_result_cv(gsea_list, 
                                       ensemble_sig_list[[index]],
                                       combineSignatureName = "combine", 
                                       training_method = "LASSO", batchCorrection = FALSE,
                                       leave_one_cv = FALSE, scale_input = TRUE, 
                                       times = 300, size = floor(0.2 * length(gsea_list)),
                                       norm_input = FALSE, runModel = TRUE)
if (!dir.exists(outfile_path)) {
  dir.create(outfile_path)
}
saveRDS(merge_result, 
        file.path(outfile_path, paste0(names(ensemble_sig_list)[index], ".RDS")))

