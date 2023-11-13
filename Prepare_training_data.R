library(curatedTBData)
library(dplyr)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TBSignatureProfiler)
# setwd("~/Desktop/practice/EnsembleSignatureAnalyze/scripts/")
# set working directory to your own path
source("~/Desktop/practice/curatedTBData_paper_results/Scripts/combine_objects.R")
#### Load data from curatedTBData ####
objects_list <- curatedTBData(study_name = "", dry.run = FALSE)

#### Edit samples from studies ####
# Subset studies remove repeated measurement
# Exclude GSE31348 (only PTB) GSE36238 (only PTB)
# Merge "GSE19435", "GSE19439", "GSE19442", "GSE19444", "GSE22098" to GSE19491
# Merge "GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832" to GSE42834

# Remove progression datasets:
geo_progression <- c("GSE112104", "GSE79362", "GSE94438", "GSE107993", 
                     "GSE107994", "GSETornheim")
# Remove studies with only TB
geo_TB_only <- c("GSE31348", "GSE36238")

# Remove patients with HIV 
geo_HIV <- c("GSE69581", "GSE83892", "GSE50834", "GSE107104")

objects_list_sub <- objects_list[-which(names(objects_list) %in% 
                                           c(geo_progression, geo_HIV, geo_TB_only))]
# length(objects_list_sub) = 37
# Expect total = 37 - 5 + 1 - 6 + 1 = 28
# Get datasets contains HIV positive subjects
geo_need_change_HIV <- lapply(objects_list_sub, function(x) {
    col_info <- colData(x) |> 
        as.data.frame() |> 
        dplyr::filter(HIVStatus == "Positive")
    if (nrow(col_info) == 0) {
        return(NULL)
    }
    return(x)
}) |> 
    plyr::compact() |> 
    names()
geo_need_change_diabetes <- lapply(objects_list_sub, function(x) {
    col_info <- colData(x) |> 
        as.data.frame()
    index <- which(colnames(col_info) == "DiabetesStatus")
    
    if (length(index) == 0) {
        return(NULL)
    }
    return(x)
}) |> 
    plyr::compact() |> 
    names()

GSE42834_geo <- c("GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", 
                  "GSE42832")
GSE19491_geo <- c("GSE19435", "GSE19439", "GSE19442", "GSE19444", "GSE22098")

geo_need_change <- c(geo_need_change_HIV, geo_need_change_diabetes, 
                         GSE42834_geo, GSE19491_geo) |> 
    unique()
geo_nochange_ref <- names(objects_list_sub)[-match(geo_need_change, 
                                                   names(objects_list_sub))]
# Remove 
##### Study without change #####
geo_nochange <- c("GSE101705", "GSE107731", "GSE107991", "GSE107992","GSE25534", 
                  "GSE28623", "GSE29536", "GSE34608", "GSE41055", "GSE6112", 
                  "GSE62525", "GSE74092")

object_edit <- lapply(geo_nochange, function(x) {
    Mobject <- objects_list_sub[[x]]
    dat <- Mobject[["assay_curated"]]
    if (is.null(dat)) {
        stop("check your assay_curated")
    } 
    SummarizedExperiment(assays = list(assay_curated = dat), 
                         colData = colData(Mobject))
})
names(object_edit) <- geo_nochange

##### Combine GSE19491 #####
GSE19491_geo <- c("GSE19435", "GSE19439", "GSE19442", "GSE19444", "GSE22098")

GSE19491_combine <- combine_objects(objects_list_sub[GSE19491_geo],
                                    experiment_name = "assay_curated")

GSE19435_baseline <- objects_list_sub$GSE19435[, objects_list_sub$GSE19435$MeasurementTime
                                           == "0_months"]
GSE19491_edit <- GSE19491_combine[, c(colnames(GSE19435_baseline)[["assay_curated"]],
                                      colnames(objects_list_sub$GSE19439)[["assay_curated"]],
                                      colnames(objects_list_sub$GSE19442)[["assay_curated"]],
                                      colnames(objects_list_sub$GSE19444)[["assay_curated"]],
                                      colnames(objects_list_sub$GSE22098)[["assay_curated"]])]
object_edit$GSE19491_edit <- GSE19491_edit

##### Combine GSE42834 #####
GSE42834_geo <- c("GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832")
GSE42834_combine <- combine_objects(objects_list_sub[GSE42834_geo],
                                    experiment_name = "assay_curated")
GSE42832_sobject_WB <- objects_list_sub$GSE42832[, objects_list_sub$GSE42832$Tissue ==
                                                 "Whole Blood"]
GSE42834_edit <- GSE42834_combine[,c(colnames(objects_list_sub$GSE42825)[["assay_curated"]],
                                     colnames(objects_list_sub$GSE42826)[["assay_curated"]],
                                     colnames(objects_list_sub$GSE42827)[["assay_curated"]],
                                     colnames(objects_list_sub$GSE42830)[["assay_curated"]],
                                     colnames(objects_list_sub$GSE42831)[["assay_curated"]],
                                     colnames(GSE42832_sobject_WB)[["assay_curated"]])]
object_edit$GSE42834_edit <- GSE42834_edit

##### Subset GSE19443 #####
# Subset samples with cell type: Neutrophils
object_edit$GSE19443_edit <- objects_list_sub$GSE19443[, objects_list_sub$GSE19443$Tissue == "Neutrophils"]

##### Subset GSE37250 #####
object_edit$GSE37250_edit <- objects_list_sub$GSE37250[, objects_list_sub$GSE37250$HIVStatus == "Negative"]

##### Subset GSE39939 #####
object_edit$GSE39939_edit <- objects_list_sub$GSE39939[, objects_list_sub$GSE39939$HIVStatus == "Negative"]

##### Subset GSE39940 #####
object_edit$GSE39940_edit <- objects_list_sub$GSE39940[, objects_list_sub$GSE39940$HIVStatus == "Negative"]

##### Subset GSE40553 #####
# Subset samples collected at baseline, or when it is measured at the first time
sample_baseline_GSE40553 <- colData(objects_list_sub$GSE40553)[, c("PatientID", "MeasurementTime")] |> 
    data.frame() |> 
    dplyr::mutate(sample_name = row.names(colData(objects_list_sub$GSE40553))) |> 
    dplyr::arrange(MeasurementTime, PatientID) |> 
    dplyr::group_by(PatientID) |> 
    dplyr::mutate(first = dplyr::first(sample_name))
GSE40553_baseline <- objects_list_sub$GSE40553[, unique(sample_baseline_GSE40553$first)]
object_edit$GSE40553_edit <- SummarizedExperiment(
    list(counts = GSE40553_baseline[["assay_curated"]]),
    colData = colData(GSE40553_baseline))

##### Subset GSE56153 #####
# Get PTB measurement at Baseline and Controls
GSE56153_baseline <- objects_list_sub$GSE56153[, objects_list_sub$GSE56153$MeasurementTime
                                           %in% c("0 weeks", NA)]
object_edit$GSE56153_edit <- SummarizedExperiment(
    list(counts = GSE56153_baseline[["assay_curated"]]),
    colData = colData(GSE56153_baseline))

##### Subset GSE54992 #####
# Only include samples at baseline
GSE54992_baseline <- objects_list_sub$GSE54992[,objects_list_sub$GSE54992$MeasurementTime
                                           == "Baseline"]
object_edit$GSE54992_edit <- SummarizedExperiment(
    list(counts = GSE54992_baseline[["assay_curated"]]),
    colData = colData(GSE54992_baseline))

##### Subset GSE62147 #####
# Only include samples prior to treatment
GSE62147_pre_treatment <- objects_list_sub$GSE62147[, objects_list_sub$GSE62147$MeasurementTime
                                                == "recruit"]
object_edit$GSE62147_edit <- SummarizedExperiment(
    list(counts = GSE62147_pre_treatment[["assay_curated"]]),
    colData = colData(GSE62147_pre_treatment))


##### Subset GSE73408 #####
object_edit$GSE73408_edit <- objects_list_sub$GSE73408[, objects_list_sub$GSE73408$DiabetesStatus
                                                       == "Negative"]
##### Subset GSE81746 #####
# Exclude the pooled PTB (Male)
GSE81746_sub <- objects_list_sub$GSE81746[, objects_list_sub$GSE81746$Gender %in% "Male"]
object_edit$GSE81746_edit <- SummarizedExperiment(
    list(counts = GSE81746_sub[["assay_curated"]]),
    colData = colData(GSE81746_sub))

##### Subset GSE83456 #####
# Exclude patients with EPTB 
GSE83456_sub <- objects_list_sub$GSE83456[, objects_list_sub$GSE83456$TBStatus != "EPTB"]
object_edit$GSE83456_edit <- SummarizedExperiment(
    list(counts = GSE83456_sub[["assay_curated"]]),
    colData = colData(GSE83456_sub))

##### Subset GSE84076 #####
# Take BCG vaccinated controls and LTBIs
# Take PTB before treatment results
GSE84076_BCG <- objects_list_sub$GSE84076[, objects_list_sub$GSE84076$BcgVaccinated %in% "Yes"]
GSE84076_beforeTreat <- objects_list_sub$GSE84076[, objects_list_sub$GSE84076$TreatmentStatus == "Treatment-naive"]
GSE84076_sub1 <- objects_list_sub$GSE84076[, c(colnames(GSE84076_BCG)[["assay_curated"]],
                                          colnames(GSE84076_beforeTreat)[["assay_curated"]])]
object_edit$GSE84076_edit <- SummarizedExperiment(list(counts = GSE84076_sub1[["assay_curated"]]),
                                                  colData = colData(GSE84076_sub1))
##### Subset GSE89403 #####
# Only include samples collected at baseline, or when it is measured at first time
# Remove samples with TBStatus NA
sample_baseline_GSE89403 <- colData(objects_list_sub$GSE89403)[, c("PatientID","MeasurementTime")] |> 
    data.frame() |> 
    dplyr::mutate(sample_name = row.names(colData(objects_list_sub$GSE89403))) |> 
    dplyr::arrange(MeasurementTime, PatientID) |> 
    dplyr::group_by(PatientID) |> 
    dplyr::mutate(first = dplyr::first(sample_name))
GS89403_baseline <- objects_list_sub$GSE89403[, unique(sample_baseline_GSE89403$first)]
GS89403_baseline_noNA <- GS89403_baseline[, !is.na(GS89403_baseline$TBStatus)]
object_edit$GSE89403_edit <- SummarizedExperiment(
    list(counts = GS89403_baseline_noNA[["assay_curated"]]),
    colData = colData(GS89403_baseline_noNA))

##### Subset GSEBruno #####
object_edit$GSEBruno_edit <- objects_list_sub$GSEBruno[, objects_list_sub$GSEBruno$DiabetesStatus == "Negative"]

##### Save training data #####
saveRDS(object_edit,
        "~/Desktop/practice/curatedTBData_paper_results/data/object_edit.RDS")

