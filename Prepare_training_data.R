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
source("combine_objects.R")
#### Load data from curatedTBData ####
objects_list <- curatedTBData(study_name = "", dry.run = FALSE)

#### Edit samples from studies ####
# Subset studies remove repeated measurement
# Exclude GSE31348 (only PTB) GSE3628 (only PTB) GSE74092(RT-PCR)
# Collapse "GSE19435", "GSE19439", "GSE19442", "GSE19444", "GSE22098" to GSE19491
# Collapse "GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832" to GSE42834
# Expect 49 - 2 - 5 + 1 - 6 + 1 = 38

##### Study without change #####
geo_nochange <- c("GSE28623", "GSE29536", "GSE34608", "GSE37250", "GSE39939",
                  "GSE39940", "GSE41055", "GSE50834", "GSE62525", "GSE73408",
                  "GSE83892", "GSE101705", "GSE107731", "GSE107991", "GSE107992",
                  "GSE107104", "GSE112104", "GSEBruno", "GSE25534", "GSE6112", "GSE74092")
object_edit <- lapply(geo_nochange, function(x) {
    Mobject <- objects_list[[x]]
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

GSE19491_combine <- combine_objects(objects_list[GSE19491_geo],
                                    experiment_name = "assay_curated")

GSE19435_baseline <- objects_list$GSE19435[, objects_list$GSE19435$MeasurementTime
                                           == "0_months"]
GSE19491_edit <- GSE19491_combine[, c(colnames(GSE19435_baseline)[["assay_curated"]],
                                      colnames(objects_list$GSE19439)[["assay_curated"]],
                                      colnames(objects_list$GSE19442)[["assay_curated"]],
                                      colnames(objects_list$GSE19444)[["assay_curated"]],
                                      colnames(objects_list$GSE22098)[["assay_curated"]])]
object_edit$GSE19491_edit <- GSE19491_edit
##### Combine GSE42834 #####
GSE42834_geo <- c("GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832")
GSE42834_combine <- combine_objects(objects_list[GSE42834_geo],
                                    experiment_name = "assay_curated")
GSE42832_sobject_WB <- objects_list$GSE42832[, objects_list$GSE42832$Tissue ==
                                                 "Whole Blood"]
GSE42834_edit <- GSE42834_combine[,c(colnames(objects_list$GSE42825)[["assay_curated"]],
                                     colnames(objects_list$GSE42826)[["assay_curated"]],
                                     colnames(objects_list$GSE42827)[["assay_curated"]],
                                     colnames(objects_list$GSE42830)[["assay_curated"]],
                                     colnames(objects_list$GSE42831)[["assay_curated"]],
                                     colnames(GSE42832_sobject_WB)[["assay_curated"]])]
object_edit$GSE42834_edit <- GSE42834_edit

##### Subset GSE19443 #####
# Subset samples with cell type: Neutrophils
object_edit$GSE19443_edit <- objects_list$GSE19443[, objects_list$GSE19443$Tissue == "Neutrophils"]

##### Subset GSE56153 #####
# Get PTB measurement at Baseline and Controls
GSE56153_baseline <- objects_list$GSE56153[, objects_list$GSE56153$MeasurementTime
                                           %in% c("0 weeks", NA)]
object_edit$GSE56153_edit <- SummarizedExperiment(
    list(counts = GSE56153_baseline[["assay_curated"]]),
    colData = colData(GSE56153_baseline))

##### Subset GSE54992 #####
# Only include samples at baseline
GSE54992_baseline <- objects_list$GSE54992[,objects_list$GSE54992$MeasurementTime
                                           == "Baseline"]
object_edit$GSE54992_edit <- SummarizedExperiment(
    list(counts = GSE54992_baseline[["assay_curated"]]),
    colData = colData(GSE54992_baseline))

##### Subset GSE62147 #####
# Only include samples prior to treatment
GSE62147_pre_treatment <- objects_list$GSE62147[, objects_list$GSE62147$MeasurementTime
                                                == "recruit"]
object_edit$GSE62147_edit <- SummarizedExperiment(
    list(counts = GSE62147_pre_treatment[["assay_curated"]]),
    colData = colData(GSE62147_pre_treatment))

##### Subset GSE69581 #####
# Exclude 10 subclinical samples
GSE69581_PTB_Latent <- objects_list$GSE69581[, objects_list$GSE69581$TBStatus
                                             %in% c("PTB", "LTBI")]
object_edit$GSE69581_edit <- SummarizedExperiment(
    list(counts = GSE69581_PTB_Latent[["assay_curated"]]),
    colData = colData(GSE69581_PTB_Latent))

##### Subset GSE79362 ####
# Use the reprocessed RNA-seq, only inlude progressors and non-progressors at baseline
counts.africa.baseline <- objects_list$GSE79362[["assay_reprocess_hg19"]]
# Max 5 filter
MaxFilter <- function(df, max.value = 10){
    df.filtered <- df[which(apply(df,1,max) >= max.value),]
    return(df.filtered)
}
counts.africa.baseline.filtered <- MaxFilter(counts.africa.baseline, 5)
# Normalization
counts.africa.baseline.norm <- TBSignatureProfiler::deseq2_norm_rle(counts.africa.baseline.filtered)
GSE79362_train_full <- SummarizedExperiment(list(counts=counts.africa.baseline.norm),
                                            colData = colData(objects_list$GSE79362))
sample_baseline_GSE79362 <- colData(GSE79362_train_full)[, c("PatientID","MeasurementTime")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(GSE79362_train_full))) %>%
    dplyr::arrange(MeasurementTime, PatientID) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))
GSE79362_baseline <- GSE79362_train_full[,unique(sample_baseline_GSE79362$first)]
object_edit$GSE79362_edit <- GSE79362_baseline ## validated with the Brazil data

##### Subset GSE81746 #####
# Exclude the pooled PTB (Male)
GSE81746_sub <- objects_list$GSE81746[,objects_list$GSE81746$Gender %in% "Male"]
object_edit$GSE81746_edit <- SummarizedExperiment(
    list(counts = GSE81746_sub[["assay_curated"]]),
    colData = colData(GSE81746_sub))

##### Subset GSE83456 #####
# Exclude patients with EPTB 
GSE83456_sub <- objects_list$GSE83456[,objects_list$GSE83456$TBStatus != "EPTB"]
object_edit$GSE83456_edit <- SummarizedExperiment(
    list(counts = GSE83456_sub[["assay_curated"]]),
    colData = colData(GSE83456_sub))

##### Subset GSE84076 #####
# Take BCG vaccinated controls and LTBIs
# Take PTB before treatment results
GSE84076_BCG <- objects_list$GSE84076[, objects_list$GSE84076$BcgVaccinated %in% "Yes"]
GSE84076_beforeTreat <- objects_list$GSE84076[, objects_list$GSE84076$TreatmentStatus == "Treatment-naive"]
GSE84076_sub1 <- objects_list$GSE84076[,c(colnames(GSE84076_BCG)[["assay_curated"]],
                                          colnames(GSE84076_beforeTreat)[["assay_curated"]])]
object_edit$GSE84076_edit <- SummarizedExperiment(list(counts = GSE84076_sub1[["assay_curated"]]),
                                                  colData = colData(GSE84076_sub1))
##### Subset GSE107994 #####
# Use the reprocessed RNA-seq, only include samples at baseline
counts.gse107994.baseline <- objects_list$GSE107994[["assay_reprocess_hg38"]]
# Max 5 filter
counts.gse107994.baseline.filtered <- MaxFilter(counts.gse107994.baseline, 5)
# Normalization
counts.gse107994.baseline.norm <- TBSignatureProfiler::deseq2_norm_rle(counts.gse107994.baseline.filtered)

GSE107994_test_full <- SummarizedExperiment(list(counts=counts.gse107994.baseline.norm),
                                            colData = colData(objects_list$GSE107994))
# index <- which(is.na(GSE107994_test_full$PatientID))
sample_baseline_GSE107994 <- colData(GSE107994_test_full)[, c("PatientID","Progression")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(GSE107994_test_full))) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))

# Patient_087 does not have baseline measurement
GSE107994_baseline <- GSE107994_test_full[, unique(sample_baseline_GSE107994$first)]

object_edit$GSE107994_edit <- GSE107994_baseline

##### Subset GSE94438 #####
# Use reprocessed RNA-seq, exclude samples with TBStatus NA
counts.gse94438.baseline <- objects_list$GSE94438[["assay_reprocess_hg19"]]
# Max 5 filter
counts.gse94438.baseline.filtered <- MaxFilter(counts.gse94438.baseline, 5)
# Normalization
counts.gse94438.baseline.norm <- TBSignatureProfiler::deseq2_norm_rle(counts.gse94438.baseline.filtered)

GSE94438_test_full <- SummarizedExperiment(list(counts=counts.gse94438.baseline.norm),
                                           colData = colData(objects_list$GSE94438))
attributes(row.names(GSE94438_test_full)) <- NULL
GSE94438_test_full_NoNA <- GSE94438_test_full[,GSE94438_test_full$Progression %in% c("Positive","Negative")]
object_edit$GSE94438_edit <- GSE94438_test_full_NoNA

##### Subset GSE89403 #####
# Only include samples collected at baseline, or when it is measured at first time
# Remove samples with TBStatus NA
sample_baseline_GS89403 <- colData(objects_list$GSE89403)[, c("PatientID","MeasurementTime")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(objects_list$GSE89403))) %>%
    dplyr::arrange(MeasurementTime, PatientID) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))
GS89403_baseline <- objects_list$GSE89403[,unique(sample_baseline_GS89403$first)]
GS89403_baseline_noNA <- GS89403_baseline[, !is.na(GS89403_baseline$TBStatus)]
object_edit$GSE89403_edit <- SummarizedExperiment(
    list(counts = GS89403_baseline_noNA[["assay_curated"]]),
    colData = colData(GS89403_baseline_noNA))

##### Subset GSE107993 #####
# Only include samples collected at baseline, or when it is measured at first time
sample_baseline_GSE107993 <- colData(objects_list$GSE107993)[, c("PatientID","MeasurementTime")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(objects_list$GSE107993))) %>%
    dplyr::arrange(MeasurementTime, PatientID) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))
GSE107993_baseline <- objects_list$GSE107993[, unique(sample_baseline_GSE107993$first)]

object_edit$GSE107993_edit <- SummarizedExperiment(
    list(counts = GSE107993_baseline[["assay_curated"]]),
    colData = colData(GSE107993_baseline)[colnames(GSE107993_baseline[["assay_curated"]]),])

##### Subset GSETornheim #####
# Subset samples collected at baseline, or when it is measured at the first time
# Remove patients with EPTB
sample_baseline_GSETornheim <- colData(objects_list$GSETornheim)[, c("PatientID", "MeasurementTime")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(objects_list$GSETornheim))) %>%
    dplyr::arrange(MeasurementTime) %>% 
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))
sample_baseline_GSEToenheim <- objects_list$GSETornheim[, unique(sample_baseline_GSETornheim$first)]
sample_baseline_GSEToenheimSub <- sample_baseline_GSEToenheim[, sample_baseline_GSEToenheim$TBStatus %in% c("Control", "PTB")]
object_edit$GSETornheim_edit <- SummarizedExperiment(
    list(counts = sample_baseline_GSEToenheimSub[["assay_curated"]]),
    colData = colData(sample_baseline_GSEToenheimSub))

##### Subset GSE40553 #####
# Subset samples collected at baseline, or when it is measured at the first time
sample_baseline_GSE40553 <- colData(objects_list$GSE40553)[, c("PatientID","MeasurementTime")] %>%
    data.frame() %>%
    dplyr::mutate(sample_name = row.names(colData(objects_list$GSE40553))) %>%
    dplyr::arrange(MeasurementTime, PatientID) %>%
    dplyr::group_by(PatientID) %>%
    dplyr::mutate(first = dplyr::first(sample_name))
GSE40553_baseline <- objects_list$GSE40553[, unique(sample_baseline_GSE40553$first)]
object_edit$GSE40553_edit <- SummarizedExperiment(
    list(counts = GSE40553_baseline[["assay_curated"]]),
    colData = colData(GSE40553_baseline))

##### Save training data #####
saveRDS(object_edit, 
        "~/Desktop/practice/curatedTBData_paper_results/data/object_edit.RDS")
