# This is script contains code to generate figures and tables in the curatedTBData paper
#### Import functions ####
source("~/Desktop/practice/curatedTBData_paper_results/Scripts/functions.R")
source("~/Desktop/Packages/curatedTBData/data-raw/UtilityFunctionForCuration.R")
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(SummarizedExperiment)
library(TBSignatureProfiler)
wd <- "~/Desktop/practice/curatedTBData_paper_results"

#### Modify datasets for curatedTBData ####
# See ~/Desktop/practice/curatedTBData_paper_results/scripts/prepare_training_data.R for details

#### Analysis ####
object_edit <- readRDS("~/Desktop/practice/curatedTBData_paper_results/data/object_edit.RDS")
object_edit_update <- lapply(object_edit, function(x) {
    dat <- assay(x)
    row_names <- row.names(dat) |> 
        update_genenames()
    row.names(dat) <- row_names
    SummarizedExperiment(assays = list(dat), colData = colData(x))
})

sigatures_names <- c("Anderson_42", "Anderson_OD_51", "Berry_393", "Berry_OD_86", 
                     "Blankley_380", "Blankley_5", "Bloom_OD_144", "Bloom_RES_268", 
                     "Bloom_RES_558", "Chen_HIV_4", "Darboe_RISK_11",
                     "Dawany_HIV_251", "Duffy_23", "Esmail_203", "Esmail_82",
                     "Esmail_OD_893", "Estevez_133", "Estevez_259",
                     "Gjoen_10", "Gjoen_7", "Gliddon_2_OD_4", "Gliddon_HIV_3",
                     "Gliddon_OD_3", "Gliddon_OD_4", "Gong_OD_4", "Heycken_FAIL_22",
                     "Hoang_OD_13", "Hoang_OD_20", "Hoang_OD_3", "Huang_OD_13",
                     "Jacobsen_3", "Jenum_8", "Kaforou_27", "Kaforou_OD_44",
                     "Kaforou_OD_53", "Kaul_3", "Kwan_186",
                     "LauxdaCosta_OD_3", "Lee_4", "Leong_24", "Leong_RISK_29",
                     "Long_RES_10", "Maertzdorf_15", "Maertzdorf_4", "Maertzdorf_OD_100",
                     "Natarajan_7", "PennNich_RISK_6", "Qian_OD_17", "Rajan_HIV_5",
                     "Roe_3", "Roe_OD_4", "Sambarey_HIV_10", "Singhania_OD_20",
                     "Sivakumaran_11", "Suliman_4", "Suliman_RISK_4", "Sweeney_OD_3", 
                     "Tabone_OD_11", "Tabone_RES_25", "Tabone_RES_27", "Thompson_9", 
                     "Thompson_FAIL_13", "Thompson_RES_5", "Tornheim_71", "Tornheim_RES_25", 
                     "Verhagen_10", "Walter_51", "Walter_PNA_119", "Walter_PNA_47", 
                     "Zak_RISK_16", "Zhao_NANO_6", "Zimmer_RES_3")
signatures_list <- TBSignatureProfiler::TBsignatures[sigatures_names] |> 
    lapply(update_genenames)

signatures_list <- signatures_list[order(sigatures_names)]

##### ssGSEA PTB vs. Control #####
# ssgsea_PTB_Control_out <- lapply(object_edit_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                          annotationCondition = c("PTB", "Control"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "ssGSEA", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(ssgsea_PTB_Control_out, file.path(wd, "data/ssgsea_PTB_Control_out.RDS"))

# ssgsea_PTB_Control_out <- readRDS(file.path(wd, "data/ssgsea_PTB_Control_out.RDS"))
# ssgsea_PTB_Control_out_combine <- combine_auc(ssgsea_PTB_Control_out,
#                                               annotationColName = "TBStatus",
#                                               signatureColNames = names(signatures_list),
#                                               num.boot = 1000, percent = 0.95)
# saveRDS(ssgsea_PTB_Control_out_combine,
#         file.path(wd, "data/ssgsea_PTB_Control_out_combine.RDS"))

ssgsea_PTB_Control_out <- readRDS(file.path(wd, "data/ssgsea_PTB_Control_out.RDS"))

study_PTB_Control <- lapply(1:length(ssgsea_PTB_Control_out), function(i) {
    sobject <- ssgsea_PTB_Control_out[[i]]
    sutdy_name <- gsub("_edit", "", names(ssgsea_PTB_Control_out)[i])
    data.frame(Study = sutdy_name, Size = ncol(sobject))
}) |> 
    dplyr::bind_rows()

ssgsea_PTB_Control_out_combine <- file.path(wd, "data/ssgsea_PTB_Control_out_combine.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_Control, "Study")

##### PLAGE PTB vs. Control #####
# plage_PTB_Control_out <- lapply(object_edit_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                                 annotationCondition = c("PTB", "Control"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "PLAGE", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(plage_PTB_Control_out, file.path(wd, "data/plage_PTB_Control_out.RDS"))

# plage_PTB_Control_out_combine <- combine_auc(plage_PTB_Control_out,
#                                              annotationColName = "TBStatus",
#                                              signatureColNames = names(signatures_list),
#                                              num.boot = 1000, percent = 0.95)
# saveRDS(plage_PTB_Control_out_combine, file.path(wd, "data/plage_PTB_Control_out_combine.RDS"))

plage_PTB_Control_out <- readRDS(file.path(wd, "data/plage_PTB_Control_out.RDS"))

plage_PTB_Control_out_combine <- file.path(wd, "data/plage_PTB_Control_out_combine.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_Control, "Study")

##### ssGSEA PTB vs. LTBI #####

# ssgsea_PTB_LTBI_out <- lapply(object_edit_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                                 annotationCondition = c("PTB", "LTBI"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "ssGSEA", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(ssgsea_PTB_LTBI_out, file.path(wd, "data/ssgsea_PTB_LTBI_out.RDS"))

# ssgsea_PTB_LTBI_out_combine <- combine_auc(ssgsea_PTB_LTBI_out,
#                                               annotationColName = "TBStatus",
#                                               signatureColNames = names(signatures_list),
#                                               num.boot = 1000, percent = 0.95)
# saveRDS(ssgsea_PTB_LTBI_out_combine, file.path(wd, "data/ssgsea_PTB_LTBI_out_combine.RDS"))

ssgsea_PTB_LTBI_out <- readRDS(file.path(wd, "data/ssgsea_PTB_LTBI_out.RDS"))

study_PTB_LTBI <- lapply(1:length(ssgsea_PTB_LTBI_out), function(i) {
    sobject <- ssgsea_PTB_LTBI_out[[i]]
    sutdy_name <- gsub("_edit", "", names(ssgsea_PTB_LTBI_out)[i])
    data.frame(Study = sutdy_name, Size = ncol(sobject))
}) |> 
    dplyr::bind_rows()

ssgsea_PTB_LTBI_out_combine <- file.path(wd, "data/ssgsea_PTB_LTBI_out_combine.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_LTBI, "Study")

##### PLAGE PTB vs. LTBI #####
# plage_PTB_LTBI_out <- lapply(object_edit_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                                 annotationCondition = c("PTB", "LTBI"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "PLAGE", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(plage_PTB_LTBI_out, file.path(wd, "data/plage_PTB_LTBI_out.RDS"))

# plage_PTB_LTBI_out_combine <- combine_auc(plage_PTB_LTBI_out,
#                                            annotationColName = "TBStatus",
#                                            signatureColNames = names(signatures_list),
#                                            num.boot = 1000, percent = 0.95)
# saveRDS(plage_PTB_LTBI_out_combine, file.path(wd, "data/plage_PTB_LTBI_out_combine.RDS"))

plage_PTB_LTBI_out <- readRDS(file.path(wd, "data/plage_PTB_LTBI_out.RDS"))

plage_PTB_LTBI_out_combine <- file.path(wd, "data/plage_PTB_LTBI_out_combine.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_LTBI, "Study")

#### Table S1 Datasets in the curatedTBData ####
data("DataSummary")
study_summary_sub <- DataSummary |> 
    dplyr::select(-Notes, -GeneralType)
# write.table(study_summary_sub, quote = FALSE, row.names = FALSE, sep = " & ",
#             file = "~/Desktop/practice/curatedTBData_paper_results/Figures_and_tables/study_summary_sub.txt") 
                                   
# Run the following command in terminal (Add '\\\hline' by the end of each line)
# awk '{print $0, " \\\\\\\hline"}' study_summary_sub.txt > study_summary_sub_for_latex.txt

#### Table S2 Clinical annotation ####
clincial_annotation <- readxl::read_excel("Desktop/practice/curatedTBData_paper_results/clinicalDataAnnotation.xlsx")
# write.table(clincial_annotation, 
#             quote = FALSE, row.names = FALSE, sep = " & ",
#             file = file.path(out_file_path, "clincial_annotation.txt")) 
# Run the following command in terminal (Add '\\\hline' by the end of each line)
# awk '{print $0, " \\\\\\\hline"}' clincial_annotation.txt > clincial_annotation_for_latex.txt
#### Extract AUC, sens, spec####
ssgsea_PTB_Control_high <- ssgsea_PTB_Control_out_combine |> 
    dplyr::filter(Signature %in% c("Blankley_5","Kaul_3", 
                                   "Tabone_RES_27"))
ssgsea_PTB_Control_high |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(mean_AUC = weighted.mean(AUC, Size), 
                     AUC_low = weighted.mean(`AUC CI lower.2.5%`, Size),
                     AUC_high = weighted.mean(`AUC CI upper.97.5%`, Size),
                     mean_sens = weighted.mean(Sensitivity, Size),
                     sens_low = weighted.mean(`Sens CI lower.2.5%`, Size),
                     sens_high = weighted.mean(`Sens CI upper.97.5%`, Size),
                     mean_spec = weighted.mean(Specificity, Size),
                     spec_low = weighted.mean(`Spec CI lower.2.5%`, Size),
                     spec_high = weighted.mean(`Spec CI upper.97.5%`, Size))

plage_PTB_Control_high <- plage_PTB_Control_out_combine |> 
    dplyr::filter(Signature %in% c("Maertzdorf_15","Suliman_4", 
                                   "Tabone_OD_11"))
plage_PTB_Control_high |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(mean_AUC = weighted.mean(AUC, Size), 
                     AUC_low = weighted.mean(`AUC CI lower.2.5%`, Size),
                     AUC_high = weighted.mean(`AUC CI upper.97.5%`, Size),
                     mean_sens = weighted.mean(Sensitivity, Size),
                     sens_low = weighted.mean(`Sens CI lower.2.5%`, Size),
                     sens_high = weighted.mean(`Sens CI upper.97.5%`, Size),
                     mean_spec = weighted.mean(Specificity, Size),
                     spec_low = weighted.mean(`Spec CI lower.2.5%`, Size),
                     spec_high = weighted.mean(`Spec CI upper.97.5%`, Size))

ssgsea_PTB_LTBI_high <- ssgsea_PTB_LTBI_out_combine |> 
    dplyr::filter(Signature %in% c("Bloom_RES_268","Estevez_133", 
                                   "Tabone_RES_27"))
ssgsea_PTB_LTBI_high |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(mean_AUC = weighted.mean(AUC, Size), 
                     AUC_low = weighted.mean(`AUC CI lower.2.5%`, Size),
                     AUC_high = weighted.mean(`AUC CI upper.97.5%`, Size),
                     mean_sens = weighted.mean(Sensitivity, Size),
                     sens_low = weighted.mean(`Sens CI lower.2.5%`, Size),
                     sens_high = weighted.mean(`Sens CI upper.97.5%`, Size),
                     mean_spec = weighted.mean(Specificity, Size),
                     spec_low = weighted.mean(`Spec CI lower.2.5%`, Size),
                     spec_high = weighted.mean(`Spec CI upper.97.5%`, Size))

plage_PTB_LTBI_high <- plage_PTB_LTBI_out_combine |> 
    dplyr::filter(Signature %in% c("Kaforou_27"))
plage_PTB_LTBI_high |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(mean_AUC = weighted.mean(AUC, Size), 
                     AUC_low = weighted.mean(`AUC CI lower.2.5%`, Size),
                     AUC_high = weighted.mean(`AUC CI upper.97.5%`, Size),
                     mean_sens = weighted.mean(Sensitivity, Size),
                     sens_low = weighted.mean(`Sens CI lower.2.5%`, Size),
                     sens_high = weighted.mean(`Sens CI upper.97.5%`, Size),
                     mean_spec = weighted.mean(Specificity, Size),
                     spec_low = weighted.mean(`Spec CI lower.2.5%`, Size),
                     spec_high = weighted.mean(`Spec CI upper.97.5%`, Size))

#### Table: Weighted AUCs (95% CI) for PTB vs. Control and PTB vs. LTBI ####
ssgsea_PTB_Control_CI <- extract_CI(ssgsea_PTB_Control_out_combine)
ssgsea_PTB_Control_CI_sub <- ssgsea_PTB_Control_CI |> 
    dplyr::mutate(ssgsea_PTB_Control = AUC) |> 
    dplyr::select(Signature, ssgsea_PTB_Control)

ssgsea_PTB_LTBI_CI <- extract_CI(ssgsea_PTB_LTBI_out_combine)
ssgsea_PTB_LTBI_CI_sub <- ssgsea_PTB_LTBI_CI |> 
    dplyr::mutate(ssgsea_PTB_LTBI = AUC) |> 
    dplyr::select(Signature, ssgsea_PTB_LTBI)

plage_PTB_Control_CI <- extract_CI(plage_PTB_Control_out_combine)
plage_PTB_Control_CI_sub <- plage_PTB_Control_CI |> 
    dplyr::mutate(plage_PTB_Control = AUC) |> 
    dplyr::select(Signature, plage_PTB_Control)

plage_PTB_LTBI_CI <- extract_CI(plage_PTB_LTBI_out_combine)
plage_PTB_LTBI_CI_sub <- plage_PTB_LTBI_CI |> 
    dplyr::mutate(plage_PTB_LTBI = AUC) |> 
    dplyr::select(Signature, plage_PTB_LTBI)

table_mean_AUC <- ssgsea_PTB_Control_CI_sub |> 
    dplyr::inner_join(plage_PTB_Control_CI_sub, by = "Signature") |> 
    dplyr::inner_join(ssgsea_PTB_LTBI_CI_sub, by = "Signature") |> 
    dplyr::inner_join(plage_PTB_LTBI_CI_sub, by = "Signature")
table_mean_AUC_sub <- table_mean_AUC |> 
    dplyr::filter(!Signature == "Hoang_OD_3")
    
table_mean_AUC_sub$Signature <- gsub("_", "\\\\_", table_mean_AUC_sub$Signature)
write.table(table_mean_AUC_sub, quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = " & ",
            file = file.path(wd, "Figures_and_tables/table_mean_AUC.txt")) 

# Run the following command in terminal (Add '\\\hline' by the end of each line)
# awk '{print $0, " \\\\\\\hline"}' table_mean_AUC.txt > table_mean_AUC_for_latex.txt

ssgsea_PTB_Control_AUCs <- gsub(" \\(.*", "",table_mean_AUC$ssgsea_PTB_Control) |> 
    as.numeric()
plage_PTB_Control_AUCs <- gsub(" \\(.*", "",table_mean_AUC$plage_PTB_Control) |> 
    as.numeric()
ssgsea_PTB_LTBI_AUCs <- gsub(" \\(.*", "",table_mean_AUC$ssgsea_PTB_LTBI) |> 
    as.numeric()
plage_PTB_LTBI_AUCs <- gsub(" \\(.*", "",table_mean_AUC$plage_PTB_LTBI) |> 
    as.numeric()

wilcox.test(ssgsea_PTB_Control_AUCs, ssgsea_PTB_LTBI_AUCs, paired = TRUE)
wilcox.test(plage_PTB_Control_AUCs, plage_PTB_LTBI_AUCs, paired = TRUE)
wilcox.test(ssgsea_PTB_Control_AUCs, plage_PTB_Control_AUCs, paired = TRUE)
wilcox.test(ssgsea_PTB_LTBI_AUCs, plage_PTB_LTBI_AUCs, paired = TRUE)

#### New analysis 20250722 in response to editors ####
library(metafor)
get_meta_random_effect <- function(df_out_combine, sigs, col_name) {
    df_out <- lapply(sigs, function(sig) {
        df_random <- df_out_combine |> 
            dplyr::filter(Signature %in% sig)
        df_random$SE <- (df_random$`AUC CI upper.97.5%` - df_random$`AUC CI lower.2.5%`) / (2 * 1.96)
        # df_random <- df_random |> dplyr::filter(!SE==0)
        res <- rma(yi = AUC, sei = df_random[,col_name], data = df_random, method = "REML")
        data.frame(Signature = sig, AUC = res$beta[,1], 
                   lower = res$ci.lb, upper = res$ci.ub) |> 
            dplyr::mutate(random_effect_AUC = sprintf("%.2f (%.2f - %.2f)", 
                                  AUC, lower, upper))
    }) |> 
        dplyr::bind_rows()
    return(df_out)
}

# PTB vs. Control
ssgsea_PTB_Control_CI_hi <- ssgsea_PTB_Control_CI |> 
    dplyr::slice(grep("^0.9", ssgsea_PTB_Control_CI$AUC)) 

plage_PTB_Control_CI_hi <- plage_PTB_Control_CI |> 
    dplyr::slice(grep("^0.9", plage_PTB_Control_CI$AUC)) 

PTB_Control_sig_hi <- intersect(ssgsea_PTB_Control_CI_hi$Signature, 
                               plage_PTB_Control_CI_hi$Signature)
PTB_Control_sig_hi <- PTB_Control_sig_hi[PTB_Control_sig_hi!='Hoang_OD_3']
# Results from bootstrap SE is similar using the formula in get_meta_random_effect function
# ssgsea_PTB_Control_SE_combine <- combine_auc(ssgsea_PTB_Control_out,
#                                               annotationColName = "TBStatus",
#                                               signatureColNames = PTB_Control_sig_hi,
#                                               num.boot = 1000, percent = 0.95)

ssgsea_PTB_Control_RE <- get_meta_random_effect(ssgsea_PTB_Control_out_combine, 
                                                PTB_Control_sig_hi, col_name = 'SE')
ssgsea_PTB_Control_compare <- ssgsea_PTB_Control_CI |> 
    dplyr::filter(Signature %in% PTB_Control_sig_hi) |> 
    dplyr::mutate(weighted_mean_AUC_ssgsea = AUC) |> 
    dplyr::select(Signature, weighted_mean_AUC_ssgsea) |> 
    dplyr::inner_join(ssgsea_PTB_Control_RE |> 
                          dplyr::select(Signature, random_effect_AUC)) |> 
    dplyr::rename(random_effect_AUC_ssgsea = random_effect_AUC)

plage_PTB_Control_RE <- get_meta_random_effect(plage_PTB_Control_out_combine, 
                                                PTB_Control_sig_hi, col_name = 'SE')

plage_PTB_Control_compare <- plage_PTB_Control_CI |> 
    dplyr::filter(Signature %in% PTB_Control_sig_hi) |> 
    dplyr::mutate(weighted_mean_AUC_plage = AUC) |> 
    dplyr::select(Signature, weighted_mean_AUC_plage) |> 
    dplyr::inner_join(plage_PTB_Control_RE |> 
                          dplyr::select(Signature, random_effect_AUC)) |> 
    dplyr::rename(random_effect_AUC_plage = random_effect_AUC)

PTB_Control_compare <- ssgsea_PTB_Control_compare |> 
    dplyr::inner_join(plage_PTB_Control_compare)

# PTB vs. LTBI
ssgsea_PTB_LTBI_CI_hi <- ssgsea_PTB_LTBI_CI |> 
    dplyr::slice(grep("^0.9", ssgsea_PTB_LTBI_CI$AUC))

plage_PTB_LTBI_CI_hi <- plage_PTB_LTBI_CI |> 
    dplyr::slice(grep("^0.9", plage_PTB_LTBI_CI$AUC))

PTB_LTBI_sig_hi <- intersect(ssgsea_PTB_LTBI_CI_hi$Signature, 
                             plage_PTB_LTBI_CI_hi$Signature)
ssgsea_PTB_LTBI_RE <- get_meta_random_effect(ssgsea_PTB_LTBI_out_combine, 
                                                PTB_LTBI_sig_hi, col_name = 'SE')
ssgsea_PTB_LTBI_compare <- ssgsea_PTB_LTBI_CI |> 
    dplyr::filter(Signature %in% PTB_LTBI_sig_hi) |> 
    dplyr::mutate(weighted_mean_AUC_ssgsea = AUC) |> 
    dplyr::select(Signature, weighted_mean_AUC_ssgsea) |> 
    dplyr::inner_join(ssgsea_PTB_LTBI_RE |> 
                          dplyr::select(Signature, random_effect_AUC)) |> 
    dplyr::rename(random_effect_AUC_ssgsea = random_effect_AUC)

plage_PTB_LTBI_RE <- get_meta_random_effect(plage_PTB_LTBI_out_combine, 
                                             PTB_LTBI_sig_hi, col_name = 'SE')

plage_PTB_LTBI_compare <- plage_PTB_LTBI_CI |> 
    dplyr::filter(Signature %in% PTB_LTBI_sig_hi) |> 
    dplyr::mutate(weighted_mean_AUC_plage = AUC) |> 
    dplyr::select(Signature, weighted_mean_AUC_plage) |> 
    dplyr::inner_join(plage_PTB_LTBI_RE |> 
                          dplyr::select(Signature, random_effect_AUC)) |> 
    dplyr::rename(random_effect_AUC_plage = random_effect_AUC)

PTB_LTBI_compare <- ssgsea_PTB_LTBI_compare |> 
    dplyr::inner_join(plage_PTB_LTBI_compare)

#### Table: Weighted Sensitivity (95% CI) for PTB vs. Control and PTB vs. LTBI ####
ssgsea_PTB_Control_CI_sens <- ssgsea_PTB_Control_CI |> 
    dplyr::mutate(ssgsea_PTB_Control = Sensitivity) |> 
    dplyr::select(Signature, ssgsea_PTB_Control)

ssgsea_PTB_LTBI_CI_sens <- ssgsea_PTB_LTBI_CI |> 
    dplyr::mutate(ssgsea_PTB_LTBI = Sensitivity) |> 
    dplyr::select(Signature, ssgsea_PTB_LTBI)

plage_PTB_Control_CI_sens <- plage_PTB_Control_CI |> 
    dplyr::mutate(plage_PTB_Control = Sensitivity) |> 
    dplyr::select(Signature, plage_PTB_Control)

plage_PTB_LTBI_CI_sens <- plage_PTB_LTBI_CI |> 
    dplyr::mutate(plage_PTB_LTBI = Sensitivity) |> 
    dplyr::select(Signature, plage_PTB_LTBI)

table_mean_sens <- ssgsea_PTB_Control_CI_sens |> 
    dplyr::inner_join(plage_PTB_Control_CI_sens, by = "Signature") |> 
    dplyr::inner_join(ssgsea_PTB_LTBI_CI_sens, by = "Signature") |> 
    dplyr::inner_join(plage_PTB_LTBI_CI_sens, by = "Signature")
table_mean_sens <- table_mean_sens |> 
    dplyr::filter(!Signature == "Hoang_OD_3")

table_mean_sens$Signature <- gsub("_", "\\\\_", table_mean_sens$Signature)
write.table(table_mean_sens, quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = " & ",
            file = file.path(wd, "Figures_and_tables/table_mean_sens.txt")) 
# awk '{print $0, " \\\\\\\hline"}' table_mean_AUC.txt > table_mean_AUC_for_latex.txt

#### Table: Weighted Specificity (95% CI) for PTB vs. Control and PTB vs. LTBI ####
ssgsea_PTB_Control_CI_spec <- ssgsea_PTB_Control_CI |> 
    dplyr::mutate(ssgsea_PTB_Control = Specificity) |> 
    dplyr::select(Signature, ssgsea_PTB_Control)

ssgsea_PTB_LTBI_CI_spec <- ssgsea_PTB_LTBI_CI |> 
    dplyr::mutate(ssgsea_PTB_LTBI = Specificity) |> 
    dplyr::select(Signature, ssgsea_PTB_LTBI)

plage_PTB_Control_CI_spec <- plage_PTB_Control_CI |> 
    dplyr::mutate(plage_PTB_Control = Specificity) |> 
    dplyr::select(Signature, plage_PTB_Control)

plage_PTB_LTBI_CI_spec <- plage_PTB_LTBI_CI |> 
    dplyr::mutate(plage_PTB_LTBI = Specificity) |> 
    dplyr::select(Signature, plage_PTB_LTBI)

table_mean_spec <- ssgsea_PTB_Control_CI_spec |> 
    dplyr::inner_join(plage_PTB_Control_CI_spec, by = "Signature") |> 
    dplyr::inner_join(ssgsea_PTB_LTBI_CI_spec, by = "Signature") |> 
    dplyr::inner_join(plage_PTB_LTBI_CI_spec, by = "Signature")
table_mean_spec <- table_mean_spec |> 
    dplyr::filter(!Signature == "Hoang_OD_3")

table_mean_spec$Signature <- gsub("_", "\\\\_", table_mean_spec$Signature)
write.table(table_mean_spec, quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = " & ",
            file = file.path(wd, "Figures_and_tables/table_mean_spec.txt")) 

#### Figure 1: Flowchart of curatedTBData #### 
##### Enseml ID for RNA-seq ####
geo <- "GSE107991"
sequencePlatform <- "GPL20301"
gse <- GEOquery::getGEO(geo, GSEMatrix = FALSE)
data_list <- readDataGSE107995(geo)
probeInfo <- c("Genes", "Gene_name", "Gene_biotype")
data_list_processed <- lapply(data_list, function(x) {
    x_counts <- x %>%
        dplyr::select(-probeInfo) %>% as.matrix()
    colnames(x_counts) <- names(GEOquery::GSMList(gse))
    row.names(x_counts) <- x$Genes
    x_counts
})

data_Non_normalized_counts <- data_list_processed$data_Non_normalized
data_normalized_counts <- data_list_processed$data_normalized %>% 
    data.frame()
data_normalized_counts$SYMBOL <- data_list$data_normalized$Gene_name

gene_freq <- data_normalized_counts$SYMBOL |> 
    table() |> 
    unlist() |> 
    sort(decreasing = T) 
data_normalized_counts |> 
    dplyr::filter(SYMBOL == names(gene_freq[2])) |> 
    dplyr::select(SYMBOL, 1:2)
# Remove temporary files
unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
# probeset for microarray 
# geo <- "GSE36238"
# sequencePlatform <- "GPL570"
# GSE36238_normalized_rma <- readRawData(geo, sequencePlatform)
# colnames(GSE36238_normalized_rma) <- gsub("\\..*", "", colnames(GSE36238_normalized_rma))
# GSE36238_normalized_data <- GSE36238_normalized_rma
# row_data <- map_gene_symbol(GSE36238_normalized_rma, sequencePlatform)
# row_data |> 
#     as.data.frame() |> 
#     dplyr::filter(SYMBOL == "SEPTIN4") |> 
#     dplyr::select(ID_REF, SYMBOL)
##### Clinical annotation ####
# GSM2475317 from GSE94438

#### Figure 2 Venn diagram for overlapping ####
mean_AUC_combine <- ssgsea_PTB_Control_CI_sub |> 
    dplyr::inner_join(plage_PTB_Control_CI_sub, by = "Signature") |> 
    dplyr::inner_join(ssgsea_PTB_LTBI_CI_sub, by = "Signature") |> 
    dplyr::inner_join(plage_PTB_LTBI_CI_sub, by = "Signature")

df <- lapply(2:ncol(mean_AUC_combine), function(x) {
    cur <- mean_AUC_combine[, x]
    vapply(strsplit(cur, split = " "), function(x) 
        as.numeric(x[1]), numeric(1))
}) |> 
    dplyr::bind_cols() |> 
    as.data.frame()
colnames(df) <- colnames(mean_AUC_combine)[-1]
out <- cbind(Signature = mean_AUC_combine$Signature, df) |> 
    dplyr::filter(!Signature == "Hoang_OD_3")
    

venn_input_0.8 <-list('PTB vs. Control (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_Control >= 0.8]),
                  'PTB vs. Control (PLAGE)' = na.omit(out$Signature[out$plage_PTB_Control >= 0.8]),
                  'PTB vs. LTBI (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_LTBI >= 0.8]),
                  'PTB vs. LTBI (PLAGE)' = na.omit(out$Signature[out$plage_PTB_LTBI >= 0.8]))

venn_input_0.9 <-list('PTB vs. Control (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_Control >= 0.9]),
                  'PTB vs. Control (PLAGE)' = na.omit(out$Signature[out$plage_PTB_Control >= 0.9]),
                  'PTB vs. LTBI (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_LTBI >= 0.9]),
                  'PTB vs. LTBI (PLAGE)' = na.omit(out$Signature[out$plage_PTB_LTBI >= 0.9]))

venn_input_low <-list('PTB vs. Control (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_Control < 0.8]),
                      'PTB vs. Control (PLAGE)' = na.omit(out$Signature[out$plage_PTB_Control < 0.8]),
                      'PTB vs. LTBI (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_LTBI < 0.8]),
                      'PTB vs. LTBI (PLAGE)' = na.omit(out$Signature[out$plage_PTB_LTBI < 0.8]))

# create venn diagram and display all the sets
figure2A_0.8 <- ggvenn::ggvenn(venn_input_0.8, show_percentage= T, show_elements = FALSE, set_name_size = 4,
               fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
               stroke_size = 0.1)
figure2B_0.9 <- ggvenn::ggvenn(venn_input_0.9, show_percentage= T, show_elements = FALSE, set_name_size = 4,
               fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
               stroke_size = 0.1)


ggsave(figure2A_0.8,
       file = file.path(wd, "Figures_and_tables/Figure2A_0.8_venn_diagram.pdf"), 
       device = "pdf", height = 10, width = 10)
ggsave(figure2B_0.9,
       file = file.path(wd, "Figures_and_tables/Figure2B_0.9_venn_diagram.pdf"), 
       device = "pdf", height = 10, width = 10)

#### Figure S1 PTB vs. Control ####
##### Figure S1A heatmap AUC ssGSEA PTB vs. Control #####

heatmap_ssgsea_PTB_Control <- ssgsea_PTB_Control_out_combine |> 
    dplyr::filter(!Signature == "Hoang_OD_3") |> 
    heatmap_auc(GSE_sig = NULL, facet = TRUE, clustering = FALSE) + 
    theme(legend.position = "bottom")

##### Figure S1B heatmap AUC PLAGE PTB vs. Control #####
heatmap_plage_PTB_Control <- plage_PTB_Control_out_combine |> 
    dplyr::filter(!Signature == "Hoang_OD_3") |> 
    heatmap_auc(GSE_sig = NULL, facet = TRUE, clustering = FALSE) +
    theme(legend.position="bottom")


##### Figure S1 #####
figureS1 <- cowplot::plot_grid(heatmap_ssgsea_PTB_Control, 
                               heatmap_plage_PTB_Control,
                               align = "h")
ggsave(figureS1,
       file = file.path(wd, "Figures_and_tables/heatmap_PTB_Control.pdf"),
       device = "pdf", height = 18, width = 20)

#### Figure S2 PTB vs. LTBI ####
##### Figure S2A heatmap AUC ssGSEA PTB vs. LTBI #####
heatmap_ssgsea_PTB_LTBI <- ssgsea_PTB_LTBI_out_combine |> 
    dplyr::filter(!Signature == "Hoang_OD_3") |> 
    heatmap_auc(GSE_sig = NULL, facet = TRUE, clustering = FALSE) +
    theme(legend.position = "bottom")

##### Figure S2B heatmap AUC PLAGE PTB vs. LTBI #####
heatmap_plage_PTB_LTBI <- plage_PTB_LTBI_out_combine |> 
    dplyr::filter(!Signature == "Hoang_OD_3") |> 
    heatmap_auc(GSE_sig = NULL, facet = TRUE, clustering = FALSE) +
    theme(legend.position="bottom")

##### Figure S2 #####
figureS2 <- cowplot::plot_grid(heatmap_ssgsea_PTB_LTBI, 
                               heatmap_plage_PTB_LTBI,
                               align = "h")
ggsave(figureS2,
       file = file.path(wd, "Figures_and_tables/heatmap_PTB_LTBI.pdf"),
       device = "pdf", height = 18, width = 18)

#### Ensemble learning ####
sigatures_names_sub <- c("Maertzdorf_4", "Maertzdorf_15", "LauxdaCosta_OD_3",
                     "Verhagen_10", "Jacobsen_3", "Sambarey_HIV_10",
                     "Leong_24", "Berry_OD_86", "Berry_393", "Bloom_OD_144",
                     "Suliman_RISK_4", "Zak_RISK_16", "Leong_RISK_29",
                     "Anderson_42", "Anderson_OD_51", "Kaforou_27",
                     "Kaforou_OD_44", "Kaforou_OD_53", "Sweeney_OD_3",
                     "Blankley_5", "Singhania_OD_20", "Thompson_9", "Esmail_82",
                     "Thompson_FAIL_13", "Tornheim_71", "Darboe_RISK_11", 
                     "Hoang_OD_20", "Roe_3", "Gong_OD_4",
                     "Gjoen_10", "Duffy_23", "Jenum_8", "Qian_OD_17", "Estevez_133")
ensemble_sig_list <- make_ensembleSignaturesList(sigatures_names_sub, n = 30, 
                                                 size = 5)

# Remove Hoang_OD_3, GSE6112, GSE74092, and study with sample size < filter_size for PTB vs. Control
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

# Remove Hoang_OD_3, GSE6112, and study with sample size < filter_size for PTB vs. LTBI 
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
# See ensemble_results.R and extract_ensemble_results.R for obtaining results for ensemble
##### ssGSEA PTB vs. Control #####
ssgsea_PTB_Control_out_edit <- file.path(wd, "data/ssgsea_PTB_Control_out_edit.RDS") |> 
    readRDS()
# ssgsea_PTB_Control_out_edit <- PTB_Control_edit(ssgsea_PTB_Control_out, 
#                                                 filter_size = 15,
#                                                 out_file_path = file.path(wd, "data/ssgsea_PTB_Control_out_edit.RDS"),
#                                                 save_file = TRUE)
ssgsea_PTB_Control_edit_AUC <- ssgsea_PTB_Control_out_edit |> 
    extract_signature_score(sigatures_names_sub) |> 
    get_AUC_for_each_signature(combineSignatureName = NULL)

ssgsea_PTB_Control_ensl <- read.delim(file.path(wd, "data/ssgsea_PTB_Control.txt"))

ssgsea_PTB_Control_ensl_combine <- find_sig_with_max_AUC(
    gsea_AUC = ssgsea_PTB_Control_edit_AUC, 
    gsea_ensl = ssgsea_PTB_Control_ensl,
    ensemble_sig_list = ensemble_sig_list, 
    gsea_method = "ssGSEA")

##### PLAGE PTB vs. Control #####
plage_PTB_Control_out_edit <- file.path(wd, "data/plage_PTB_Control_out_edit.RDS") |> 
    readRDS()
# plage_PTB_Control_out_edit <- PTB_Control_edit(plage_PTB_Control_out,
#                                                filter_size = 15,
#                                                out_file_path = file.path(wd, "data/plage_PTB_Control_out_edit.RDS"),
#                                                save_file = TRUE)

plage_PTB_Control_edit_AUC <- plage_PTB_Control_out_edit |> 
    extract_signature_score(sigatures_names_sub) |> 
    get_AUC_for_each_signature(combineSignatureName = NULL)

plage_PTB_Control_ensl <- read.delim(file.path(wd, "data/plage_PTB_Control.txt")) |> 
    dplyr::mutate(method = "PLAGE")

plage_PTB_Control_ensl_combine <- find_sig_with_max_AUC(
    gsea_AUC = plage_PTB_Control_edit_AUC, 
    gsea_ensl = plage_PTB_Control_ensl,
    ensemble_sig_list = ensemble_sig_list, 
    gsea_method = "PLAGE")

##### Figure S3A Plot: PTB vs. Control ensemble #####
PTB_Control_ensl <- list()
PTB_Control_ensl$AUC_final <- rbind(ssgsea_PTB_Control_ensl_combine$AUC_final, 
                                    plage_PTB_Control_ensl_combine$AUC_final)
PTB_Control_ensl$signature_with_max_AUC <- rbind(ssgsea_PTB_Control_ensl_combine$signature_with_max_AUC,
                                                 plage_PTB_Control_ensl_combine$signature_with_max_AUC)

p_PTB_Control_ensemble <- ggplot(PTB_Control_ensl$AUC_final,
                                 aes(Method, AUC, fill = Signature_type)) +
        labs(fill = "Signature type") +
        geom_split_violin() +
        facet_wrap(~Set) +
        geom_text(data = PTB_Control_ensl$signature_with_max_AUC,
                  aes(x = Method, y = 0.88, label = Signature), size = 3, angle = 10) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 8, face = "bold"))

# ggsave(file.path(wd, "Figures_and_tables/FigureS3_PTB_Control_ensemble.pdf"),
#        p_PTB_Control_ensemble, device = "pdf", width = 14, height = 10)

##### ssGSEA PTB vs. LTBI #####
# ssgsea_PTB_LTBI_out_edit <- PTB_LTBI_edit(ssgsea_PTB_LTBI_out,
#                                           filter_size = 15,
#                                           out_file_path = file.path(wd, "data/ssgsea_PTB_LTBI_out_edit.RDS"),
#                                           save_file = TRUE)
ssgsea_PTB_LTBI_out_edit <- file.path(wd, "data/ssgsea_PTB_LTBI_out_edit.RDS") |> 
    readRDS()
ssgsea_PTB_LTBI_edit_AUC <- ssgsea_PTB_LTBI_out_edit |> 
    extract_signature_score(sigatures_names_sub) |> 
    get_AUC_for_each_signature(combineSignatureName = NULL)

ssgsea_PTB_LTBI_ensl <- read.delim(file.path(wd, "data/ssgsea_PTB_LTBI.txt"))

ssgsea_PTB_LTBI_ensl_combine <- find_sig_with_max_AUC(
    gsea_AUC = ssgsea_PTB_LTBI_edit_AUC, 
    gsea_ensl = ssgsea_PTB_LTBI_ensl,
    ensemble_sig_list = ensemble_sig_list, 
    gsea_method = "ssGSEA")

##### PLAGE PTB vs. LTBI #####
# plage_PTB_LTBI_out_edit <- PTB_LTBI_edit(plage_PTB_LTBI_out,
#                                          filter_size = 15,
#                                          out_file_path = file.path(wd, "data/plage_PTB_LTBI_out_edit.RDS"),
#                                          save_file = TRUE)
# plage_PTB_LTBI_out_edit <- PTB_LTBI_edit(plage_PTB_LTBI_out,
#                                          filter_size = 15,
#                                          out_file_path = NULL,
#                                          save_file = FALSE)
plage_PTB_LTBI_out_edit <- file.path(wd, "data/plage_PTB_LTBI_out_edit.RDS") |> 
    readRDS()
plage_PTB_LTBI_edit_AUC <- plage_PTB_LTBI_out_edit |> 
    extract_signature_score(sigatures_names_sub) |> 
    get_AUC_for_each_signature(combineSignatureName = NULL)

plage_PTB_LTBI_ensl <- read.delim(file.path(wd, "data/plage_PTB_LTBI.txt"))

plage_PTB_LTBI_ensl_combine <- find_sig_with_max_AUC(
    gsea_AUC = plage_PTB_LTBI_edit_AUC, 
    gsea_ensl = plage_PTB_LTBI_ensl,
    ensemble_sig_list = ensemble_sig_list, 
    gsea_method = "PLAGE")

##### Figure S3B: PTB vs. LTBI ensemble #####
PTB_LTBI_ensl <- list()
PTB_LTBI_ensl$AUC_final <- rbind(ssgsea_PTB_LTBI_ensl_combine$AUC_final, 
                                 plage_PTB_LTBI_ensl_combine$AUC_final)
PTB_LTBI_ensl$signature_with_max_AUC <- rbind(ssgsea_PTB_LTBI_ensl_combine$signature_with_max_AUC,
                                              plage_PTB_LTBI_ensl_combine$signature_with_max_AUC)

p_PTB_LTBI_ensemble <- ggplot(PTB_LTBI_ensl$AUC_final,
                              aes(Method, AUC, fill = Signature_type)) +
    labs(fill = "Signature type") +
    geom_split_violin() +
    facet_wrap(~Set) +
    geom_text(data = PTB_LTBI_ensl$signature_with_max_AUC,
              aes(x = Method, y = 0.84, label = Signature), size = 3, angle = 10) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, face = "bold"))

# ggsave(file.path(wd, "Figures_and_tables/FigureS3_PTB_LTBI_ensemble.pdf"),
#        p_PTB_LTBI_ensemble, device = "pdf", width = 14, height = 10)

##### Figure 3: Sub Plot for ensemble distribution #####
my_set_viollin <- c("Set 6", "Set 21", "Set 22")

AUC_final_PTB_Control <- PTB_Control_ensl$AUC_final |> 
    dplyr::filter(Set %in% my_set_viollin) |> 
    dplyr::mutate(Comparison = "PTB vs. Control")
signature_with_max_AUC_PTB_Control <- PTB_Control_ensl$signature_with_max_AUC |> 
    dplyr::filter(Set %in% my_set_viollin) |> 
    dplyr::mutate(Comparison = "PTB vs. Control")
AUC_final_PTB_LTBI <- PTB_LTBI_ensl$AUC_final |> 
    dplyr::filter(Set %in% my_set_viollin) |> 
    dplyr::mutate(Comparison = "PTB vs. LTBI")
signature_with_max_AUC_PTB_LTBI <- PTB_LTBI_ensl$signature_with_max_AUC |> 
    dplyr::filter(Set %in% my_set_viollin) |> 
    dplyr::mutate(Comparison = "PTB vs. LTBI")

AUC_final_sub <- rbind(AUC_final_PTB_Control, AUC_final_PTB_LTBI)

signature_with_max_AUC_sub <- rbind(signature_with_max_AUC_PTB_Control,
                                    signature_with_max_AUC_PTB_LTBI)

p_ensl_sub_PTB_Control <- AUC_final_sub |> 
    dplyr::filter(Comparison == "PTB vs. Control") |> 
    ggplot(aes(Method, AUC, fill = Signature_type)) +
    labs(fill = "Signature type") +
    geom_split_violin() + 
    facet_grid(cols = vars(Set)) +
    geom_text(data = signature_with_max_AUC_PTB_Control,
              aes(x = Method, y = 0.90, label = Signature), 
              size = 3, angle = 10) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, face = "bold"))
# ggsave(file.path(wd, "Figures_and_tables/p_ensl_sub_PTB_Control.pdf"),
#        p_ensl_sub_PTB_Control, device = "pdf", width = 10, height = 3)

p_ensl_sub_PTB_LTBI <- AUC_final_sub |> 
    dplyr::filter(Comparison == "PTB vs. LTBI") |> 
    ggplot(aes(Method, AUC, fill = Signature_type)) +
    labs(fill = "Signature type") +
    geom_split_violin() + 
    facet_grid(cols = vars(Set)) +
    geom_text(data = signature_with_max_AUC_PTB_LTBI,
              aes(x = Method, y = 0.88, label = Signature), 
              size = 3, angle = 10) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, face = "bold"))

# ggsave(file.path(wd, "Figures_and_tables/p_ensl_sub_PTB_LTBI.pdf"),
#        p_ensl_sub_PTB_LTBI, device = "pdf", width = 10, height = 3)
##### p values for viollin plot comparisons #####
compare_single_ensl <- function(df_ensl, set_name, df_single, single_name) {
    AUC_set_ensl <- df_ensl |> 
        dplyr::filter(Set == set_name)
    AUC_set_single <- df_single |> 
        dplyr::filter(Signature == single_name)
    out <- AUC_set_ensl |> 
        dplyr::inner_join(AUC_set_single, by = "Study")
    
    re <- weights::wtd.t.test(x = out$AUC.x, 
                              y = out$AUC.y, 
                              weight = out$sample_size.x)    
    data.frame(Set = set_name, p_value = re$coefficients["p.value"],
               row.names = NULL)
}

get_AUC_CI_ensl_single <- function(df, set,type, 
                                   comparison, method) {
    df |> 
        dplyr::filter(Set == set,
                      Comparison == comparison,
                      Signature_type == type,
                      Method == method) |> 
        dplyr::pull(AUC) |> 
        stats::quantile(c(0.025, 0.5, 0.975), na.rm = TRUE)
}

get_pvalue_ensl <- function(sig_with_max_AUC, method, df_ensl, df_single) {
    df <- sig_with_max_AUC |> 
        dplyr::filter(Method == method)
    lapply(1:nrow(df), function(i) {
        df_one <- df[i,,drop=FALSE]
        set_name <- gsub(" ", "", df_one$Set) |> 
            as.character()
        single_name <- df_one$Signature
        # print(c(i, set_name, single_name))
        compare_single_ensl(df_ensl, set_name, df_single, single_name) 
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::mutate(Method = method)
}

p_value_plage_PTB_Control <- PTB_Control_ensl$signature_with_max_AUC |> 
    get_pvalue_ensl("PLAGE", plage_PTB_Control_ensl, plage_PTB_Control_edit_AUC)

p_value_ssgsea_PTB_Control <- PTB_Control_ensl$signature_with_max_AUC |> 
    get_pvalue_ensl("ssGSEA", ssgsea_PTB_Control_ensl, ssgsea_PTB_Control_edit_AUC)

p_value_plage_PTB_LTBI <- PTB_LTBI_ensl$signature_with_max_AUC |> 
    get_pvalue_ensl("PLAGE", plage_PTB_LTBI_ensl, plage_PTB_LTBI_edit_AUC)

p_value_ssgsea_PTB_LTBI <- PTB_LTBI_ensl$signature_with_max_AUC |> 
    get_pvalue_ensl("ssGSEA", ssgsea_PTB_LTBI_ensl, ssgsea_PTB_LTBI_edit_AUC)

###### Set 6 ######
# PTB vs. Control
p_value_plage_PTB_Control |> 
    dplyr::filter(Set == "Set6") # 8.996961e-03
p_value_ssgsea_PTB_Control |> 
    dplyr::filter(Set == "Set6") # 0.0001411727
## PLAGE
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 6", type = "Ensemble",
                       comparison = "PTB vs. Control", method = "PLAGE")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 6", type = "Single",
                       comparison = "PTB vs. Control", method = "PLAGE")
## ssGSEA
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 6", type = "Ensemble",
                       comparison = "PTB vs. Control", method = "ssGSEA")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 6", type = "Single",
                       comparison = "PTB vs. Control", method = "ssGSEA")
# PTB vs. LTBI
p_value_plage_PTB_LTBI |> 
    dplyr::filter(Set == "Set6") # 0.01043963
p_value_ssgsea_PTB_LTBI |> 
    dplyr::filter(Set == "Set6") # 9.470883e

###### Set 21 ######
# PTB vs. Control
p_value_plage_PTB_Control |> 
    dplyr::filter(Set == "Set21") # 0.7820349
p_value_ssgsea_PTB_Control |> 
    dplyr::filter(Set == "Set21") # 1.294722e-06

get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Ensemble",
                       comparison = "PTB vs. Control", method = "PLAGE")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Single",
                       comparison = "PTB vs. Control", method = "PLAGE")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Ensemble",
                       comparison = "PTB vs. Control", method = "ssGSEA")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Single",
                       comparison = "PTB vs. Control", method = "ssGSEA")


p_value_plage_PTB_LTBI |> 
    dplyr::filter(Set == "Set21") # 2.105078e-09
p_value_ssgsea_PTB_LTBI |> 
    dplyr::filter(Set == "Set21") # 9.907839e-06


get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Ensemble",
                       comparison = "PTB vs. LTBI", method = "PLAGE")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Single",
                       comparison = "PTB vs. LTBI", method = "PLAGE")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Ensemble",
                       comparison = "PTB vs. LTBI", method = "ssGSEA")
get_AUC_CI_ensl_single(AUC_final_sub, set = "Set 21", type = "Single",
                       comparison = "PTB vs. LTBI", method = "ssGSEA")
###### Set 22 ######
p_value_plage_PTB_Control |> 
    dplyr::filter(Set == "Set22") # 0.1178086
p_value_ssgsea_PTB_Control |> 
    dplyr::filter(Set == "Set22") # 0.9434711

p_value_plage_PTB_LTBI |> 
    dplyr::filter(Set == "Set22") # 0.01167012
p_value_ssgsea_PTB_LTBI |> 
    dplyr::filter(Set == "Set22") # 1.267017e-06

##### Figure: Sub Ridge plot for PTB vs. Control #####
my_set_ridge <- c("Set6", "Set21", "Set22")
get_ridge_plot <- function(df_ensl1, df_single1, method1,
                           df_ensl2, df_single2, method2, 
                           set_name, threshold = 0.8) {
    df_ridge1 <- get_results_for_ridge(df_ensl1, df_single1, set_name, threshold) |> 
        dplyr::mutate(Method = method1)
    df_ridge2 <- get_results_for_ridge(df_ensl2, df_single2, set_name, threshold) |> 
        dplyr::mutate(Method = method2)
    df_ridge <- rbind(df_ridge1, df_ridge2)
    
    ggplot(df_ridge, aes(x = AUC, y = Signature)) +
        geom_density_ridges() +
        facet_wrap(~ Method) +
        ylab(NULL) + xlab(NULL) +
        ggtitle(paste("Set", gsub("Set", "", set_name))) +
        theme_bw() +
        theme(axis.text.y=element_text(angle=30))
}

p_ridge_Control_sub_list <- lapply(my_set_ridge, function(set_name) {
    get_ridge_plot(ssgsea_PTB_Control_ensl, ssgsea_PTB_Control_edit_AUC, "ssGSEA",
                   plage_PTB_Control_ensl, plage_PTB_Control_edit_AUC, "PLAGE",
                   set_name, threshold = 0.5)
})
p_ridge_sub_PTB_Control <- do.call("grid.arrange", 
                               c(p_ridge_Control_sub_list, ncol= 3))
# ggsave(file.path(wd, "Figures_and_tables/p_ridge_sub_PTB_Control.pdf"),
#        p_ridge_sub_PTB_Control, device = "pdf", width = 16, height = 8)
##### Figure: Sub Ridge plot for PTB vs. LTBI #####
p_ridge_LTBI_sub_list <- lapply(my_set_ridge, function(set_name) {
    get_ridge_plot(ssgsea_PTB_LTBI_ensl, ssgsea_PTB_LTBI_edit_AUC, "ssGSEA",
                   plage_PTB_LTBI_ensl, plage_PTB_LTBI_edit_AUC, "PLAGE",
                   set_name, threshold = 0.5)
})
p_ridge_sub_PTB_LTBI <- do.call("grid.arrange", 
                            c(p_ridge_LTBI_sub_list, ncol= 3))
# ggsave(file.path(wd, "Figures_and_tables/p_ridge_sub_PTB_LTBI.pdf"),
#        p_ridge_sub_PTB_LTBI, device = "pdf", width = 16, height = 8)

##### Weighted SD for ridge plot ####
get_wtd_var <- function(x, wt) {
    wt <- wt / sum(wt)
    n <- length(x)
    xm <- weighted.mean(x,  wt)
    sum(wt * (x - xm) ^ 2) / (n - 1)
}

get_sd_for_ridge <- function(df_ensl, df_single, set_name, single_name, 
                             threshold) {
    df_ridge <- get_results_for_ridge(df_ensl, df_single, set_name, threshold)
    df_ridge |> 
        dplyr::group_by(Signature) |> 
        dplyr::summarise(standard_sd = sd(AUC),
                         wtd_var = get_wtd_var(AUC, sample_size)) |> 
        dplyr::mutate(wtd_sd = sqrt(wtd_var)) |> 
        dplyr::filter(Signature %in% c(set_name, single_name))
}
get_sd_ensl <- function(sig_with_max_AUC, method, 
                        df_ensl, df_single, threshold = 0.5) {
    df <- sig_with_max_AUC |> 
        dplyr::filter(Method == method)
    lapply(1:nrow(df), function(i) {
        df_one <- df[i,,drop=FALSE]
        set_name <- gsub(" ", "", df_one$Set) |> 
            as.character()
        single_name <- df_one$Signature
        get_sd_for_ridge(df_ensl, df_single, set_name, single_name, 
                         threshold) |> 
            dplyr::mutate(Set = set_name)
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::mutate(Method = method)
}
get_sd_compare <- function(sd_method_re) {
    sd_method_re |> 
        dplyr::group_by(Set) |> 
        dplyr::mutate(diff_wtd_sd = wtd_sd - wtd_sd[!Signature == Set],
                      diff_std_sd = standard_sd - standard_sd[!Signature == Set]) |> 
        dplyr::filter(!diff_wtd_sd == 0, !diff_std_sd == 0)
}

###### Compute all set ######
threshold <- 0.5
# PLAGE PTB vs. Control
sd_plage_PTB_Control_sig <- plage_PTB_Control_edit_AUC |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(standard_sd = sd(AUC),
                     wtd_var = get_wtd_var(AUC, sample_size),
                     wtd_std = sqrt(wtd_var))
sd_plage_PTB_Control <- get_sd_ensl(sig_with_max_AUC = PTB_Control_ensl$signature_with_max_AUC, 
                                    method = "PLAGE", 
                                    df_ensl = plage_PTB_Control_ensl, 
                                    df_single = plage_PTB_Control_edit_AUC,
                                    threshold = threshold)
sd_plage_PTB_Control_compare <- get_sd_compare(sd_plage_PTB_Control)
table(sd_plage_PTB_Control_compare$diff_wtd_sd <= 0)
# ssGSEA PTB vs. Control
sd_ssgsea_PTB_Control_sig <- ssgsea_PTB_Control_edit_AUC |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(standard_sd = sd(AUC),
                     wtd_var = get_wtd_var(AUC, sample_size),
                     wtd_std = sqrt(wtd_var))
sd_ssgsea_PTB_Control <- get_sd_ensl(PTB_Control_ensl$signature_with_max_AUC, "ssGSEA", 
                               ssgsea_PTB_Control_ensl, ssgsea_PTB_Control_edit_AUC,
                               threshold = threshold)
sd_ssgsea_PTB_Control_compare <- get_sd_compare(sd_ssgsea_PTB_Control)
table(sd_ssgsea_PTB_Control_compare$diff_wtd_sd <= 0)
# PLAGE PTB vs. LTBI
sd_plage_PTB_LTBI_sig <- plage_PTB_LTBI_edit_AUC |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(standard_sd = sd(AUC),
                     wtd_var = get_wtd_var(AUC, sample_size),
                     wtd_std = sqrt(wtd_var))
sd_plage_PTB_LTBI <- get_sd_ensl(PTB_LTBI_ensl$signature_with_max_AUC, "PLAGE", 
                            plage_PTB_LTBI_ensl, plage_PTB_LTBI_edit_AUC,
                            threshold = threshold)
sd_plage_PTB_LTBI_compare <- get_sd_compare(sd_plage_PTB_LTBI)
table(sd_plage_PTB_LTBI_compare$diff_wtd_sd <= 0)
# ssGSEA PTB vs. LTBI
sd_ssgsea_PTB_LTBI_sig <- ssgsea_PTB_LTBI_edit_AUC |> 
    dplyr::group_by(Signature) |> 
    dplyr::summarise(standard_sd = sd(AUC),
                     wtd_var = get_wtd_var(AUC, sample_size),
                     wtd_std = sqrt(wtd_var))
sd_ssgsea_PTB_LTBI <- get_sd_ensl(PTB_LTBI_ensl$signature_with_max_AUC, "ssGSEA", 
                             ssgsea_PTB_LTBI_ensl, ssgsea_PTB_LTBI_edit_AUC,
                             threshold = threshold)
sd_ssgsea_PTB_LTBI_compare <- get_sd_compare(sd_ssgsea_PTB_LTBI)
table(sd_ssgsea_PTB_LTBI_compare$diff_wtd_sd <= 0)


###### Set 6 ###### 
# PTB vs. Control
## PLAGE
sd_plage_PTB_Control |> 
    dplyr::filter(Set == "Set6")
sd_plage_PTB_Control_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set6"]]) |> 
    dplyr::arrange(wtd_std)
## ssGSEA
sd_ssgsea_PTB_Control |> 
    dplyr::filter(Set == "Set6")
sd_ssgsea_PTB_Control_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set6"]]) |> 
    dplyr::arrange(wtd_std)

# PTB vs. LTBI
## PLAGE
sd_plage_PTB_LTBI |> 
    dplyr::filter(Set == "Set6")
sd_plage_PTB_LTBI_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set6"]]) |> 
    dplyr::arrange(wtd_std)
## ssGSEA
sd_ssgsea_PTB_LTBI |> 
    dplyr::filter(Set == "Set6")
sd_ssgsea_PTB_LTBI_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set6"]]) |> 
    dplyr::arrange(wtd_std)

###### Set 21 ###### 
# PTB vs Control
## PLAGE
sd_plage_PTB_Control |> 
    dplyr::filter(Set == "Set21")
sd_plage_PTB_Control_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set21"]]) |> 
    dplyr::arrange(wtd_std)
## ssGSEA
sd_ssgsea_PTB_Control_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set21"]]) |> 
    dplyr::arrange(wtd_std)
sd_ssgsea_PTB_Control |> 
    dplyr::filter(Set == "Set21")

# PTB vs. LTBI
## PLAGE
sd_plage_PTB_LTBI |> 
    dplyr::filter(Set == "Set21")
sd_plage_PTB_LTBI_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set21"]]) |> 
    dplyr::arrange(wtd_std)
## ssGSEA
sd_ssgsea_PTB_LTBI |> 
    dplyr::filter(Set == "Set21")
sd_ssgsea_PTB_LTBI_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set21"]]) |> 
    dplyr::arrange(wtd_std)

###### Set 22 ###### 
# PTB vs Control
## PLAGE
sd_plage_PTB_Control |> 
    dplyr::filter(Set == "Set22")
sd_plage_PTB_Control_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set22"]]) |> 
    dplyr::arrange(wtd_std)
## ssGSEA
sd_ssgsea_PTB_Control_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set22"]]) |> 
    dplyr::arrange(wtd_std)
sd_ssgsea_PTB_Control |> 
    dplyr::filter(Set == "Set22")
# PTB vs. LTBI
## PLAGE
sd_plage_PTB_LTBI |> 
    dplyr::filter(Set == "Set22")
sd_plage_PTB_LTBI_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set22"]]) |> 
    dplyr::arrange(wtd_std)
## ssGSEA
sd_ssgsea_PTB_LTBI_sig |> 
    dplyr::filter(Signature %in% ensemble_sig_list[["Set22"]]) |> 
    dplyr::arrange(wtd_std)
sd_ssgsea_PTB_LTBI |> 
    dplyr::filter(Set == "Set22")

##### Figure S4A Coefficient PTB vs. Control #####
ssgsea_PTB_Control_coef <- wd |> 
    file.path("data/coef_output/ssgsea_PTB_Control_coef.RDS") |> 
    readRDS()
plage_PTB_Control_coef <- wd |> 
    file.path("data/coef_output/plage_PTB_Control_coef.RDS") |> 
    readRDS()
PTB_Control_coef <- mapply(function(x, y) {
    x <- x |> 
        dplyr::mutate(method = "ssGSEA")
    y <- y |> 
        dplyr::mutate(method = "PLAGE")
    rbind(x, y)
}, ssgsea_PTB_Control_coef, plage_PTB_Control_coef, SIMPLIFY = FALSE)

p_PTB_Control_coef_list <- lapply(1:length(PTB_Control_coef), function(i) {
    coef_df <- PTB_Control_coef[[i]]
    ggplot(coef_df, aes(x = Signature, y = Coefficient)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
        geom_boxplot() +
        scale_y_continuous(limits = c(-2, 2), 
                           breaks = seq(-2, 2, by = 1)) +
        ggtitle(names(PTB_Control_coef)[i]) + xlab(NULL) + ylab(NULL) +
        coord_flip() + 
        theme_bw() +
        facet_wrap(~method, ncol = 2)
})
p_PTB_Control_coef <- do.call("grid.arrange", 
                              c(p_PTB_Control_coef_list, ncol= 6))
# ggsave(file.path(wd, "Figures_and_tables/FigureS4_PTB_Control_coef.pdf"),
#        p_PTB_Control_coef, device = "pdf", width = 20, height = 10)

##### Figure S4B Coefficient PTB vs. LTBI #####
ssgsea_PTB_LTBI_coef <- file.path(wd, "data/coef_output/ssgsea_PTB_LTBI_coef.RDS") |> 
    readRDS()
plage_PTB_LTBI_coef <- file.path(wd, "data/coef_output/plage_PTB_LTBI_coef.RDS") |> 
    readRDS()

PTB_LTBI_coef <- mapply(function(x, y) {
    x <- x |> 
        dplyr::mutate(method = "ssGSEA")
    y <- y |> 
        dplyr::mutate(method = "PLAGE")
    rbind(x, y)
}, ssgsea_PTB_LTBI_coef, plage_PTB_LTBI_coef, SIMPLIFY = FALSE)

p_PTB_LTBI_coef_list <- lapply(1:length(PTB_LTBI_coef), function(i) {
    coef_df <- PTB_LTBI_coef[[i]]
    ggplot(coef_df, aes(x = Signature, y = Coefficient)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
        geom_boxplot() +
        scale_y_continuous(limits = c(-2, 2), 
                           breaks = seq(-2, 2, by = 1)) +
        ggtitle(names(PTB_LTBI_coef)[i]) + xlab(NULL) + ylab(NULL) +
        coord_flip() + 
        theme_bw() +
        facet_wrap(~method, ncol = 2)
})
p_PTB_LTBI_coef <- do.call("grid.arrange", 
                           c(p_PTB_LTBI_coef_list, ncol= 6))
# ggsave(file.path(wd, "Figures_and_tables/FigureS4_PTB_LTBI_coef.pdf"),
#        p_PTB_LTBI_coef, device = "pdf", width = 20, height = 10)

##### Figure S5 #####
sigatures_names_ensl <- c("Maertzdorf_4", "Maertzdorf_15", "LauxdaCosta_OD_3",
                     "Verhagen_10", "Jacobsen_3", "Sambarey_HIV_10",
                     "Leong_24", "Berry_OD_86", "Berry_393", "Bloom_OD_144",
                     "Suliman_RISK_4", "Zak_RISK_16", "Leong_RISK_29",
                     "Anderson_42", "Anderson_OD_51", "Kaforou_27",
                     "Kaforou_OD_44", "Kaforou_OD_53", "Sweeney_OD_3",
                     "Blankley_5", "Singhania_OD_20", "Thompson_9", "Esmail_82",
                     "Thompson_FAIL_13", "Tornheim_71", "Darboe_RISK_11", 
                     "Hoang_OD_20", "Roe_3", "Gong_OD_4",
                     "Gjoen_10", "Duffy_23", "Jenum_8", "Qian_OD_17", "Estevez_133")
ensemble_sig_list <- make_ensembleSignaturesList(sigatures_names_ensl, n = 30, 
                                                 size = 5)

get_score_distribution <- function(ensemble_sig_list, gsea_sobject) {
    gsea_scores <- lapply(1:length(ensemble_sig_list), function(i) {
        ensemble_sig <- ensemble_sig_list[[i]]
        lapply(1:length(ssgsea_PTB_Control_out_edit), function(j) {
            sobject <- ssgsea_PTB_Control_out_edit[[j]]
            colData(sobject) |> 
                as.data.frame() |> 
                dplyr::select(all_of(ensemble_sig), TBStatus) |> 
                dplyr::mutate(Study = names(ssgsea_PTB_Control_out_edit)[j])
        }) |> 
            dplyr::bind_rows() |> 
            reshape2::melt() |> 
            dplyr::mutate(Set_name = names(ensemble_sig_list)[i])
    }) 
    p_list <- lapply(gsea_scores, function(df) {
        ggplot(df, aes(x = variable, y = value, color = TBStatus)) +
            geom_boxplot() +
            xlab(NULL) + ylab(NULL) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_text(size = 6)) 
    })
    return(p_list)
    
}

get_score_all <- function(gsea_sobject) {
    gsea_scores <- lapply(gsea_sobject, function(x) {
        colData(x) |> 
            as.data.frame() |> 
            dplyr::select(all_of(sigatures_names_ensl), TBStatus)
    }) |> 
        dplyr::bind_rows() |> 
        reshape2::melt()
    return(gsea_scores)
    
}

get_score_mean_split <- function(df, compare) {
    df |> 
        dplyr::group_by(TBStatus, variable) |> 
        dplyr::summarize(value = mean(value)) |> 
        reshape2::dcast(variable ~ TBStatus) |> 
        dplyr::mutate(diff = PTB - .data[[compare]])
}
ssgsea_PTB_Control_all <- get_score_all(ssgsea_PTB_Control_out_edit) |> 
    dplyr::mutate(Method = "ssGSEA")
ssgsea_PTB_LTBI_all <- get_score_all(ssgsea_PTB_LTBI_out_edit) |> 
    dplyr::mutate(Method = "ssGSEA")


plage_PTB_Control_all <- get_score_all(plage_PTB_Control_out_edit) |> 
    dplyr::mutate(Method = "PLAGE") 
plage_PTB_LTBI_all <- get_score_all(plage_PTB_LTBI_out_edit) |> 
    dplyr::mutate(Method = "PLAGE")



PTB_Control_all <- rbind(ssgsea_PTB_Control_all, plage_PTB_Control_all)
PTB_Control_all$Method <- factor(PTB_Control_all$Method, 
                                 levels = c("ssGSEA", "PLAGE"))

PTB_Control_all |> 
    dplyr::group_by(Method, TBStatus, variable) |> 
    dplyr::summarise(mean_val = mean(value)) |> 
    dplyr::group_by(TBStatus, Method) |> 
    dplyr::summarise(sd_val = sd(mean_val))
p_PTB_Control <- ggplot(PTB_Control_all, 
                        aes(x = variable, y = value, color = TBStatus)) +
    geom_boxplot() +
    coord_flip() +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~ Method) +
    theme_bw()

# ggsave(filename = file.path(wd, "Figures_and_tables/FigureS5_PTB_Control_scores.pdf"),
#        p_PTB_Control, width = 10, height = 10, device = "pdf")

PTB_LTBI_all <- rbind(ssgsea_PTB_LTBI_all, plage_PTB_LTBI_all)
PTB_LTBI_all$Method <- factor(PTB_LTBI_all$Method, 
                              levels = c("ssGSEA", "PLAGE"))
PTB_LTBI_all |> 
    dplyr::group_by(Method, TBStatus, variable) |> 
    dplyr::summarise(mean_val = mean(value)) |> 
    dplyr::group_by(TBStatus, Method) |> 
    dplyr::summarise(sd_val = sd(mean_val))
p_PTB_LTBI <- ggplot(PTB_LTBI_all, 
                     aes(x = variable, y = value, color = TBStatus)) +
    geom_boxplot() +
    coord_flip() +
    xlab(NULL) + ylab(NULL) +
    facet_wrap(~ Method) +
    theme_bw()

# ggsave(filename = file.path(wd, "Figures_and_tables/FigureS5_PTB_LTBI_scores.pdf"),
#        p_PTB_LTBI, width = 10, height = 10, device = "pdf")


##### Figure S6 Rideplot for ensemble ####
n_col <- length(ensemble_sig_list) |> 
    sqrt() |> 
    floor()

p_ridge_Control_list <- lapply(names(ensemble_sig_list), function(set_name) {
    get_ridge_plot(ssgsea_PTB_Control_ensl, ssgsea_PTB_Control_edit_AUC, "ssGSEA",
                   plage_PTB_Control_ensl, plage_PTB_Control_edit_AUC, "PLAGE",
                   set_name, threshold = 0.5)
})
p_ridge_Control <- do.call("grid.arrange", 
                           c(p_ridge_Control_list, ncol= n_col))
# ggsave(filename = file.path(wd, "Figures_and_tables/FigureS6_PTB_Control_ridge.pdf"),
#        p_ridge_Control, width = 18, height = 20, device = "pdf")

p_ridge_LTBI_list <- lapply(names(ensemble_sig_list), function(set_name) {
    get_ridge_plot(ssgsea_PTB_LTBI_ensl, ssgsea_PTB_LTBI_edit_AUC, "ssGSEA",
                   plage_PTB_LTBI_ensl, plage_PTB_LTBI_edit_AUC, "PLAGE",
                   set_name, threshold = 0.5)
})
p_ridge_LTBI <- do.call("grid.arrange", 
                        c(p_ridge_LTBI_list, ncol= n_col))
# ggsave(filename = file.path(wd, "Figures_and_tables/FigureS6_PTB_LTBI_ridge.pdf"),
#        p_ridge_LTBI, width = 18, height = 20, device = "pdf")
##### Tables with metrics for ensemble ####
get_group_output <- function(df_ensl, df_ensl_combine, df_out_combine) {
    df_out_combine_group <- df_out_combine |> 
        dplyr::group_by(Signature) |> 
        dplyr::summarise(AUC = weighted.mean(AUC, Size),
                         Sensitivity = weighted.mean(Sensitivity, Size),
                         Specificity = weighted.mean(Specificity, Size))
    df_sig_re <- df_ensl_combine$signature_with_max_AUC |> 
        dplyr::inner_join(df_out_combine_group, by = "Signature") |> 
        dplyr::mutate(Set = gsub(" ", "", Set, fixed = TRUE))
    df_ensl_group <- df_ensl |> 
        dplyr::group_by(Set) |> 
        dplyr::summarise(AUC = weighted.mean(AUC, sample_size),
                         Sensitivity = weighted.mean(Sensitivity, sample_size),
                         Specificity = weighted.mean(Specificity, sample_size)) |> 
        dplyr::mutate(Signature = Set, Signature_type = "Ensemble", 
                      Method = unique(df_sig_re$Method))
    index_col <- match(colnames(df_sig_re), colnames(df_ensl_group))
    df_ensl_group <- df_ensl_group[, index_col]
    df_out <- rbind(df_ensl_group, df_sig_re)
    return(df_out)
}

make_group_tab <- function(df_out) {
    set_names <- unique(df_out$Set) |> 
        gtools::mixedsort()
    df_tab_reorder <- lapply(set_names, function(set_name) {
        df_out |> 
            dplyr::filter(Set == set_name) |> 
            dplyr::select(-c("Set", "Method")) |>
            dplyr::mutate(AUC = round(AUC, 2), 
                          Sensitivity = round(Sensitivity, 2),
                          Specificity = round(Specificity, 2)) |> 
            rbind("")
    }) |> 
        dplyr::bind_rows()
    df_tab <- df_tab_reorder |> 
        ggtexttable(rows = NULL, theme = ttheme("blank")) |> 
        tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) |> 
        tab_add_hline(at.row = nrow(df_tab_reorder), 
                      row.side = "bottom", linewidth = 2)
    return(df_tab)
}
# Make numeric table for ensemble results
ssgsea_PTB_Control_group <- ssgsea_PTB_Control_ensl |> 
    get_group_output(ssgsea_PTB_Control_ensl_combine, 
                     ssgsea_PTB_Control_out_combine)
ssgsea_PTB_Control_tab <- make_group_tab(ssgsea_PTB_Control_group)
# ggsave("~/Desktop/practice/curatedTBData_paper_results/Figures_and_tables/ssgsea_PTB_Control_ensl.pdf", 
#        ssgsea_PTB_Control_tab, width = 8, height = 25)

plage_PTB_Control_group <- plage_PTB_Control_ensl |> 
    get_group_output(plage_PTB_Control_ensl_combine, 
                     plage_PTB_Control_out_combine)
plage_PTB_Control_tab <- make_group_tab(plage_PTB_Control_group)
# ggsave("~/Desktop/practice/curatedTBData_paper_results/Figures_and_tables/plage_PTB_Control_ensl.pdf", 
#        plage_PTB_Control_tab, width = 8, height = 25)

ssgsea_PTB_LTBI_group <- ssgsea_PTB_LTBI_ensl |> 
    get_group_output(ssgsea_PTB_LTBI_ensl_combine, 
                     ssgsea_PTB_LTBI_out_combine)
ssgsea_PTB_LTBI_tab <- make_group_tab(ssgsea_PTB_LTBI_group)
# ggsave("~/Desktop/practice/curatedTBData_paper_results/Figures_and_tables/ssgsea_PTB_LTBI_ensl.pdf", 
#        ssgsea_PTB_LTBI_tab, width = 8, height = 25)

plage_PTB_LTBI_group <- plage_PTB_LTBI_ensl |> 
    get_group_output(plage_PTB_LTBI_ensl_combine, 
                     plage_PTB_LTBI_out_combine)
plage_PTB_LTBI_tab <- make_group_tab(plage_PTB_LTBI_group)
# ggsave("~/Desktop/practice/curatedTBData_paper_results/Figures_and_tables/plage_PTB_LTBI_ensl.pdf", 
#        plage_PTB_LTBI_tab, width = 8, height = 25)

##### Analysis in discussion ####
get_score_mean_split(ssgsea_PTB_Control_all, "Control") |> 
    dplyr::filter(diff < 0)
get_score_mean_split(ssgsea_PTB_LTBI_all, "LTBI") |> 
    dplyr::filter(diff < 0)
get_score_mean_split(plage_PTB_Control_all, "Control") |> 
    dplyr::filter(diff > 0)
get_score_mean_split(plage_PTB_LTBI_all, "LTBI") |> 
    dplyr::filter(diff > 0)

PTB_Control_ensl_AUC <- PTB_Control_ensl$AUC_final
PTB_Control_ensl_AUC_mean <- PTB_Control_ensl_AUC |> 
    dplyr::group_by(Method, Signature, Set, Signature_type) |> 
    dplyr::summarize(AUC = mean(AUC)) |> 
    reshape2::dcast(Set + Method ~ Signature_type) |> 
    dplyr::mutate(Diff_AUC = Ensemble - Single, Set_new = gsub(" ", "", Set))

plage_PTB_Control_compare <- PTB_Control_ensl_AUC_mean |> 
    dplyr::filter(Method == "PLAGE") |> 
    dplyr::inner_join(sd_plage_PTB_Control_compare, by = c("Set_new" = "Signature")) |> 
    dplyr::inner_join(p_value_plage_PTB_Control, by = c("Set_new" = "Set")) |> 
    dplyr::select(-c(Set.x,  Method.x))
plage_PTB_Control_compare |> filter(Diff_AUC > 0)
plage_PTB_Control_compare |> filter(Diff_AUC > 0, p_value < 0.05)
ifelse(plage_PTB_Control_compare |> 
           filter(Diff_AUC > 0) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()
ifelse(plage_PTB_Control_compare |> 
           filter(Diff_AUC > 0, p_value < 0.05) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()

ssgsea_PTB_Control_compare <- PTB_Control_ensl_AUC_mean |> 
    dplyr::filter(Method == "ssGSEA") |> 
    dplyr::inner_join(sd_ssgsea_PTB_Control_compare, by = c("Set_new" = "Signature")) |> 
    dplyr::inner_join(p_value_ssgsea_PTB_Control, by = c("Set_new" = "Set")) |> 
    dplyr::select(-c(Set.x,  Method.x))
ssgsea_PTB_Control_compare |> filter(Diff_AUC > 0)
ssgsea_PTB_Control_compare |> filter(Diff_AUC > 0, p_value < 0.05)
ifelse(ssgsea_PTB_Control_compare |> 
           filter(Diff_AUC > 0) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()
ifelse(ssgsea_PTB_Control_compare |> 
           filter(Diff_AUC > 0, p_value < 0.05) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()
ssgsea_PTB_Control_compare |> filter(abs(Diff_AUC) <= 1e-2) |> 
    dplyr::pull(diff_std_sd)


# compare_single_ensl(plage_PTB_Control_ensl, "Set6", 
#                     plage_PTB_Control_edit_AUC, "Roe_3")

PTB_LTBI_ensl_AUC <- PTB_LTBI_ensl$AUC_final
PTB_LTBI_ensl_AUC_mean <- PTB_LTBI_ensl_AUC |> 
    dplyr::group_by(Method, Signature, Set, Signature_type) |> 
    dplyr::summarize(AUC = mean(AUC)) |> 
    reshape2::dcast(Set + Method ~ Signature_type) |> 
    dplyr::mutate(Diff_AUC = Ensemble - Single, Set_new = gsub(" ", "", Set))


plage_PTB_LTBI_compare <- PTB_LTBI_ensl_AUC_mean |> 
    dplyr::filter(Method == "PLAGE") |> 
    dplyr::inner_join(sd_plage_PTB_LTBI_compare, by = c("Set_new" = "Signature")) |> 
    dplyr::inner_join(p_value_plage_PTB_LTBI, by = c("Set_new" = "Set")) |> 
    dplyr::select(-c(Set.x,  Method.x))
# plage_PTB_LTBI_compare |> filter(Diff_AUC > 0, p_value < 0.05/30)
ifelse(plage_PTB_LTBI_compare |> 
           filter(Diff_AUC > 0) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()
ifelse(plage_PTB_LTBI_compare |> 
           filter(Diff_AUC < 0) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()
ifelse(plage_PTB_LTBI_compare |> 
    filter(Diff_AUC > 0, p_value < 0.05) |> 
    dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()

ssgsea_PTB_LTBI_compare <- PTB_LTBI_ensl_AUC_mean |> 
    dplyr::filter(Method == "ssGSEA") |> 
    dplyr::inner_join(sd_ssgsea_PTB_LTBI_compare, by = c("Set_new" = "Signature")) |> 
    dplyr::inner_join(p_value_ssgsea_PTB_LTBI, by = c("Set_new" = "Set")) |> 
    dplyr::select(-c(Set.x,  Method.x))
# ssgsea_PTB_LTBI_compare |> filter(Diff_AUC > 0, p_value < 0.05/30)
ifelse(ssgsea_PTB_LTBI_compare |> 
           filter(Diff_AUC > 0) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()

ifelse(ssgsea_PTB_LTBI_compare |> 
           filter(Diff_AUC > 0, p_value < 0.05) |> 
           dplyr::pull(diff_wtd_sd) < 0 , "small", "large") |> 
    table()

union(plage_PTB_LTBI_compare |> 
          dplyr::filter(Diff_AUC > 0) |> 
          dplyr::pull(Set_new),
      ssgsea_PTB_LTBI_compare |> 
          dplyr::filter(Diff_AUC > 0) |> 
          dplyr::pull(Set_new))

#### Additional analysis for HIV positive ####
# # Find studies with HIV patients
# # All patients are HIV co-infected
# geo_HIV_all <- c("GSE69581", "GSE83892", "GSE50834", "GSE107104")
# # Some patients are HIV co-infected
# geo_HIV_some <- c("GSE37250", "GSE39939", "GSE39940")

# geo_HIV <- c(geo_HIV_all, geo_HIV_some)
# objects_list <- curatedTBData(study_name = geo_HIV, dry.run = FALSE)
objects_list_HIV <- lapply(geo_HIV_some, function(x) {
    obj <- objects_list[[x]]
    obj[, obj$HIVStatus == "Positive"]
    # obj[, (obj$HIVStatus == "Positive" & obj$TBStatus == "PTB") |
    #         (obj$HIVStatus == "Negative" & obj$TBStatus == "LTBI")]
})
names(objects_list_HIV) <- geo_HIV
objects_list_HIV_update <- lapply(objects_list_HIV, function(x) {
    dat <- x[["assay_curated"]]
    row_names <- row.names(dat) |>
        update_genenames()
    row.names(dat) <- row_names
    SummarizedExperiment(assays = list(dat), colData = colData(x))
})
##### ssGSEA PTB vs. Control #####
# ssgsea_PTB_Control_out_HIV <- lapply(objects_list_HIV_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                                 annotationCondition = c("PTB", "Control"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "ssGSEA", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(ssgsea_PTB_Control_out_HIV,
#         file.path(wd, "data/ssgsea_PTB_Control_out_HIV.RDS"))
# 
# ssgsea_PTB_Control_out_combine_HIV <- combine_auc(ssgsea_PTB_Control_out_HIV,
#                                                   annotationColName = "TBStatus",
#                                                   signatureColNames = names(signatures_list),
#                                                   num.boot = 1000, percent = 0.95)
# saveRDS(ssgsea_PTB_Control_out_combine_HIV,
#         file.path(wd, "data/ssgsea_PTB_Control_out_combine_HIV.RDS"))

ssgsea_PTB_Control_out_HIV <- file.path(wd, "data/ssgsea_PTB_Control_out_HIV.RDS") |> 
    readRDS()
study_PTB_Control_HIV <- lapply(1:length(ssgsea_PTB_Control_out_HIV), function(i) {
    sobject <- ssgsea_PTB_Control_out_HIV[[i]]
    sutdy_name <- gsub("_edit", "", names(ssgsea_PTB_Control_out_HIV)[i])
    data.frame(Study = sutdy_name, Size = ncol(sobject))
}) |> 
    dplyr::bind_rows()

ssgsea_PTB_Control_out_combine_HIV <- file.path(wd, "data/ssgsea_PTB_Control_out_combine_HIV.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_Control_HIV, "Study")

# ssgsea_PTB_Control_CI_HIV <- extract_CI(ssgsea_PTB_Control_out_combine_HIV)
# ssgsea_PTB_Control_CI_sub_HIV <- ssgsea_PTB_Control_CI_HIV |> 
#     dplyr::mutate(ssgsea_PTB_Control = AUC) |> 
#     dplyr::select(Signature, ssgsea_PTB_Control)

##### PLAGE PTB vs. Control #####
# plage_PTB_Control_out_HIV <- lapply(objects_list_HIV_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                                 annotationCondition = c("PTB", "Control"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "PLAGE", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(plage_PTB_Control_out_HIV,
#         file.path(wd, "data/plage_PTB_Control_out_HIV.RDS"))
# plage_PTB_Control_out_combine_HIV <- combine_auc(plage_PTB_Control_out_HIV,
#                                                   annotationColName = "TBStatus",
#                                                   signatureColNames = names(signatures_list),
#                                                   num.boot = 1000, percent = 0.95)
# saveRDS(plage_PTB_Control_out_combine_HIV,
#         file.path(wd, "data/plage_PTB_Control_out_combine_HIV.RDS"))

plage_PTB_Control_out_HIV <- file.path(wd, "data/plage_PTB_Control_out_HIV.RDS") |> 
    readRDS()
plage_PTB_Control_out_combine_HIV <- file.path(wd, "data/plage_PTB_Control_out_combine_HIV.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_Control_HIV, "Study")

# plage_PTB_Control_CI_HIV <- extract_CI(plage_PTB_Control_out_combine_HIV)
# plage_PTB_Control_CI_sub_HIV <- plage_PTB_Control_CI_HIV |> 
#     dplyr::mutate(plage_PTB_Control = AUC) |> 
#     dplyr::select(Signature, plage_PTB_Control)
##### ssGSEA PTB vs. LTBI #####
# ssgsea_PTB_LTBI_out_HIV <- lapply(objects_list_HIV_update, function(x) {
#     out <- subset_curatedTBData(x, annotationColName = "TBStatus",
#                                 annotationCondition = c("PTB", "LTBI"))
#     if (is.null(out)) {
#         return(NULL)
#     }
#     runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
#                      algorithm = "ssGSEA", update_genes = FALSE,
#                      combineSigAndAlgorithm = FALSE)
# }) |>
#     plyr::compact()
# saveRDS(ssgsea_PTB_LTBI_out_HIV,
#         file.path(wd, "data/ssgsea_PTB_LTBI_out_HIV.RDS"))

# ssgsea_PTB_LTBI_out_combine_HIV <- combine_auc(ssgsea_PTB_LTBI_out_HIV,
#                                                   annotationColName = "TBStatus",
#                                                   signatureColNames = names(signatures_list),
#                                                   num.boot = 1000, percent = 0.95)
# saveRDS(ssgsea_PTB_LTBI_out_combine_HIV,
#         file.path(wd, "data/ssgsea_PTB_LTBI_out_combine_HIV.RDS"))

ssgsea_PTB_LTBI_out_HIV <- file.path(wd, "data/ssgsea_PTB_LTBI_out_HIV.RDS") |> 
    readRDS()
study_PTB_LTBI_HIV <- lapply(1:length(ssgsea_PTB_LTBI_out_HIV), function(i) {
    sobject <- ssgsea_PTB_LTBI_out_HIV[[i]]
    sutdy_name <- gsub("_edit", "", names(ssgsea_PTB_LTBI_out_HIV)[i])
    data.frame(Study = sutdy_name, Size = ncol(sobject))
}) |> 
    dplyr::bind_rows()

ssgsea_PTB_LTBI_out_combine_HIV <- file.path(wd, "data/ssgsea_PTB_LTBI_out_combine_HIV.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_LTBI_HIV, "Study")

ssgsea_PTB_LTBI_CI_HIV <- extract_CI(ssgsea_PTB_LTBI_out_combine_HIV)
ssgsea_PTB_LTBI_CI_sub_HIV <- ssgsea_PTB_LTBI_CI_HIV |> 
    dplyr::mutate(ssgsea_PTB_LTBI = AUC) |> 
    dplyr::select(Signature, ssgsea_PTB_LTBI)

##### PLAGE PTB vs. LTBI #####
plage_PTB_LTBI_out_HIV <- lapply(objects_list_HIV_update, function(x) {
    out <- subset_curatedTBData(x, annotationColName = "TBStatus",
                                annotationCondition = c("PTB", "LTBI"))
    if (is.null(out)) {
        return(NULL)
    }
    runTBsigProfiler(input = out, useAssay = 1, signatures = signatures_list,
                     algorithm = "PLAGE", update_genes = FALSE,
                     combineSigAndAlgorithm = FALSE)
}) |>
    plyr::compact()
saveRDS(plage_PTB_LTBI_out_HIV,
        file.path(wd, "data/plage_PTB_LTBI_out_HIV.RDS"))

plage_PTB_LTBI_out_combine_HIV <- combine_auc(plage_PTB_LTBI_out_HIV,
                                               annotationColName = "TBStatus",
                                               signatureColNames = names(signatures_list),
                                               num.boot = 1000, percent = 0.95)
saveRDS(plage_PTB_LTBI_out_combine_HIV,
        file.path(wd, "data/plage_PTB_LTBI_out_combine_HIV.RDS"))

plage_PTB_LTBI_out_HIV <- file.path(wd, "data/plage_PTB_LTBI_out_HIV.RDS") |> 
    readRDS()
study_PTB_LTBI_HIV <- lapply(1:length(plage_PTB_LTBI_out_HIV), function(i) {
    sobject <- plage_PTB_LTBI_out_HIV[[i]]
    sutdy_name <- gsub("_edit", "", names(plage_PTB_LTBI_out_HIV)[i])
    data.frame(Study = sutdy_name, Size = ncol(sobject))
}) |> 
    dplyr::bind_rows()

plage_PTB_LTBI_out_combine_HIV <- file.path(wd, "data/plage_PTB_LTBI_out_combine_HIV.RDS") |> 
    readRDS() |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_LTBI_HIV, "Study")


plage_PTB_LTBI_CI_HIV <- extract_CI(plage_PTB_LTBI_out_combine_HIV)
plage_PTB_LTBI_CI_sub_HIV <- plage_PTB_LTBI_CI_HIV |> 
    dplyr::mutate(plage_PTB_LTBI = AUC) |> 
    dplyr::select(Signature, plage_PTB_LTBI)

##### Run analysis #####
# selected_gene_sets <- c("Gliddon_HIV_3", "Chen_HIV_4", "Rajan_HIV_5",
#                    "Sambarey_HIV_10", "Dawany_HIV_251")
# selected_gene_sets <- c("Blankley_5", "Kaforou_27", "Kaul_3", "Tabone_RES_27",
#                         "Dawany_HIV_251", "Lee_4", "Thompson_FAIL_13", "Verhagen_10", "Walter_PNA_47")

# ssgsea_PTB_Control_HIV <- ssgsea_PTB_Control_CI |> 
#     dplyr::filter(Signature %in% selected_gene_sets) |> 
#     dplyr::inner_join(ssgsea_PTB_Control_CI_HIV |> 
#                           dplyr::filter(Signature %in% selected_gene_sets), 
#                       by = "Signature")
# index_match <- match(selected_gene_sets, ssgsea_PTB_Control_HIV$Signature)
# 
# write.table(ssgsea_PTB_Control_HIV[index_match, ], 
#             quote = FALSE, row.names = FALSE, sep = " & ",
#             file = file.path(wd, "data/ssgsea_PTB_Control_HIV.txt")) 

# plage_PTB_Control_HIV <- plage_PTB_Control_CI |> 
#     dplyr::filter(Signature %in% selected_gene_sets) |> 
#     dplyr::inner_join(plage_PTB_Control_CI_HIV |> 
#                           dplyr::filter(Signature %in% selected_gene_sets), 
#                       by = "Signature")
# write.table(plage_PTB_Control_HIV[index_match, ], 
#             quote = FALSE, row.names = FALSE, sep = " & ",
#             file = file.path(wd, "data/plage_PTB_Control_HIV.txt")) 

ssgsea_PTB_Control_HIV_all <- ssgsea_PTB_Control_out_combine_HIV |> 
    dplyr::mutate(Group = "HIV subjects") |> 
    rbind(ssgsea_PTB_Control_out_combine |> 
              dplyr::mutate(Group = "Non-HIV subjects")) |>
    dplyr::mutate(Method = "ssGSEA")
plage_PTB_Control_HIV_all <- plage_PTB_Control_out_combine_HIV |> 
    dplyr::mutate(Group = "HIV subjects") |> 
    rbind(plage_PTB_Control_out_combine |> 
              dplyr::mutate(Group = "Non-HIV subjects")) |>
    dplyr::mutate(Method = "PLAGE")

#    dplyr::filter(Signature %in% selected_gene_sets) |> 
p_PTB_Control_HIV <- rbind(ssgsea_PTB_Control_HIV_all, plage_PTB_Control_HIV_all) |> 
    ggplot(aes(x = AUC, y = Signature, fill = Group)) +
        geom_density_ridges(alpha=0.6) +
        theme_ridges() + ylab(NULL) +
        facet_wrap(Method ~. ) +
        theme(legend.title = element_blank()
        )
library(ggridges)
library(ggplot2)
ggsave(p_PTB_Control_HIV, 
       file = file.path(wd, "Figures_and_tables/PTB_Control_HIV.pdf"),
       device = "pdf", height = 18, width = 8)
# ssgsea_PTB_LTBI_HIV <- ssgsea_PTB_LTBI_CI |> 
#     dplyr::filter(Signature %in% selected_gene_sets) |> 
#     dplyr::inner_join(ssgsea_PTB_LTBI_CI_HIV |> 
#                           dplyr::filter(Signature %in% selected_gene_sets), 
#                       by = "Signature")
# write.table(ssgsea_PTB_LTBI_HIV[index_match, ], 
#             quote = FALSE, row.names = FALSE, sep = " & ",
#             file = file.path(wd, "data/ssgsea_PTB_LTBI_HIV.txt")) 
# 
# plage_PTB_LTBI_HIV <- plage_PTB_LTBI_CI |> 
#     dplyr::filter(Signature %in% selected_gene_sets) |> 
#     dplyr::inner_join(plage_PTB_LTBI_CI_HIV |> 
#                           dplyr::filter(Signature %in% selected_gene_sets), 
#                       by = "Signature")
# write.table(plage_PTB_LTBI_HIV[index_match, ], 
#             quote = FALSE, row.names = FALSE, sep = " & ",
#             file = file.path(wd, "data/plage_PTB_LTBI_HIV.txt"))

#### Additional analysis for diabetes ####
# objects_list_dia <- lapply(objects_list, function(x) {
#     if ("DiabetesStatus" %in% colnames(colData(x))) {
#         x_sub <- x[, x$DiabetesStatus == "Positive"]
#         if (ncol(x_sub[["assay_curated"]]) > 0) {
#             return(x_sub)
#         }
#     } 
#     return(NULL)
# }) |> 
#     plyr::compact()
objects_list_dia <- curatedTBData(c("GSEBruno", "GSE73408"), dry.run = FALSE)
GSEBruno <- SummarizedExperiment(list(objects_list_dia$GSEBruno[["assay_curated"]]),
                                 colData = colData(objects_list_dia$GSEBruno))

pp <- paste0(GSEBruno$TBStatus, "_",GSEBruno$DiabetesStatus)
for (i in 1:length(pp)) {
    p <- pp[i]
    if (p == "Control_Negative") {
        pp[i] <- "HCs"
    } else if (p == "OD_Positive") {
        pp[i] <- "Diabetes"
    } else if (p == "PTB_Negative") {
        pp[i] <- "PTB"
    } else {
        pp[i] <- "PTB+Diabetes"
    }
}
GSEBruno$pp <- pp

ssgsea_GSEBruno <- runTBsigProfiler(input = GSEBruno, useAssay = 1, 
                             signatures = signatures_list,
                     algorithm = "ssGSEA", update_genes = FALSE,
                     combineSigAndAlgorithm = FALSE)
common_sigs <- intersect(names(signatures_list), 
                         colnames(colData(ssgsea_GSEBruno)))
# p_dia_ssgsea <- signatureBoxplot_edit(inputData = ssgsea_GSEBruno,
#                  name = NULL,
#                  signatureColNames = common_sigs,
#                  annotationColName = "pp") +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust=0.95, vjust=1),
#           legend.title = element_blank())
# ggsave(p_dia_ssgsea, 
#        file = file.path(wd, "Figures_and_tables/p_dia_ssgsea.pdf"),
#        device = "pdf", height = 10, width = 17)

plage_GSEBruno <- runTBsigProfiler(input = GSEBruno, useAssay = 1, 
                                    signatures = signatures_list,
                                    algorithm = "PLAGE", update_genes = FALSE,
                                    combineSigAndAlgorithm = FALSE)
# Choose 9 as an example
sigs_example <- c("Bloom_OD_144", "Chen_HIV_4", "Darboe_RISK_11", 
                  "Hoang_OD_20", "Sweeney_OD_3", "Suliman_4",
                  "Zak_RISK_16", "Walter_51", "Zimmer_RES_3")
ssgsea_GSEBruno_long <- colData(ssgsea_GSEBruno) |> 
    as.data.frame() |> 
    dplyr::select(sigs_example, pp) |> 
    reshape2::melt() |> 
    dplyr::mutate(Method = "ssGSEA")
plage_GSEBruno_long <- colData(plage_GSEBruno) |> 
    as.data.frame() |> 
    dplyr::select(sigs_example, pp) |> 
    reshape2::melt() |> 
    dplyr::mutate(Method = "PLAGE")

GSEBruno_long <- rbind(ssgsea_GSEBruno_long, plage_GSEBruno_long)
GSEBruno_long$pp <- factor(GSEBruno_long$pp, 
                           levels = c("HCs", "PTB", "Diabetes", "PTB+Diabetes"))
p_dia_example <- ggplot(GSEBruno_long, aes(x = pp, y = value, fill = Method)) +
    geom_boxplot() +
    facet_wrap(~variable, scale = "free_y") +
    ylab("Scores") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust=0.95, vjust=1))
ggsave(p_dia_example,
       file = file.path(wd, "Figures_and_tables/p_dia_example.pdf"),
       device = "pdf", height = 6, width = 8)
# common_sigs <- intersect(names(signatures_list), colnames(colData(GSEBruno)))
# p_dia_plage <- signatureBoxplot_edit(inputData = plage_GSEBruno,
#                                  name = NULL,
#                                  signatureColNames = common_sigs,
#                                  annotationColName = "pp") +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust=0.95, vjust=1),
#           legend.title = element_blank())
# ggsave(p_dia_plage, 
#        file = file.path(wd, "Figures_and_tables/p_dia_plage.pdf"),
#        device = "pdf", height = 10, width = 17)

#### Additional analysis for treatment results ####
# objects_list_trt <- curatedTBData(c("GSE36238", "GSE42831", "GSE56153",
#                                     "GSE84076", "GSE89403"), dry.run = FALSE)

# objects_list_trt <- lapply(objects_list, function(x) {
#     if ("TreatmentResult" %in% colnames(colData(x))) {
#         return(x)
#         }
#     return(NULL)
# }) |> 
#     plyr::compact()
#### Table for gene signatures list ####
signatures_list <- TBSignatureProfiler::TBsignatures[sigatures_names] |> 
    lapply(update_genenames)
table_signatures <- lapply(names(signatures_list), function(sig_name) {
    sig_name_split <- strsplit(sig_name, "_") |> 
        unlist()
    n <- length(sig_name_split)
    if (n == 2) {
        compare <- "PTB vs. (HCs and/or LTBI)"
    } else if (n == 3) {
        type <- sig_name_split[2]
        if (type == "FAIL") {
            compare <- "Failure of TB treatment"
        } else if (type == "HIV") {
            compare <- "HIV co-infection for PTB vs. (LTBI and/or OD)"
        } else if (type == "NANO") {
            compare <- "PTB vs. (HCs and/or LTBI) (NanoString
nCounter Platform)"
        } else if (type == "OD") {
            compare <- "PTB vs. OD"
        } else if (type == "PNA") {
            compare <- "PTB vs. Pneumonia"
        } else if (type == "RES") {
            compare <- "Response to TB Treatment"
        } else if (type == "RISK") {
            compare <- "Progressor vs. Non-progressor"
        }
    } else if (n == 4) { # "Gliddon_2_OD_4" 
        compare <- "PTB vs. OD"
    }
    n_genes <- sig_name_split[length(sig_name_split)]
    genes <- paste0(signatures_list[[sig_name]], collapse = ", ")
    data.frame(Signatures = sig_name, 
               Comparison = compare,
               `Number of genes` = n_genes,
               `Genes` = genes)
}) |> 
    dplyr::bind_rows()

library(xlsx)

write.xlsx(table_signatures, 
           file.path(wd, "Supplemental_files/Supplemental_file_1_gene_lists.xlsx"),
           row.names = FALSE)

