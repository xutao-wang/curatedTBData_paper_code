# This is script contains code to generate figures and tables in the curatedTBData paper
source("~/Desktop/practice/curatedTBData_paper_results/Scripts/functions.R")
source("~/Desktop/Packages/curatedTBData/data-raw/UtilityFunctionForCuration.R")
library(cowplot)
library(dplyr)
library(ggplot2)
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
                     "Kaforou_OD_53", "Kaul_3", "Kulkarni_HIV_2", "Kwan_186",
                     "LauxdaCosta_OD_3", "Lee_4", "Leong_24", "Leong_RISK_29",
                     "Long_RES_10", "Maertzdorf_15", "Maertzdorf_4", "Maertzdorf_OD_100",
                     "Natarajan_7", "PennNich_RISK_6", "Qian_OD_17", "Rajan_HIV_5",
                     "Roe_3", "Roe_OD_4", "Sambarey_HIV_10", "Singhania_OD_20",
                     "Sivakumaran_11", "Suliman_4", "Suliman_RISK_4", "Sweeney_OD_3", 
                     "Tabone_OD_11", "Tabone_RES_25", "Tabone_RES_27", "Thompson_9", 
                     "Thompson_FAIL_13", "Thompson_RES_5", "Tornheim_71", "Tornheim_RES_25", 
                     "Verhagen_10", "Walter_51", "Walter_PNA_119", "Walter_PNA_47", 
                     "Zak_RISK_16", "Zhao_NANO_6", "Zimmer_RES_3")
signatures_list <- TBSignatureProfiler::TBsignatures[sigatures_names]
signatures_list <- lapply(signatures_list, update_genenames)
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
# 
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

ssgsea_PTB_Control_out_combine <- readRDS(file.path(wd, "data/ssgsea_PTB_Control_out_combine.RDS")) |> 
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

plage_PTB_Control_out_combine <- readRDS(file.path(wd, "data/plage_PTB_Control_out_combine.RDS")) |> 
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
# 
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

ssgsea_PTB_LTBI_out_combine <- readRDS(file.path(wd, "data/ssgsea_PTB_LTBI_out_combine.RDS")) |> 
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
# 
# plage_PTB_LTBI_out_combine <- combine_auc(plage_PTB_LTBI_out,
#                                            annotationColName = "TBStatus",
#                                            signatureColNames = names(signatures_list),
#                                            num.boot = 1000, percent = 0.95)
# saveRDS(plage_PTB_LTBI_out_combine, file.path(wd, "data/plage_PTB_LTBI_out_combine.RDS"))

plage_PTB_LTBI_out <- readRDS(file.path(wd, "data/plage_PTB_LTBI_out.RDS"))

plage_PTB_LTBI_out_combine <- readRDS(file.path(wd, "data/plage_PTB_LTBI_out_combine.RDS")) |> 
    dplyr::mutate(Study = gsub("_edit", "", Study)) |> 
    dplyr::inner_join(study_PTB_LTBI, "Study")

#### Table S1 Datasets in the curatedTBData ####
data("DataSummary")
study_summary_sub <- DataSummary |> 
    dplyr::select(-Notes, -GeneralType)
write.table(study_summary_sub, quote = FALSE, row.names = FALSE, sep = " & ",
            file = "~/Desktop/practice/curatedTBData_paper_results/Figures_and_tables/study_summary_sub.txt") 
                                   
# Run the following command in terminal (Add '\\\hline' by the end of each line)
# awk '{print $0, " \\\\\\\hline"}' study_summary_sub.txt > study_summary_sub_for_latex.txt

#### Table S2 Clinical annotation ####
clincial_annotation <- readxl::read_excel("Desktop/practice/curatedTBData_paper_results/clinicalDataAnnotation.xlsx")
write.table(clincial_annotation, 
            quote = FALSE, row.names = FALSE, sep = " & ",
            file = file.path(out_file_path, "clincial_annotation.txt")) 
# Run the following command in terminal (Add '\\\hline' by the end of each line)
# awk '{print $0, " \\\\\\\hline"}' clincial_annotation.txt > clincial_annotation_for_latex.txt

#### Table S3: Weighted Mean AUC and 95% confidence interval for PTB vs. Control or PTB vs. LTBI ####
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

tableS3_mean_AUC <- ssgsea_PTB_Control_CI_sub |> 
    dplyr::inner_join(plage_PTB_Control_CI_sub, by = "Signature") |> 
    dplyr::inner_join(ssgsea_PTB_LTBI_CI_sub, by = "Signature") |> 
    dplyr::inner_join(plage_PTB_LTBI_CI_sub, by = "Signature")

tableS3_mean_AUC$Signature <- gsub("_", "\\\\_", tableS3_mean_AUC$Signature)
write.table(tableS3_mean_AUC, quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = " & ",
            file = file.path(wd, "Figures_and_tables/tableS3_mean_AUC.txt")) 

# Run the following command in terminal (Add '\\\hline' by the end of each line)
# awk '{print $0, " \\\\\\\hline"}' tableS3_mean_AUC.txt > tableS3_mean_AUC_for_latex.txt

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
out <- cbind(Signature = mean_AUC_combine$Signature, df)

venn_input <-list('PTB vs. Control (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_Control >= 0.8]),
                  'PTB vs. Control (PLAGE)' = na.omit(out$Signature[out$plage_PTB_Control >= 0.8]),
                  'PTB vs. LTBI (ssGSEA)' = na.omit(out$Signature[out$ssgsea_PTB_LTBI >= 0.8]),
                  'PTB vs. LTBI (PLAGE)' = na.omit(out$Signature[out$plage_PTB_LTBI >= 0.8]))

# create venn diagram and display all the sets
figure2 <- ggvenn::ggvenn(venn_input, show_percentage= T, show_elements = FALSE, set_name_size = 4,
               fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
               stroke_size = 0.1)
# ggsave(file.path(wd, "Figures_and_tables/figure2_venn_diagram.pdf"), device = "pdf",
#        height = 10, width = 10)
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
# See ensemble_results.R and extract_ensemble_results.R for data preparation
##### ssGSEA PTB vs. Control #####
ssgsea_PTB_Control_out_edit <- PTB_Control_edit(ssgsea_PTB_Control_out, 
                                                filter_size = 15,
                                                save_file = FALSE)
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
plage_PTB_Control_out_edit <- PTB_Control_edit(plage_PTB_Control_out,
                                               filter_size = 15,
                                               save_file = FALSE)
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
                  aes(x = Method, y = 0.82, label = Signature), size = 3, angle = 10) +
        theme_bw() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(size = 8, face = "bold"))

# ggsave(file.path(wd, "Figures_and_tables/FigureS3_PTB_Control_ensemble.pdf"),
#        p_PTB_Control_ensemble, device = "pdf", width = 14, height = 10)

##### ssGSEA PTB vs. LTBI #####
ssgsea_PTB_LTBI_out_edit <- PTB_LTBI_edit(ssgsea_PTB_LTBI_out,
                                          filter_size = 15,
                                          save_file = FALSE)
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
plage_PTB_LTBI_out_edit <- PTB_LTBI_edit(plage_PTB_LTBI_out,
                                         filter_size = 15,
                                         save_file = FALSE)
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

##### Figure 3 Plot: Sub ensemble #####
my_set <- c("Set 6", "Set 9", "Set 21", "Set 22")

AUC_final_PTB_Control <- PTB_Control_ensl$AUC_final |> 
    dplyr::filter(Set %in% my_set) |> 
    dplyr::mutate(Comparison = "PTB vs. Control")
signature_with_max_AUC_PTB_Control <- PTB_Control_ensl$signature_with_max_AUC |> 
    dplyr::filter(Set %in% my_set) |> 
    dplyr::mutate(Comparison = "PTB vs. Control")

AUC_final_PTB_LTBI <- PTB_LTBI_ensl$AUC_final |> 
    dplyr::filter(Set %in% my_set) |> 
    dplyr::mutate(Comparison = "PTB vs. LTBI")
signature_with_max_AUC_PTB_LTBI <- PTB_LTBI_ensl$signature_with_max_AUC |> 
    dplyr::filter(Set %in% my_set) |> 
    dplyr::mutate(Comparison = "PTB vs. LTBI")

AUC_final_sub <- rbind(AUC_final_PTB_Control, AUC_final_PTB_LTBI)
signature_with_max_AUC_sub <- rbind(signature_with_max_AUC_PTB_Control,
                                    signature_with_max_AUC_PTB_LTBI)

p_ensl_sub <- ggplot(AUC_final_sub, 
       aes(Method, AUC, fill = Signature_type)) +
    labs(fill = "Signature type") +
    geom_split_violin() + 
    facet_grid(rows = vars(Comparison), cols = vars(Set), switch = "y") +
    geom_text(data = signature_with_max_AUC_sub,
              aes(x = Method, y = 0.83, label = Signature), size = 3, angle = 10) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 8, face = "bold"))

# ggsave(file.path(wd, "Figures_and_tables/Figure3_ensemble_sub.pdf"),
#        p_ensl_sub, device = "pdf", width = 14, height = 6)

##### p values for comparisons #####
###### Set 6 ######
# ssGSEA
AUC_set6_PTB_Control_ssgsea <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set6", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
mean(AUC_set6_PTB_Control_ssgsea$AUC)
stats::quantile(AUC_set6_PTB_Control_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
stats::quantile(AUC_set6_PTB_Control_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

AUC_set6_PTB_Control_ssgsea_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Bloom_OD_144", 
                  Set == "Set 6",
                  Method == "ssGSEA")

wilcox.test(AUC_set6_PTB_Control_ssgsea$AUC, 
            AUC_set6_PTB_Control_ssgsea_single$AUC)

AUC_set6_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set6", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set6_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Bloom_OD_144", 
                  Set == "Set 6",
                  Method == "ssGSEA")
mean(AUC_set6_PTB_LTBI_ssgsea$AUC)
stats::quantile(AUC_set6_PTB_LTBI_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
stats::quantile(AUC_set6_PTB_LTBI_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

wilcox.test(AUC_set6_PTB_LTBI_ssgsea$AUC, 
            AUC_set6_PTB_LTBI_ssgsea_single$AUC)

# PLAGE
AUC_set6_PTB_Control_PLAGE <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set6", 
                  Signature_type == "Ensemble",
                  Method == "PLAGE")
# mean(AUC_set6_PTB_Control_PLAGE$AUC)
# stats::quantile(AUC_set6_PTB_Control_PLAGE$AUC, prob = 0.025, na.rm = TRUE)
# stats::quantile(AUC_set6_PTB_Control_PLAGE$AUC, prob = 0.975, na.rm = TRUE)

AUC_set6_PTB_Control_PLAGE_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Roe_3", 
                  Set == "Set 6",
                  Method == "PLAGE")

wilcox.test(AUC_set6_PTB_Control_PLAGE$AUC, 
            AUC_set6_PTB_Control_PLAGE_single$AUC)

AUC_set6_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set6", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set6_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Bloom_OD_144", 
                  Set == "Set 6",
                  Method == "ssGSEA")
mean(AUC_set6_PTB_LTBI_ssgsea$AUC)
stats::quantile(AUC_set6_PTB_LTBI_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
stats::quantile(AUC_set6_PTB_LTBI_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

wilcox.test(AUC_set6_PTB_LTBI_ssgsea$AUC, 
            AUC_set6_PTB_LTBI_ssgsea_single$AUC)

###### Set 9 ######
# ssGSEA
AUC_set9_PTB_Control_ssgsea <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set9", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
# mean(AUC_set9_PTB_Control_ssgsea$AUC)
# stats::quantile(AUC_set9_PTB_Control_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
# stats::quantile(AUC_set9_PTB_Control_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

AUC_set9_PTB_Control_ssgsea_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 9",
                  Method == "ssGSEA")

wilcox.test(AUC_set9_PTB_Control_ssgsea$AUC, 
            AUC_set9_PTB_Control_ssgsea_single$AUC) # p value 0.1187

AUC_set9_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set9", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set9_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 9",
                  Method == "ssGSEA")
# mean(AUC_set9_PTB_LTBI_ssgsea$AUC)
# stats::quantile(AUC_set6_PTB_LTBI_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
# stats::quantile(AUC_set6_PTB_LTBI_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

wilcox.test(AUC_set9_PTB_LTBI_ssgsea$AUC, 
            AUC_set9_PTB_LTBI_ssgsea_single$AUC)

# PLAGE
AUC_set9_PTB_Control_PLAGE <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set9", 
                  Signature_type == "Ensemble",
                  Method == "PLAGE")
# mean(AUC_set6_PTB_Control_PLAGE$AUC)
# stats::quantile(AUC_set6_PTB_Control_PLAGE$AUC, prob = 0.025, na.rm = TRUE)
# stats::quantile(AUC_set6_PTB_Control_PLAGE$AUC, prob = 0.975, na.rm = TRUE)

AUC_set9_PTB_Control_PLAGE_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 9",
                  Method == "PLAGE")

wilcox.test(AUC_set9_PTB_Control_PLAGE$AUC, 
            AUC_set9_PTB_Control_PLAGE_single$AUC)

AUC_set9_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set9", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set9_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 9",
                  Method == "ssGSEA")
mean(AUC_set9_PTB_LTBI_ssgsea$AUC)
stats::quantile(AUC_set9_PTB_LTBI_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
stats::quantile(AUC_set9_PTB_LTBI_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

wilcox.test(AUC_set9_PTB_LTBI_ssgsea$AUC, 
            AUC_set9_PTB_LTBI_ssgsea_single$AUC)


###### Set 21 ######
# ssGSEA
AUC_set21_PTB_Control_ssgsea <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set21", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")

AUC_set21_PTB_Control_ssgsea_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Tornheim_71", 
                  Set == "Set 21",
                  Method == "ssGSEA")

wilcox.test(AUC_set21_PTB_Control_ssgsea$AUC, 
            AUC_set21_PTB_Control_ssgsea_single$AUC) # p value 0.1187

AUC_set21_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set21", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set21_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Bloom_OD_144", 
                  Set == "Set 21",
                  Method == "ssGSEA")

wilcox.test(AUC_set21_PTB_LTBI_ssgsea$AUC, 
            AUC_set21_PTB_LTBI_ssgsea_single$AUC)

# PLAGE
AUC_set21_PTB_Control_PLAGE <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set21", 
                  Signature_type == "Ensemble",
                  Method == "PLAGE")

AUC_set21_PTB_Control_PLAGE_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "LauxdaCosta_OD_3", 
                  Set == "Set 21",
                  Method == "PLAGE")

wilcox.test(AUC_set21_PTB_Control_PLAGE$AUC, 
            AUC_set21_PTB_Control_PLAGE_single$AUC)

AUC_set21_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set21", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set21_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Anderson_OD_51", 
                  Set == "Set 21",
                  Method == "ssGSEA")
mean(AUC_set21_PTB_LTBI_ssgsea$AUC)
stats::quantile(AUC_set21_PTB_LTBI_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
stats::quantile(AUC_set21_PTB_LTBI_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

wilcox.test(AUC_set21_PTB_LTBI_ssgsea$AUC, 
            AUC_set21_PTB_LTBI_ssgsea_single$AUC)


###### Set 22 ######
# ssGSEA
AUC_set22_PTB_Control_ssgsea <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set22", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")

AUC_set22_PTB_Control_ssgsea_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 22",
                  Method == "ssGSEA")

wilcox.test(AUC_set22_PTB_Control_ssgsea$AUC, 
            AUC_set22_PTB_Control_ssgsea_single$AUC) # p value 0.1187

AUC_set22_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set22", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set22_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 22",
                  Method == "ssGSEA")

wilcox.test(AUC_set22_PTB_LTBI_ssgsea$AUC, 
            AUC_set22_PTB_LTBI_ssgsea_single$AUC)

# PLAGE
AUC_set22_PTB_Control_PLAGE <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Set22", 
                  Signature_type == "Ensemble",
                  Method == "PLAGE")

AUC_set22_PTB_Control_PLAGE_single <- AUC_final_PTB_Control |> 
    dplyr::filter(Signature == "Zak_RISK_16", 
                  Set == "Set 22",
                  Method == "PLAGE")

wilcox.test(AUC_set22_PTB_Control_PLAGE$AUC, 
            AUC_set22_PTB_Control_PLAGE_single$AUC)

AUC_set22_PTB_LTBI_ssgsea <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Set22", 
                  Signature_type == "Ensemble",
                  Method == "ssGSEA")
AUC_set22_PTB_LTBI_ssgsea_single <- AUC_final_PTB_LTBI |> 
    dplyr::filter(Signature == "Berry_393", 
                  Set == "Set 22",
                  Method == "ssGSEA")
mean(AUC_set22_PTB_LTBI_ssgsea$AUC)
stats::quantile(AUC_set22_PTB_LTBI_ssgsea$AUC, prob = 0.025, na.rm = TRUE)
stats::quantile(AUC_set22_PTB_LTBI_ssgsea$AUC, prob = 0.975, na.rm = TRUE)

wilcox.test(AUC_set22_PTB_LTBI_ssgsea$AUC, 
            AUC_set22_PTB_LTBI_ssgsea_single$AUC)


#### Figure S1 PTB vs. Control ####
##### Figure S1A heatmap AUC ssGSEA PTB vs. Control #####
heatmap_ssgsea_PTB_Control <- heatmap_auc(ssgsea_PTB_Control_out_combine, 
                                          GSE_sig = NULL, facet = TRUE, 
                                          clustering = FALSE) + 
    theme(legend.position="bottom")
# ggsave(heatmap_ssgsea_PTB_Control, 
#        file = file.path(wd, "Figures_and_tables/figS1A_heatmap_ssgsea_PTB_Control.pdf"),
#        device = "pdf", height = 18, width = 12)


##### Figure S1B heatmap AUC PLAGE PTB vs. Control #####
heatmap_plage_PTB_Control <- heatmap_auc(plage_PTB_Control_out_combine, GSE_sig = NULL, 
                                         facet = TRUE, clustering = FALSE) +
    theme(legend.position="bottom")
# ggsave(heatmap_plage_PTB_Control, 
#        file = file.path(wd, "Figures_and_tables/figS1B_heatmap_plage_PTB_Control.pdf"),
#        device = "pdf", height = 18, width = 12)

##### Figure S1 #####
figureS1 <- cowplot::plot_grid(heatmap_ssgsea_PTB_Control, 
                               heatmap_plage_PTB_Control,
                               align = "h")
ggsave(figureS1, 
       file = file.path(wd, "Figures_and_tables/FigureS1_PTB_Control.pdf"),
       device = "pdf", height = 18, width = 24)

#### Figure S2 PTB vs. LTBI ####
##### Figure S2A heatmap AUC ssGSEA PTB vs. LTBI #####
heatmap_ssgsea_PTB_LTBI <- heatmap_auc(ssgsea_PTB_LTBI_out_combine, 
                                       GSE_sig = NULL, facet = TRUE, 
                                       clustering = FALSE) +
    theme(legend.position="bottom")
# ggsave(heatmap_ssgsea_PTB_LTBI, 
#        file = file.path(wd, "Figures_and_tables/figS2A_heatmap_ssgsea_PTB_LTBI.pdf"),
#        device = "pdf", height = 18, width = 10)

##### Figure S2B heatmap AUC PLAGE PTB vs. LTBI #####
heatmap_plage_PTB_LTBI <- heatmap_auc(plage_PTB_LTBI_out_combine, GSE_sig = NULL, 
                                      facet = TRUE, clustering = FALSE) +
    theme(legend.position="bottom")
# ggsave(heatmap_plage_PTB_LTBI, 
#        file = file.path(wd, "Figures_and_tables/figS2B_heatmap_plage_PTB_LTBI.pdf"),
#        device = "pdf", height = 18, width = 10)
##### Figure S2 #####
figureS2 <- cowplot::plot_grid(heatmap_ssgsea_PTB_LTBI, 
                               heatmap_plage_PTB_LTBI,
                               align = "h")
ggsave(figureS2, 
       file = file.path(wd, "Figures_and_tables/FigureS2_PTB_LTBI.pdf"),
       device = "pdf", height = 18, width = 20)



