rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/results/8.0.cluster_degs/")
dir.create(od)
pacman::p_load(tidyverse,GSVA,enrichplot,SummarizedExperiment)

set.seed(1110)
sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))
load(str_glue("{run_home}/data/expr_data.RData"))

x <- "Cluster4"
pdata <- sample_cluster %>% dplyr::rename(sample = "Sample") %>% 
            mutate(group = ifelse(Group == x,x,"Others")) %>% dplyr::select(sample,group)
degs <- limma_deg(DEG_exp = expr_dat, DEG_pdata = pdata, od = str_glue("{od}/{x}/"),
                    controlLabel = "Others",caseLabel = x,color_fun = c("#80253F","#015694", "#e6e1e1"),
                    DEG_FC = 0.585, DEG_P = 0.05, prefix = x)
save(degs,file = str_glue("{od}/{x}/degs.Rdata"))

load(str_glue("{run_home}/data/raw_expr.Rdata"))
rowdata <- rowData(all_dat)
mrna_dat <- all_dat[rowdata$gene_type == "protein_coding",]
symbol_mrna <- rowData(mrna_dat)$gene_name
expr_dat <- assay(mrna_dat, i = "fpkm_unstrand") %>% as.data.frame()
expr_dat <- cbind(data.frame(symbol_mrna),
                  as.data.frame(expr_dat))
expr_dat <- expr_dat %>% as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before=2) %>% 
  arrange(desc(meanrow)) %>% 
  distinct(symbol_mrna,.keep_all=T) %>% 
  select(-meanrow) %>% 
  column_to_rownames(var = "symbol_mrna") %>% 
  as.data.frame()                  
expr_dat <- log2(expr_dat + 1) %>% as.data.frame() %>% dplyr::select(grep("^.{13}01A.*$|^.{13}11A.*$", colnames(.)))

colnames(expr_dat) <- lapply(colnames(expr_dat), function(x) {substring(x, 1, 16)}) %>% unlist()
expr_dat <- expr_dat %>% select(-colnames(expr_dat)[duplicated(colnames(expr_dat))])
save(expr_dat,file = str_glue("{od}/all_expr_data.RData"))

pdata <- data.frame(sample = colnames(expr_dat),group = ifelse(substr(colnames(expr_dat),14,16) == "01A","Tumor","Control"))
all_degs <- limma_deg(DEG_exp = expr_dat, DEG_pdata = pdata, od = str_glue("{od}/tumor_control/"),
                    controlLabel = "Control",caseLabel = "Tumor",color_fun = c("#80253F","#015694", "#e6e1e1"),
                    DEG_FC = 0.585, DEG_P = 0.05, prefix = "ALL")
save(all_degs,file = str_glue("{od}/tumor_control/degs.Rdata"))

load(str_glue("{od}/Cluster4/degs.Rdata"))
load(str_glue("{od}/tumor_control/degs.Rdata"))
venn.plot <- VennDiagram::venn.diagram(
    x = list(
        subtype_deg = degs$DEGs,tumor_deg = all_degs$DEGs
    ),
    cat.pos = c(-135, 135),
    filename = NULL,
    col = "white",scaled =F,
    fill = c( "#2590CF","#F2737F"),
    #lty = "blank",
    alpha = 0.50,
    cat.cex = 1,cat.dist = 0.1,
    cat.fontface = "bold",
    margin = 0.1
)
pdf(file = str_glue("{od}/degs_venn.pdf"), height = 4, width = 4)
grid::grid.draw(venn.plot)
dev.off()


model_gene <- intersect(degs$DEGs,all_degs$DEGs)
write.csv(model_gene,file = str_glue("{od}/model_gene.csv"))

load(str_glue("{run_home}/data/expr_data.RData"))
data_survival <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_survival.csv"))
identical(colnames(expr_dat),data_survival$sample)
cox_res <- signature_cox(
    signaturelist = model_gene, exp = expr_dat, clin = data_survival, coxp = 0.01,
    timecol = "OS.time", statuscol = "OS",
    bygroup = TRUE, xlab = str_c("OS days"), savekmplot = TRUE
)
save(cox_res, file = str_glue("{od}/cox_res.RData"))
write.csv(cox_res$coxResult, file = str_glue("{od}/SupplementaryTable_coxResult.csv"), quote = FALSE, row.names = T)
write.csv(cox_res$sigcoxResult, file = str_glue("{od}/SupplementaryTable_sigcoxResult.csv"), quote = FALSE, row.names = T)
write.csv(cox_res$cox_signature, file = str_glue("{od}/SupplementaryTable_cox_signature.csv"), quote = FALSE, row.names = FALSE)

forest_plot(
    od = od, input = cox_res$sigcoxResult %>% rownames_to_column("gene") %>% arrange(pvalue), plotN = 11,
    dataset = "TCGA", h = 3.5, w = 8,
    signaturecol = "gene", pvaluecol = "pvalue", HRcol = "HR", lower95col = "Low 95%CI",
    upper95col = "High 95%CI"
)

signature_coef <- lasso_model(
        signaturelist = cox_res$cox_signature,
        od = od,timecol = "OS.time", statuscol = "OS",
        exp = expr_dat, clin = data_survival,
        signaturetype = "cox_gene",seed = 1110)
save(signature_coef,file = str_glue("{od}/lasso_signature_coef.RData"))

rm(list=ls())
run_home <- "./"
options(bitmapType = "cairo")
od <- file.path(run_home, "/results/8.0.cluster_degs/")
dir.create(od)
pacman::p_load(tidyverse)
load(str_glue("{od}/lasso_signature_coef.RData"))
model_res <- datacenter_validation_v2(
    signature_coef = signature_coef, od = str_glue("{od}/lasso_verify"), cancer = "LUSC",
    outcome_list = c("OS"),
    version = "v1", time_points = c(1, 3, 5), best_cut = T, prognostic_independence = FALSE,
    factors = c("Stage", "Age", "Gender"), model = "lasso", cohort_list = NULL, minprop = 0.2,
    save_exp_clinical_stat = F, surtime_unit = 365, alias_substitution = TRUE, coxp = 0.05,
    color_fun = c("#2590CF","#F2737F"), data_center = "/Pub/Users/database/", other_plot = T
)

rm(list=ls())
run_home <- "./"
od <- file.path(run_home, "/results/8.0.cluster_degs/")
dir.create(od)
pacman::p_load(tidyverse,survival,survminer)
load(str_glue("{od}/cox_res.RData"))
load(str_glue("{run_home}/data/expr_data.RData"))
data_survival <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_survival.csv"))
infor <- expr_dat[cox_res$cox_signature, ] %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(., "sample") %>%
        merge(.,data_survival, by = "sample") %>% dplyr::rename(time = "OS.time", status = "OS")
multi_cox_model <- as.formula(paste0("Surv(time, status) ~", str_c("`", cox_res$cox_signature, "`", collapse = "+"))) %>% coxph(data = infor)
signature_coef <- summary(multi_cox_model)$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("gene")
save(signature_coef,file = str_glue("{od}/multi_cox_signature_coef.RData"))
write.csv(signature_coef,file = str_glue("{od}/multi_cox_signature_coef.csv"))


train_model_res <- modelscore_km_roc_v2(
    signature_coef = signature_coef,
    od = str_glue("{od}/train/"),
    time_points = c(1, 3, 5),
    no_roc_pheatmap = F,
    best_cut = F,
    minprop = 0.4,
    exp = infor[,1:12] %>% column_to_rownames(var = "sample") %>% t() %>% as.data.frame(),
    clin = infor[,c(1,38,39)],timecol = "time",statuscol = "status",
    dataset = "TCGA",
    xlab = str_c("OS", " days")
)
train_model_res$infor
save(train_model_res, file = str_glue("{od}/train_model_res.RData"))

rm(list=ls())
run_home <- "./"
libSources <- list.files("/Pub/Users/cuiye/RCodes/ProjectCode/", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources, source))
options(bitmapType = "cairo")
od <- file.path(run_home, "/results/8.0.cluster_degs/")
dir.create(od)
pacman::p_load(tidyverse)
load(str_glue("{od}/multi_cox_signature_coef.RData"))
model_res <- datacenter_validation_v2(
    signature_coef = signature_coef, od = str_glue("{od}/multi_cox_verify"), cancer = "LUSC",
    outcome_list = c("OS"),
    version = "v1", time_points = c(1, 3, 5), best_cut = T, prognostic_independence = FALSE,
    factors = c("Stage", "Age", "Gender"), model = "multi_cox", cohort_list = NULL, minprop = 0.2,
    save_exp_clinical_stat = F, surtime_unit = 365, alias_substitution = TRUE, coxp = 0.05,
    color_fun = c("#2590CF","#F2737F"), data_center = "/Pub/Users/database/", other_plot = T
)

