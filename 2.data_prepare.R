rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/results/1.0.data_prepare/")
dir.create(od)
pacman::p_load(tidyverse,GSVA,ComplexHeatmap,RColorBrewer,IOBR,factoextra,FactoMineR,circlize)

set.seed(1110)
load(str_glue("{run_home}/data/expr_data.RData"))

data_survival <- readxl::read_excel(str_glue("{run_home}/data/mmc1.xlsx"),1) %>% .[,-1] %>% filter(type == "LUSC") %>% 
                  filter(bcr_patient_barcode %in% colnames(expr_dat))
data_survival <- clinical_pro_v2(input_clin = data_survival, samplecol = "bcr_patient_barcode", timecol = NULL, statuscol = NULL, 
                        stagecol = "ajcc_pathologic_tumor_stage", agecol = "age_at_initial_pathologic_diagnosis",
                            gendercol = "gender", gradecol = "histological_type")
load(str_glue("{run_home}/data/data_clinical.RData"))

write.csv(data_survival,file = str_glue("{od}/lusc_survival.csv"))
write.csv(expr_dat,file = str_glue("{od}/lusc_expr.csv"))
write.table(expr_dat, file = str_glue("{od}/lusc_expr.txt"), sep = "\t",quote = F, col.names = TRUE, row.names = T)

expr_tide <- sweep(expr_dat, 1, apply(expr_dat,1,mean,na.rm = T))
write.table(expr_tide, file = str_glue("{od}/lusc_expr_tide.txt"), sep = "\t",quote = F, col.names = TRUE, row.names = T)

fges_info <- readxl::read_xlsx(str_glue("{run_home}/data/29fges_Cancer_cell.xlsx"),sheet = 1,skip = 1) %>% 
            arrange(ID,"Process_Signature")
ssgsea_geneSets <- list()
ssgsea_geneSets <- map(1:nrow(fges_info),function(x){
    fges_info[x,4:ncol(fges_info)] %>% c() %>% unlist() %>% na.omit() %>% as.character()
})
names(ssgsea_geneSets) <- fges_info$`Process_Signature` %>% as.character()

fges_res <- gsva(expr = as.matrix(expr_dat), gset.idx.list = ssgsea_geneSets, kcdf = "Gaussian", method = "ssgsea",
                abs.ranking = FALSE, min.sz = 2, ssgsea.norm = TRUE) %>% t() %>%
                as.data.frame() %>% rownames_to_column(var = "sample")
write.csv(fges_res,file = str_glue("{od}/lusc_fges.csv"))

geneSets <- GSEABase::getGmt(str_glue("{run_home}/data/h.all.v2024.1.Hs.symbols.gmt"))
hallmark_res <- gsva(expr = as.matrix(expr_dat), gset.idx.list = geneSets, kcdf = "Gaussian", method = "ssgsea",
                abs.ranking = FALSE, min.sz = 2, ssgsea.norm = TRUE) %>% t() %>%
                as.data.frame() %>% rownames_to_column(var = "sample")
write.csv(hallmark_res,file = str_glue("{od}/lusc_hallmark.csv"))

estimate_res <- deconvo_tme(eset = as.data.frame(expr_dat), method = "estimate") %>% dplyr::rename(sample = ID)
write.csv(estimate_res,file = str_glue("{od}/lusc_estimate.csv"))

mcpcounter_res <- deconvo_tme(eset = as.data.frame(expr_dat), method = "mcpcounter") %>% dplyr::rename(sample = ID)
write.csv(mcpcounter_res,file = str_glue("{od}/lusc_mcpcounter.csv"))

cibersort_res <- deconvo_tme(eset = as.data.frame(expr_dat), method = "cibersort", arrays = F, perm = 1000) %>%
                dplyr::rename(sample = ID) %>% dplyr::select(sample:Neutrophils_CIBERSORT)
write.csv(cibersort_res,file = str_glue("{od}/lusc_cibersort.csv"))

data_ic50 <- oncoPredict_calcPhenotype(
    valid_exp = expr_dat, database = "GDSC2",
    drug_list = NULL,
    train_exp = NULL, train_ptype = NULL, od = str_glue("{od}/")
)
