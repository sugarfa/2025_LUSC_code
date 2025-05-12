rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/results/2.0.cluster_km/")
dir.create(od)
pacman::p_load(tidyverse,GSVA,ComplexHeatmap,RColorBrewer,IOBR,factoextra,FactoMineR,circlize)

set.seed(1110)

data_survival <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_survival.csv"))
fges_res <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_fges.csv")) %>% select(-1)
fges_nmf <- fges_res %>% column_to_rownames("sample") %>% + 0.5
nmf_res <- get_nmf_res(
            input = fges_nmf %>% t(), od = str_glue("{od}/"),
            rank_space = 2:5, thread = 150, cluster_res = FALSE, nrun = 10, method = "brunet", seed = 123456)
save(nmf_res,file = str_glue("{od}/nmf_res.RData"))

load(str_glue("{od}/nmf_res.RData"))
nmf_res$rank_opt_coph
table(nmf_res$group_opt$Group)
nmf_res$group_opt <- nmf_res$group_opt %>% mutate(Group = str_c("Cluster",Group))
write.csv(nmf_res$group_opt, file = str_glue("{od}/nmf_fges_cluster.csv"), quote = FALSE, row.names = FALSE)

dat_km <- nmf_res$group_opt %>% inner_join(data_survival %>% dplyr::select(Sample = "bcr_patient_barcode",status = "OS",time ="OS.time"))
kmplot <- survfit(Surv(time, status) ~ Group, data = dat_km)
p <- survminer::ggsurvplot(kmplot,
                data = dat_km, conf.int = FALSE, pval = TRUE, conf.int.style = "step",
                risk.table = "absolute",
                pval.size = 5, palette = color_four,
                fontsize = 4,tables.height = 0.2,
                risk.table.y.text = FALSE, ncensor.plot = FALSE,
                xlab = str_c("OS", " days"))
pdf(file = str_glue("{od}/fges_nmf_k4.pdf"), height = 4, width = 4)
print(p)
dev.off()

