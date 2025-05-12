rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/results/3.0.cluster_tme/")
dir.create(od)
pacman::p_load(tidyverse,GSVA,ComplexHeatmap,RColorBrewer,IOBR,factoextra,FactoMineR,circlize)
source(str_glue("{run_home}/src/function.r"))

set.seed(1110)

sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))

fges_info <- readxl::read_xlsx(str_glue("{run_home}/data/29fges_Cancer_cell.xlsx"),sheet = 1,skip = 1) %>% 
            arrange(ID,"Process/Signature")
fges_res <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_fges.csv"),check.names = F) %>% select(-1)
characteristics_plot_by_group(
    characteristics_score = fges_res,
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/fges"), type = "fges",
    color_fun = color_four,#brewer.pal(3, "Set1")
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
plot_data <- fges_res %>% column_to_rownames("sample") %>% 
            dplyr::select(fges_info$`Process_Signature`) %>% 
            scale() %>% t() %>% as.data.frame() %>% select(sample_cluster$Sample)

row_ann <- rowAnnotation("fges type" = fges_info$Type,
                col = list("fges type" = setNames(brewer.pal(4, "Set1"), unique(fges_info$Type))),
                simple_anno_size = unit(0.15, "cm"),
                annotation_legend_param = list("fges type" = list(title = "fges type",
                                                labels_gp = gpar(fontsize = 5),
                                                title_gp = gpar(fontsize = 8))))
col_ann <- HeatmapAnnotation("TME subtype" = sample_cluster$Group,
                col = list("TME subtype" = setNames(color_four, sort(unique(sample_cluster$Group)))),
                simple_anno_size = unit(0.15, "cm"),
                annotation_legend_param = list("TME subtype" = list(title = "TME subtype",
                                                labels_gp = gpar(fontsize = 5),
                                                labels = c("Cluster1", "Cluster2", "Cluster3","Cluster4"),
                                                title_gp = gpar(fontsize = 8))))

pdf(file = str_glue("{od}/fges/fges_group_heatmap.pdf"), height = 5, width = 8)
Heatmap(as.matrix(plot_data),
        column_order = sample_cluster$Sample,name = "fges",
        #col = colorRamp2(c(-5, 0, 5), c("green", "white", "red")),
        show_column_names = F,cluster_rows = F,cluster_column_slices = F,
        row_names_gp = gpar(fontsize = 8),
        column_title_gp = gpar(fontsize = 5),row_title_gp = gpar(fontsize = 5),
        column_split = sample_cluster$Group,row_split = fges_info$Type,
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 5),title_gp = gpar(fontsize = 8)),
        left_annotation = row_ann,top_annotation = col_ann)
dev.off()

pca_dat <- fges_res %>% slice(match(sample_cluster$Sample,fges_res$sample)) %>% column_to_rownames("sample")
identical(rownames(pca_dat),sample_cluster$Sample)
pca_res <- PCA(pca_dat,scale.unit = F, ncp = 2, graph = F)
p <- fviz_pca_ind(
    X = pca_res, 
    axes = 1:2,
    geom = "point",
    habillage = as.factor(sample_cluster$Group),
    legend.title = "Cluster",
    # palette = "lancet",
    palette = color_four,
    addEllipses = T, 
    ellipse.level = 0.95, 
    title = "PCA Plot",
    mean.point = T 
) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 20), panel.grid = element_blank())
ggsave(filename = str_glue("{od}/fges_pca.pdf"),p,height = 3,width = 3.5)

mcpcounter_dat <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_mcpcounter.csv"),check.names = F) %>% select(-1)
colnames(mcpcounter_dat) <- gsub("_MCPcounter","",colnames(mcpcounter_dat))
identical(mcpcounter_dat$sample,sample_cluster$Sample)
characteristics_plot_by_group(
    characteristics_score = mcpcounter_dat,
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/mcpcounter"), type = "mcpcounter",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

cibersort_dat <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_cibersort.csv"),check.names = F) %>% select(-1)
colnames(cibersort_dat) <- gsub("_CIBERSORT","",colnames(cibersort_dat))
identical(cibersort_dat$sample,sample_cluster$Sample)
characteristics_plot_by_group(
    characteristics_score = cibersort_dat,
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/cibersort"), type = "cibersort",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
ciber_dat <- cibersort_dat %>% column_to_rownames("sample") %>% colMeans() %>% 
                as.data.frame() %>% rownames_to_column("cell_type") %>% rename(proportion = ".") %>% 
                mutate(proportion = 100 * proportion) %>% arrange(proportion)
ciber_dat$cell_type <- factor(ciber_dat$cell_type,levels = ciber_dat$cell_type)
p <- ggplot(ciber_dat %>% filter(proportion > 0.5), aes(x = cell_type,y = proportion)) +
  geom_bar(aes(fill = cell_type),stat = "identity") +
  theme_classic() +
  labs(x = "", y = "Proportion of immune cells(%)") + coord_flip() +
  theme(legend.position = "none",axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 6),axis.title.x = element_text(size = 8))
ggsave(filename = str_glue("{od}/ciber_cell_proportion.pdf"), plot = p, width = 4, height = 3.5)

identical(cibersort_dat$sample,rownames(fges_res))
fges_res <- fges_res %>% column_to_rownames("sample") %>% .[match(rownames(.),cibersort_dat$sample),]
cor <- psych::corr.test(cibersort_dat %>% column_to_rownames("sample"), fges_res, method = "spearman", adjust = "none")
r <- cor$r %>% reshape2::melt(varnames = c("immune_cell", "fges"), value.name = "cor")
p <- cor$p %>% reshape2::melt(varnames = c("immune_cell", "fges"), value.name = "pvalue")
res <- merge(r, p, by = c("immune_cell", "fges"))
write.table(res, file = str_glue("{od}/ciber_fges_cor_pval.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = T)
res <- res %>% inner_join(fges_info[,c("Type","ID")],by = join_by("fges" == "ID")) %>% arrange(Type)
res$fges <- factor(res$fges,levels = unique(res$fges))
p <- res %>%
    mutate(pvalue = ifelse(pvalue < 0.05, pvalue, NA)) %>%
    ggplot(data = ., mapping = aes(x = fges, y = immune_cell)) +
    geom_point(mapping = aes(colour = cor, size = -log(pvalue))) +
    theme_bw(10) +
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 5), axis.text.y = element_text(size = 6),panel.grid = element_blank()) +
    scale_size(range = c(0.2, 2.5)) +
    scale_color_gradient2(low = "#377EB8", high = "firebrick3", space = "Lab") +
    guides(colour = guide_colorbar(title = "Correlation"))
ggsave(filename = str_glue("{od}/ciber_fges_cor_corpoint.pdf"), plot = p, width = 5.5, height = 4.2)

estimate_dat <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_estimate.csv"),check.names = F) %>% select(-1)
colnames(estimate_dat) <- gsub("_estimate","",colnames(estimate_dat))
identical(estimate_dat$sample,sample_cluster$Sample)
characteristics_plot_by_group(
    characteristics_score = estimate_dat,
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/estimate"), type = "estimate",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

hallmark_dat <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_hallmark.csv"),check.names = F) %>% select(-1)
identical(hallmark_dat$sample,sample_cluster$Sample)
characteristics_plot_by_group(
    characteristics_score = hallmark_dat,
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/hallmark"), type = "hallmark",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
pdf(file = str_glue("{od}/fig_hallmark_heatmap.pdf"), height = 4, width = 5.5)
p <- Heatmap_manul_v2(data_input = hallmark_dat %>% column_to_rownames("sample") %>% t(), 
                  key_color = group_color,
                  group_infor = sample_cluster %>% rename(group ="Group") %>% column_to_rownames("Sample"),
                  col_in_heat = c("#2590CF", "white", "#F2737F"),
                  Colored_OtherInfor = F,
                  w = 3, h = 6, heatmap_name = " ",
                  show_rownames = T, DoWilcox.test = T, rownames_fontsize = 3,
                  cluster_rows = T, cluster_columns = T)
print(p$plot)
dev.off()

ic50_dat <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/calcPhenotype_Output/DrugPredictions.csv")) %>% rename(sample =X)
identical(ic50_dat$sample,sample_cluster$Sample)
res_sign <- characteristics_plot_by_group(
    characteristics_score = ic50_dat,
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/ic50"), type = "ic50",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

drug_p <- res_sign %>% rownames_to_column("drug") %>% filter(!p.signif == " ") %>% pull(drug)
drug_info <- read.csv(file = str_glue("{run_home}/data/PANCANCER_ANOVA_Tue_Oct_29_06_43_35_2024.csv"), check.names = FALSE) %>% 
                    mutate(drug = paste(`Drug name`,`Drug ID`,sep ="_")) %>% dplyr::select(drug,pathway = "Target Pathway") %>% unique()
drug_info_p <- data.frame(drug = drug_p) %>% left_join(drug_info,by = "drug")
identical(drug_info_p$drug,rownames(ic50))
pdf(file = str_glue("{od}/fig_ic50_heatmap.pdf"), height = 10, width = 5.5)
p <- Heatmap_manul_v2(data_input = ic50_dat %>% select("sample",drug_p) %>% column_to_rownames("sample") %>% t(), 
                  key_color = group_color,
                  group_infor = sample_cluster %>% rename(group ="Group") %>% column_to_rownames("Sample"),
                  col_in_heat = c("#2590CF", "white", "#F2737F"),
                  row_type = data.frame(obs = drug_info_p$drug,Type = drug_info_p$pathway),
                  Colored_OtherInfor = F,
                  w = 3, h = 6, heatmap_name = " ",
                  show_rownames = T, DoWilcox.test = T, rownames_fontsize = 3,
                  cluster_rows = T, cluster_columns = T)
print(p$plot)
dev.off()

characteristics_plot_by_group(
    characteristics_score = ic50_dat[,c("sample",res_sign %>% rownames_to_column("name") %>% filter(pvalue< 0.05) %>% pull(name))],
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/ic50"), type = "ic50_sign",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

load(str_glue("{run_home}/data/expr_data.RData"))
chemokines <- read.csv(str_glue("{run_home}/data/immunomodulator.txt"),sep = "\t") %>% 
            filter(type == "chemokine") %>% pull(gene)
expr_dat <- expr_dat %>% t() %>% as.data.frame() %>% rownames_to_column("sample") 
characteristics_plot_by_group(
    characteristics_score = expr_dat %>% select(c("sample",chemokines)),
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/chemokine"), type = "chemokine",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
identical(expr_dat$sample,sample_cluster$Sample)
pdf(file = str_glue("{od}/fig_chemokines_heatmap.pdf"), height = 4, width = 5.5)
p <- Heatmap_manul_v2(data_input = expr_dat %>% select(c("sample",chemokines)) %>% column_to_rownames("sample") %>% t(), 
                  key_color = group_color,
                  group_infor = sample_cluster %>% rename(group ="Group") %>% column_to_rownames("Sample"),
                  col_in_heat = c("#2590CF", "white", "#F2737F"),
                  Colored_OtherInfor = F,
                  w = 3, h = 6, heatmap_name = " ",
                  show_rownames = T, DoWilcox.test = T, rownames_fontsize = 3,
                  cluster_rows = T, cluster_columns = T)
print(p$plot)
dev.off()

chemokines_dat <- merge(expr_dat %>% select(c("sample",chemokines)),sample_cluster %>% rename(sample = Sample))
dat <- chemokines_dat %>% pivot_longer(cols = c(2:42)) %>% filter(!name=="CCL27")
sort(unique(dat$name))
p <- ggplot(data = dat,aes(name,value, fill = Group)) +
        geom_boxplot(aes(col = Group),outlier.size = .2, outlier.shape = NA,
          width = .75,position = position_dodge(width = .9)) + xlab(NULL) +
        ggthemes::theme_few() + ylim(0,10)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
          legend.position = "right",
          legend.title = element_blank()
        ) +  ggpubr::stat_compare_means(
          aes(group = Group,
            label = after_stat(p.signif))) +
        ggplot2::coord_cartesian(clip = "on") + scale_fill_manual(values = ggplot2::alpha(c("#d2985f","#799c54","#5c6fdb","#ae4fa2"), 0.8)) +
        scale_color_manual(values = c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))
ggsave(p,filename = str_glue("{od}/chemokines_cluster.pdf"),width = 15,height = 4)    

walk(unique(dat$name),function(x){
    dat_p <- dat %>% filter(name == x)
    p <- Box_in_One_p(dat,x_ = "name",y_ = "value",group_ = "Group",style = "box2",
        color_in_p = c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))
    ggsave(p,filename = str_glue("{od}/chemokines_cluster.pdf"),width = 10,height = 4)    
})


checkpoints <- read.csv(file = str_glue("{run_home}/data/checkpoints_InAcMarker_extend.csv")) %>%
    mutate(Type = case_when(
        Role.with.Immunity %in% c("Activate", "Active") ~ "Activate",
        Role.with.Immunity %in% c("Inhibit") ~ "Inhibit",
        Role.with.Immunity %in% c("TwoSide") ~ "TwoSide"
    )) %>%
    arrange(Type, Symbol)

checkpoints <- checkpoints %>% filter(Symbol %in% rownames(expr_dat))
data_group <- sample_cluster %>% rename(group ="Group") %>% column_to_rownames("Sample")
pdf(file = str_glue("{od}/fig_checkpoints_heatmap.pdf"), height = 4, width = 5.5)
p <- Heatmap_manul_v2(data_input = expr_dat[checkpoints$Symbol,], 
                  key_color = group_color, group_infor = data_group,
                  col_in_heat = c("#2590CF", "white", "#F2737F"),
                  row_type = data.frame(obs = checkpoints$Symbol,Type = checkpoints$Type),
                  Colored_OtherInfor = F,
                  w = 3, h = 6, heatmap_name = " ",
                  show_rownames = T, DoWilcox.test = T, rownames_fontsize = 3,
                  cluster_rows = T, cluster_columns = T)
print(p$plot)
dev.off()

characteristics_plot_by_group(
    characteristics_score = expr_dat[,colnames(expr_dat) %in% c("sample",checkpoints$Symbol)],
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/checkpoints"), type = "all_chemokine",
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)

walk(unique(checkpoints$Type), function(x) {
    genelist <- checkpoints %>%
        filter(Type == x) %>%
        pull(Symbol)
    characteristics_plot_by_group(
    characteristics_score = expr_dat[,colnames(expr_dat) %in% c("sample",genelist)],
    Group = sample_cluster %>% column_to_rownames("Sample"),
    od = str_glue("{od}/checkpoints"), type = x,
    color_fun = color_four,
    heatplot_by_scale = TRUE, cluster_rows = TRUE
)
})

sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))
tidy_res <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_tidy_res.csv")) %>% 
           rename(Sample = "Patient")
dir.create(str_glue("{od}/tide"))
tidy_res <- merge(sample_cluster ,tidy_res) %>% arrange(desc(TIDE)) 

label <- paste("fisher.test p value =",round(fisher.test(table(tidy_res$Group,tidy_res$Responder))$p.value,3))
p1 <- ggplot(tidy_res, aes(x = 1:nrow(tidy_res), 
                            y = TIDE, 
                            fill = Responder)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values =  c("#f87669","#2fa1dd"))+
  xlab("patient")+
  annotate("text", x = 250, y = 2.5, label = label,size = 5) +  
  theme_minimal()
fisher.test(dat$n)
dat <- count(tidy_res,Group,Responder)
dat <- dat %>% group_by(Group) %>% 
  summarise(Responder = Responder,n = n/sum(n))
fisher.test(table(tidy_res[,c("Group","Responder")]))
# p-value = 0.07681
dat$Responder <- factor(dat$Responder,levels = c("False","True"))

p2 <- ggplot(data = dat)+
  geom_bar(aes(x = Group,y = n,
               fill = Responder),
           stat = "identity") +
  scale_fill_manual(values =  c("#f87669","#2fa1dd"))+
  geom_label(aes(x = Group,y = n,
                label = scales::percent(n),
                fill = Responder),
            color = "white",
            size = 3,label.size = 0,
            show.legend = F,
            position = position_fill(vjust = 0.5))+
  ggthemes::theme_few()+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),legend.position = "top")+
  guides(label = "none")
ggsave(p2,filename = str_glue("{od}/tide/tide_lusc_responder.pdf"),width = 3,height = 3)
p <- p1 + p2 + patchwork::plot_layout(widths = c(3,2),guides = "collect")
ggsave(p,filename = str_glue("{od}/tide/tide_lusc_tide.pdf"),width = 10,height = 5)

dat <- tidy_res[,c(2,4:6,8:10,12:16)] %>% pivot_longer(cols = c(3:12))

walk(unique(dat$name),function(x){
    dat_p <- dat %>% filter(name == x)
    p <- Box_in_One_p(dat_p,x_ = "name",y_ = "value",group_ = "Responder",style = "box2",
        color_in_p = c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))
    ggsave(p,filename = str_glue("{od}/tide/lusc_responde_{x}.pdf"),width = 3.5,height = 3)    
})

walk(unique(dat$name),function(x){
    dat_p <- dat %>% filter(name == x)
    p <- Box_in_One_p(dat_p,x_ = "name",y_ = "value",group_ = "Group",style = "box2",
        color_in_p = c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))
    ggsave(p,filename = str_glue("{od}/tide/lusc_cluster_{x}.pdf"),width = 3.5,height = 3)    
})

sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))
immunecellai_res <- read.csv(str_glue("{run_home}/results/1.0.data_prepare/lusc_immunecell_response.txt"),sep = "\t") %>% 
           rename(Sample = "X")
dir.create(str_glue("{od}/immunecellai"))
immunecellai_res <- merge(sample_cluster,immunecellai_res) %>% arrange(desc(Score)) 

label <- paste("fisher.test p value =",round(fisher.test(table(immunecellai_res$Group,immunecellai_res$Response))$p.value,3))
p1 <- ggplot(immunecellai_res, aes(x = 1:nrow(immunecellai_res), 
                            y = Score#,fill = Response 
                            )) +
  geom_bar(stat = "identity") +
  xlab("patient")+
  annotate("text", x = 250, y = 2.5, label = label,size = 5) +  
  theme_minimal()

dat <- count(immunecellai_res,Group,Response)
dat <- dat %>% group_by(Group) %>% 
  summarise(Responder = Response,n = n/sum(n))
dat$Responder <- factor(dat$Responder,levels = c("0","1"))

p2 <- ggplot(data = dat)+
  geom_bar(aes(x = Group,y = n,
               fill = Responder),
           stat = "identity") +
  scale_fill_manual(values =  c("#f87669","#2fa1dd"))+
  geom_label(aes(x = Group,y = n,
                label = scales::percent(n),
                fill = Responder),
            color = "white",
            size = 4,label.size = 0,
            show.legend = F,
            position = position_fill(vjust = 0.5))+
  theme_minimal()+
  guides(label = "none")
p <- p1 + p2 + patchwork::plot_layout(widths = c(3,2),guides = "collect")
ggsave(p,filename = str_glue("{od}/immunecellai/lusc_immunecellai.pdf"),width = 10,height = 5)

dat <- immunecellai_res[,c(2,3,5:29)] %>% rename(Responder = "Response") %>% pivot_longer(cols = c(3:27))

walk(unique(dat$name),function(x){
    dat_p <- dat %>% filter(name == x)
    p <- Box_in_One_p(dat_p,x_ = "name",y_ = "value",group_ = "Responder",style = "box2",
        color_in_p = c("#F2737F","#2590CF"))
    ggsave(p,filename = str_glue("{od}/immunecellai/lusc_group_{x}.pdf"),width = 3.5,height = 3)    
})

walk(unique(dat$name),function(x){
    dat_p <- dat %>% filter(name == x)
    p <- Box_in_One_p(dat_p,x_ = "name",y_ = "value",group_ = "Group",style = "box2",
        color_in_p = c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))
    ggsave(p,filename = str_glue("{od}/immunecellai/lusc_cluster_{x}.pdf"),width = 3.5,height = 3)    
})
