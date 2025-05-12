rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/results/4.0.cluster_deg/")
dir.create(od)
pacman::p_load(tidyverse,GSVA,enrichplot)

set.seed(1110)
drug_info <- read.csv(file = str_glue("{run_home}/data/PANCANCER_ANOVA_Tue_Oct_29_06_43_35_2024.csv"), check.names = FALSE) %>% 
                    mutate(drug = paste(`Drug name`,`Drug ID`,sep ="_")) %>% dplyr::select(drug,pathway = "Target Pathway") %>% unique()
sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))
load(str_glue("{run_home}/data/expr_data.RData"))
walk(sort(unique(sample_cluster$Group)),function(x){
    # x <- "Cluster2"
    pdata <- sample_cluster %>% dplyr::rename(sample = "Sample") %>% 
            mutate(group = ifelse(Group == x,x,"Others")) %>% dplyr::select(sample,group)
    degs <- limma_deg(DEG_exp = expr_dat, DEG_pdata = pdata, od = str_glue("{od}/{x}/"),
                    controlLabel = "Others",caseLabel = x,color_fun = c("#80253F","#015694", "#e6e1e1"),
                    DEG_FC = 0.585, DEG_P = 0.05, prefix = x)
    save(degs,file = str_glue("{od}/{x}/degs.Rdata"))
    load(str_glue("{od}/{x}/degs.Rdata"))
    gene_go_kegg <- enrich(
        genetype = x, genelist = degs$DEGs, od = str_glue("{od}/{x}/go_kegg"),
        pAdjustMethod = "none", color_fun = c("#2590CF", "#B9E5F8", "#F2737F"),w=6,h=4
    )
    gene_enrich <- enrich_Reactome(
        genetype = x, genelist = degs$DEGs, od = str_glue("{od}/{x}/reactome"),w=4,h=4,
        pAdjustMethod = "none",color_fun = c("#2590CF", "#B9E5F8", "#F2737F")
    )
    data <- degs[["allDEG"]] %>% as.data.frame() %>% rownames_to_column("gene") %>% 
            filter(gene %in% degs$DEGs)
    logFC <- data$logFC
    names(logFC) <- data$gene
    p <- cnetplot(gene_enrich$Reactome,
        showCategory = 5,
        circular = F, # 是否环形展示
        foldChange = logFC,
        layout = "kk",
        node_label = "category",
        cex_category = 0.8,
        cex_gene = 0.5,
        cex_label_category = 0.8,
        cex_label_gene = 0.5,
        color_category = "#C06CAB",
        color_gene = "#8A9FD1"
    ) + labs(title = str_c(x, " Reactome")) + theme(plot.title = element_text(hjust = 0.6, vjust = 1, size = 10))
    ggsave(filename = str_glue("{od}/{x}/Figure_{x}_Reactome_cnetplot.pdf"), plot = p, width = 5, height = 4)
    p <- emapplot(pairwise_termsim(gene_enrich$Reactome),
        showCategory = 10,
        color = "pvalue",
        cex_category = 0.8,
        cex_label_category = 0.8,
        min_edge = 0.1,
        repel = TRUE
    ) + labs(title = str_c(x, " Reactome")) + theme(plot.title = element_text(hjust = 0.7, vjust = 1, size = 10))
    ggsave(filename = str_glue("{od}/{x}/Figure_{x}_Reactome_emapplot.pdf"), plot = p, width = 5.5, height = 5)

    load(str_glue("{od}/{x}/degs.Rdata"))
    dir.create(str_glue("{od}/{x}/ic50"))
    data_ic50 <- read.csv(file = str_glue("{run_home}/results/1.0.data_prepare/calcPhenotype_Output/DrugP#F2737Fictions.csv"), row.names = 1, check.names = FALSE)
    genelist <- degs$DEGs
    # identical(colnames(expr_dat),rownames(data_ic50))
    exprs <- t(expr_dat) %>% as.data.frame()
    cor <- psych::corr.test(exprs[,genelist[genelist %in% colnames(exprs)]], data_ic50, method = "spearman", adjust = "none")
    r <- cor$r %>% reshape2::melt(varnames = c("gene", "ic50"), value.name = "cor")
    p <- cor$p %>% reshape2::melt(varnames = c("gene", "ic50"), value.name = "pvalue")
    res <- merge(r, p, by = c("ic50", "gene"))
    res <- res %>% inner_join(drug_info,by = join_by("ic50" == "drug"))
    write.table(res, file = str_glue("{od}/{x}/ic50/gene_ic50_cor_pval.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = T)
    
    res <- data.table::fread(str_glue("{od}/{x}/ic50/gene_ic50_cor_pval.txt")) %>% mutate(drug = str_split_i(ic50,"_",1))
    res <- res %>% filter(pvalue < 0.05,pathway != "",pathway != "Other") %>% 
            mutate(p_or_n_cor = ifelse(cor > 0,"positive","negative"))
    show_feature <- res %>%
        group_by(ic50) %>%
        summarise(cor_gene_num = n()) %>%
        arrange(desc(cor_gene_num)) %>% head(30) %>% pull(ic50)

    plot_feature <- res %>% filter(ic50 %in% show_feature) %>% 
        group_by(ic50,p_or_n_cor) %>%
        summarise(cor_gene_num = n()) %>%
        arrange(desc(cor_gene_num)) %>% 
        mutate(drug = str_split_i(ic50,"_",1))

    p <- ggplot(plot_feature, aes(x = cor_gene_num,y = drug)) +
                geom_bar(stat = "identity",aes(fill = p_or_n_cor), position = position_dodge(),width = 0.8) +
                theme_classic() + scale_fill_manual(values = c("#2590CF","#F2737F"))+
                labs(x = str_glue("Number of associated genes"), y = "") +
                theme(legend.position = "top",axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 5)) #+
                #ggrepel::geom_text_repel(aes(x = cor_gene_num,y = drug,label = cor_gene_num),size = 2)
    ggsave(filename = str_glue("{od}/{x}/ic50/gene_drug_cor.pdf"), plot = p, width = 4, height = 5)

    p <- ggplot(plot_feature %>% filter(p_or_n_cor == "positive"), aes(x = cor_gene_num,y = drug)) +
                geom_bar(stat = "identity", width = 0.8,fill = "#F2737F") +
                theme_classic() + 
                labs(x = str_glue("Number of positive associated genes"), y = "") +
                theme(legend.position = "none",axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 6)) +
                geom_text(aes(x = cor_gene_num,y = drug,label = cor_gene_num),size = 1.5)
    ggsave(filename = str_glue("{od}/{x}/ic50/gene_drug_positive.pdf"), plot = p, width = 3, height = 4)

    p <- ggplot(plot_feature %>% filter(p_or_n_cor == "negative"), aes(x = cor_gene_num,y = drug)) +
                geom_bar(stat = "identity", width = 0.8,fill = "#2590CF") +
                theme_classic() + 
                labs(x = str_glue("Number of negative associated genes"), y = "") +
                theme(legend.position = "none",axis.text.x = element_text(size = 5),axis.text.y = element_text(size = 6)) +
                geom_text(aes(x = cor_gene_num,y = drug,label = cor_gene_num),size = 1.5)
    ggsave(filename = str_glue("{od}/{x}/ic50/gene_drug_negative.pdf"), plot = p, width = 3, height = 4)

    pathway_feature <- res %>% filter(ic50 %in% show_feature) %>% 
        group_by(pathway) %>%
        summarise(cor_gene_num = n()) %>%
        arrange(desc(cor_gene_num)) %>% head(1) %>% pull(pathway)
    print(pathway_feature)

    show_genes <- res %>% filter(pathway == pathway_feature,ic50 %in% show_feature) %>% 
        group_by(gene) %>%
        summarise(cor_gene_num = n()) %>%
        arrange(desc(cor_gene_num)) %>% head(30) %>% pull(gene)

    show_drug <- res %>% filter(pathway == pathway_feature,ic50 %in% show_feature) %>% pull(drug) %>% unique()
    p <- res %>%
        filter(gene %in% show_genes,drug %in% show_drug) %>%
        ggplot(data = ., mapping = aes(x = drug, y = gene)) +
        geom_point(mapping = aes(colour = cor, size = -log(pvalue))) +
        theme_bw(10) + 
        labs(x = "", y = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 5), axis.text.y = element_text(size = 5),panel.grid = element_blank()) +
        scale_size(range = c(0.5, 4)) +
        scale_color_gradient2(low = "#377EB8", high = "firebrick3", space = "Lab") +
        guides(colour = guide_colorbar(title = "Correlation"))
    ggsave(filename = str_glue("{od}/{x}/ic50/drug_gene_corpoint.pdf"), plot = p, width = 4, height = 5)
})
