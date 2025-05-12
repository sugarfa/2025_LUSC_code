rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/results/5.0.cluster_mut/")
dir.create(od)
pacman::p_load(tidyverse,maftools,RColorBrewer,IOBR,factoextra,FactoMineR,circlize,enrichplot)
source(str_glue("{run_home}/src/function.r"))

set.seed(1110)

sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))
load(str_glue("{run_home}/data/data_snv.RData"))

data_snv <- data_snv %>% mutate(sample = substr(Tumor_Sample_Barcode,1,12))
sample_cluster <- read.csv(str_glue("{run_home}/results/2.0.cluster_km/nmf_fges_cluster.csv"))

get_TMB <- function(snv_dat) {
  # 需要用到的列
  use_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", 
    "HGVSc", "t_depth", "t_alt_count"
  )
  # 删除这些突变类型
  mut_type <- c(
    "5'UTR", "3'UTR", "3'Flank", "5'Flank", "Intron", "IGR",
    "Silent", "RNA", "Splice_Region"
  )
  # 读取文件
  df <- snv_dat %>% select(use_cols)
  data <- df %>% subset(!Variant_Classification %in% mut_type) %>%
    # 计算 VAF
    mutate(vaf = t_alt_count / t_depth) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(mut_num = n(), TMB = mut_num / 30, MaxVAF = max(vaf))
  return(data)
}
tmb_data <- get_TMB(data_snv) %>% filter(substr(Tumor_Sample_Barcode,14,15)=="01") %>% mutate(Sample = substr(Tumor_Sample_Barcode,1,12))
tmb_data <- merge(sample_cluster,tmb_data)

p <- ggplot(tmb_data,aes(x = Group,y = TMB,fill = Group))+
  geom_boxplot()+ ggthemes::theme_few()+ ylim(0,50)+
  ggpubr::stat_compare_means(aes(group = Group,label = after_stat(p.signif)))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))+
  scale_fill_manual(values =  c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))
ggsave(p,filename = str_glue("{od}/tmb_boxplot.pdf"),width = 4,height = 3)    

wilcox.test(tmb_data$TMB,g = tmb_data$Group)

kruskal.test(tmb_data$TMB,g = tmb_data$Group)
snv_dat <- read.maf(data_snv)

vc_cols = c("#1F6CB6","red3","#70BC6B","#F4B9C5","#784194","#B81A7B","#A65628","#9E1F63")  # RdYlBu
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Del',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
col1 = c("black","gray50","#CCCCCC");
names(col1) = c('MSI-H','MSI-L','MSS')
col2 = c("#784193","#CCCCCC","darkgreen","skyblue4","white");
names(col2) = c('Stage I','Stage II',"Stage III","Stage IV","NA")


pdf(file = str_glue("{od}/lusc_mutation_summary.pdf"), height = 5, width = 8)
plotmafSummary(maf = snv_dat, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf(file = str_glue("{od}/lusc_mutation_oncoplot.pdf"), height = 5, width = 5)
oncoplot(maf = snv_dat,top = 20,colors = vc_cols,annotationColor = list(MSI_STATUS=col1,TUMOR_STAGE_2009=col2))
dev.off()

walk(unique(sample_cluster$Group),function(x){
    dir.create(str_glue("{od}/{x}"))
    cluster_sample <- sample_cluster %>% filter(Group == x) %>% pull(Sample)
    snv_dat_sub <- read.maf(data_snv %>% filter(sample %in% cluster_sample))
    pdf(file = str_glue("{od}/{x}/lusc_mutation_snp_summary.pdf"), height = 5, width = 8)
    plotmafSummary(maf = snv_dat_sub, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    pdf(file = str_glue("{od}/{x}/lusc_mutation_snp_oncoplot.pdf"), height = 5, width = 6)
    oncoplot(maf = snv_dat_sub,,colors = vc_cols,annotationColor = list(MSI_STATUS=col1,TUMOR_STAGE_2009=col2))
    dev.off()

    walk(snv_dat_sub@gene.summary$Hugo_Symbol %>% head(5),function(y){
        pdf(file = str_glue("{od}/{x}/lollipop_{y}.pdf"), height = 5, width = 4)
        lollipopPlot(maf = snv_dat, gene = y, showDomainLabel = FALSE)
        dev.off()
    })
})

gene_per <- map_df(sort(unique(sample_cluster$Group)),function(x){
    #dir.create(str_glue("{od}/{x}"))
    cluster_sample <- sample_cluster %>% filter(Group == x) %>% pull(Sample)
    snv_dat_sub <- read.maf(data_snv %>% filter(sample %in% cluster_sample))
    gene_per <- snv_dat_sub@gene.summary %>% as.data.frame() %>% 
                slice(1:20) %>% select(Hugo_Symbol,total,AlteredSamples) %>% 
                mutate(group = x, mut = AlteredSamples,non_mut = total - AlteredSamples,mut_p = round(AlteredSamples/total,3),non_mut_p = 1 - mut_p) %>% 
                select(-total,-AlteredSamples)
})
gene <- sort(table(gene_per$Hugo_Symbol),decreasing = T) %>% head(10) %>% names()
walk(gene,function(g){
   p <- gene_per %>% filter(Hugo_Symbol %in% g) %>% select(mut,non_mut) %>% chisq.test()
   print(c(g,p$p.value))
})

walk(c("MUC16","TP53","USH2A"),function(g){
  p <- gene_per %>% filter(Hugo_Symbol %in% g) %>% select(mut,non_mut) %>% chisq.test()
  label <- paste("chisq.test p value =",round(p$p.value,3))
  plot_dat <- gene_per %>% filter(Hugo_Symbol %in% g) %>% pivot_longer(cols = c(5,6)) 
  
  # p_plot <- ggplot(plot_dat,aes(x= group,y = percent,fill = group)) + geom_bar(stat = "identity") +
  #   scale_fill_manual(values =  c("#d2985f","#799c54","#5c6fdb","#ae4fa2"))+
  #   ggthemes::theme_few() + xlab(label = label)+ ylim(c(0,1))+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),legend.position = "top")
  p2 <- ggplot(data = plot_dat)+
  geom_bar(aes(x = group,y = value,
               fill = name),
           stat = "identity") +
  scale_fill_manual(values =  c("#f87669","#2fa1dd"))+
  geom_label(aes(x = group,y = value,
                label = scales::percent(value),
                fill = name),
            color = "black",
            size = 3,label.size = 0,
            show.legend = F,
            position = position_fill(vjust = 0.5))+
  ggthemes::theme_few()+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),legend.position = "top")+
  guides(label = "none") + xlab(label = label)

  ggsave(p2 ,filename = str_glue("{od}/lusc_{g}_mut_per.pdf"),width = 3,height = 4)
})

snv_dat <- read.maf(data_snv)

vc_cols = c("#1F6CB6","red3","#70BC6B","#F4B9C5","#784194","#B81A7B","#A65628","#9E1F63")  # RdYlBu
names(vc_cols) = c(
  'Missense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Del',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
col1 = c("black","gray50","#CCCCCC");
names(col1) = c('MSI-H','MSI-L','MSS')
col2 = c("#784193","#CCCCCC","darkgreen","skyblue4","white");
names(col2) = c('Stage I','Stage II',"Stage III","Stage IV","NA")


pdf(file = str_glue("{od}/lusc_mutation_summary.pdf"), height = 5, width = 8)
plotmafSummary(maf = snv_dat, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf(file = str_glue("{od}/lusc_mutation_oncoplot.pdf"), height = 5, width = 5)
oncoplot(maf = snv_dat,top = 20,colors = vc_cols,annotationColor = list(MSI_STATUS=col1,TUMOR_STAGE_2009=col2))
dev.off()

walk(unique(sample_cluster$Group),function(x){
    dir.create(str_glue("{od}/{x}"))
    cluster_sample <- sample_cluster %>% filter(Group == x) %>% pull(Sample)
    snv_dat_sub <- read.maf(data_snv %>% filter(sample %in% cluster_sample))
    pdf(file = str_glue("{od}/{x}/lusc_mutation_snp_summary.pdf"), height = 5, width = 8)
    plotmafSummary(maf = snv_dat_sub, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
    dev.off()
    pdf(file = str_glue("{od}/{x}/lusc_mutation_snp_oncoplot.pdf"), height = 5, width = 6)
    oncoplot(maf = snv_dat_sub,,colors = vc_cols,annotationColor = list(MSI_STATUS=col1,TUMOR_STAGE_2009=col2))
    dev.off()

    walk(snv_dat_sub@gene.summary$Hugo_Symbol %>% head(5),function(y){
        pdf(file = str_glue("{od}/{x}/lollipop_{y}.pdf"), height = 5, width = 4)
        lollipopPlot(maf = snv_dat, gene = y, showDomainLabel = FALSE)
        dev.off()
    })
})
