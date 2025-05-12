rm(list=ls())
run_home <- "./"

od <- file.path(run_home, "/data")
dir.create(od)
pacman::p_load(TCGAbiolinks,SummarizedExperiment,tidyverse)

query.exp <- GDCquery(
                project = paste0("TCGA-","LUSC"),
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                experimental.strategy = "RNA-Seq",
                workflow.type = "STAR - Counts")
GDCdownload(query.exp)
all_dat <- GDCprepare(query = query.exp)
save(all_dat,file = str_glue("{od}/raw_expr.Rdata"))
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
expr_dat <- log2(expr_dat + 1) %>% as.data.frame() %>% dplyr::select(grep("^.{13}01A.*$", colnames(.)))
colnames(expr_dat) <- lapply(colnames(expr_dat), function(x) {substring(x, 1, 12)}) %>% unlist()
expr_dat <- expr_dat %>% select(-colnames(expr_dat)[duplicated(colnames(expr_dat))])
save(expr_dat,file = str_glue("{od}/expr_data.RData"))

query <- GDCquery(
    project = paste0("TCGA-","LUSC"), 
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "bcr xml"
)
GDCdownload(query)
data_clinical <- GDCprepare_clinic(query, clinical.info = "patient")
save(data_clinical, file = str_glue("{od}/data_clinical.RData"))

query <- GDCquery(
    project = paste0("TCGA-","LUSC"), 
    data.category = "Simple Nucleotide Variation", 
    access = "open",
    data.type = "Masked Somatic Mutation", 
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
data_snv <- GDCprepare(query)
save(data_snv, file = str_glue("{od}/data_snv.RData"))

query <- GDCquery(project = paste0("TCGA-","LUSC"), 
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment")
GDCdownload(query)
data_cnv <- GDCprepare(query)
save(data_cnv, file = str_glue("{od}/data_cnv.RData"))
