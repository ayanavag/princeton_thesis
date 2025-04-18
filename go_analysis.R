library(AnnotationDbi)
library(org.Hs.eg.db)
library(magrittr)
library(dplyr)
library(tidyr)
library(gtools)
library(minfi) 
library(missMethyl)
library(readxl)
library(sva)
library(qqman)
library(stringr)
library(IlluminaHumanMethylation450kmanifest)
library(FDRestimation)
library(bacon)
library(ENmix)
library(limma)
library(ewastools)
library(DMRcate)
library(annotatr)
library(clusterProfiler)


reference <- as.data.frame(readRDS("/drives/drive1/geneticData/referenceDataSets/Methylation/Zhou/EPIC.hg19.manifest.rds"))

# read in the 3 matches found from Aim 2
dmr_match_r <- readRDS("/drives/drive1/ff_chaya/Ayanava/dmr_match_r.rds")

dmrs <- dmr_match_r[1:3, ]
dmrs <- dmrs %>% select(c("chr", "start", "end"))
dmrs$gene <- c(0, 0, 0)
for (i in 1:3) {
  dmrs$gene[i] <- reference %>% filter(reference$seqnames %in% dmrs$chr[i]) %>%
    filter(start == dmrs$start[i]) %>% select(gene)
}

genes <- (unlist(strsplit(unlist(dmrs$gene), ";")))

# add in the genes identified from Aim 1
gene <- "SOCS3"
genes <- c(genes, gene)
gene <- "CPT1A"
genes <- c(genes, gene)

genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# GO analysis for biological processes
go_results_bio <- enrichGO(genes_entrez$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", 
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_results_bio_df <- as.data.frame(go_results_bio@result)


# GO analysis for the CpGs from the FF MR 
mr_results1 <- readRDS("/drives/drive1/ff_chaya/Ayanava/mr_ff_results_20kb.rds")
reference_cpgs1 <- reference[rownames(reference) %in% mr_results1$id.exposure, ]
reference_cpgs1$id.exposure <- rownames(reference_cpgs1)
mr_genes <- merge(mr_results1, reference_cpgs1[, c("id.exposure", "gene")], by = "id.exposure", all.X = TRUE)
mr_genes$adj.pval <- p.adjust(mr_genes$pval, method = "BH")
mr_genes <- mr_genes %>% filter(adj.pval < 0.05)
mr_genes <- na.omit(mr_genes)

# checking the number of unique genes and the CpGs corresponding to it
unique_mr_genes <- mr_genes %>%
  group_by(id.exposure) %>% 
  group_by(gene) %>%
  summarize(count = n_distinct(id.exposure), CpGs = paste(unique(id.exposure), collapse = ",")) %>% 
  ungroup()

saveRDS(unique_mr_genes, "/drives/drive1/ff_chaya/Ayanava/mr_genes_ff.rds")


genes_mr_list <- (unlist(strsplit(unlist(mr_genes$gene), ";")))

genes_entrez_mr <- bitr(genes_mr_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# creating background genes from EPIC array 
background_gene <- reference %>% select(gene)
background_gene <- na.omit(background_gene)
bg_gene_list <- (unlist(strsplit(unlist(background_gene$gene), ";")))
bg_gene_list <- unique(bg_gene_list)

bg_gene_entrez <- bitr(bg_gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO analysis for biological processes for FFCWS MR
go_results_bio_mr <- enrichGO(genes_entrez_mr$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", 
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_results_bio_mr_df <- as.data.frame(go_results_bio_mr@result)

# GO analysis with the background gene set 
go_results_bio_mr_bg <- enrichGO(genes_entrez_mr$SYMBOL, universe = bg_gene_entrez$SYMBOL, OrgDb = org.Hs.eg.db, 
                              keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_results_bio_mr_bg_df <- as.data.frame(go_results_bio_mr@result)



# keep the genes that have at least 3 CpGs associated with them
sig_genes_mr <- unique_mr_genes %>% filter(count > 3)
sig_genes_mr_list <- (unlist(strsplit(unlist(sig_genes_mr$gene), ";")))

sig_genes_entrez_mr <- bitr(sig_genes_mr_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

sig_go_results_bio_mr <- enrichGO(sig_genes_entrez_mr$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", 
                              ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
sig_go_results_bio_mr_df <- as.data.frame(sig_go_results_bio_mr@result)




# doing the same GO Term Analysis, but for the UK GWAS
mr_results2 <- readRDS("/drives/drive1/ff_chaya/Ayanava/mr_uk_results_20kb.RDS")
reference_cpgs2 <- reference[rownames(reference) %in% mr_results2$id.exposure, ]
reference_cpgs2$id.exposure <- rownames(reference_cpgs2)
mr_genes2 <- merge(mr_results2, reference_cpgs2[, c("id.exposure", "gene")], by = "id.exposure", all.X = TRUE)
mr_genes2$adj.pval <- p.adjust(mr_genes2$pval, method = "BH")
mr_genes2 <- mr_genes2 %>% filter(adj.pval < 0.05)
mr_genes2 <- na.omit(mr_genes2)

# checking the number of unique genes and the CpGs corresponding to it
unique_mr_genes2 <- mr_genes2 %>%
  group_by(id.exposure) %>% 
  group_by(gene) %>%
  summarize(count = n_distinct(id.exposure), CpGs = paste(unique(id.exposure), collapse = ",")) %>% 
  ungroup()

# finding the intersection between genes of FFCWS and UK Biobank
common_mr_genes <- intersect(unique_mr_genes$gene, unique_mr_genes2$gene)
com_mr1 <- unique_mr_genes[unique_mr_genes$gene %in% common_mr_genes, ]
com_mr2 <- unique_mr_genes2[unique_mr_genes2$gene %in% common_mr_genes, ]
com_mr <- merge(com_mr1, com_mr2, by = "gene")


# running GO term on all genes from UK Biobank 
genes_mr_list2 <- (unlist(strsplit(unlist(mr_genes2$gene), ";")))

genes_entrez_mr2 <- bitr(genes_mr_list2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# GO analysis for biological processes for UK Biobank MR 
go_results_bio_mr2 <- enrichGO(genes_entrez_mr2$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", 
                              ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_results_bio_mr_df2 <- as.data.frame(go_results_bio_mr2@result)


# GO term analysis on the genes that overlap between the two cohorts
genes_mr_list3 <- (unlist(strsplit(unlist(common_mr_genes), ";")))

genes_entrez_mr3 <- bitr(genes_mr_list3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO analysis for biological processes
go_results_bio_mr3 <- enrichGO(genes_entrez_mr3$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", 
                               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_results_bio_mr_df3 <- as.data.frame(go_results_bio_mr3@result)



# GO Term analysis on the CpGs that are significant on both the UK and FF cohorts
common_cpgs <- intersect(mr_results1$id.exposure, mr_results2$id.exposure)
common_mr1 <- mr_results1[mr_results1$id.exposure %in% common_cpgs, ]
common_mr1$adj.pval <- p.adjust(common_mr1$pval, method = "BH")
common_mr2 <- mr_results2[mr_results2$id.exposure %in% common_cpgs, ]
common_mr2$adj.pval <- p.adjust(common_mr2$pval, method = "BH")

# how many are the same at p<0.05 
sig_mr1 <- common_mr1 %>% filter(adj.pval < 0.05)
sig_mr2 <- common_mr2 %>% filter(adj.pval < 0.05)
common_sig_cpgs <- intersect(sig_mr1$id.exposure, sig_mr2$id.exposure) 
reference_cpgs3 <- reference[rownames(reference) %in% common_sig_cpgs, ]
reference_cpgs3$id.exposure <- rownames(reference_cpgs3)
mr_genes3 <- reference_cpgs3 %>% select("id.exposure", "gene")
mr_genes3 <- na.omit(mr_genes3)

# checking the number of unique genes and the CpGs corresponding to it
unique_mr_genes3 <- mr_genes3 %>%
  group_by(id.exposure) %>% 
  group_by(gene) %>%
  summarize(count = n_distinct(id.exposure), CpGs = paste(unique(id.exposure), collapse = ",")) %>% 
  ungroup()


genes_mr_list4 <- (unlist(strsplit(unlist(mr_genes3$gene), ";")))

genes_entrez_mr4 <- bitr(genes_mr_list4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# GO analysis for biological processes for the genes with overlapping significant CpGs
go_results_bio_mr4 <- enrichGO(genes_entrez_mr4$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", 
                              ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_results_bio_mr_df4 <- as.data.frame(go_results_bio_mr4@result)



