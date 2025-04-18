# Mendelian randomization
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
library(data.table)
library(GenomicRanges)
library(snpStats)
library(MatrixEQTL)
library(TwoSampleMR)

# Setting up GWAS to generate SNP-CIMT associations
new_genetic <- read.plink("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/ff_psych_cvid")
fam <- new_genetic$fam
cimt_le8 <- readRDS("/drives/drive1/ff_chaya/Ayanava/cimt_le8.RDS")
colnames(cimt_le8)[colnames(cimt_le8) == "age_at_clinic"] <- "Age_k7"
colnames(cimt_le8)[colnames(cimt_le8) == "cvid"] <- "CVID"
cimtdataset <- readRDS("/drives/drive1/ff_chaya/R_versions/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2_datafile.RDS")
biomarker <- readRDS("/drives/drive1/ff_chaya/R_versions/FF_biomarker_CVID_5.24/FF_biomarker_CVID_5.24_datafile.RDS")
ancestry <- read.csv("/drives/drive1/ff_chaya/Ada/PsychChip2b.ancestryPCs.CVID.csv")
cimt <- cimtdataset %>% select("CVID", "Mean_Mean_CIMT")
surveydata <- readRDS("/drives/drive1/ff_chaya/R_versions/FF_allwaves_wY22_CVID_3.24/FF_allwaves_wY22_CVID_3.24_datafile.RDS")
cimt_AND_survey <- surveydata %>% filter(surveydata$CVID %in% cimt$CVID)
cimt_AND_survey <- merge(cimt_AND_survey, cimt, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, biomarker, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, ancestry, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, cimt_le8[, c("CVID", "Age_k7")], by = "CVID", all = TRUE)
cimt_AND_survey <- cimt_AND_survey[!duplicated(as.list(cimt_AND_survey))]
cimt_AND_survey <- cimt_AND_survey[!duplicated(cimt_AND_survey$CVID, )]
cimt_AND_survey$log_Mean_Mean_CIMT <- log(cimt_AND_survey$Mean_Mean_CIMT)
cimt_AND_survey$k7inches <- ((cimt_AND_survey$k7i40)*12 + cimt_AND_survey$k7i41)

cimt_mr <- cimt_AND_survey[!duplicated(cimt_AND_survey$CVID), ] # it is worth checking up on what is being lost
cimt_mr <- cimt_mr %>% filter(CVID %in% surveydata$CVID) %>% filter(CVID %in% cimt$CVID)


rownames(cimt_mr) <- cimt_mr$CVID
cimt_mr$ck7ethrace <- as.factor(cimt_mr$ck7ethrace)
cimt_mr$cm1bsex <- droplevels(cimt_mr$cm1bsex)


sex_dummy <- model.matrix(~ cm1bsex, data = cimt_mr)
sex_dummy <- as.data.frame(sex_dummy[, 2])
colnames(sex_dummy) <- c("cm1bsex2")

cimt_mr <- cbind(cimt_mr, sex_dummy)

# set up covariates, just using Sex and Age
covariates2 <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + cm1bsex, data = cimt_mr)
covariates2 <- t(covariates2)

# set up Model, just using Sex and Age
mod2 <- model.matrix(~ Mean_Mean_CIMT  + PC1 + PC2 + PC3 + PC4 + cm1bsex2, data = cimt_mr)

rows_mr <- rownames(mod2)
covariates2 <- data.frame(covariates2)
colnames(covariates2) <- sub("^.", "", colnames(covariates2))
covariates2 <- covariates2 %>% select(rows_mr)

# read in SNP genotyping data 
snp2_meqtl_fixed <- readRDS("/drives/drive1/ff_chaya/Ayanava/snp2_meqtl_fixed.rds")
snp_meqtl_mr <- snp2_meqtl_fixed

# run the same imputation protocol with TopMed, except with a different sample size than in the meQTLs
num_NAs <- as.data.frame(colSums(is.na(t(snp_meqtl_mr))))
colnames(num_NAs) <- "num"
good_snps <- num_NAs %>% filter(num < 0.05*ncol(snp_meqtl_mr))
snp_ids <- row.names(good_snps)
# filter out SNPs that have more than 5% missing sample values
snp_meqtl_mr <- snp_meqtl_mr[row.names(snp_meqtl_mr) %in% snp_ids, ]
num_Samples_NAs <- as.data.frame(colSums(is.na(snp_meqtl_mr)))
colnames(num_Samples_NAs) <- "num"
good_samples <- num_Samples_NAs %>% filter(num < 0.02*nrow(snp_meqtl_mr))
samp_ids <- row.names(good_samples)
# filter out samples that are missing more than 2% of SNP values 
snp_meqtl_mr <- snp_meqtl_mr[, colnames(snp_meqtl_mr) %in% samp_ids]
# remove the _0 that is at the end of all the SNP names
rownames(snp_meqtl_mr) <- sub("_0$", "", rownames(snp_meqtl_mr)) 

# filter out mitochondria, sex-linked, non-mapping SNPs 
bim <- as.data.frame(new_genetic$map)
# seems as though there are chromosomes that are numbered 23, 24, and 26- I'm removing all those
bim2 <- bim %>% filter(chromosome < 23)
common_snps <- intersect(rownames(bim2), rownames(snp_meqtl_mr))
snp_meqtl_mr_filtered <- snp_meqtl_mr[common_snps, , drop = FALSE]

snp_meqtl_mr_filtered <- snp_meqtl_mr_filtered %>% select(rows_mr)

# only keep the SNPs used in the meQTL analysis 
snp_meqtl_w7_imputed <- readRDS("/drives/drive1/ff_chaya/Ayanava/snp_with_imputed_values.rds")
snp_meqtl_w7_nonas <- na.omit(snp_meqtl_w7_imputed)
snp_meqtl_mr_filtered <- snp_meqtl_mr_filtered[rownames(snp_meqtl_mr_filtered) %in% rownames(snp_meqtl_w7_nonas), ]


# # impute all of the remaining NA values that are present in the SNPs
# imputed_genetic_bim <- fread("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/topmedimputed/genofiltered.bim")
original_genetic_bim <- new_genetic$map
snp_imputed <- fread("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/topmedimputed/cimt.meqtl.CVID.topmed.GRCh38.raw", header = TRUE, stringsAsFactors = FALSE)

# converting GrCh38 imputed data into GrCh37 positions
library(rtracklayer)
chain <- import.chain(("/drives/drive1/ff_chaya/Ayanava/hg38ToHg19.over.chain"))
colnames(imputed_genetic_bim) <- c("chr", "id", "gen_dist", "pos", "allele1", "allele2")
imputed_genetic_bim$chr <- paste0("chr", imputed_genetic_bim$chr)
gr <- GRanges(seqnames = imputed_genetic_bim$chr, ranges = IRanges(start = imputed_genetic_bim$pos, end = imputed_genetic_bim$pos))
lifted <- liftOver(gr, chain)
imputed_genetic_bim$pos_grch37 <- sapply(lifted, function(x) if (length(x) >0) start(x)[1] else NA)
imputed_genetic_bim_grch37 <- na.omit(imputed_genetic_bim)


# # see if the NAs can be imputed in 
snps_with_na_mr <- rownames(snp_meqtl_mr_filtered)[rowSums(is.na(snp_meqtl_mr_filtered)) > 0]
bim_nas_mr <- original_genetic_bim %>% filter(snp.name %in% snps_with_na_mr)
colnames(bim_nas_mr) <- c("chr", "snp.name", "cM", "pos_grch37", "allele1", "allele2")
bim_nas_mr$chr <- paste0("chr", bim_nas_mr$chr)
imputed_na_bim_mr <-semi_join(imputed_genetic_bim_grch37, bim_nas_mr, by = c("chr", "pos_grch37", "allele1"))

snp_imputed <- as.data.frame(snp_imputed)
snp_imputed$IID <- substr(snp_imputed$IID, 1, nchar(snp_imputed$IID) - 4)

# # remove the dups 
snp_imputed_nodups <- snp_imputed[!duplicated(snp_imputed$IID), ]

for (i in 1:nrow(snp_imputed_nodups)) {
  rownames(snp_imputed_nodups)[i] <- snp_imputed_nodups[i,2]
}

snp_impute <- snp_imputed_nodups %>% select(-c(1:6))

# # generated when matching by allele
imputed_merge <- merge(imputed_na_bim_mr, bim_nas_mr, by = c("chr", "pos_grch37"), all.x = TRUE)
imputed_merge1 <- imputed_merge[!duplicated(imputed_merge$snp.name), ]
imputed_merge1 <- imputed_merge1[!duplicated(imputed_merge1$id), ]

# # fixing up the column names of the imputed snp matrix
snp_impute1 <- snp_impute
colnames(snp_impute1) <- sub("..$", "", colnames(snp_impute))
colnames(snp_impute1) <- gsub("\\.", ":", colnames(snp_impute1))

t  <- intersect(colnames(snp_impute1), imputed_merge1$id)

snp_impute_filter <- snp_impute1 %>% select(t)
snp_impute_filter <- data.frame(t(snp_impute_filter))
colnames(snp_impute_filter) <- sub("^.", "", colnames(snp_impute_filter))

# # how many samples are in the imputed and also the original SNP datasets 
matched_ids <- intersect(rownames(snp_impute), colnames(snp_meqtl_mr_filtered)) # 1058 matched CVIDS
snp_impute_filtered <- snp_impute_filter %>% select(matched_ids)

snp_meqtl_mr_impute <- snp_meqtl_mr_filtered %>% select(matched_ids)

# 
# # ensuring that the bim file for the imputed has the same number of SNPs 
imputed_merge1 <- imputed_merge1 %>% filter(imputed_merge1$id %in% t)
# # converting the rownames of the imputed SNPs to the actual SNP names
rownames(snp_impute_filtered) <- imputed_merge1$snp.name[match(rownames(snp_impute_filtered), imputed_merge1$id)]
# 
# swapping out the NAs in the original with the values in the imputed
match_snps <- intersect(rownames(snp_impute_filtered), rownames(snp_meqtl_mr_impute))
original_subset <- snp_meqtl_mr_impute[match_snps, ]
impute_subset <- snp_impute_filtered[match_snps, ]

impute_subset <- impute_subset[match(rownames(original_subset), rownames(impute_subset)), ]
impute_subset <- impute_subset[, match(colnames(original_subset), colnames(impute_subset))]

# swaps the NAs with the imputed values
original_subset[is.na(original_subset)] <- impute_subset[is.na(original_subset)]

# swapping out the imputed rows into the original matrix
original_subset <- original_subset[, match(colnames(original_subset), colnames(snp_meqtl_mr_impute))]

snp_meqtl_mr_imputed <- snp_meqtl_mr_filtered

for (i in 1:nrow(original_subset)) {
  name <- rownames(original_subset)[i]
  index <- which(rownames(snp_meqtl_mr_imputed) == name)
  snp_meqtl_mr_imputed[index, ] <- original_subset[i, ]
}

common_rows <- intersect(rownames(snp_meqtl_mr_imputed), rownames(original_subset))
common_samples <- intersect(colnames(snp_meqtl_mr_imputed), colnames(original_subset))
df1_common <- snp_meqtl_mr_imputed[common_rows, common_samples] 
df2_common <- original_subset[common_rows, common_samples]

df1_common[is.na(df1_common)] <- df2_common[is.na(df1_common)]

snp_meqtl_mr_imputed[common_rows, common_samples] <- df1_common

# remove all remaining SNPs with NAs
snp_meqtl_mr_nonas <- na.omit(snp_meqtl_mr_imputed)

covariates2 <- covariates2[, colnames(covariates2) %in% colnames(snp_meqtl_mr_nonas)]

# preparation for PLINK GWAS

# generating phenotype file
cimt_phenotype <- fam %>% select(c("pedigree", "member"))
colnames(cimt_phenotype) <- c("CVID", "member")
cimt_phenotype <- cimt_phenotype[rownames(cimt_phenotype) %in% colnames(covariates2), ]
cimt_phenotype <- left_join(cimt_phenotype, cimt_AND_survey[, c("CVID", "Mean_Mean_CIMT")], by = "CVID")
colnames(cimt_phenotype) <- c("pedigree", "member", "Mean_Mean_CIMT")
cimt_phenotype <- cimt_phenotype[!duplicated(cimt_phenotype$pedigree, ),]
rownames(cimt_phenotype) <- cimt_phenotype$pedigree

write.table(cimt_phenotype, "/drives/drive1/ff_chaya/Ayanava/cimt_phenotypes.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

s <- read.table("/drives/drive1/ff_chaya/Ayanava/cimt_phenotypes.txt")

# generating the covariates file 
t_covariates2 <- data.frame(t(covariates2))
cov_file <- fam %>% select(c("pedigree", "member"))
cov_file <- merge(cov_file, t_covariates2, by = "row.names", all = TRUE)
rownames(cov_file) <- cov_file$Row.names
cov_file <- cov_file[, -1]
cov_file <- na.omit(cov_file)
cov_file$X.Intercept. <- NULL
cov_file_order <- cov_file[order(factor(cov_file$pedigree, levels = s$V1)), ]
write.table(cov_file_order, "/drives/drive1/ff_chaya/Ayanava/covariates_file.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

covs <-  read.table("/drives/drive1/ff_chaya/Ayanava/covariates_file.txt")


# to make the bim file
bim <- new_genetic$map
bim_new  <- bim[rownames(bim) %in% rownames(snp_meqtl_mr_nonas), ]
bim_new$cM <- 0

write.table(bim_new, file = "/drives/drive1/ff_chaya/Ayanava/ff_psych_cvid_subset_final.bim", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# to make the fam file
fam <- new_genetic$fam
fam_new  <- fam[rownames(fam) %in% colnames(snp_meqtl_mr_nonas), ]
cimt_pheno <- read.table("/drives/drive1/ff_chaya/Ayanava/cimt_phenotypes.txt", header = FALSE)
colnames(cimt_pheno) <- c("pedigree", "member", "Mean_Mean_CIMT")
fam_merged <- merge(fam_new, cimt_pheno, by = "pedigree", all.x = TRUE)
fam_merged$member.y <- NULL
colnames(fam_merged)[2] <- "member"
fam_merged$affected <- fam_merged$Mean_Mean_CIMT
fam_merged$Mean_Mean_CIMT <- NULL

fam_merged_order <- fam_merged[order(factor(fam_merged$pedigree, levels = s$V1)), ]

write.table(fam_merged_order, file = "/drives/drive1/ff_chaya/Ayanava/ff_psych_cvid_subset_final.fam", 
            quote = FALSE, col.names = FALSE, row.names = FALSE)


# make the set of samples and SNPs to keep for filtering the genotyping .bed file data (done in PLINK)
samples <- fam_new[, 1:2]
write.table(samples, "/drives/drive1/ff_chaya/Ayanava/keep_samples.txt", 
            quote = FALSE, row.names = F, col.names = F, sep = "\t")

snp_ids <- bim_new[, 2]
write.table(snp_ids, "/drives/drive1/ff_chaya/Ayanava/keep_snps.txt", 
            quote = FALSE, row.names = F, col.names = F, sep = "\t")



# ensuring that they all have the same matches in the sample names 
pheno_geno <- intersect((cimt_phenotype$member), (fam_new$member))
cov_geno <- intersect((cov_file$pedigree), (fam_new$pedigree))



# trying to use the UK Biobank GWAS summary stats reference 
meqtls_20kb <- readRDS("/drives/drive1/ff_chaya/Ayanava/sig_cis_meqtls_20kb.rds")
meqtls_20kb$se <- abs((meqtls_20kb$beta) / meqtls_20kb$statistic)
colnames(meqtls_20kb) <- c("SNP", "CpG", "statistic", "pval", "FDR", "beta", "se")

gwas_reference <- fread("/drives/drive1/ff_chaya/Ayanava/47gd92pbn9-1/cIMTmean_UKBiobank_MAFandINFOfiltered_Yeung_et_al")
snps_common <- intersect(meqtls_20kb$SNP, gwas_reference$SNP)

# manhattan plot of UK GWAS
library(qqman)
png("/drives/drive1/ff_chaya/Ayanava/Plots/gwas_manhattan_uk.png", width = 6.5, height = 6, units = "in", res = 300)
manhattan(gwas_reference, main = "UK Biobank GWAS", p = "pval")
dev.off()

# only testing the SNPs in common between UK Biobank and FFCWS
gwas_reference_common <- gwas_reference %>% filter(SNP %in% snps_common)
gwas_reference_common <- gwas_reference_common %>% select(SNP, BETA, SE, ALLELE1, ALLELE0, A1FREQ, pval)
colnames(gwas_reference_common) <- c("SNP", "beta.outcome", "se.outcome", 
                                                          "effect_allele.outcome", "other_allele.outcome", 
                                                          "eaf.outcome", "outcome")

meqtls_reference <- meqtls_20kb %>% filter(SNP %in% snps_common)
new_genetic <- read.plink("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/ff_psych_cvid")
bim_mr <- new_genetic$map
colnames(bim_mr) <- c("chr", "SNP", "cM", "position", "effect_allele.exposure", "other_allele.exposure")
meqtls_reference <- merge(meqtls_reference, bim_mr[, c("SNP", "effect_allele.exposure", "other_allele.exposure")], by = "SNP", all.x = TRUE)

# minor allele frequencies generated in PLINK
eaf_plink <- read.table("/drives/drive1/ff_chaya/Ayanava/subset_freq_check.frq", header = TRUE)

meqtls_reference <- merge(meqtls_reference, eaf_plink, by = "SNP", all.x = TRUE)

meqtls_reference_common <- meqtls_reference %>% select(SNP, CpG, beta, se, effect_allele.exposure, other_allele.exposure, 
                                                       MAF, pval)
colnames(meqtls_reference_common) <- c("SNP", "id.exposure", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                                       "other_allele.exposure", "eaf.exposure", "exposure")

gwas_reference_common$id.outcome <- "CIMT"

# harmonizing the GWAS and meQTL data
harmonized_data <- harmonise_data(exposure_dat = meqtls_reference_common, outcome_dat = gwas_reference_common)

# running MR
mr_results <- mr(harmonized_data)


# running the Mendelian Randomization with FFCWS GWAS

# GWAS run using PLINK
gwas_results <- read.table("/drives/drive1/ff_chaya/Ayanava/nobatch_gwas_results.assoc.linear", header = TRUE)

gwas_nonas <- gwas_results[!is.na(gwas_results$BETA) & !is.na(gwas_results$P), ]
gwas_summarized <- gwas_nonas %>% filter(TEST == "ADD")

gwas_ff_man <- gwas_summarized

# Manhattan plot
png("/drives/drive1/ff_chaya/Ayanava/Plots/gwas_manhattan_ff.png", width = 6.5, height = 6, units = "in", res = 300)
manhattan(gwas_ff_man, ylim = c(0, 8), main = "FFCWS GWAS")
dev.off()

# seeing how many significant SNPs in both cohorts
sig_gwas_ff <- gwas_summarized %>% filter(P < 0.00000005)
sig_gwas_uk <- gwas_reference %>% filter(pval < 0.00000005)


common_snps <- intersect(gwas_summarized$SNP, meqtls_20kb$SNP)

gwas_common <- gwas_summarized %>% filter(SNP %in% common_snps)
meqtls_common <- meqtls_20kb %>% filter(SNP %in% common_snps)

# using Future of Families Data
meqtls_common_FF <- merge(meqtls_common, bim_mr[, c("SNP", "effect_allele.exposure", "other_allele.exposure")], by = "SNP", all.x = TRUE)
meqtls_common_FF <- merge(meqtls_common_FF, eaf_plink, by = "SNP", all.x = TRUE)
meqtls_common_FF$statistic <- NULL
meqtls_common_FF$FDR <- NULL
meqtls_common_FF <- meqtls_common_FF %>% select(SNP, CpG, pval, beta, se, effect_allele.exposure, other_allele.exposure, MAF)
colnames(meqtls_common_FF) <- c("SNP", "id.exposure", "exposure", "beta.exposure", "se.exposure", "effect_allele.exposure", 
                                "other_allele.exposure", "eaf.exposure")

gwas_common$id.outcome <- "CIMT"
gwas_common_FF <- merge(gwas_common, bim_mr[, c("SNP", "effect_allele.exposure", "other_allele.exposure")], by = "SNP", all.x = TRUE)
gwas_common_FF <- merge(gwas_common_FF, eaf_plink[, c("SNP", "MAF")], by = "SNP", all.x = TRUE)
gwas_common_FF$CHR <- NULL
gwas_common_FF$BP <- NULL
gwas_common_FF$A1 <- NULL
gwas_common_FF$TEST <- NULL
gwas_common_FF$NMISS <- NULL
gwas_common_FF$se.outcome <- abs((gwas_common_FF$BETA) / gwas_common_FF$STAT)
gwas_common_FF$STAT <- NULL
colnames(gwas_common_FF) <- c("SNP", "beta.outcome","outcome", "id.outcome",  
                                     "effect_allele.outcome", "other_allele.outcome", 
                                     "eaf.outcome", "se.outcome")


# running harmonizing and MR protocol
harmonized_data2 <- harmonise_data(exposure_dat = meqtls_common_FF, outcome_dat = gwas_common_FF)
saveRDS(harmonized_data2, "/drives/drive1/ff_chaya/Ayanava/harmonized_data_ff_gwas.rds")

mr_results2 <- mr(harmonized_data2)




# comparing the results of the two GWASes 
mr_results1 <- readRDS("/drives/drive1/ff_chaya/Ayanava/mr_ff_results_20kb.rds")
mr_results2 <- readRDS("/drives/drive1/ff_chaya/Ayanava/mr_uk_results_20kb.RDS")


common_cpgs <- intersect(mr_results1$id.exposure, mr_results2$id.exposure)
common_mr1 <- mr_results1[mr_results1$id.exposure %in% common_cpgs, ]
common_mr1$adj.pval <- p.adjust(common_mr1$pval, method = "BH")
common_mr2 <- mr_results2[mr_results2$id.exposure %in% common_cpgs, ]
common_mr2$adj.pval <- p.adjust(common_mr2$pval, method = "BH")

# how many are the same at adj.p<0.05 
sig_mr1 <- common_mr1 %>% filter(adj.pval < 0.05)
sig_mr2 <- common_mr2 %>% filter(adj.pval < 0.05)
common_sig_cpgs <- intersect(sig_mr1$id.exposure, sig_mr2$id.exposure) 

common_sig_mr1 <- sig_mr1[sig_mr1$id.exposure %in% common_sig_cpgs, ]
common_sig_mr2 <- sig_mr2[sig_mr2$id.exposure %in% common_sig_cpgs, ]

# identify location of the significant CpGs relative to an island
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno)
anno <- anno %>% select(Name, chr, pos, UCSC_RefGene_Name, Relation_to_Island)
anno_cpgs_mr1 <- common_sig_mr1 %>% left_join(anno, by = c("id.exposure" ="Name"))
anno_cpgs_mr1 <- anno_cpgs_mr1[!duplicated(anno_cpgs_mr1$id.exposure), ]
table(anno_cpgs_mr1$Relation_to_Island)

anno_cpgs_mr2 <- sig_mr1 %>% left_join(anno, by = c("id.exposure" ="Name"))
anno_cpgs_mr2 <- anno_cpgs_mr2[!duplicated(anno_cpgs_mr2$id.exposure), ]
table(anno_cpgs_mr2$Relation_to_Island)

anno_cpgs_mr3 <- sig_mr2 %>% left_join(anno, by = c("id.exposure" ="Name"))
anno_cpgs_mr3 <- anno_cpgs_mr3[!duplicated(anno_cpgs_mr3$id.exposure), ]
table(anno_cpgs_mr3$Relation_to_Island)


# analysis of important genes identified from MR
target_genes <- c("ZFP57", "HOXA3", "HOXA-AS2", "MGST3", "HCG27","TNXB", "FADS1", "FADS2", "C6orf10")
pattern <- paste(target_genes, collapse = "|")

mr_genes_ff <- readRDS("/drives/drive1/ff_chaya/Ayanava/mr_genes_ff.rds")
mr_genes_ff <- as.data.frame(mr_genes_ff)
mr_genes_ff <- mr_genes_ff %>% filter(grepl(pattern, .data$gene))

# unnesting the genes, and merging CpGs 
mr_genes_ff <- mr_genes_ff %>% separate_rows(gene, sep = ";")
mr_genes_ff <- mr_genes_ff %>% separate_rows(CpGs, sep = ";")
mr_genes_ff <- mr_genes_ff %>% group_by(gene) %>% 
  summarise(CpGs = paste(unique(CpGs), collapse  = ","), count = sum(count), .groups = "drop")

mr_genes_ff <- mr_genes_ff %>% filter(grepl(pattern, .data$gene))


mr_results_ff <- mr_results1
mr_results_ff$adj.pval <- p.adjust(mr_results_ff$pval, method = "BH")
mr_results_ff <- mr_results_ff %>% filter(adj.pval < 0.05)

reference_cpgs1 <- reference[rownames(reference) %in% mr_results1$id.exposure, ]
reference_cpgs1$id.exposure <- rownames(reference_cpgs1)

# function to generate forest plots for the CpG effect sizes of each gene of interest
make_forest <- function(gene_index, gene_list, reference, mr_results, title) {
  cpgs_gene <- unlist(strsplit(gene_list$CpGs[gene_index], split = ","))
  cpg_gene_mr <- mr_results[mr_results$id.exposure %in% cpgs_gene, ]
  cpg_gene_mr <- merge(cpg_gene_mr, reference[, c("seqnames", "start", "id.exposure")], by = "id.exposure")
  cpg_gene_mr <- cpg_gene_mr %>% group_by(id.exposure) %>%
    filter(adj.pval == min(adj.pval)) %>% 
    ungroup()
  cpg_gene_mr <- cpg_gene_mr %>% arrange(start)
  cpg_gene_mr$id.exposure <- factor(cpg_gene_mr$id.exposure, levels = rev(unique(cpg_gene_mr$id.exposure)))
  p <- ggplot(cpg_gene_mr, aes(x = b, y = id.exposure)) + 
    geom_point(aes(color = ifelse(b < 0, "below_zero", "above_zero")), size = 4, shape = 16) + 
    geom_errorbarh(aes(xmin = b - 1.96*se, 
                       xmax = b + 1.96*se), height = 0.2) + 
    labs(x = "Beta estimate", y = "CpG Sites", title = title) + 
    theme(axis.text.x = element_text(size = 12, face = "bold"), 
          axis.text.y = element_text(size = 12, face = "bold"), 
          axis.title.x = element_text(size = 14, face = "bold"), 
          axis.title.y = element_text(size = 14, face = "bold"), 
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5), 
          plot.margin = margin(r = 10)) + 
    geom_vline(xintercept = 0, size = 0.5) + 
    scale_color_manual(values = c("below_zero" = "blue", "above_zero" = "red"), guide = "none")
  print(p)
  title_part1 <- "/drives/drive1/ff_chaya/Ayanava/Plots/forest_"
  title_part2 <- ".png"
  graph_title <- paste0(title_part1, title, title_part2, sep = "")
  ggsave(graph_title, dpi = 600, width = 6.5, height = 6)
  return(cpg_gene_mr)
}

# make forest plots of each gene, and save the locations of each CpG
gene <- "TSBP1"
tsbp1 <- make_forest(1, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "FADS1_FADS2"
fads12 <- make_forest(2, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "HCG27"
hcg27 <- make_forest(4, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "HOXA-AS2"
hoxa_as2 <- make_forest(5, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "HOXA3"
hoxa3 <- make_forest(6, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "MGST3"
mgst3 <- make_forest(7, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "TNXB"
tnxb <- make_forest(8, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)
gene <- "ZFP57"
zfp57 <- make_forest(9, mr_genes_ff, reference_cpgs1, mr_results_ff, gene)


# check the logFC values of these CpGs from the EWAS
test_meqtl <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_meqtl_results.rds")
tsbp1_ewas <- test_meqtl[rownames(test_meqtl) %in% tsbp1$id.exposure, ]
fads12_ewas <- test_meqtl[rownames(test_meqtl) %in% fads12$id.exposure, ]
hcg27_ewas <- test_meqtl[rownames(test_meqtl) %in% hcg27$id.exposure, ]
hoxa_as2_ewas <- test_meqtl[rownames(test_meqtl) %in% hoxa_as2$id.exposure, ]
hoxa3_ewas <- test_meqtl[rownames(test_meqtl) %in% hoxa3$id.exposure, ]
mgst3_ewas <- test_meqtl[rownames(test_meqtl) %in% mgst3$id.exposure, ]
tnxb_ewas <-test_meqtl[rownames(test_meqtl) %in% tnxb$id.exposure, ]
zfp57_ewas <- test_meqtl[rownames(test_meqtl) %in% zfp57$id.exposure, ]




