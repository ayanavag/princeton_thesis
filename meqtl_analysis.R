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

library(data.table)
library(GenomicRanges)
library(snpStats)
library(MatrixEQTL)

# reading in the methylation data for Y22
ids_convert3 <- read.csv("/drives/drive1/ff_chaya/SampleSheets/Combined_EPIC.SS.correctedSex.csv")
methyl_y22 <- readRDS("/drives/drive1/ff_chaya/beta_CVH.blood.RDS")
methyl_y22 <- as.data.frame(methyl_y22)
badsamples <- read.table("/drives/drive1/ff_chaya/SampleSheets/y22.blood_and_saliva_badsamples_list.Sep24.txt")
badsamples <- data.frame(badsamples[-1, ])
colnames(badsamples) <- "Array_ID"
dupsamples <- read.csv("/drives/drive1/ff_chaya/SampleSheets/y22.duplicates_to_remove.csv")

# REMOVE BAD SAMPLES
for (i in 1:nrow(badsamples)) {
  id <- badsamples[i, ]
  if (id %in% colnames(methyl_y22)) {
    methyl_y22 <- methyl_y22 %>% select(-c(id))
  }
}

for (i in 1:nrow(dupsamples)) {
  id <- dupsamples[i, 2]
  if (id %in% colnames(methyl_y22)) {
    methyl_y22 <- methyl_y22 %>% select(-c(id))
  }
}

# read in the batch variables for blood
batch_blood <- readRDS("/drives/drive1/ff_chaya/Ayanava/ctrlsva_epic_blood.RDS")
batch_blood <- as.data.frame(batch_blood)
blood_int <- intersect(rownames(batch_blood), colnames(methyl_y22))
batch_blood <- batch_blood %>% filter(rownames(batch_blood) %in% blood_int)


for (i in 1:nrow(batch_blood)) {
  match <- ids_convert3 %>% 
    filter(Collection == "PaxGeneDNA") %>% 
    filter(Array_ID == rownames(batch_blood)[i])
  newname <- match[1,5]
  rownames(batch_blood)[i] <- newname
}

for (i in 1:length(methyl_y22)) {
  match <- ids_convert3 %>% 
    filter(Collection == "PaxGeneDNA") %>% 
    filter(Array_ID == colnames(methyl_y22)[i])
  newname <- match[1,5]
  colnames(methyl_y22)[i] <- newname
}

# convert to M-value
methyl_w7_epic <- log2((methyl_y22) / (1 - methyl_y22))


# exclude non-autosomal chromosomes, but keep masked CpGs
reference <- as.data.frame(readRDS("/drives/drive1/geneticData/referenceDataSets/Methylation/Zhou/EPIC.hg19.manifest.rds"))
reference$ID <- rownames(reference)
refSub_mr <- c(reference[reference$seqnames == "chrX","ID"],
            reference[reference$seqnames == "chrY","ID"],
            reference[reference$seqnames == "chrM","ID"])
methyl_w7_epic_mr <- methyl_w7_epic[!rownames(methyl_w7_epic) %in% refSub_mr,]

# read in PLINK file of genotyping data
new_genetic <- read.plink("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/ff_psych_cvid")

# setting up the dataframe from the PLINK file
library(data.table)
snp2 <- fread("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/ff_psych_cvid_recode.raw", header = TRUE, stringsAsFactors = FALSE)

snp2 <- as.data.frame(snp2)
for (i in 1:nrow(snp2)) {
  rownames(snp2)[i] <- snp2[i,1]
}

snp2_NAs <- data.frame(colSums(is.na(snp2)))
snp2_meqtl <- snp2[, -c(1:6)]

snp2_meqtl_fixed <- t(snp2_meqtl)
snp2_meqtl_fixed <- as.data.frame(snp2_meqtl_fixed)


# removing all of the SNPs that are diff names but same loci
bim <- new_genetic$map
dup_b <- bim %>% group_by(chromosome, position) %>%
  filter(n() > 1) %>% 
  ungroup
dup_b <- dup_b %>% arrange(chromosome, position)

dup_groups <- split(dup_b$snp.name, list(dup_b$chromosome, dup_b$position))
dup_groups <- Filter(function(x) length(x) > 0, dup_groups)
select_snps <- sapply(dup_groups, function(snp_group) {
  na.counts = sapply(snp_group, function(snp) sum(is.na(snp2_meqtl_fixed[[snp]])))
  snp_group[which.min(na.counts)]
})
select_snps <- unique(select_snps)

# remove the _0 that is at the end of all the SNP names
rownames(snp2_meqtl_fixed) <- sub("_0$", "", rownames(snp2_meqtl_fixed)) 

allsnp <- rownames(snp2_meqtl_fixed)
dupsnps <- dup_b$snp.name
nodupssnps <- setdiff(allsnp, dupsnps)

snps_keep <- c(nodupssnps, select_snps)

snp2_meqtl_fixed <- snp2_meqtl_fixed[rownames(snp2_meqtl_fixed) %in% snps_keep, ]

# filter for samples that are present in both the methylation and genotyping datasets
common_samples_w7 <- intersect(colnames(snp2_meqtl_fixed), colnames(methyl_w7_epic_mr))
snp2_meqtl_w7 <- snp2_meqtl_fixed %>% select(common_samples_w7)
methyl_meqtl_w7 <- methyl_w7_epic_mr[, common_samples_w7]



# generating covariates
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

# adding in the cell proportions
cp_blood <- estimateLC(methyl_w7_epic_mr, ref = "Salas")
cp_blood$CVID <- colnames(methyl_w7_epic_mr)
cols <- colnames(methyl_w7_epic_mr)
cimt_AND_survey <- cimt_AND_survey %>% filter(CVID %in% cols)

cimt_AND_survey <- cimt_AND_survey[!duplicated(cimt_AND_survey$CVID), ] # it is worth checking up on what is being lost
cimt_AND_survey <- cimt_AND_survey[match(cols, cimt_AND_survey$CVID), ]
cimt_AND_survey <- merge(cimt_AND_survey, cp, by = "CVID", all = TRUE)

cimt_methyl <- cimt_AND_survey
batch2 <- batch_blood
colnames(batch2) <- c("sv1", "sv2", "sv3", "sv4")
batch2$CVID <- rownames(batch2)
cimt_methyl <- merge(cimt_methyl, batch2, by = "CVID", all = TRUE)


cimt_methyl <- cimt_methyl %>% filter(CVID %in% common_samples_w7)

# remove the duplicates
cimt_methyl <- cimt_methyl[!duplicated(cimt_methyl$CVID), ]
rownames(cimt_methyl) <- cimt_methyl$CVID
cimt_methyl$ck7ethrace <- as.factor(cimt_methyl$ck7ethrace)
cimt_methyl$cm1bsex <- droplevels(cimt_methyl$cm1bsex)


# dummy variables for sex
sex_dummy <- model.matrix(~ cm1bsex, data = cimt_methyl)
sex_dummy <- as.data.frame(sex_dummy[, 2])
colnames(sex_dummy) <- c("cm1bsex2")
# cimt_methyl <- cbind(cimt_methyl, race_dummy)
cimt_methyl <- cbind(cimt_methyl, sex_dummy)

# generate covariates for Model 2 for meQTLs
covariates <- model.matrix(~ Age_k7 + sv1 + sv2 + sv3 + sv4+ PC1 + PC2 + PC3 + PC4 + GR + CD4+ CD8+ cm1bsex2, data = cimt_methyl)
covariates <- t(covariates)
covariates <- as.data.frame(covariates)


# generate Model 2 matrix
mod <- model.matrix(~ Mean_Mean_CIMT + Age_k7 + sv1 + sv2 + sv3 + sv4+PC1 + PC2+PC3+PC4 + GR + CD4+ CD8+ cm1bsex2, data = cimt_methyl)
rows <- rownames(mod)

# convert the methyl and the SNP data to only have the same number of samples as the covariates
methyl_meqtl_w7 <- methyl_meqtl_w7 %>% select(rows)
snp2_meqtl_w7 <- snp2_meqtl_w7 %>% select(rows)
covariates <- covariates %>% select(rows)

saveRDS(methyl_meqtl_w7, "/drives/drive1/ff_chaya/Ayanava/methyl_meqtl_w7.rds")
saveRDS(snp2_meqtl_w7, "/drives/drive1/ff_chaya/Ayanava/snp2_meqtl_w7.rds")
saveRDS(covariates, "/drives/drive1/ff_chaya/Ayanava/covariates_ctrlsva.rds")

# to remove the intercept, since it is collinear
covariates <- as.data.frame(covariates[-1, ])

# setting up last steps of methylation data
Meth_data <- SlicedData$new()
Meth_data$CreateFromMatrix(as.matrix(methyl_meqtl_w7))


# setting up SNP data 
num_NAs <- as.data.frame(colSums(is.na(t(snp2_meqtl_w7))))
colnames(num_NAs) <- "num"
good_snps <- num_NAs %>% filter(num < 0.05*ncol(snp2_meqtl_w7))
snp_ids <- row.names(good_snps)
# filter out SNPs that have more than 5% missing sample values
snp2_meqtl_w7 <- snp2_meqtl_w7[row.names(snp2_meqtl_w7) %in% snp_ids, ]
num_Samples_NAs <- as.data.frame(colSums(is.na(snp2_meqtl_w7)))
colnames(num_Samples_NAs) <- "num"
good_samples <- num_Samples_NAs %>% filter(num < 0.02*nrow(snp2_meqtl_w7))
samp_ids <- row.names(good_samples)
# filter out samples that are missing more than 2% of SNP values 
snp2_meqtl_w7 <- snp2_meqtl_w7[, colnames(snp2_meqtl_w7) %in% samp_ids]

# filter out mitochondria, sex-linked, non-mapping SNPs 
bim <- as.data.frame(new_genetic$map)
# remove all SNPs on non-autosomal chromosomes
bim2 <- bim %>% filter(chromosome < 23)
common_snps <- intersect(rownames(bim2), rownames(snp2_meqtl_w7))
snp2_meqtl_w7_filtered <- snp2_meqtl_w7[common_snps, , drop = FALSE]

# impute all of the remaining NA values that are present in the SNPs from the TopMed imputation
imputed_genetic_bim <- fread("/drives/drive1/ff_chaya/PsychChipRelease2b_052520_cvid/topmedimputed/genofiltered.bim")
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

saveRDS(imputed_genetic_bim_grch37, "/drives/drive1/ff_chaya/Ayanava/imputed_hg19_snps.bim")

# see if the NAs can be imputed in 
snps_with_na <- rownames(snp2_meqtl_w7_filtered)[apply(snp2_meqtl_w7_filtered, 1, function(x) any(is.na(x)))]
bim_nas <- original_genetic_bim %>% filter(snp.name %in% snps_with_na)
colnames(bim_nas) <- c("chr", "snp.name", "cM", "pos_grch37", "allele1", "allele2")
bim_nas$chr <- paste0("chr", bim_nas$chr)
# match by allele, chr, and position
imputed_na_bim <-semi_join(imputed_genetic_bim_grch37, bim_nas, by = c("chr", "pos_grch37", "allele1"))

snp_imputed <- data.frame(snp_imputed)
snp_imputed$IID <- substr(snp_imputed$IID, 1, nchar(snp_imputed$IID) - 4)

# remove the dups 
snp_imputed_nodups <- snp_imputed[!duplicated(snp_imputed$IID), ]

for (i in 1:nrow(snp_imputed_nodups)) {
  rownames(snp_imputed_nodups)[i] <- snp_imputed_nodups[i,2]
}

snp_impute <- snp_imputed_nodups %>% select(-c(1:6))


imputed_merge <- merge(imputed_na_bim, bim_nas, by = c("chr", "pos_grch37"), all.x = TRUE)
imputed_merge1 <- imputed_merge[!duplicated(imputed_merge$snp.name), ]
imputed_merge1 <- imputed_merge1[!duplicated(imputed_merge1$id), ]

# fixing up the column names of the imputed snp matrix
snp_impute1 <- snp_impute
colnames(snp_impute1) <- sub("..$", "", colnames(snp_impute))
colnames(snp_impute1) <- gsub("\\.", ":", colnames(snp_impute1))

# find intersection of samples between the genotyping and the imputed datasets
t  <- intersect(colnames(snp_impute1), imputed_merge1$id)

# contains 58,187 SNPs, and the 1058 samples
snp_impute_filter <- snp_impute1 %>% select(t)
snp_impute_filter <- data.frame(t(snp_impute_filter))
colnames(snp_impute_filter) <- sub("^.", "", colnames(snp_impute_filter))

# how many samples are in the imputed and also the original SNP datasets 
matched_ids <- intersect(rownames(snp_impute), colnames(snp2_meqtl_w7_filtered)) # 1058 matched CVIDS
snp_impute_filtered <- snp_impute_filter %>% select(matched_ids)

snp2_meqtl_w7_filtered <- snp2_meqtl_w7_filtered %>% select(matched_ids)


# ensuring that the bim file for the imputed has the same number of SNPs 
imputed_merge1 <- imputed_merge1 %>% filter(imputed_merge1$id %in% t)
# converting the rownames of the imputed SNPs to the actual SNP names
rownames(snp_impute_filtered) <- imputed_merge1$snp.name[match(rownames(snp_impute_filtered), imputed_merge1$id)]

# swapping out the NAs in the original with the values in the imputed
match_snps <- intersect(rownames(snp_impute_filtered), rownames(snp2_meqtl_w7_filtered))
original_subset <- snp2_meqtl_w7_filtered[match_snps, ]
impute_subset <- snp_impute_filtered[match_snps, ]

impute_subset <- impute_subset[match(rownames(original_subset), rownames(impute_subset)), ]
impute_subset <- impute_subset[, match(colnames(original_subset), colnames(impute_subset))]

# swaps the NAs with the imputed values
original_subset[is.na(original_subset)] <- impute_subset[is.na(original_subset)]

# swapping out the imputed rows into the original matrix 
original_subset <- original_subset[, match(colnames(original_subset), colnames(snp2_meqtl_w7_filtered))]

snp2_meqtl_imputed <- snp2_meqtl_w7_filtered

for (i in 1:nrow(original_subset)) {
  name <- rownames(original_subset)[i]
  index <- which(rownames(snp2_meqtl_imputed) == name)
  snp2_meqtl_imputed[index, ] <- original_subset[i, ]
}

# remove any remaining SNPs with NAs
snp_meqtl_w7_nonas <- na.omit(snp_meqtl_w7_imputed)


SNP_data <- SlicedData$new()
SNP_data$CreateFromMatrix(as.matrix(snp_meqtl_w7_nonas)) 

# generate the SNP position dataframe for cis-meQTLs
snp_pos <- new_genetic$map
snp_pos <- snp_pos[rownames(snp_pos) %in% rownames(snp_meqtl_w7_nonas), ]
snp_pos <- snp_pos %>% select(chromosome, position)
snp_pos$snpid <- rownames(snp_pos)
snp_pos <- snp_pos[, c("snpid", "chromosome", "position")]
colnames(snp_pos) <- c("snpid", "chr", "pos")

# generate the CpG position dataframe for cis-meQTLs (labelled gene for the ease of inputting into the meQTL package)
gene_pos <- reference[rownames(reference) %in% rownames(methyl_meqtl_w7), ]
gene_pos <- gene_pos %>% select(seqnames, start, end)
gene_pos$geneid <- rownames(gene_pos)
gene_pos <- gene_pos[, c("geneid", "seqnames", "start", "end")]
colnames(gene_pos) <- c("geneid", "chr", "left", "right")
gene_pos$chr <- as.numeric(gsub("chr", "", gene_pos$chr))

cis_methyl_meqtl_w7 <- methyl_meqtl_w7[rownames(methyl_meqtl_w7) %in% rownames(gene_pos), ]



# last clean up steps to ensure equal number of samples
cvids <- colnames(snp_meqtl_w7_nonas)
methyl_meqtl_w7 <- methyl_meqtl_w7 %>% select(cvids)
cis_methyl_meqtl_w7 <- cis_methyl_meqtl_w7 %>% select(cvids)
covariates <- covariates %>% select(cvids)
covars <- SlicedData$new()
covars$CreateFromMatrix(as.matrix(covariates))


# whole meQTL analysis controlling for Age + Sex + Race + CP + Batch
for (i in seq(1, nrow(methyl_meqtl_w7), by = 10000)) {
  subset_meth <- methyl_meqtl_w7[i:(min(i+9999, nrow(methyl_meqtl_w7))), ]
  Meth <- SlicedData$new()
  Meth$CreateFromMatrix(as.matrix(subset_meth)) 
  results <- Matrix_eQTL_main(snps = SNP_data, gene = Meth, cvrt = covars, 
                              output_file_name = NULL, 
                              pvOutputThreshold = 1e-5, useModel = modelLINEAR, errorCovariance = numeric(), 
                              verbose = TRUE)
  saveRDS(results, paste0("/drives/drive1/ff_chaya/Ayanava/meQTL_batch_new_", i, ".rds"))
}



# 20kb threshold cis-meQTL analysis controlling for Age + Sex + Race + CP + Batch
for (i in seq(1, nrow(cis_methyl_meqtl_w7), by = 10000)) {
  subset_meth <- cis_methyl_meqtl_w7[i:(min(i+9999, nrow(cis_methyl_meqtl_w7))), ]
  Meth <- SlicedData$new()
  Meth$CreateFromMatrix(as.matrix(subset_meth)) 
  results <- Matrix_eQTL_main(snps = SNP_data, gene = Meth, cvrt = covars, 
                              output_file_name = NULL, 
                              pvOutputThreshold.cis = 1e-5, snpspos = snp_pos, 
                              genepos = gene_pos, cisDist = 20000, useModel = modelLINEAR, errorCovariance = numeric(), 
                              verbose = TRUE)
  saveRDS(results, paste0("/drives/drive1/ff_chaya/Ayanava/cis_meQTL_batch_new_", i, ".rds"))
}


# cis-meQTL analysis, using a bigger window of 1 million base-pairs
for (i in seq(1, nrow(cis_methyl_meqtl_w7), by = 10000)) {
  subset_meth <- cis_methyl_meqtl_w7[i:(min(i+9999, nrow(methyl_meqtl_w7))), ]
  Meth <- SlicedData$new()
  Meth$CreateFromMatrix(as.matrix(subset_meth)) 
  results <- Matrix_eQTL_main(snps = SNP_data, gene = Meth, cvrt = covars, 
                              output_file_name = NULL, 
                              pvOutputThreshold.cis = 1e-5, snpspos = snp_pos, 
                              genepos = gene_pos, cisDist = 1000000, useModel = modelLINEAR, errorCovariance = numeric(), 
                              verbose = TRUE)
  saveRDS(results, paste0("/drives/drive1/ff_chaya/Ayanava/cis_meQTL_batch_new_1mb_", i, ".rds"))
}



# extract all of the significant cis-meQTLs into one dataframe
file_indices <- seq(1, 839001, by = 10000)
file_paths_cis <- paste0("/drives/drive1/ff_chaya/Ayanava/cis_meQTL_batch_new_", file_indices, ".rds")

df_list <- lapply(file_paths_cis, function(file) {
  obj <- readRDS(file)
  obj$cis$eqtls
})

set_cis_meqtls <- bind_rows(df_list)


# extract all of the significant cis-meQTLs 1mb distance into one dataframe
file_paths_1mb <- paste0("/drives/drive1/ff_chaya/Ayanava/cis_meQTL_batch_new_1mb_", file_indices, ".rds")

df_list_1mb <- lapply(file_paths_1mb, function(file) {
  obj <- readRDS(file)
  obj$cis$eqtls
})

set_cis_meqtls_1mb <- bind_rows(df_list_1mb)


# extract all of the significant meQTLs into one dataframe 
file_paths_trans <- paste0("/drives/drive1/ff_chaya/Ayanava/meQTL_batch_new_", file_indices, ".rds")

df_list_trans <- lapply(file_paths_trans, function(file) {
  obj <- readRDS(file)
  obj$all$eqtls
})

for (i in 1:length(file_paths_trans)) {
  r <- readRDS(file_paths_trans[i])
  eqtls <- r$all$eqtls
  df_list_trans[[i]] <- (eqtls)
}

set_trans_meqtls <- bind_rows(df_list_trans)



sig_cis_meqtls2 <- set_cis_meqtls %>% filter(pvalue < 0.00000005 & FDR < 0.05) # 333,477 cis-meQTLs

sig_cis_meqtls_1mb <- set_cis_meqtls_1mb %>% filter(pvalue < 0.00000005 & FDR < 0.05)

sig_trans_meqtls <- set_trans_meqtls %>% filter(pvalue < 0.00000005 & FDR < 0.05)

saveRDS(sig_cis_meqtls2, "/drives/drive1/ff_chaya/Ayanava/sig_cis_meqtls_20kb.rds")
saveRDS(sig_cis_meqtls_1mb, "/drives/drive1/ff_chaya/Ayanava/sig_cis_meqtls_1mb.rds")
saveRDS(sig_trans_meqtls, "/drives/drive1/ff_chaya/Ayanava/all_sig_meqtls.rds")


