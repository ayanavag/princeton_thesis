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

source("/drives/drive1/ff_chaya/Ayanava/limma_functions.R")
# load(file = "/drives/drive1/ff_chaya/Ayanava/main_environment.RData")

# read in saliva methylation data (for Y9 and Y15)
new_methyl <- readRDS("/drives/drive1/methylation_processed_data/BetaMatrix.saliva.n2627.Sep16.RDS")
new_methyl <- data.frame(new_methyl, check.names = FALSE)
ids_convert <- read.csv("/drives/drive1/ff_chaya/Ayanava/Combined_EPIC_with_y22.arrayID.CVID.csv")
ids_convert2 <- read.table("/drives/drive1/ff_chaya/SampleSheets/epic_ss_cvid_nodups.txt")
ids_convert3 <- read.csv("/drives/drive1/ff_chaya/SampleSheets/Combined_EPIC.SS.correctedSex.csv")


# read in the batch variables created by Eve using ctrlsva
batch_blood <- readRDS("/drives/drive1/ff_chaya/Ayanava/ctrlsva_epic_blood.RDS")
batch_blood <- as.data.frame(batch_blood)
batch_saliva <- readRDS("/drives/drive1/ff_chaya/Ayanava/ctrlsva_epic_saliva.RDS")
batch_saliva <- as.data.frame(batch_saliva)


# remove duplicate samples in the batch variables
dup_ids <- read.table("/drives/drive1/ff_chaya/Ayanava/dup_ids.txt")

dups_to_keep <- read.table("/drives/drive1/ff_chaya/Ayanava/dups_to_keep.txt")

dup_ids_keep <- dup_ids[dup_ids$Array_ID %in% dups_to_keep$Array_ID, ]

diff <- setdiff(ids_convert3$Array_ID, dup_ids$Array_ID)
ids_to_keep <- ids_convert3[ids_convert3$Array_ID %in% diff,]
ids_to_keep <- rbind(ids_to_keep, dup_ids_keep)


batch_saliva <- batch_saliva[rownames(batch_saliva) %in% ids_to_keep$Array_ID, ]

for (i in 1:nrow(batch_saliva)) {
  match <- ids_to_keep %>% 
    filter(Collection == "OGR600") %>% 
    filter(Array_ID == rownames(batch_saliva)[i])
  newname <- paste(match[1,5], match[1,4], sep = "")
  rownames(batch_saliva)[i] <- newname
}


# filter for the same samples between the saliva data and saliva batch 
new_methyl <- new_methyl[, colnames(new_methyl) %in% ids_to_keep$Array_ID]

for (i in 1:length(new_methyl)) {
  match <- ids_convert3 %>% 
    filter(Collection == "OGR600") %>% 
    filter(Array_ID == colnames(new_methyl)[i])
  newname <- paste(match[1,5], match[1,4], sep = "")
  colnames(new_methyl)[i] <- newname
}



# read in the blood methylation data 
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

# convert saliva methylation to M-value 
mval_meth <- log2((new_methyl) / (1 - new_methyl))

# filter for w5 (Y9) and w6 (Y15)
methyl_w5_epic <- mval_meth %>% select(ends_with("FC09"))
methyl_w6_epic <- mval_meth %>% select(ends_with("FC15"))

# w5_cols <- grep("FC09", names(mval_meth), value = FALSE)
# methyl_w5_epic <- mval_meth[, w5_cols, drop = FALSE] 
# w6_cols <- grep("FC15", names(mval_meth), value = FALSE)
# methyl_w6_epic <- mval_meth[, w6_cols, drop = FALSE] 

# convert blood methylation to M-value (w7 : Y22)
methyl_w7_epic <- log2((methyl_y22) / (1 - methyl_y22))

for ( i in 1:ncol(methyl_w5_epic)){
  colnames(methyl_w5_epic)[i] <-  sub("FC09", "", colnames(methyl_w5_epic)[i])
}
for ( i in 1:ncol(methyl_w6_epic)){
  colnames(methyl_w6_epic)[i] <-  sub("FC15", "", colnames(methyl_w6_epic)[i])
}

# cleaning up the batch variables for saliva
batch_saliva <- as.data.frame(t(batch_saliva))
batch_saliva_w5 <- batch_saliva %>% select(ends_with("FC09"))
batch_saliva_w6 <- batch_saliva %>% select(ends_with("FC15"))

for ( i in 1:ncol(batch_saliva_w5)){
  colnames(batch_saliva_w5)[i] <-  sub("FC09", "", colnames(batch_saliva_w5)[i])
}
for ( i in 1:ncol(methyl_w6_epic)){
  colnames(batch_saliva_w6)[i] <-  sub("FC15", "", colnames(batch_saliva_w6)[i])
}

batch_saliva_w5 <- as.data.frame(t(batch_saliva_w5))
batch_saliva_w6 <- as.data.frame(t(batch_saliva_w6))


# exclude bad probes and sex chromosomes
reference <- as.data.frame(readRDS("/drives/drive1/geneticData/referenceDataSets/Methylation/Zhou/EPIC.hg19.manifest.rds"))
reference$ID <- rownames(reference)
refSub <- c(reference[reference$MASK_general==TRUE,"ID"],
            reference[reference$seqnames == "chrX","ID"],
            reference[reference$seqnames == "chrY","ID"],
            reference[reference$seqnames == "chrM","ID"])
methyl_w5_epic <- methyl_w5_epic[!rownames(methyl_w5_epic) %in% refSub,]
methyl_w6_epic <- methyl_w6_epic[!rownames(methyl_w6_epic) %in% refSub,]
methyl_w7_epic_e <- methyl_w7_epic[!rownames(methyl_w7_epic) %in% refSub,]


# read in the CIMT, candidate CpG list, and the covariate data
cimt_le8 <- readRDS("/drives/drive1/ff_chaya/Ayanava/cimt_le8.RDS")
colnames(cimt_le8)[colnames(cimt_le8) == "age_at_clinic"] <- "Age_k7"
colnames(cimt_le8)[colnames(cimt_le8) == "cvid"] <- "CVID"


cpgsites <- read.csv("/drives/drive1/ff_chaya/Ayanava/cpgsites_ayanava_revised_new.csv")
cimtdataset <- readRDS("/drives/drive1/ff_chaya/R_versions/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2_datafile.RDS")
biomarker <- readRDS("/drives/drive1/ff_chaya/R_versions/FF_biomarker_CVID_5.24/FF_biomarker_CVID_5.24_datafile.RDS")
ancestry <- read.csv("/drives/drive1/ff_chaya/Ada/PsychChip2b.ancestryPCs.CVID.csv")
# y22_age <- read.csv("/drives/drive1/ff_chaya/Ayanava/CVID_y22_age.csv")

cpgsites <- cpgsites %>% mutate("slope", "pval", "rsq")
colnames(cpgsites) = c("CpGs", "slope", "pval", "rsq")
cpgsites <- cpgsites %>% distinct(CpGs, .keep_all = TRUE)
threshold_cpgs <- read.csv("/drives/drive1/ff_chaya/Ayanava/savedcpgs_post_saliva_blood.csv")
cpgsites <- cpgsites %>% filter(cpgsites$CpGs %in% threshold_cpgs$CpGs)


# checking where these CpGs are located
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- as.data.frame(anno)
anno <- anno %>% select(Name, chr, pos, UCSC_RefGene_Name, Relation_to_Island)
anno_cpgs <- cpgsites %>% left_join(anno, by = c("CpGs" ="Name"))
table(anno_cpgs$Relation_to_Island)


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

# estimated cell proportions added to dataset
cp_blood <- estimateLC(methyl_w7_epic_e, ref = "Salas")
cp_blood$CVID <- colnames(methyl_w7_epic_e)
cp_saliva_w6 <- estimateLC(methyl_w6_epic, ref = "salivaEPIC")
cp_saliva_w6$CVID <- colnames(methyl_w6_epic)
colnames(cp_saliva_w6) <- c("w6_leu", "w6_epi", "CVID")
cp_saliva_w5 <- estimateLC(methyl_w5_epic, ref = "salivaEPIC")
cp_saliva_w5$CVID <- colnames(methyl_w5_epic)
colnames(cp_saliva_w5) <- c("w5_leu", "w5_epi", "CVID")


cimt_AND_survey <- merge(cimt_AND_survey, cp_blood, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, cp_saliva_w5, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, cp_saliva_w6, by = "CVID", all = TRUE)


# create covariate dataframe and candidate CpG methylation dataframe for each year
cpgs_w5_epic <- clean_methyl(methyl_w5_epic, cpgsites, cimt_AND_survey)
cpgs_w6_epic <- clean_methyl(methyl_w6_epic, cpgsites, cimt_AND_survey)
cpgs_w7_epic <- clean_methyl(methyl_w7_epic_e, cpgsites, cimt_AND_survey)

cimt_w5_epic <- data.frame(cpgs_w5_epic[2])
cpgs_w5_epic <- data.frame(cpgs_w5_epic[1])

cimt_w6_epic <- data.frame(cpgs_w6_epic[2])
cpgs_w6_epic <- data.frame(cpgs_w6_epic[1])

cimt_w7_epic <- data.frame(cpgs_w7_epic[2])
cpgs_w7_epic <- data.frame(cpgs_w7_epic[1])

# create the matrices that will be used in the EWAS for each year
allcovs_5 <- c("Mean_Mean_CIMT", "log_Mean_Mean_CIMT", "ck7ethrace", "cm1bsex", "k5me_age", 
               "ch5cbmi", "ch5chtcm", "w5_leu", "w5_epi", 
               "PC1", "PC2", "PC3", "PC4")
mod5_epic <- cov_matrix(allcovs_5, cimt_w5_epic, cpgs_w5_epic)
mod5_epic$ck7ethrace <- as.factor(mod5_epic$ck7ethrace)
mod5_epic$cm1bsex <- droplevels(mod5_epic$cm1bsex)

allcovs_6 <- c("Mean_Mean_CIMT", "log_Mean_Mean_CIMT", "ck7ethrace", "cm1bsex", "k6me_age", 
               "ck6cbmi", "ck6chtcm", "w6_epi", "w6_leu", "PC1", 
               "PC2", "PC3", "PC4")
mod6_epic <- cov_matrix(allcovs_6, cimt_w6_epic, cpgs_w6_epic)
mod6_epic$ck7ethrace <- as.factor(mod6_epic$ck7ethrace)
mod6_epic$cm1bsex <- droplevels(mod6_epic$cm1bsex)

allcovs_7 <- c("Mean_Mean_CIMT", "log_Mean_Mean_CIMT", "ck7ethrace", "cm1bsex", "Age_k7", 
               "ck7bmi", "k7i41", "GR", "MO", "NK", "B", "CD4", "CD8", "PC1", "PC2", 
               "PC3", "PC4")
mod7_epic <- cov_matrix(allcovs_7, cimt_w7_epic, cpgs_w7_epic)
mod7_epic$ck7ethrace <- as.factor(mod7_epic$ck7ethrace)
mod7_epic$cm1bsex <- droplevels(mod7_epic$cm1bsex)

# remove all missing values
mod5_epic_nona <- na.omit(mod5_epic)
mod6_epic_nona <- na.omit(mod6_epic)
mod7_epic_nona <- na.omit(mod7_epic)
cpgs_w5_epic_nona <- cpgs_w5_epic %>% filter(rownames(cpgs_w5_epic) %in% rownames(mod5_epic_nona))
cpgs_w6_epic_nona <- cpgs_w6_epic %>% filter(rownames(cpgs_w6_epic) %in% rownames(mod6_epic_nona))
cpgs_w7_epic_nona <- cpgs_w7_epic %>% filter(rownames(cpgs_w7_epic) %in% rownames(mod7_epic_nona))

# Model 1 on Candidate Set of CpGs for Aim 1
design_arscp <- c("ck7ethrace", "k5me_age", "cm1bsex", "w5_epi", "w5_leu") 
design_arscp6 <- c("ck7ethrace", "k6me_age", "cm1bsex", "w6_epi", "w6_leu")
design_arscp7 <- c("ck7ethrace", "Age_k7", "cm1bsex", "GR", "CD4")

test5e_arscp <- run_ewas(mod5_epic_nona, design_arscp, cpgs_w5_epic_nona, log = FALSE, batch = batch_saliva_w5)

test6e_arscp <- run_ewas(mod6_epic_nona, design_arscp6, cpgs_w6_epic_nona, log = FALSE, batch = batch_saliva_w6)

test7e_arscp <- run_ewas(mod7_epic_nona, design_arscp7, cpgs_w7_epic_nona, log = FALSE, batch = batch_blood) # 10 hits

# Model 2 on Candidate Set of CpGs for Aim 1
design_aancscp <- c("k5me_age", "cm1bsex", "w5_epi", "w5_leu", "PC1", "PC2", "PC3", "PC4") 
design_aancscp6 <- c("k6me_age", "cm1bsex", "w6_epi", "w6_leu", "PC1", "PC2", "PC3", "PC4")
design_aancscp7 <- c("Age_k7", "cm1bsex", "GR", "CD4", "CD8", "PC1", "PC2", "PC3", "PC4")

test5e_aancscp <- run_ewas(mod5_epic_nona, design_aancscp, cpgs_w5_epic_nona, log = FALSE, batch = batch_saliva_w5)

test6e_aancscp <- run_ewas(mod6_epic_nona, design_aancscp6, cpgs_w6_epic_nona, log = FALSE, batch = batch_saliva_w6)

test7e_aancscp <- run_ewas(mod7_epic_nona, design_aancscp7, cpgs_w7_epic_nona, log = FALSE, batch = batch_blood) # 10 hits


title <- "Y22 Model 1"
p5 <- make_volcano(test7e_arscp, title)
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y22_r_volcano.png", p5, dpi = 300, width = 6.5, height = 6)

title <- "Y15 Model 1"
p9 <- make_volcano(test6e_arscp, title)
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y15_r_volcano.png", p9, dpi = 300, width = 6.5, height = 6)

title <- "Y9 Model 1"
p10 <- make_volcano(test5e_arscp, title)
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y9_r_volcano.png", p10, dpi = 300, width = 6.5, height = 6)



title <- "Y22 Model 2"
p6 <- make_volcano(test7e_aancscp, title)
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y22_anc_volcano.png", p6, dpi = 300, width = 6.5, height = 6)

title <- "Y15 Model 2"
p7 <- make_volcano(test6e_aancscp, title)
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y15_anc_volcano.png", p7, dpi = 300, width = 6.5, height = 6)

title <- "Y9 Model 2"
p8 <- make_volcano(test5e_aancscp, title)
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y9_anc_volcano.png", p8, dpi = 300, width = 6.5, height = 6)



# seeing what the trend CpGs are for the two models
trends_arscp <- get_trendcpgs(test5e_arscp, test6e_arscp, test7e_arscp) 
trends_aancscp <- get_trendcpgs(test5e_aancscp, test6e_aancscp, test7e_aancscp) 


# ancestry trends graph 
graph_trend_cpgs_anc <- trends_aancscp %>% select(logFC_5, logFC_6, logFC_7)
graph_trend_cpgs_anc <- as.data.frame(t(graph_trend_cpgs_anc))
graph_trend_cpgs_anc$Age <- c(9, 15, 22)

p11 <- ggplot(data = graph_trend_cpgs_anc, aes(x = Age)) + 
  geom_line(aes(y = graph_trend_cpgs_anc[, 1], color = colnames(graph_trend_cpgs_anc[1]))) + 
  geom_line(aes(y = graph_trend_cpgs_anc[, 2], color = colnames(graph_trend_cpgs_anc[2]))) + 
  geom_line(aes(y = graph_trend_cpgs_anc[, 3], color = colnames(graph_trend_cpgs_anc[3]))) + 
  geom_line(aes(y = graph_trend_cpgs_anc[, 4], color = colnames(graph_trend_cpgs_anc[4]))) +
  geom_line(aes(y = graph_trend_cpgs_anc[, 5], color = colnames(graph_trend_cpgs_anc[5]))) +
  geom_line(aes(y = graph_trend_cpgs_anc[, 6], color = colnames(graph_trend_cpgs_anc[6]))) +
  ylab("logFC") + 
  scale_x_continuous(breaks = c(9, 15, 22)) + 
  geom_abline(slope = 0, intercept = 0) + 
  scale_color_manual("", breaks = colnames(graph_trend_cpgs_anc), values = c("red", "cyan", "green", "blue", "orange", "purple")) + 
  labs(title = "Trend CpGs of Model 2") + 
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),  
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  scale_y_continuous(breaks = seq(-1.5, 1, by = 0.5), limits = c(-1.7, 0.6))

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/trends_anc.png", p11, dpi = 300, width =6.5, height = 6)

# the race trends graph
graph_trend_cpgs_r <- trends_arscp %>% select(logFC_5, logFC_6, logFC_7)
graph_trend_cpgs_r <- as.data.frame(t(graph_trend_cpgs_r))
graph_trend_cpgs_r$Age <- c(9, 15, 22)

p12 <- ggplot(data = graph_trend_cpgs_r, aes(x = Age)) + 
  geom_line(aes(y = graph_trend_cpgs_r[, 1], color = colnames(graph_trend_cpgs_r[1]))) + 
  geom_line(aes(y = graph_trend_cpgs_r[, 2], color = colnames(graph_trend_cpgs_r[2]))) + 
  geom_line(aes(y = graph_trend_cpgs_r[, 3], color = colnames(graph_trend_cpgs_r[3]))) + 
  geom_line(aes(y = graph_trend_cpgs_r[, 4], color = colnames(graph_trend_cpgs_r[4]))) +
  geom_line(aes(y = graph_trend_cpgs_r[, 5], color = colnames(graph_trend_cpgs_r[5]))) +
  ylab("logFC") + 
  scale_x_continuous(breaks = c(9, 15, 22)) + 
  geom_abline(slope = 0, intercept = 0) + 
  scale_color_manual("", breaks = colnames(graph_trend_cpgs_r), values = c("red", "green", "blue", "orange", "purple")) + 
  labs(title = "Trend CpGs of Model 1") + 
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) + 
  scale_y_continuous(breaks = seq(-1.5, 1, by = 0.5), limits = c(-1.7, 0.6))

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/trends_r.png", p12, dpi = 300, width =6.5, height = 6)


# # EWAS runs of Aim 1 

# Setting up the methylation dataframes and covariate dataframes for EWAS
methyl_w5_ewas <- methyl_w5_epic %>% select(cimt_w5_epic$CVID)
methyl_w5_ewas <- data.frame(t(methyl_w5_ewas))
ewas5_epic <- cov_matrix(allcovs_5, cimt_w5_epic, methyl_w5_ewas)
ewas5_epic$ck7ethrace <- as.factor(ewas5_epic$ck7ethrace)
ewas5_epic$cm1bsex <- droplevels(ewas5_epic$cm1bsex)
# 
methyl_w6_ewas <- methyl_w6_epic %>% select(cimt_w6_epic$CVID)
methyl_w6_ewas <- data.frame(t(methyl_w6_ewas))
ewas6_epic <- cov_matrix(allcovs_6, cimt_w6_epic, methyl_w6_ewas)
ewas6_epic$ck7ethrace <- as.factor(ewas6_epic$ck7ethrace)
ewas6_epic$cm1bsex <- droplevels(ewas6_epic$cm1bsex)

methyl_w7_ewas <- methyl_w7_epic_e %>% select(cimt_w7_epic$CVID)
methyl_w7_ewas <- data.frame(t(methyl_w7_ewas))
ewas7_epic <- cov_matrix(allcovs_7, cimt_w7_epic, methyl_w7_ewas)
ewas7_epic$ck7ethrace <- as.factor(ewas7_epic$ck7ethrace)
ewas7_epic$cm1bsex <- droplevels(ewas7_epic$cm1bsex)

# remove NAs
ewas5_epic_nona <- ewas5_epic %>% filter(rownames(ewas5_epic) %in% rownames(mod5_epic_nona))
ewas6_epic_nona <- ewas6_epic %>% filter(rownames(ewas6_epic) %in% rownames(mod6_epic_nona))
ewas7_epic_nona <- ewas7_epic %>% filter(rownames(ewas7_epic) %in% rownames(mod7_epic_nona))
methyl_w5_ewas_nona <- methyl_w5_ewas %>% filter(rownames(methyl_w5_ewas) %in% rownames(mod5_epic_nona))
methyl_w6_ewas_nona <- methyl_w6_ewas %>% filter(rownames(methyl_w6_ewas) %in% rownames(mod6_epic_nona))
methyl_w7_ewas_nona <- methyl_w7_ewas %>% filter(rownames(methyl_w7_ewas) %in% rownames(mod7_epic_nona))


# EWAS with Model 1
test5ewas_arscp <- run_ewas(ewas5_epic_nona, design_arscp, methyl_w5_ewas_nona, log = FALSE, batch = batch_saliva_w5)
test6ewas_arscp <- run_ewas(ewas6_epic_nona, design_arscp6, methyl_w6_ewas_nona, log = FALSE, batch = batch_saliva_w6)
test7ewas_arscp <- run_ewas(ewas7_epic_nona, design_arscp7, methyl_w7_ewas_nona, log = FALSE, batch = batch_blood)

# EWAS with Model 2
test5ewas_aancscp <- run_ewas(ewas5_epic_nona, design_aancscp, methyl_w5_ewas_nona, log = FALSE, batch = batch_saliva_w5)
test6ewas_aancscp <- run_ewas(ewas6_epic_nona, design_aancscp6, methyl_w6_ewas_nona, log = FALSE, batch = batch_saliva_w6)
test7ewas_aancscp <- run_ewas(ewas7_epic_nona, design_aancscp7, methyl_w7_ewas_nona, log = FALSE, batch = batch_blood)




saveRDS(test5ewas_aancscp, "/drives/drive1/ff_chaya/Ayanava/ewas_results_w5_aancscp.R")
saveRDS(test6ewas_aancscp, "/drives/drive1/ff_chaya/Ayanava/ewas_results_w6_aancscp.R")
saveRDS(test7ewas_aancscp, "/drives/drive1/ff_chaya/Ayanava/ewas_results_w7_aancscp.R")
saveRDS(test5ewas_arscp, "/drives/drive1/ff_chaya/Ayanava/ewas_results_w5_arscp.R")
saveRDS(test6ewas_arscp, "/drives/drive1/ff_chaya/Ayanava/ewas_results_w6_arscp.R")
saveRDS(test7ewas_arscp, "/drives/drive1/ff_chaya/Ayanava/ewas_results_w7_arscp.R")

# checking for potential qualitative DMRs from the 6 trend CpGs from Model 2
main_cpgs_anc <- as.data.frame(colnames(graph_trend_cpgs_anc))
main_cpgs_anc <- as.data.frame(main_cpgs_anc[1:6,])
colnames(main_cpgs_anc) <- "cpgs"
anno_cpgs_anc <- reference %>% filter(rownames(reference) %in% main_cpgs_anc$cpgs)

write.csv(anno_cpgs_anc, "/drives/drive1/ff_chaya/Ayanava/trend_anc_w_anno.csv")
write.csv(trends_aancscp, "/drives/drive1/ff_chaya/Ayanava/trend_anc_w_pvalratio_logfc.csv")

# time to identify the CpGs in and around these 6 CpGs - let's establish a +/- threshold of 10,000 bp 

cpg1 <- reference %>% filter(reference$seqnames %in% anno_cpgs_anc$seqnames[1])
cpg1 <- cpg1 %>% filter((cpg1$start > (anno_cpgs_anc$start[1] - 10000)) & (cpg1$end < (anno_cpgs_anc$end[1] + 10000)))

cpg2 <- reference %>% filter(reference$seqnames %in% anno_cpgs_anc$seqnames[2])
cpg2 <- cpg2 %>% filter((cpg2$start > (anno_cpgs_anc$start[2] - 10000)) & (cpg2$end < (anno_cpgs_anc$end[2] + 10000)))

cpg3 <- reference %>% filter(reference$seqnames %in% anno_cpgs_anc$seqnames[3])
cpg3 <- cpg3 %>% filter((cpg3$start > (anno_cpgs_anc$start[3] - 10000)) & (cpg3$end < (anno_cpgs_anc$end[3] + 10000)))

cpg4 <- reference %>% filter(reference$seqnames %in% anno_cpgs_anc$seqnames[4])
cpg4 <- cpg4 %>% filter((cpg4$start > (anno_cpgs_anc$start[4] - 10000)) & (cpg4$end < (anno_cpgs_anc$end[4] + 10000)))

cpg5 <- reference %>% filter(reference$seqnames %in% anno_cpgs_anc$seqnames[5])
cpg5 <- cpg5 %>% filter((cpg5$start > (anno_cpgs_anc$start[5] - 10000)) & (cpg5$end < (anno_cpgs_anc$end[5] + 10000)))

cpg6 <- reference %>% filter(reference$seqnames %in% anno_cpgs_anc$seqnames[6])
cpg6 <- cpg6 %>% filter((cpg6$start > (anno_cpgs_anc$start[6] - 10000)) & (cpg6$end < (anno_cpgs_anc$end[6] + 10000)))

# function to run the EWAS model on a subset of CpGs, like defined for the above 6
cpg_modelrun <- function(cpg, methyl, cimt, covs, design, log, wave, batch) {
  cpg_sites <- as.data.frame(rownames(cpg))
  colnames(cpg_sites) <- "CpGs"
  cpg_sites <- cpg_sites %>% mutate("slope", "pval", "rsq")
  colnames(cpg_sites) = c("CpGs", "slope", "pval", "rsq")
  cpg_w <- clean_methyl(methyl, cpg_sites, cimt)
  cimt_w <- data.frame(cpg_w[2])
  cpg_w <- data.frame(cpg_w[1])
  mod_cpg <- cov_matrix(covs, cimt_w, cpg_w)
  mod_cpg$ck7ethrace <- as.factor(mod_cpg$ck7ethrace)
  mod_cpg$cm1bsex <- droplevels(mod_cpg$cm1bsex)
  mod_cpg_nona <- na.omit(mod_cpg)
  cpg_w_nona <- cpg_w %>% filter(rownames(cpg_w) %in% rownames(mod_cpg_nona))
  trial_cpg <- run_ewas(mod_cpg_nona, design, cpg_w_nona, log = log, batch = batch)
  return (trial_cpg)
}

# this is for the EWAS needed for Aim 3 (with the meQTL and MR results)
methyl_meqtl_w7 <- readRDS("/drives/drive1/ff_chaya/Ayanava/methyl_meqtl_w7.rds")
test_meqtl <- cpg_modelrun(cpg = methyl_meqtl_w7, methyl = methyl_meqtl_w7, cimt = cimt_AND_survey, 
                     covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)
saveRDS(test_meqtl, "/drives/drive1/ff_chaya/Ayanava/ewas_meqtl_results.rds")


cpg1_results <- cpg_modelrun(cpg1, methyl = methyl_w7_epic_e, cimt = cimt_AND_survey, 
                             covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)

cpg2_results <- cpg_modelrun(cpg2, methyl = methyl_w7_epic_e, cimt = cimt_AND_survey, 
                             covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)

cpg3_results <- cpg_modelrun(cpg3, methyl = methyl_w7_epic_e, cimt = cimt_AND_survey, 
                             covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)

cpg4_results <- cpg_modelrun(cpg4, methyl = methyl_w7_epic_e, cimt = cimt_AND_survey, 
                             covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)

cpg5_results <- cpg_modelrun(cpg5, methyl = methyl_w7_epic_e, cimt = cimt_AND_survey, 
                             covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)

cpg6_results <- cpg_modelrun(cpg6, methyl = methyl_w7_epic_e, cimt = cimt_AND_survey, 
                             covs = allcovs_7, design = design_aancscp7, log = FALSE, wave = 7, batch = batch_blood)


# graphing the p-values against the CpGs' relative position

# the chr11 triplet
cpg2_filter <- cpg2 %>% filter(rownames(cpg2) %in% rownames(cpg2_results))
cpg2_filter <- cpg2_filter %>% select(c("seqnames", "start", "end"))
cpg2_filter$P.Value <- cpg2_results$P.Value[match(rownames(cpg2_filter), rownames(cpg2_results))]
cpg2_filter$logFC <- cpg2_results$logFC[match(rownames(cpg2_filter), rownames(cpg2_results))]
cpg2_filter$Target <- FALSE
cpg2_filter$Target[rownames(cpg2_filter) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg2_filter, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg2_filter, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg2_filter, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), name = "-log10(p-value)", limits = c(-1.75, 6.25), 
                     sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y22")


ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y22_chr11.png", dpi = 600, width = 7, height = 6)

# plotting the chr9 CpG region's logFC and p-value 
cpg1_filter <- cpg1 %>% filter(rownames(cpg1) %in% rownames(cpg1_results))
cpg1_filter <- cpg1_filter %>% select(c("seqnames", "start", "end"))
cpg1_filter$P.Value <- cpg1_results$P.Value[match(rownames(cpg1_filter), rownames(cpg1_results))]
cpg1_filter$logFC <- cpg1_results$logFC[match(rownames(cpg1_filter), rownames(cpg1_results))]
cpg1_filter$Target <- FALSE
cpg1_filter$Target[rownames(cpg1_filter) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg1_filter, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg1_filter, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg1_filter, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1.25, 3.5), name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y22")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y22_chr9.png", dpi = 600, width = 7, height = 6)

# chr17 plot
cpg5_filter <- cpg5 %>% filter(rownames(cpg5) %in% rownames(cpg5_results))
cpg5_filter <- cpg5_filter %>% select(c("seqnames", "start", "end"))
cpg5_filter$P.Value <- cpg5_results$P.Value[match(rownames(cpg5_filter), rownames(cpg5_results))]
cpg5_filter$logFC <- cpg5_results$logFC[match(rownames(cpg5_filter), rownames(cpg5_results))]
cpg5_filter$Target <- FALSE
cpg5_filter$Target[rownames(cpg5_filter) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg5_filter, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg5_filter, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg5_filter, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1),  limits = c(-1.5, 3.3),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y22")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y22_chr17.png", dpi = 600, width = 7, height = 6)




# checking with wave 5 and 6 for chr11
cpg4_results_w5 <- cpg_modelrun(cpg4, methyl = methyl_w5_epic, cimt = cimt_AND_survey, 
                             covs = allcovs_5, design = design_aancscp, log = FALSE, wave = 5, batch = batch_saliva_w5)

cpg4_results_w6 <- cpg_modelrun(cpg4, methyl = methyl_w6_epic, cimt = cimt_AND_survey, 
                                covs = allcovs_6, design = design_aancscp6, log = FALSE, wave = 6, batch = batch_saliva_w6)


cpg4_filter_w5 <- cpg4 %>% filter(rownames(cpg4) %in% rownames(cpg4_results_w5))
cpg4_filter_w5 <- cpg4_filter_w5 %>% select(c("seqnames", "start", "end"))
cpg4_filter_w5$P.Value <- cpg4_results_w5$P.Value[match(rownames(cpg4_filter_w5), rownames(cpg4_results_w5))]
cpg4_filter_w5$logFC <- cpg4_results_w5$logFC[match(rownames(cpg4_filter_w5), rownames(cpg4_results_w5))]
cpg4_filter_w5$Target <- FALSE
cpg4_filter_w5$Target[rownames(cpg4_filter_w5) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg4_filter_w5, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg4_filter_w5, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg4_filter_w5, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1),  limits = c(-1.75, 6.25),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y9")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y9_chr11.png", dpi = 600, width = 7, height = 6)



cpg4_filter_w6 <- cpg4 %>% filter(rownames(cpg4) %in% rownames(cpg4_results_w6))
cpg4_filter_w6 <- cpg4_filter_w6 %>% select(c("seqnames", "start", "end"))
cpg4_filter_w6$P.Value <- cpg4_results_w6$P.Value[match(rownames(cpg4_filter_w6), rownames(cpg4_results_w6))]
cpg4_filter_w6$logFC <- cpg4_results_w6$logFC[match(rownames(cpg4_filter_w6), rownames(cpg4_results_w6))]
cpg4_filter_w6$Target <- FALSE
cpg4_filter_w6$Target[rownames(cpg4_filter_w6) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg4_filter_w6, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg4_filter_w6, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg4_filter_w6, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1.75, 6.25), name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y15")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y15_chr11.png", dpi = 600, width = 7, height = 6)


# doing it for the chr17 
cpg5_results_w5 <- cpg_modelrun(cpg5, methyl = methyl_w5_epic, cimt = cimt_AND_survey, 
                                covs = allcovs_5, design = design_aancscp, log = FALSE, wave = 5, batch = batch_saliva_w5)

cpg5_results_w6 <- cpg_modelrun(cpg5, methyl = methyl_w6_epic, cimt = cimt_AND_survey, 
                                covs = allcovs_6, design = design_aancscp6, log = FALSE, wave = 6, batch = batch_saliva_w6)


cpg5_filter_w5 <- cpg5 %>% filter(rownames(cpg5) %in% rownames(cpg5_results_w5))
cpg5_filter_w5 <- cpg5_filter_w5 %>% select(c("seqnames", "start", "end"))
cpg5_filter_w5$P.Value <- cpg5_results_w5$P.Value[match(rownames(cpg5_filter_w5), rownames(cpg5_results_w5))]
cpg5_filter_w5$logFC <- cpg5_results_w5$logFC[match(rownames(cpg5_filter_w5), rownames(cpg5_results_w5))]
cpg5_filter_w5$Target <- FALSE
cpg5_filter_w5$Target[rownames(cpg5_filter_w5) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg5_filter_w5, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg5_filter_w5, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg5_filter_w5, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1.5, 3.3),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y9")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y9_chr17.png", dpi = 600, width = 7, height = 6)


cpg5_filter_w6 <- cpg5 %>% filter(rownames(cpg5) %in% rownames(cpg5_results_w6))
cpg5_filter_w6 <- cpg5_filter_w6 %>% select(c("seqnames", "start", "end"))
cpg5_filter_w6$P.Value <- cpg5_results_w6$P.Value[match(rownames(cpg5_filter_w6), rownames(cpg5_results_w6))]
cpg5_filter_w6$logFC <- cpg5_results_w6$logFC[match(rownames(cpg5_filter_w6), rownames(cpg5_results_w6))]
cpg5_filter_w6$Target <- FALSE
cpg5_filter_w6$Target[rownames(cpg5_filter_w6) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg5_filter_w6, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg5_filter_w6, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg5_filter_w6, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1.5, 3.3),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y15")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y15_chr17.png", dpi = 600, width = 7, height = 6)



# doing it for chr9
cpg1_results_w5 <- cpg_modelrun(cpg1, methyl = methyl_w5_epic, cimt = cimt_AND_survey, 
                                covs = allcovs_5, design = design_aancscp, log = FALSE, wave = 5, batch = batch_saliva_w5)

cpg1_results_w6 <- cpg_modelrun(cpg1, methyl = methyl_w6_epic, cimt = cimt_AND_survey, 
                                covs = allcovs_6, design = design_aancscp6, log = FALSE, wave = 6, batch = batch_saliva_w6)


cpg1_filter_w5 <- cpg1 %>% filter(rownames(cpg1) %in% rownames(cpg1_results_w5))
cpg1_filter_w5 <- cpg1_filter_w5 %>% select(c("seqnames", "start", "end"))
cpg1_filter_w5$P.Value <- cpg1_results_w5$P.Value[match(rownames(cpg1_filter_w5), rownames(cpg1_results_w5))]
cpg1_filter_w5$logFC <- cpg1_results_w5$logFC[match(rownames(cpg1_filter_w5), rownames(cpg1_results_w5))]
cpg1_filter_w5$Target <- FALSE
cpg1_filter_w5$Target[rownames(cpg1_filter_w5) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg1_filter_w5, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg1_filter_w5, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg1_filter_w5, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1.25, 3.5),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y9")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y9_chr9.png", dpi = 600, width = 7, height = 6)

cpg1_filter_w6 <- cpg1 %>% filter(rownames(cpg1) %in% rownames(cpg1_results_w6))
cpg1_filter_w6 <- cpg1_filter_w6 %>% select(c("seqnames", "start", "end"))
cpg1_filter_w6$P.Value <- cpg1_results_w6$P.Value[match(rownames(cpg1_filter_w6), rownames(cpg1_results_w6))]
cpg1_filter_w6$logFC <- cpg1_results_w6$logFC[match(rownames(cpg1_filter_w6), rownames(cpg1_results_w6))]
cpg1_filter_w6$Target <- FALSE
cpg1_filter_w6$Target[rownames(cpg1_filter_w6) %in% colnames(graph_trend_cpgs_anc)] <- TRUE
ggplot(cpg1_filter_w6, aes(x = start)) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = -log10(P.Value)), data = subset(cpg1_filter_w6, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = logFC), data = subset(cpg1_filter_w6, Target == TRUE), 
             color = "red", size = 5, shape = 18, alpha = 1) + 
  geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
  geom_point(aes(y = logFC, color = "logFC")) + 
  geom_line(aes(y = logFC, color = "logFC")) + 
  scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1.25, 3.5),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
  labs(x = "CpG Position", color = "Metric") + 
  geom_hline(yintercept = 0, color = "gray50") +
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        legend.position = "right", 
        legend.justification = "bottom",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_color_manual(values = c("-log10(p-value)" = "#0072B2", "logFC" = "#009E73")) + 
  labs(title = "Y15")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/y15_chr9.png", dpi = 600, width = 7, height = 6)



# Manhattan plot of the EWAS from Model 1 and Model 2 
library(qqman)

# function that generates Manhattan plot given EWAS results, reference manifest, and title of plot
make_manhattan <- function(ewas, reference, title) {
  ewas$CpG <- rownames(ewas)
  reference$CpG <- rownames(reference)
  ewas7_man <- merge(ewas, reference[, c("CpG", "seqnames", "start")], , by = "CpG")
  colnames(ewas7_man)[colnames(ewas7_man) == "seqnames"] <- "CHR"
  ewas7_man$CHR <- gsub("chr", "", ewas7_man$CHR)
  ewas7_man$CHR <- as.numeric(ewas7_man$CHR)
  colnames(ewas7_man)[colnames(ewas7_man) == "start"] <- "BP"
  colnames(ewas7_man)[colnames(ewas7_man) == "CpG"] <- "SNP"
  colnames(ewas7_man)[colnames(ewas7_man) == "P.Value"] <- "P"
  manhattan(ewas7_man, genomewideline = -log10(5e-8), ylim = c(0, 8), main = title)
}


test5ewas_aancscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w5_aancscp.R")
test6ewas_aancscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w6_aancscp.R")
test7ewas_aancscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w7_aancscp.R")

test5ewas_arscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w5_arscp.R")
test6ewas_arscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w6_arscp.R")
test7ewas_arscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w7_arscp.R")

title <- "Y22"
png("/drives/drive1/ff_chaya/Ayanava/Plots/ewas_manhattan_y22_anc.png", width = 6.5, height = 6, units = "in", res = 300)
make_manhattan(test7ewas_aancscp, reference, title)
dev.off()

title <- "Y15"
png("/drives/drive1/ff_chaya/Ayanava/Plots/ewas_manhattan_y15_anc.png", width = 6.5, height = 6, units = "in", res = 300)
make_manhattan(test6ewas_aancscp, reference, title)
dev.off()

title <- "Y9"
png("/drives/drive1/ff_chaya/Ayanava/Plots/ewas_manhattan_y9_anc.png", width = 6.5, height = 6, units = "in", res = 300)
make_manhattan(test5ewas_aancscp, reference, title)
dev.off()


title <- "Y22"
png("/drives/drive1/ff_chaya/Ayanava/Plots/ewas_manhattan_y22_r.png", width = 6.5, height = 6, units = "in", res = 300)
make_manhattan(test7ewas_arscp, reference, title)
dev.off()

title <- "Y15"
png("/drives/drive1/ff_chaya/Ayanava/Plots/ewas_manhattan_y15_r.png", width = 6.5, height = 6, units = "in", res = 300)
make_manhattan(test6ewas_arscp, reference, title)
dev.off()

title <- "Y9"
png("/drives/drive1/ff_chaya/Ayanava/Plots/ewas_manhattan_y9_r.png", width = 6.5, height = 6, units = "in", res = 300)
make_manhattan(test5ewas_arscp, reference, title)
dev.off()



