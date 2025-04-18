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

# read in Zhou annotation manifest
reference <- as.data.frame(readRDS("/drives/drive1/geneticData/referenceDataSets/Methylation/Zhou/EPIC.hg19.manifest.rds"))

# read in saliva methylation data
new_methyl <- readRDS("/drives/drive1/methylation_processed_data/BetaMatrix.saliva.n2627.Sep16.RDS")
new_methyl <- data.frame(new_methyl, check.names = FALSE)
ids_convert <- read.csv("/drives/drive1/ff_chaya/Ayanava/Combined_EPIC_with_y22.arrayID.CVID.csv")
ids_convert2 <- read.table("/drives/drive1/ff_chaya/SampleSheets/epic_ss_cvid_nodups.txt")
ids_convert3 <- read.csv("/drives/drive1/ff_chaya/Ayanava/Combined_EPIC.SS.correctedSex.csv")


colnames(ids_convert2) <- c("MethID", "CVID", "wave", "Sample_Plate", 
                            "cm1bsex", "cm1bethrace")
ids_convert2 <- ids_convert2[-1, ]


for (i in 1:length(new_methyl)) {
  match <- ids_convert3 %>% 
    filter(Collection == "OGR600") %>% 
    filter(Array_ID == colnames(new_methyl)[i])
  newname <- paste(match[1,5], match[1,4], sep = "")
  colnames(new_methyl)[i] <- newname
}

# read in blood methylation data
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


for (i in 1:length(methyl_y22)) {
  match <- ids_convert3 %>% 
    filter(Collection == "PaxGeneDNA") %>% 
    filter(Array_ID == colnames(methyl_y22)[i])
  newname <- match[1,5]
  colnames(methyl_y22)[i] <- newname
}

mval_meth <- log2((new_methyl) / (1 - new_methyl))


methyl_w7_saliva <- select(mval_meth, contains("FC22"))



for ( i in 1:ncol(methyl_w7_saliva)){
  colnames(methyl_w7_saliva)[i] <-  sub("FC22", "", colnames(methyl_w7_saliva)[i])
}

methyl_w7_blood <- log2((methyl_y22) / (1 - methyl_y22))


overlap <- intersect(colnames(methyl_w7_blood), colnames(methyl_w7_saliva))
methyl_w7_blood_common <- methyl_w7_blood[, overlap, drop = FALSE]
methyl_w7_saliva_common <- methyl_w7_saliva[, overlap, drop = FALSE]

overlap <- intersect(rownames(methyl_w7_blood_common), rownames(methyl_w7_saliva_common))
methyl_w7_blood_common <- methyl_w7_blood_common[overlap, , drop = FALSE]
methyl_w7_saliva_common <- methyl_w7_saliva_common[overlap, , drop = FALSE]

# remove all CpGs that are masked or on non-autosomal chromosomes
reference$ID <- rownames(reference)
refSub <- c(reference[reference$MASK_general==TRUE,"ID"],
            reference[reference$seqnames == "chrX","ID"],
            reference[reference$seqnames == "chrY","ID"],
            reference[reference$seqnames == "chrM","ID"])

methyl_w7_saliva_common <- methyl_w7_saliva_common[!rownames(methyl_w7_saliva_common) %in% refSub,]
methyl_w7_blood_common <- methyl_w7_blood_common[!rownames(methyl_w7_blood_common) %in% refSub,]

# PCA Analysis

library(stats)

methyl_w7_combined <- cbind(methyl_w7_blood_common, methyl_w7_saliva_common)
pca_result <- prcomp(t(methyl_w7_combined), scale. =TRUE)
pc_scores <- as.data.frame(pca_result$x)
pc_scores$Sample_Type <- c(rep("Blood", ncol(methyl_w7_blood_common)), rep("Saliva", ncol(methyl_w7_saliva_common)))
p4 <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = Sample_Type)) + 
  geom_point(size = 2) + 
  labs(title = "PCA of Saliva vs Blood Methylation Data", color = "Sample Type") + 
  theme(axis.title.x = element_text(size = "14", face = "bold"), 
        axis.title.y = element_text(size = "14", face = "bold"), 
        plot.title = element_text(size = "16", face = "bold", hjust = 0.5),
        legend.text = element_text(size = "12"),
        legend.title = element_text(size = "14"))

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/pc_plot.png", p4, dpi = 300, width = 7, height = 6)


# adjust for the effects of the first PC
pc_to_remove <- pca_result$x[,1]
pc_loadings <- pca_result$rotation[,1]
adjusted_methylation <- methyl_w7_combined - pc_to_remove %*% t(pc_loadings)

adjusted_methyl_w7_blood <- adjusted_methylation[, 1:ncol(methyl_w7_blood_common)]
adjusted_methyl_w7_saliva <- adjusted_methylation[, (ncol(methyl_w7_blood_common)+1):ncol(methyl_w7_combined)]

# non-adjusted methyl data 
w7_saliva_means <- data.frame(Saliva = rowMeans(methyl_w7_saliva_common, na.rm = TRUE))
w7_blood_means <- data.frame(Blood = rowMeans(methyl_w7_blood_common, na.rm = TRUE))
w7_means <- cbind(w7_blood_means, w7_saliva_means)



# using the adjusted methyl values 
adj_saliva_means <- data.frame(Saliva = rowMeans(adjusted_methyl_w7_saliva, na.rm = TRUE))
adj_blood_means <- data.frame(Blood = rowMeans(adjusted_methyl_w7_blood, na.rm = TRUE))
adj_means <- cbind(adj_blood_means, adj_saliva_means)



# highlighting where the candidate CpGs are in the scatterplot 

# visualizing with the M-values 

cpgsites <- read.csv("/drives/drive1/ff_chaya/Ayanava/cpgsites_ayanava_revised_new.csv")
cpgsites <- cpgsites %>% mutate("slope", "pval", "rsq")
colnames(cpgsites) = c("CpGs", "slope", "pval", "rsq")
cpgsites <- cpgsites %>% distinct(CpGs, .keep_all = TRUE)

adj_means$Target <- FALSE
cpgs_present  <- rownames(adj_means)
cpgsites_vect <- cpgsites$CpGs
candidates <- cpgsites_vect[cpgsites_vect %in% cpgs_present]
adj_means$Target[rownames(adj_means) %in% candidates] <- TRUE

p <- ggplot(data = adj_means, aes(x = Blood, y = Saliva)) + 
  geom_hex(bins = 100) + 
  geom_point(data = subset(adj_means, Target == TRUE), aes(x = Blood, y = Saliva), 
             color = "gold", size = 2, alpha = 0.8) + 
  geom_abline(intercept = 0, slope = 1, color = "white", linetype = "dashed", size = 1) + 
  geom_abline(intercept = -0.9, slope = 1, color = "grey", size = 1) + 
  geom_abline(intercept = 0.9, slope = 1, color = "grey", size = 1) +   
  theme(axis.title.x = element_text(size = "16", face = "bold"), 
  axis.title.y = element_text(size = "16", face = "bold"), 
  plot.title = element_text(size = "18", face = "bold", hjust = 0.5)) + 
  xlim(-10, 10) + 
  ylim(-10, 10) + 
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "log",
                      breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1000", "10000"), 
                      name = "CpG Count") + 
  labs(title = "Comparison of M-Values")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/mvalue_comp.png", p, dpi = 300, width = 6.5, height = 6)


# visualizing with the beta values 

adj_means_beta <- 2^adj_means / (1 + 2^adj_means)
adj_means_beta$Target <- FALSE
adj_means_beta$Target[rownames(adj_means_beta) %in% candidates] <- TRUE
adj_means_beta$Special <- FALSE
adj_means_beta$Special[rownames(adj_means_beta) %in% special] <- TRUE

# bounds set at +/- 15% difference in methylation, derived from previous literature
p2 <- ggplot(data = adj_means_beta, aes(x = Blood, y = Saliva)) + 
  geom_hex(bins = 100) + 
  geom_point(data = subset(adj_means_beta, Target == TRUE), aes(x = Blood, y = Saliva), 
             color = "gold", size = 2, alpha = 0.8) + 
  geom_abline(intercept = 0, slope = 1, color = "white", linetype = "dashed", size = 1) + 
  geom_abline(intercept = -0.15, slope = 1, color = "white", size = 1) + 
  geom_abline(intercept = 0.15, slope = 1, color = "white", size = 1) + 
  theme(axis.title.x = element_text(size = "16", face = "bold"), 
      axis.title.y = element_text(size = "16", face = "bold"), 
      plot.title = element_text(size = "18", face = "bold", hjust = 0.5)) + 
  xlim(0, 1) + 
  ylim(0, 1) + 
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "log",
                      breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1000", "10000"), 
                      name = "CpG Count") + 
  labs(title = "Comparison of Beta Values")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/betavalue_comp.png", p2, dpi = 300, width = 6.5, height = 6)


# isolating the target CpGs within the bounds identified from previous literature
threshold <- 0.15 
cpgs_within_threshold <- adj_means_beta %>% filter(Blood - threshold < Saliva & Saliva < Blood + threshold)
cpgs_within_threshold <- data.frame(CpGs = rownames(cpgs_within_threshold))
candidate_cpgs <- cpgs_within_threshold %>% filter(CpGs %in% candidates)

write.csv(candidate_cpgs, "/drives/drive1/ff_chaya/Ayanava/savedcpgs_post_saliva_blood.csv", row.names = FALSE)

# visualizing the saliva and blood data separately, as density plots
adj_means$CpGs <- rownames(adj_means)

ggplot(data = adj_means, aes(x = Blood)) + 
  geom_density(fill = "lightblue", alpha = 0.5) + 
  geom_point(data = (adj_means[adj_means$Special == TRUE, ]), aes(x = Blood, y = 0, color = CpGs), 
            size = 1.5) 

ggplot(data = adj_means, aes(x = Saliva)) + 
  geom_density(fill = "lightblue", alpha = 0.5) + 
  geom_point(data = (adj_means[adj_means$Special == TRUE, ]), aes(x = Saliva, y = 0, color = CpGs), 
             size = 1.5) 



