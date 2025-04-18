source("/drives/drive1/ff_chaya/Ayanava/limma_functions.R")
# source("/drives/drive1/ff_chaya/Ayanava/cpglimma_epic.R")

test5ewas_aancscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w5_aancscp.R")
test6ewas_aancscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w6_aancscp.R")
test7ewas_aancscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w7_aancscp.R")

test5ewas_arscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w5_arscp.R")
test6ewas_arscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w6_arscp.R")
test7ewas_arscp <- readRDS("/drives/drive1/ff_chaya/Ayanava/ewas_results_w7_arscp.R")

reference <- as.data.frame(readRDS("/drives/drive1/geneticData/referenceDataSets/Methylation/Zhou/EPIC.hg19.manifest.rds"))

# running the DMR comb-p on the race model 

# set up combp dataframes
combp_df_5_r <- get_combp_df(reference, test5ewas_arscp)
combp_df_6_r <- get_combp_df(reference, test6ewas_arscp)
combp_df_7_r <- get_combp_df(reference, test7ewas_arscp)

# run combp with Model 1 for each year (distance = 1000 bp)
combp(combp_df_5_r, dist.cutoff=1000,bin.size=310,seed=0.05,
      region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)
combp(combp_df_6_r, dist.cutoff=1000,bin.size=310,seed=0.05,
      region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)
combp(combp_df_7_r, dist.cutoff=1000,bin.size=310,seed=0.05,
      region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)


# running comb-p on the ancestry data 
combp_df_5_anc <- get_combp_df(reference, test5ewas_aancscp)
combp_df_6_anc <- get_combp_df(reference, test6ewas_aancscp)
combp_df_7_anc <- get_combp_df(reference, test7ewas_aancscp)

combp(combp_df_5_anc, dist.cutoff=1000,bin.size=310,seed=0.05,
      region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)
combp(combp_df_6_anc, dist.cutoff=1000,bin.size=310,seed=0.05,
      region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)
combp(combp_df_7_anc, dist.cutoff=1000,bin.size=310,seed=0.05,
      region_plot=TRUE,mht_plot=TRUE,nCores=10,verbose=TRUE)

# reading in and using the Model 1 results 
combp_results_w6 <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas6_new_r.csv")
combp_results_w5 <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas5_new_r.csv")
combp_results_w7 <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas7_new_r.csv")

combp_results_w5 <- combp_results_w5 %>% mutate("filter" = start %/% 10000)
combp_results_w6 <- combp_results_w6 %>% mutate("filter" = start %/% 10000)
combp_results_w7 <- combp_results_w7 %>% mutate("filter" = start %/% 10000)

# matches between wave 6 and 7 
dmr_match_67 <- combp_results_w6 %>% filter(combp_results_w6$filter %in% combp_results_w7$filter)
df <- combp_results_w7 %>% filter(combp_results_w7$filter %in% combp_results_w6$filter)
dmr_match_67 <- rbind(dmr_match_67, df)
dmr_match_67$wave <- c(6, 7)

# matches between wave 5 and 6 
dmr_match_56 <- combp_results_w5 %>% filter(combp_results_w5$filter %in% combp_results_w6$filter)
df <- combp_results_w6 %>% filter(combp_results_w6$filter %in% combp_results_w5$filter)
dmr_match_56 <- rbind(dmr_match_56, df)
dmr_match_56$wave <- c(5, 5, 6, 6)

# matches between wave 5 and 7 
dmr_match_57 <- combp_results_w5 %>% filter(combp_results_w5$filter %in% combp_results_w7$filter)
df <- combp_results_w7 %>% filter(combp_results_w7$filter %in% combp_results_w5$filter)
dmr_match_57 <- rbind(dmr_match_57, df)
dmr_match_57$wave <- c(5, 7)

# finds matches between all three years
dmr_match <- function(combp_w1, combp_w2, combp_w3) {
  combp_1 <- combp_w1 %>% mutate("filter" = start %/% 10000, "wave" = 5)
  combp_2 <- combp_w2 %>% mutate("filter" = start %/% 10000, "wave" = 6)
  combp_3 <- combp_w3 %>% mutate("filter" = start %/% 10000, "wave" = 7)
  dmr_match <- combp_1 %>% filter(combp_1$filter %in% combp_2$filter
                                  & combp_1$filter %in% combp_3$filter)
  df <- combp_2 %>% filter(combp_2$filter %in% combp_1$filter & combp_2$filter %in% combp_3$filter)
  df2 <- combp_3 %>% filter(combp_3$filter %in% combp_1$filter & combp_3$filter %in% combp_2$filter)
  dmr_match <- rbind(dmr_match, df)
  dmr_match <- rbind(dmr_match, df2)
  return (dmr_match)
}

combp_results_w6_anc <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas6_new_anc.csv")
combp_results_w5_anc <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas5_new_anc.csv")
combp_results_w7_anc <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas7_new_anc.csv")

# check overall matches for Model 2
dmr_match_anc <- dmr_match(combp_results_w5_anc, combp_results_w6_anc, combp_results_w7_anc)

combp_results_w6_r <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas6_new_r.csv")
combp_results_w5_r <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas5_new_r.csv")
combp_results_w7_r <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas7_new_r.csv")

# check overall matches for Model 1
dmr_match_r <- dmr_match(combp_results_w5_r, combp_results_w6_r, combp_results_w7_r)

saveRDS(dmr_match_r, "/drives/drive1/ff_chaya/Ayanava/dmr_match_r.rds")


# 
# make venn diagram of matches for Model 1
library(VennDiagram)
venn.plot <- draw.triple.venn(
  area1 = 33, 
  area2 = 40, 
  area3 = 73, 
  n12 = 10, 
  n23 = 4, 
  n13 = 5, 
  n123 = 3, 
  category = c("Y9", "Y15", "Y22"), 
  fill = c("skyblue", "pink", "lightgreen"), 
  cat.cex = 2, 
  cex = 2,
  cat.fontface = "bold"
)
png("/drives/drive1/ff_chaya/Ayanava/Plots/venn_diagram.png", width = 7, height = 7, units = "in", res = 300)
grid.draw(venn.plot)
dev.off()


# graphing the logFC values of the dmr_match_r


# analzying EWAS DMRs on logFC/pval graphs


# Wave 7 DMRs

dmr_match_r <- readRDS("/drives/drive1/ff_chaya/Ayanava/dmr_match_r.rds")


# DMR on chr1
dmr1 <- reference %>% filter(reference$seqnames %in% dmr_match_r$chr[7])
dmr1_cpgs <- as.data.frame(strsplit(dmr_match_r$probe[7], ";")) 
colnames(dmr1_cpgs) = "CpGs"

dmr1 <- dmr1 %>% filter((dmr1$start > (dmr_match_r$start[7] - 10000)) & (dmr1$end < (dmr_match_r$end[7] + 10000)))

dmr1_results <- test7ewas_arscp[rownames(test7ewas_arscp) %in% rownames(dmr1), ]
dmr1_filter <- dmr1 %>% filter(rownames(dmr1) %in% rownames(dmr1_results))
dmr1_filter <- dmr1_filter %>% select(c("seqnames", "start", "end"))
dmr1_filter$logFC <- dmr1_results$logFC[match(rownames(dmr1_filter), rownames(dmr1_results))]
dmr1_filter$Target <- FALSE
dmr1_filter$Target[rownames(dmr1_filter) %in% dmr1_cpgs$CpGs] <- TRUE
ggplot(dmr1_filter, aes(x = start, y = logFC)) + 
  geom_point() + 
  geom_line(color = "blue") + 
  geom_point(data = subset(dmr1_filter, Target == TRUE), 
             color = "red", size = 4, shape = 18, alpha = 1) + 
  geom_abline(slope = 0, intercept = 0, color = "black") + 
  scale_y_continuous(limits = c(-2, 0.5))  + 
  scale_x_continuous(limits = c(dmr_match_r$start[7] - 10000, dmr_match_r$end[7] + 10000))+ 
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 9, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10)) + 
  geom_abline(slope = 0, intercept = 0) + 
  labs(x = "CpG Position", title = "Chr1") 

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/chr1_dmr_logFC.png", dpi = 600, width = 6.5, height = 6)

# DMR on chr21
dmr2 <- reference %>% filter(reference$seqnames %in% dmr_match_r$chr[8])
dmr2 <- dmr2 %>% filter((dmr2$start > (dmr_match_r$start[8] - 10000)) & (dmr2$end < (dmr_match_r$end[8] + 10000)))
dmr2_cpgs <- as.data.frame(strsplit(dmr_match_r$probe[8], ";")) 
colnames(dmr2_cpgs) = "CpGs"

dmr2_results <- test7ewas_arscp[rownames(test7ewas_arscp) %in% rownames(dmr2), ]
dmr2_filter <- dmr2 %>% filter(rownames(dmr2) %in% rownames(dmr2_results))
dmr2_filter <- dmr2_filter %>% select(c("seqnames", "start", "end"))
dmr2_filter$logFC <- dmr2_results$logFC[match(rownames(dmr2_filter), rownames(dmr2_results))]
dmr2_filter$logFC <- dmr2_results$logFC[match(rownames(dmr2_filter), rownames(dmr2_results))]
dmr2_filter$Target <- FALSE
dmr2_filter$Target[rownames(dmr2_filter) %in% dmr2_cpgs$CpGs] <- TRUE
ggplot(dmr2_filter, aes(x = start, y = logFC)) + 
  geom_point() + 
  geom_line(color = "blue") + 
  geom_point(data = subset(dmr2_filter, Target == TRUE), 
             color = "red", size = 4, shape = 18, alpha = 1) + 
  geom_abline(slope = 0, intercept = 0, color = "black") + 
  scale_y_continuous(limits = c(-1, 4))  + 
  scale_x_continuous(limits = c(dmr_match_r$start[8] - 10000, dmr_match_r$end[8] + 10000))+ 
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10)) + 
  geom_abline(slope = 0, intercept = 0) + 
  labs(x = "CpG Position", title = "Chr11") 

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/chr11_dmr_logFC.png", dpi = 600, width = 6.5, height = 6)

# DMR on chr11
dmr3 <- reference %>% filter(reference$seqnames %in% dmr_match_r$chr[9])
dmr3 <- dmr3 %>% filter((dmr3$start > (dmr_match_r$start[9] - 10000)) & (dmr3$end < (dmr_match_r$end[9] + 10000)))
dmr3_cpgs <- as.data.frame(strsplit(dmr_match_r$probe[9], ";")) 
colnames(dmr3_cpgs) = "CpGs"


dmr3_results <- test7ewas_arscp[rownames(test7ewas_arscp) %in% rownames(dmr3), ]
dmr3_filter <- dmr3 %>% filter(rownames(dmr3) %in% rownames(dmr3_results))
dmr3_filter <- dmr3_filter %>% select(c("seqnames", "start", "end"))
dmr3_filter$logFC <- dmr3_results$logFC[match(rownames(dmr3_filter), rownames(dmr3_results))]
dmr3_filter$Target <- FALSE
dmr3_filter$Target[rownames(dmr3_filter) %in% dmr3_cpgs$CpGs] <- TRUE
ggplot(dmr3_filter, aes(x = start, y = logFC)) + 
  geom_point() + 
  geom_line(color = "blue") + 
  geom_point(data = subset(dmr3_filter, Target == TRUE), 
             color = "red", size = 4, shape = 18, alpha = 1) + 
  geom_abline(slope = 0, intercept = 0, color = "black") + 
  scale_y_continuous(limits = c(-1.25, 1.25))  + 
  scale_x_continuous(limits = c(dmr_match_r$start[9] - 10000, dmr_match_r$end[9] + 10000))+ 
  theme(axis.title.x = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"), 
        axis.text.x = element_text(size = 8, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        plot.margin = margin(t = 10, r = 30, b = 10, l = 10)) + 
  geom_abline(slope = 0, intercept = 0) + 
  labs(x = "CpG Position", title = "Chr21") 

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/chr21_dmr_logFC.png", dpi = 600, width = 6.5, height = 6)

# graph the p-values and logFCs of all CpGs within a +/- 10,000 bp range of the DMR
graph_dmr <- function(dmr_match, ewas, reference, index, title) {
  dmr <- reference %>% filter(reference$seqnames %in% dmr_match$chr[index])
  dmr <- dmr %>% filter((dmr$start > (dmr_match$start[index] - 10000)) & (dmr$end < (dmr_match$end[index] + 10000)))
  dmr_cpgs <- as.data.frame(strsplit(dmr_match$probe[index], ";")) 
  colnames(dmr_cpgs) = "CpGs"
  dmr_results <- ewas[rownames(ewas) %in% rownames(dmr),]
  dmr_filter <- dmr %>% filter(rownames(dmr) %in% rownames(dmr_results))
  dmr_filter <- dmr_filter %>% select(c("seqnames", "start", "end"))
  dmr_filter$P.Value <- dmr_results$P.Value[match(rownames(dmr_filter), rownames(dmr_results))]
  dmr_filter$logFC <- dmr_results$logFC[match(rownames(dmr_filter), rownames(dmr_results))]
  dmr_filter$Target <- FALSE
  dmr_filter$Target[rownames(dmr_filter) %in% dmr_cpgs$CpGs] <- TRUE
  ggplot(dmr_filter, aes(x = start)) + 
    geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
    geom_line(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
    geom_point(aes(y = -log10(P.Value)), data = subset(dmr_filter, Target == TRUE), 
               color = "red", size = 5, shape = 18, alpha = 1) + 
    geom_point(aes(y = logFC), data = subset(dmr_filter, Target == TRUE), 
               color = "red", size = 5, shape = 18, alpha = 1) + 
    geom_point(aes(y = -log10(P.Value), color = "-log10(p-value)")) + 
    geom_point(aes(y = logFC, color = "logFC")) + 
    geom_line(aes(y = logFC, color = "logFC")) + 
    scale_y_continuous(breaks = seq(-2, 6, by = 1), limits = c(-1, 6.5),name = "-log10(p-value)", sec.axis = sec_axis(~ ., name = "logFC", breaks = seq(-2, 6, by = 1))) + 
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
    labs(title = title)
}

# chr11
graph_dmr(dmr_match_r, test7ewas_arscp, reference, 9, "Y22")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr11_y22.png", dpi = 600, width = 6.5, height = 6)
graph_dmr(dmr_match_r, test6ewas_arscp, reference, 5, "Y15")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr11_y15.png", dpi = 600, width = 6.5, height = 6)
graph_dmr(dmr_match_r, test5ewas_arscp, reference, 3, "Y9")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr11_y9.png", dpi = 600, width = 6.5, height = 6)

#chr1
graph_dmr(dmr_match_r, test7ewas_arscp, reference, 7, "Y22")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr1_y22.png", dpi = 600, width = 6.5, height = 6)
graph_dmr(dmr_match_r, test6ewas_arscp, reference, 4, "Y15")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr1_y15.png", dpi = 600, width = 6.5, height = 6)
graph_dmr(dmr_match_r, test5ewas_arscp, reference, 1, "Y9")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr1_y9.png", dpi = 600, width = 6.5, height = 6)


#chr21
graph_dmr(dmr_match_r, test7ewas_arscp, reference, 8, "Y22")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr21_y22.png", dpi = 600, width = 6.5, height = 6)
graph_dmr(dmr_match_r, test6ewas_arscp, reference, 6, "Y15")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr21_y15.png", dpi = 600, width = 6.5, height = 6)
graph_dmr(dmr_match_r, test5ewas_arscp, reference, 2, "Y9")
ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/dual_chr21_y9.png", dpi = 600, width = 6.5, height = 6)

# annotate the genes to the DMRs found 

combp_results_w6_anc <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas6_new_anc.csv")
combp_results_w5_anc <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas5_new_anc.csv")
combp_results_w7_anc <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas7_new_anc.csv")

combp_results_w6_r <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas6_new_r.csv")
combp_results_w5_r <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas5_new_r.csv")
combp_results_w7_r <- read.csv("/drives/drive1/ff_chaya/Ayanava/resu_combp_ewas7_new_r.csv")


dmrs_w5_anc_gene <- merge(combp_results_w5_anc, reference[, c("seqnames", "start", "gene")], by.x = c("chr", "start"), by.y = c("seqnames", "start"))
dmrs_w5_anc_gene$range <- paste(paste(dmrs_w5_anc_gene$chr, dmrs_w5_anc_gene$start, sep = ":"), dmrs_w5_anc_gene$end, sep = "-")
write.csv(dmrs_w5_anc_gene, "/drives/drive1/ff_chaya/Ayanava/Supp_Materials/dmrs_w5_anc.csv")

dmrs_w6_anc_gene <- merge(combp_results_w6_anc, reference[, c("seqnames", "start", "gene")], by.x = c("chr", "start"), by.y = c("seqnames", "start"))
dmrs_w6_anc_gene$range <- paste(paste(dmrs_w6_anc_gene$chr, dmrs_w6_anc_gene$start, sep = ":"), dmrs_w6_anc_gene$end, sep = "-")
write.csv(dmrs_w6_anc_gene, "/drives/drive1/ff_chaya/Ayanava/Supp_Materials/dmrs_w6_anc.csv")

dmrs_w7_anc_gene <- merge(combp_results_w7_anc, reference[, c("seqnames", "start", "gene")], by.x = c("chr", "start"), by.y = c("seqnames", "start"))
dmrs_w7_anc_gene$range <- paste(paste(dmrs_w7_anc_gene$chr, dmrs_w7_anc_gene$start, sep = ":"), dmrs_w7_anc_gene$end, sep = "-")
write.csv(dmrs_w7_anc_gene, "/drives/drive1/ff_chaya/Ayanava/Supp_Materials/dmrs_w7_anc.csv")

dmrs_w5_r_gene <- merge(combp_results_w5_r, reference[, c("seqnames", "start", "gene")], by.x = c("chr", "start"), by.y = c("seqnames", "start"))
dmrs_w5_r_gene$range <- paste(paste(dmrs_w5_r_gene$chr, dmrs_w5_r_gene$start, sep = ":"), dmrs_w5_r_gene$end, sep = "-")
write.csv(dmrs_w5_r_gene, "/drives/drive1/ff_chaya/Ayanava/Supp_Materials/dmrs_w5_r.csv")

dmrs_w6_r_gene <- merge(combp_results_w6_r, reference[, c("seqnames", "start", "gene")], by.x = c("chr", "start"), by.y = c("seqnames", "start"))
dmrs_w6_r_gene$range <- paste(paste(dmrs_w6_r_gene$chr, dmrs_w6_r_gene$start, sep = ":"), dmrs_w6_r_gene$end, sep = "-")
write.csv(dmrs_w6_r_gene, "/drives/drive1/ff_chaya/Ayanava/Supp_Materials/dmrs_w6_r.csv")

dmrs_w7_r_gene <- merge(combp_results_w7_r, reference[, c("seqnames", "start", "gene")], by.x = c("chr", "start"), by.y = c("seqnames", "start"))
dmrs_w7_r_gene$range <- paste(paste(dmrs_w7_r_gene$chr, dmrs_w7_r_gene$start, sep = ":"), dmrs_w7_r_gene$end, sep = "-")
write.csv(dmrs_w7_r_gene, "/drives/drive1/ff_chaya/Ayanava/Supp_Materials/dmrs_w7_r.csv")



