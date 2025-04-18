# Functions to use for Limma analysis- Ayanava Ganguly

# SVA function 
run_sva <- function(cimt, design, cpgs, log) {
  if (length(design) == 0) {
    if (log == TRUE) {
      mod <- model.matrix(~log_Mean_Mean_CIMT, data = cimt)
    }
    else {
      mod <- model.matrix(~Mean_Mean_CIMT, data = cimt)
    }
    mod0 <- model.matrix(~1, data = cimt)
  }
  else {
    f <- paste(design, collapse = "+")
    if (log == TRUE) {
      f2 <- as.formula(paste("~log_Mean_Mean_CIMT+", f))
    }
    else {
      f2 <- as.formula(paste("~Mean_Mean_CIMT+", f))
    }
    f <- as.formula(paste("~", f))
    mod <- model.matrix(f2, data = cimt)
    mod0 <- model.matrix(f, data = cimt)
  }
  cpgs2 <- cpgs %>% filter((rownames(cpgs) %in% rownames(cimt)) %>% replace_na(TRUE))
  n.sv = num.sv(t(cpgs2), mod, method = "be")
  svobj = sva(t(cpgs2), mod, mod0, n.sv = n.sv)
  print("SVA analysis complete") 

  pValues = f.pvalue(t(cpgs2), mod, mod0)
  qValues = p.adjust(pValues, method = "BH")

  modSv = cbind(mod,svobj$sv)
  mod0Sv = cbind(mod0,svobj$sv)
  pValuesSv = f.pvalue(t(cpgs2),modSv,mod0Sv)
  qValuesSv <- as.data.frame(p.adjust(pValuesSv,method="BH"))
  names(qValuesSv) <- "qValue"
  qValuesSv <- cbind(qValuesSv, as.data.frame(pValuesSv))

  fit <- lmFit(t(cpgs2), modSv)
  fit <- eBayes(fit)
  if (log == TRUE){
    top <- topTable(fit, coef = "log_Mean_Mean_CIMT", number = length(cpgs2))
    
  }
  else {
    top <- topTable(fit, coef = "Mean_Mean_CIMT", number = length(cpgs2))
  }
  # top <- qValuesSv
  return(top)
}


# uses ctrlsva batch variables, instead of generating them from SVA
run_ewas <- function(cimt, design, cpgs, log, batch) {
  if (length(design) == 0) {
    if (log == TRUE) {
      mod <- model.matrix(~log_Mean_Mean_CIMT, data = cimt)
    }
    else {
      mod <- model.matrix(~Mean_Mean_CIMT, data = cimt)
    }
  }
  else {
    f <- paste(design, collapse = "+")
    if (log == TRUE) {
      f2 <- as.formula(paste("~log_Mean_Mean_CIMT+", f))
    }
    else {
      f2 <- as.formula(paste("~Mean_Mean_CIMT+", f))
    }
    f <- as.formula(paste("~", f))
    mod <- model.matrix(f2, data = cimt)
  }
  cpgs2 <- cpgs %>% filter((rownames(cpgs) %in% rownames(cimt)) %>% replace_na(TRUE))
  mod <- as.data.frame(mod)
  batch2 <- batch[rownames(batch) %in% rownames(mod), ]
  modSv <- merge(mod, batch2, by = "row.names", all = TRUE)
  rownames(modSv) <- modSv$Row.names
  modSv <- modSv[, -1]
  
  fit <- lmFit(t(cpgs2), modSv)
  fit <- eBayes(fit)
  if (log == TRUE){
    top <- topTable(fit, coef = "log_Mean_Mean_CIMT", number = length(cpgs2))
    
  }
  else {
    top <- topTable(fit, coef = "Mean_Mean_CIMT", number = length(cpgs2))
  }
  return(top)
}




# Generate covariate matrix for SVA
cov_matrix <- function(design, cimt, cpgs){
  test <- cimt %>% select(all_of(design))
  rownames(test) <- rownames(cpgs)
  test_cpg <- merge(cpgs, test, by = 0)
  test_cpg2 <- test_cpg[, -1]
  rownames(test_cpg2) <- test_cpg[, 1]
  return(test_cpg2)
}


# set up methyl/survey/CIMT data matrix 
clean_methyl <- function(methyl, cpgs, cimtsurvey) {
  cpgsite <- methyl %>% filter(rownames(methyl) %in% cpgs$CpGs)
  matched_cpgs <- cpgs %>% filter(cpgs$CpGs %in% rownames(cpgsite))
  ids <- colnames(cpgsite)
  cimt_epic <- cimtsurvey %>% filter(CVID %in% ids)
  cimt_epic <- cimt_epic[!is.na(cimt_epic$Mean_Mean_CIMT), ]
  cpgsite <- cpgsite %>% select(cimt_epic$CVID)
  cpgsite <- data.frame(t(cpgsite))
  cimt_epic <- cimt_epic[!duplicated(cimt_epic$CVID), ]
  rownames(cpgsite) <- cimt_epic$CVID
  return(list(cpgsite, cimt_epic))
}

# set up methyl/survey/CIMT data matrix for y22 blood samples (which need cell proportion estimates)
clean_methyl7 <- function(methyl, cpgs, cimtsurvey) {
  cpgsite <- methyl %>% filter(rownames(methyl) %in% cpgs$CpGs)
  matched_cpgs <- cpgs %>% filter(cpgs$CpGs %in% rownames(cpgsite))
  ids <- colnames(cpgsite)
  cimt_epic <- cimtsurvey %>% filter(CVID %in% ids)
  cimt_epic <- cimt_epic[!duplicated(cimt_epic$CVID), ]
  lc <- estimateLC(methyl, ref = "Salas")
  cimt_epic <- cbind(cimt_epic, lc)
  cimt_epic <- cimt_epic[!is.na(cimt_epic$Mean_Mean_CIMT), ]
  cpgsite <- cpgsite %>% select(cimt_epic$CVID)
  cpgsite <- data.frame(t(cpgsite))
  rownames(cpgsite) <- cimt_epic$CVID
  return(list(cpgsite, cimt_epic))
}


# identify CpGs with significant P-value, large effect size, and progressive methylation from Y9 to Y22
get_trendcpgs <- function(wave5, wave6, wave7) {
  sig_cpgs_w7 <- wave7 %>% filter((wave7$logFC > log2(1.5) | wave7$logFC < -log2(1.5)) & (wave7$adj.P.Val < 0.05))
  tracked_cpgs_w6 <- wave6 %>% filter(rownames(wave6) %in% rownames(sig_cpgs_w7))
  tracked_cpgs_w5 <- wave5 %>% filter(rownames(wave5) %in% rownames(sig_cpgs_w7))
  tracked_cpgs_w7 <- sig_cpgs_w7 %>% filter(rownames(sig_cpgs_w7) %in% rownames(tracked_cpgs_w5))
  tracked_cpgs_w7 <- tracked_cpgs_w7 %>% filter(rownames(tracked_cpgs_w7) %in% rownames(tracked_cpgs_w6))
  tracked_cpgs_w7 <- tracked_cpgs_w7[order((rownames(tracked_cpgs_w7))), ]
  tracked_cpgs_w5 <- tracked_cpgs_w5[order((rownames(tracked_cpgs_w5))), ]
  tracked_cpgs_w6 <- tracked_cpgs_w6[order((rownames(tracked_cpgs_w6))), ]
  
  all_pval_log <- tracked_cpgs_w5 %>% select("logFC_5" = logFC, "pval_5" = P.Value) %>% mutate("logFC_6" = tracked_cpgs_w6$logFC) %>%
    mutate("pval_6" = tracked_cpgs_w6$P.Value) %>% mutate("logFC_7" = tracked_cpgs_w7$logFC) %>%
    mutate("pval_7" = tracked_cpgs_w7$P.Value)
  
  pval_log_ratios <- all_pval_log
  pval_log_ratios$pval_5 <- all_pval_log$pval_7 / all_pval_log$pval_5 
  pval_log_ratios$pval_6 <- all_pval_log$pval_7 / all_pval_log$pval_6 
  pval_log_ratios$pval_7 <- all_pval_log$pval_7 / all_pval_log$pval_7
  sig_pval_ratios <- pval_log_ratios
  # sig_pval_ratios <- pval_log_ratios %>% filter(pval_log_ratios$pval_5 < pval_log_ratios$pval_6) not necessarily needed to determine
  trend_cpgs <- sig_pval_ratios %>% filter(((sig_pval_ratios$logFC_5 > sig_pval_ratios$logFC_6) & (sig_pval_ratios$logFC_6 > sig_pval_ratios$logFC_7)) | ((sig_pval_ratios$logFC_5 < sig_pval_ratios$logFC_6) & (sig_pval_ratios$logFC_6 < sig_pval_ratios$logFC_7)))
}


# create the dataframe needed to run DMR analysis using combp 
get_combp_df <- function(reference, ewas) {
  combp_df <- reference %>% select(seqnames, start, end)
  combp_df <- combp_df %>% filter(rownames(combp_df) %in% rownames(ewas))
  df <- ewas %>% filter(rownames(ewas) %in% rownames(combp_df))
  combp_df <- combp_df[order(rownames(combp_df)), ]
  df <- df[order(rownames(df)), ]
  # generates a df with chr, start, end, p, and probe 
  combp_df <- combp_df %>% mutate("p" = df$P.Value, "probe" = rownames(df))
  colnames(combp_df) <- c("chr", "start", "end", "p", "probe")
  return (combp_df)
}

# make volcano plot of candidate CpGs
make_volcano <- function(datatable, plottitle) {
  table <- datatable
  table$diffexpr <- "Not_sig"
  table$diffexpr[table$logFC > log2(1.5) & table$adj.P.Val < (0.05)] <- "Upreg"
  table$diffexpr[table$logFC < -log2(1.5) & table$adj.P.Val < (0.05)] <- "Downreg"
  group.colors <- c(Downreg = "red", Not_sig = "black", Upreg = "blue")
  pval_thres <- max(table$P.Value[table$adj.P.Val < 0.05], na.rm = TRUE)
  if (pval_thres < 0) {
    pval_thres = 0.05 / nrow(table)
  }
  ggplot(table, aes(x = logFC, y = -log10(P.Value), col = diffexpr)) + 
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), col = '#4682b4', linetype = 'dashed') +
    geom_hline(yintercept = -log10(pval_thres), col = '#4682b4', linetype = 'dashed') +
    geom_point() + 
    scale_color_manual(values = group.colors, name = "Methylation Status", 
                       labels = c(Downreg = "Hypomethylated", Not_sig = "Not significant", 
                                  Upreg = "Hypermethylated")) + 
    xlim(-2.1, 2.1) + 
    ylim(0, 6.5) + 
    labs(title = title) + 
    theme(plot.title = element_text(size = "18", face = "bold", hjust = 0.5), 
          axis.title.x = element_text(size = "14", face = "bold"), 
          axis.title.y = element_text(size = "14", face = "bold"), 
          legend.title = element_text(face = "bold"), 
          legend.text = element_text(face = "bold"))
}



