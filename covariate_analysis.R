source("/drives/drive1/ff_chaya/Ayanava/limma_functions.R")
library(ggplot2)

# read in the CIMT and covariate data
cimtdataset <- readRDS("/drives/drive1/ff_chaya/R_versions/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2_datafile.RDS")
biomarker <- readRDS("/drives/drive1/ff_chaya/R_versions/FF_biomarker_CVID_5.24/FF_biomarker_CVID_5.24_datafile.RDS")
ancestry <- read.csv("/drives/drive1/ff_chaya/Ada/PsychChip2b.ancestryPCs.CVID.csv")
y22_age <- read.csv("/drives/drive1/ff_chaya/Ayanava/CVID_y22_age.csv")

cimt <- cimtdataset %>% select("CVID", "Mean_Mean_CIMT")
surveydata <- readRDS("/drives/drive1/ff_chaya/R_versions/FF_allwaves_wY22_CVID_3.24/FF_allwaves_wY22_CVID_3.24_datafile.RDS")
cimt_AND_survey <- surveydata %>% filter(surveydata$CVID %in% cimt$CVID)
cimt_AND_survey <- merge(cimt_AND_survey, cimt, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, biomarker, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, ancestry, by = "CVID", all = TRUE)
cimt_AND_survey <- merge(cimt_AND_survey, y22_age, by = "CVID", all = TRUE)
cimt_AND_survey <- cimt_AND_survey[!duplicated(as.list(cimt_AND_survey))]
cimt_AND_survey <- cimt_AND_survey[!duplicated(cimt_AND_survey$CVID, )]
cimt_AND_survey$log_Mean_Mean_CIMT <- log(cimt_AND_survey$Mean_Mean_CIMT)
cimt_AND_survey$k7inches <- ((cimt_AND_survey$k7i40)*12 + cimt_AND_survey$k7i41)


# testing the associations between different possible covariates and the CIMT outcome variable
test_cov <- function(data, covariate, outcome, factor) {
  if (factor == TRUE) {
    data[[covariate]] <- as.factor(data[[covariate]])
    formula <- as.formula(paste(outcome, "~", covariate))
    model_cov <- lm(formula, data = data)
    plot <- ggplot(data = data, aes_string(x = covariate, y = outcome)) + 
      geom_boxplot()
  }
  else {
    data_clean <- data[!is.na(data[[covariate]]), ]
    data_clean <- data_clean %>% filter(data_clean[[covariate]] > -1)
    formula <- as.formula(paste(outcome, "~", covariate))
    model_cov <- lm(formula, data = data_clean)
    plot <- ggplot(data = data_clean, aes_string(x = covariate, y = outcome)) + 
      geom_point() + 
      geom_smooth(method = "lm")
  }
  pval <- summary(model_cov)$coefficients[, 4]
  beta <- summary(model_cov)$coefficients[2, 1]
  pval <- pval[-1]
  print(plot)
  if (any(pval < 0.05)) {
    return (print(paste0("Keep this covariate! P-val = ", pval, "and beta= ", beta)))
  }
  else {
    return (print(paste0("Do not keep this covariate P-val = ", pval, "and beta= ", beta)))
  }
}

# set of tests for age
test_cov(cimt_AND_survey, "Age_k7", "Mean_Mean_CIMT", FALSE) 
test_cov(cimt_AND_survey, "k6me_age", "Mean_Mean_CIMT", FALSE) 
test_cov(cimt_AND_survey, "k5me_age", "Mean_Mean_CIMT", FALSE) 

# set of tests for race
test_cov(cimt_AND_survey, "ck7ethrace", "Mean_Mean_CIMT", TRUE) 

# test for poverty 
test_cov(cimt_AND_survey, "ck7povca", "Mean_Mean_CIMT", TRUE) 
test_cov(cimt_AND_survey, "cp6povca", "Mean_Mean_CIMT", TRUE) 
test_cov(cimt_AND_survey, "cm5povca", "Mean_Mean_CIMT", TRUE) 

# test for sex 
test_cov(cimt_AND_survey, "cm1bsex", "Mean_Mean_CIMT", TRUE) 

# test for height
test_cov(cimt_AND_survey, "ch5chtcm", "Mean_Mean_CIMT", FALSE) 
test_cov(cimt_AND_survey, "ck6chtcm", "Mean_Mean_CIMT", FALSE) 
test_cov(cimt_AND_survey, "k7inches", "Mean_Mean_CIMT", FALSE) 

# test for bmi 
test_cov(cimt_AND_survey, "ch5cbmi", "Mean_Mean_CIMT", FALSE) 
test_cov(cimt_AND_survey, "ck6cbmi", "Mean_Mean_CIMT", FALSE)
test_cov(cimt_AND_survey, "ck7bmi", "Mean_Mean_CIMT", FALSE)

# test for ancestry 
test_cov(cimt_AND_survey, "PC1", "Mean_Mean_CIMT", FALSE)  
test_cov(cimt_AND_survey, "PC2", "Mean_Mean_CIMT", FALSE)  
test_cov(cimt_AND_survey, "PC3", "Mean_Mean_CIMT", FALSE)  
test_cov(cimt_AND_survey, "PC4", "Mean_Mean_CIMT", FALSE) 
test_cov(cimt_AND_survey, "PC5", "Mean_Mean_CIMT", FALSE) 


# testing collinearity between age and BMI
test_cov(cimt_AND_survey, "ch5cbmi", "k5me_age", FALSE) 
test_cov(cimt_AND_survey, "ck6cbmi", "k6me_age", FALSE) 
test_cov(cimt_AND_survey, "ck7bmi", "Age_k7", FALSE) 


# testing collinearity between age and height
test_cov(cimt_AND_survey, "ch5chtcm", "k5me_age", FALSE) 
test_cov(cimt_AND_survey, "ck6chtcm", "k6me_age", FALSE) 
test_cov(cimt_AND_survey, "k7inches", "Age_k7", FALSE) 

# testing collinearity between sex and height
test_cov(cimt_AND_survey, "cm1bsex", "ch5chtcm", TRUE) 
test_cov(cimt_AND_survey, "cm1bsex", "ck6chtcm", TRUE)
test_cov(cimt_AND_survey, "cm1bsex", "k7inches", TRUE) 

# testing collinearity between age and sex
test_cov(cimt_AND_survey, "cm1bsex", "k5me_age", TRUE) 
test_cov(cimt_AND_survey, "cm1bsex", "k6me_age", TRUE)
test_cov(cimt_AND_survey, "cm1bsex", "Age_k7", TRUE) 

# testing collinearity between sex and BMI
test_cov(cimt_AND_survey, "cm1bsex", "ch5cbmi", TRUE) 
test_cov(cimt_AND_survey, "cm1bsex", "ck6cbmi", TRUE)
test_cov(cimt_AND_survey, "cm1bsex", "ck7bmi", TRUE) 

# testing collinearity between PC1 and BMI
test_cov(cimt_AND_survey, "ch5cbmi", "PC1", FALSE) 
test_cov(cimt_AND_survey, "ck6cbmi", "PC1", FALSE) 
test_cov(cimt_AND_survey, "ck7bmi", "PC1", FALSE) 

test_cov(cimt_AND_survey, "ch5cbmi", "ch5chtcm", FALSE) 
test_cov(cimt_AND_survey, "ck6cbmi", "ck6chtcm", FALSE) 
test_cov(cimt_AND_survey, "ck7bmi", "k7inches", FALSE) 

# testing collinearity between race and height
test_cov(cimt_AND_survey, "ck7ethrace","cm5povca", TRUE) 
test_cov(cimt_AND_survey, "ck7ethrace","cp6povca", TRUE) 
test_cov(cimt_AND_survey, "ck7ethrace","ck7povca", TRUE) 

# testing collinearity between PC1 and Height
test_cov(cimt_AND_survey, "ch5chtcm", "PC1", FALSE) 
test_cov(cimt_AND_survey, "ck6chtcm", "PC1", FALSE) 
test_cov(cimt_AND_survey, "k7inches", "PC1", FALSE) 

# testing collinearity between PC1 and Race
test_cov(cimt_AND_survey, "ck7ethrace", "PC1", TRUE) 

# testing collinearity between PC1 and Sex 
test_cov(cimt_AND_survey, "cm1bsex", "PC1", TRUE)

#testing collinearity between age and PC1 
test_cov(cimt_AND_survey, "k5me_age", "PC1", FALSE) 
test_cov(cimt_AND_survey, "k6me_age", "PC1", FALSE) 
test_cov(cimt_AND_survey, "Age_k7", "PC1", FALSE) 

test_cov(cimt_AND_survey, "k5me_age", "ck7ethrace", FALSE) 
test_cov(cimt_AND_survey, "k6me_age", "ck7ethrace", FALSE) 
test_cov(cimt_AND_survey, "Age_k7", "ck7ethrace", FALSE) 

# testing collinearity between cell proportions 
test_cov(mod5_epic_nona, "w5_leu", "w5_epi", FALSE) 
test_cov(batch_blood, "PC2", "PC3", FALSE)
test_cov(batch_blood, "PC1", "PC4", FALSE)
test_cov(mod7_epic_nona, "GR", "B", FALSE)



library(ggpubr)

# race
ggplot(cimt_AND_survey, aes(x = factor(ck7ethrace), y = Mean_Mean_CIMT, fill = factor(ck7ethrace))) + 
  geom_boxplot(outlier.shape = 21, outlier.fill = "red") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)), 
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)), 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) + 
  scale_fill_brewer(palette ="Blues") + 
  labs(title = "CIMT Measures Compared by Race/Ancestry", x = "Self-reported Race/Ethnicity", y = "Mean-Mean CIMT (mm)") + 
  scale_x_discrete(labels = c("1" = "White", "2" = "Black", "3" = "Hispanic", "4" = "Other", "5" = "Multiracial")) + 
  stat_compare_means(method = "anova", label = "p.signif", color = "blue", size = 5, label.y.npc = 0.8) +
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"),
                                        c("1", "4"), c("1", "5")),
                     method = "t.test",
                     label = "p.signif", size = 5, 
                     p.adjust.method = "BH")

ggsave("/drives/drive1/ff_chaya/Ayanava/Plots/cimt_race.png", p3, dpi = 300, width = 6.5, height = 6)

# sample sizes
cimt_AND_survey %>% filter(!is.na(Mean_Mean_CIMT)) %>% 
  group_by(ck7ethrace) %>% 
  summarise(n = n())

# sex 
s1 <- ggplot(cimt_AND_survey, aes(x = droplevels(factor(cm1bsex)), y = Mean_Mean_CIMT, fill = factor(cm1bsex))) + 
  geom_boxplot(outlier.shape = 21, outlier.fill = "red") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)), 
        axis.title.x = element_text(hjust = 0.35, size = 14, face = "bold", margin = margin(t = 10)), 
        plot.title = element_text(hjust = 0.02, size = 16, face = "bold")) + 
  scale_fill_brewer(palette ="Blues") + 
  labs(title = "CIMT Measures Compared by Sex", x = "Sex", y = "Mean-Mean CIMT (mm)") + 
  scale_x_discrete(labels = c("1" = "Male", "2" = "Female")) + 
  stat_compare_means(method = "anova", label = "p.signif", color = "blue", size = 5, label.y.npc = 1.0) +
  stat_compare_means(comparisons = list(c("1", "2")),
                     method = "t.test",
                     label = "p.signif", size = 5, 
                     p.adjust.method = "BH") + 
  ylim(0.3, 0.8)

ggsave("/drives/drive1/ff_chaya/Ayanava/Supp_Materials/sex_cov.png", s1, dpi = 300, width = 6.5, height = 6)


# sample sizes
cimt_AND_survey %>% filter(!is.na(Mean_Mean_CIMT)) %>% 
  group_by(cm1bsex) %>% 
  summarise(n = n())


# poverty 

cimt_pov <- cimt_AND_survey %>% filter(!is.na(ck7povca), ck7povca >= 1)

s2 <- ggplot(cimt_pov, aes(x = factor(ck7povca), y = Mean_Mean_CIMT, fill = factor(ck7povca))) + 
  geom_boxplot(outlier.shape = 21, outlier.fill = "red") + 
  theme_minimal() + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = 12, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)), 
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)), 
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) + 
  scale_fill_brewer(palette ="Blues") + 
  labs(title = "CIMT Measures Compared by Poverty Level", x = "% of Federal Poverty Line", y = "Mean-Mean CIMT (mm)") + 
  scale_x_discrete(labels = c("1" = "<50%", "2" = "50-99%", "3" = "100-199%", "4" = "200-299%", "5" = ">300%")) + 
  stat_compare_means(method = "anova", label = "p.signif", color = "blue", size = 5, label.y.npc = 0.93) +
  stat_compare_means(comparisons = list(c("1", "2"), c("1", "3"),
                                        c("1", "4"), c("1", "5")),
                     method = "t.test",
                     label = "p.signif", size = 5, 
                     p.adjust.method = "BH")

ggsave("/drives/drive1/ff_chaya/Ayanava/Supp_Materials/pov_cov.png", s2, dpi = 300, width = 6.5, height = 6)


# sample sizes
cimt_AND_survey %>% filter(!is.na(Mean_Mean_CIMT)) %>% 
  group_by(cp6povca) %>% 
  summarise(n = n())



# age 
library(broom)


make_cov_plot <- function(data, covariate, outcome, title) {
  data_clean <- data[!is.na(data[[covariate]]), ]
  data_clean <- data_clean %>% filter(data_clean[[covariate]] > -1)
  formula <- as.formula(paste(outcome, "~", covariate))
  model_cov <- lm(formula, data = data_clean)
  model_summ <- glance(model_cov)
  r2 <- round(model_summ$r.squared, 3)
  pval <- model_summ$p.value
  pstar <- case_when(
    pval < 0.0001 ~ "p < 0.0001",
    pval < 0.001 ~ "p < 0.001", 
    pval < 0.01 ~ "p < 0.01", 
    pval < 0.05 ~ "p < 0.05", 
    TRUE ~ "ns"
  )
  label_text <- paste0("R-squared = ", r2, ", ", pstar)
  x_pos <- max(data_clean[[covariate]], na.rm = T) - 0.22 * diff(range(data_clean[[covariate]], na.rm = T))
  y_pos <- max(data_clean[[outcome]], na.rm = T) - 0.05 * diff(range(data_clean[[outcome]], na.rm = T))
  plot <- ggplot(data = data_clean, aes_string(x = covariate, y = outcome)) + 
    geom_point() + 
    geom_smooth(method = "lm") + 
    coord_cartesian(xlim = c(min(data_clean[[covariate]], na.rm = T), max(data_clean[[covariate]], na.rm = T)), 
                    ylim = c(min(data_clean[[outcome]], na.rm = T), max(data_clean[[outcome]], na.rm = T))) + 
    annotate("text", label = label_text, x = x_pos, y = y_pos, fontface = "bold") + 
    theme(axis.text.x = element_text(size = 12, face = "bold"), 
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)), 
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)), 
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) + 
    labs(title = paste0("CIMT ~ ", title), x = title)
  print(plot)
  file <- paste0("/drives/drive1/ff_chaya/Ayanava/Supp_Materials/cov_", covariate, ".png")
  ggsave(file, plot, dpi = 300, width = 6.5, height = 6)
}

make_cov_plot(cimt_AND_survey, "Age_k7", "Mean_Mean_CIMT", "Age at Y22")
make_cov_plot(cimt_AND_survey, "k5me_age", "Mean_Mean_CIMT", "Age at Y9")
make_cov_plot(cimt_AND_survey, "k6me_age", "Mean_Mean_CIMT", "Age at Y15")
make_cov_plot(cimt_AND_survey, "ch5cbmi", "Mean_Mean_CIMT", "BMI at Y9")
make_cov_plot(cimt_AND_survey, "ck6cbmi", "Mean_Mean_CIMT", "BMI at Y15")
make_cov_plot(cimt_AND_survey, "ck7bmi", "Mean_Mean_CIMT", "BMI at Y22")
make_cov_plot(cimt_AND_survey, "ch5chtcm", "Mean_Mean_CIMT", "Height at Y9")
make_cov_plot(cimt_AND_survey, "ck6chtcm", "Mean_Mean_CIMT", "Height at Y15")
make_cov_plot(cimt_AND_survey, "k7inches", "Mean_Mean_CIMT", "Height at Y22")

make_cov_plot(cimt_AND_survey, "PC1", "Mean_Mean_CIMT", "Ancestry PC1")
make_cov_plot(cimt_AND_survey, "PC2", "Mean_Mean_CIMT", "Ancestry PC2")
make_cov_plot(cimt_AND_survey, "PC3", "Mean_Mean_CIMT", "Ancestry PC3")
make_cov_plot(cimt_AND_survey, "PC4", "Mean_Mean_CIMT", "Ancestry PC4")
make_cov_plot(cimt_AND_survey, "PC5", "Mean_Mean_CIMT", "Ancestry PC5")





