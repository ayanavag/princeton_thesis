library(magrittr)
library(dplyr)
library(tidyr)
library(TwoSampleMR)

# cat("Loading data...\n")
harmon_dat <- readRDS("/drives/drive1/ff_chaya/Ayanava/harmonized_data_ff_gwas.rds")
cat("Pleio 1...\n")
pleio <- mr_pleiotropy_test(harmon_dat)
saveRDS(pleio, "/drives/drive1/ff_chaya/Ayanava/mr_ff_pleio.rds")
cat("Heterogen 1...\n")
hetero <- mr_heterogeneity(harmon_dat)
saveRDS(hetero, "/drives/drive1/ff_chaya/Ayanava/mr_ff_heterogen.rds")


harmonized_data <- readRDS("/drives/drive1/ff_chaya/Ayanava/harmonized_data_uk_gwas.RDS")
cat("Pleio 2...\n")
pleio2 <- mr_pleiotropy_test(harmonized_data)
saveRDS(pleio2, "/drives/drive1/ff_chaya/Ayanava/mr_uk_pleio.rds")
cat("Heterogen 2...\n")
hetero2 <- mr_heterogeneity(harmonized_data)
saveRDS(hetero2, "/drives/drive1/ff_chaya/Ayanava/mr_uk_heterogen.rds")


