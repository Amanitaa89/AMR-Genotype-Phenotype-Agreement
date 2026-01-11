library(tidyverse)

# 1) Load CGE report – بمسارك الصحيح
cge_path <- "/Users/aat/Documents/Enteritidis/WGS_CGE_Report.csv"
cge <- read_csv(cge_path, show_col_types = FALSE)

# 2) Standardize text
cge2 <- cge %>%
  mutate(
    CGE_Pheno = tolower(`CGE Predicted Phenotype`),
    CGE_Pheno = str_replace_all(CGE_Pheno, ",", " , ")
  )

# 3) Define all antibiotics you want to extract (حسب حاجتك)
all_drugs <- c(
  "ciprofloxacin",
  "nalidixic acid",
  "tetracycline",
  "doxycycline",
  "chloramphenicol",
  "ampicillin",
  "amoxicillin",
  "piperacillin",
  "ticarcillin"
)

# 4) Convert text → 0/1 per drug
for (drug in all_drugs) {
  colname <- paste0(str_replace_all(drug, " ", "_"), "_CGE")
  cge2[[colname]] <- ifelse(str_detect(cge2$CGE_Pheno, drug), 1, 0)
}

# 5) Save expanded file
write_csv(cge2, "/Users/aat/Documents/Enteritidis/WGS_CGE_Report_expanded.csv")

message("✅ CGE phenotype parsing complete → WGS_CGE_Report_expanded.csv")
