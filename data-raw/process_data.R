#' _____________________________________________________________________________

#' Code to generate .Rd data files for the package so that when it is built the
#' datasets are immediately available and do not have to be read in and assigned

#' _____________________________________________________________________________

# clear workspace, load in dplyr
rm(list =ls())
library(dplyr)
library(readr)

#' _____________________________________________________________________________


# 'raw' data to be saved the the data folder as .Rd files

dwelltime_bystypestage <- read.csv(file = "~/MCTsymptomatic/data-raw/dwell_times.csv")
sens_bystypestage <- read.csv(file = "~/MCTsymptomatic/data-raw/sensitivity.csv")
waiting_bystage_symp <- read.csv(file = "~/MCTsymptomatic/data-raw/waiting_bystage_symp.csv")
netsurv_bystage <-read.csv(file = "~/MCTsymptomatic/data-raw/netsurv_5year.csv")
symp_cluster_prev <- read.csv(file = "~/MCTsymptomatic/data-raw/symp_cluster_prev.csv")

## Process CPRD data first
CPRD_master <- read.csv(file = "~/MCTsymptomatic/data-raw/CPRD_master.csv", na.string = "")

CPRD_revised <- CPRD_master %>%
  fill(time_frame, .direction = "down") %>%
  mutate(num_staged_cancer = case_when(stage == "Overall" ~ num_staged_cancer,
                                       stage == "1" ~ n_by_stage,
                                       stage == "2" ~ n_by_stage,
                                       stage == "3" ~ n_by_stage,
                                       stage == "4" ~ n_by_stage,
                                       .default = NA)) %>%
  mutate(stage = case_when(stage == "Overall" ~ "Overall",
                                       stage == "1" ~ "I",
                                       stage == "2" ~ "II",
                                       stage == "3" ~ "III",
                                       stage == "4" ~ "IV",
                                       .default = NA)) %>%
  select(-n_by_stage) %>%
  mutate(cancer_type = case_when( cancer_site == "Bone and Connective Tissue" ~ "Other",
                                  cancer_site == "Bowel" ~ "Colorectal",
                                  cancer_site == "Breast" ~ "Breast",
                                  cancer_site == "CNS"~ "Other",
                                  cancer_site == "Gastro-Oesophageal"~ "Oesphageal",
                                  cancer_site == "Head and neck" ~ "Head and neck",
                                  cancer_site == "Hepatobiliary" ~"Liver",
                                  cancer_site == "Leukaemia" ~ NA,
                                  cancer_site == "Lung" ~ "Lung",
                                  cancer_site == "Lymphoma"  ~ "Lymphoma",
                                  cancer_site == "Melanoma of skin"    ~ NA,
                                  cancer_site == "Multiple sites" ~ "Other",
                                  cancer_site == "Myeloma" ~ NA,
                                  cancer_site == "Other"  ~ "Other" ,
                                  cancer_site == "Ovarian" ~ "Ovarian" ,
                                  cancer_site == "Pancreas"  ~ "Pancreatic",
                                  cancer_site == "Prostate" ~ "Prostate" ,
                                  cancer_site == "Renal Tract" ~ "Bladder" ,
                                  cancer_site == "Testis" ~ "Other",
                                  cancer_site == "Uterine" ~ "Uterine",
                                  cancer_site == "Vulvovaginal" ~ "Other" ,
                                  .default = NA), .after = cancer_site)

# For a cohort with "any non-specific symptoms"
CPRD_revised_any <- CPRD_revised %>%
  filter(symptom == "HTA") %>%
  filter(!is.na(cancer_type)) %>%
  mutate(num_cancer = as.numeric(parse_number(word(num_cancer, start = 1)))) %>%
mutate(num_staged_cancer = as.numeric(parse_number(word(num_staged_cancer, start = 1)))) %>%
  mutate(inflator = num_cancer/num_staged_cancer, .after = num_cancer) %>%
  fill(inflator) %>%
  mutate(num_staged_cancer_interpolated = inflator*num_staged_cancer, .after= num_staged_cancer) %>%
  filter(stage != "Overall") %>%
  filter(time_frame == "Cancer diagnosis within 6 months of the symptom ") %>%
  select(cancer_type, stage, num_staged_cancer_interpolated) %>%
  group_by(cancer_type, stage) %>%
  summarise(num_staged_cancer_interpolated_bytype = sum(num_staged_cancer_interpolated, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(HTA = round((num_staged_cancer_interpolated_bytype / sum(num_staged_cancer_interpolated_bytype)) * 100, 2)) %>%
  ungroup()

# Add a made up column for now - this part wiill change when there's actualy data"

CPRD_revised_any <- CPRD_revised_any %>%
  select(cancer_type, stage , HTA) %>%
  mutate(madeup = HTA )

# transfer to data folder ----

usethis::use_data(CPRD_revised_any,
                  dwelltime_bystypestage,
                  sens_bystypestage,
                  waiting_bystage_symp,
                  netsurv_bystage,
                  symp_cluster_prev,
                  overwrite = TRUE)


