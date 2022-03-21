############################
# Forest Plot: Self-report #
############################

library(readr)
library(metafor)
library(purrr)
library(dplyr)

setwd("C:\\Users\\Person\\Desktop\\MR\\cause\\")

selfSED_LBM <- read_csv("SELF_SED_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "SED")
selfSED_BMI <- read_csv("SELF_SED_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "SED")
selfSED_WC <- read_csv("SELF_SED_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "SED")
selfSED_WHR <- read_csv("SELF_SED_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "SED")
selfSED_BFP <- read_csv("SELF_SED_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "SED")
selfSED_VAT <- read_csv("SELF_SED_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "SED")

selfMODR_LBM <- read_csv("SELF_MODR_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "MODR")
selfMODR_BMI <- read_csv("SELF_MODR_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "MODR")
selfMODR_WC <- read_csv("SELF_MODR_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "MODR")
selfMODR_WHR <- read_csv("SELF_MODR_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "MODR")
selfMODR_BFP <- read_csv("SELF_MODR_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "MODR")
selfMODR_VAT <- read_csv("SELF_MODR_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "MODR")

selfMVPA_LBM <- read_csv("SELF_MVPA_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "MVPA")
selfMVPA_BMI <- read_csv("SELF_MVPA_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "MVPA")
selfMVPA_WC <- read_csv("SELF_MVPA_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "MVPA")
selfMVPA_WHR <- read_csv("SELF_MVPA_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "MVPA")
selfMVPA_BFP <- read_csv("SELF_MVPA_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "MVPA")
selfMVPA_VAT <- read_csv("SELF_MVPA_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "MVPA")

selfVIGR_LBM <- read_csv("SELF_VIGR_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "VIGR")
selfVIGR_BMI <- read_csv("SELF_VIGR_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "VIGR")
selfVIGR_WC <- read_csv("SELF_VIGR_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "VIGR")
selfVIGR_WHR <- read_csv("SELF_VIGR_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "VIGR")
selfVIGR_BFP <- read_csv("SELF_VIGR_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "VIGR")
selfVIGR_VAT <- read_csv("SELF_VIGR_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "VIGR")

selfSSOE_LBM <- read_csv("SELF_SSOE_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "SSOE")
selfSSOE_BMI <- read_csv("SELF_SSOE_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "SSOE")
selfSSOE_WC <- read_csv("SELF_SSOE_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "SSOE")
selfSSOE_WHR <- read_csv("SELF_SSOE_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "SSOE")
selfSSOE_BFP <- read_csv("SELF_SSOE_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "SSOE")
selfSSOE_VAT <- read_csv("SELF_SSOE_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "SSOE")

library(data.table)
library(stringr)

df <- rbind(
  selfSED_LBM,
  selfSED_BMI,
  selfSED_WC,
  selfSED_WHR,
  selfSED_BFP,
  selfSED_VAT,
  
  selfMODR_LBM,
  selfMODR_BMI,
  selfMODR_WC,
  selfMODR_WHR,
  selfMODR_BFP,
  selfMODR_VAT,
  
  selfMVPA_LBM,
  selfMVPA_BMI,
  selfMVPA_WC,
  selfMVPA_WHR,
  selfMVPA_BFP,
  selfMVPA_VAT,
  
  selfVIGR_LBM,
  selfVIGR_BMI,
  selfVIGR_WC,
  selfVIGR_WHR,
  selfVIGR_BFP,
  selfVIGR_VAT,
  
  selfSSOE_LBM,
  selfSSOE_BMI,
  selfSSOE_WC,
  selfSSOE_WHR,
  selfSSOE_BFP,
  selfSSOE_VAT
) %>%
  mutate(., p = nafill(p, type = "locf")) %>%
  filter(., model == "Causal") %>%
  mutate(b = word(gamma, 1, sep = "\\(")) %>%
  mutate(ci = word(gamma, 2, sep = "\\(")) %>%
  mutate(ci = str_sub(ci, 1, nchar(ci)-1) ) %>%
  mutate(sig = ifelse(p < 0.05, "*",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.001, "***", "")))) %>%
  mutate(ci = paste0(ci,sig)) %>%
  dplyr::select(., c("outcome","exp","b", "ci"))

write.csv(df, "C:\\Users\\Person\\OneDrive\\z1_Research\\THESIS\\manuscript\\final\\revision\\t1.csv")


#########################
#Accelerometer measures #
#########################

setwd("C:\\Users\\Person\\Desktop\\MR\\cause\\")

ACCLSED_LBM <- read_csv("ACCL_SED_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "SED")
ACCLSED_BMI <- read_csv("ACCL_SED_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "SED")
ACCLSED_WC <- read_csv("ACCL_SED_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "SED")
ACCLSED_WHR <- read_csv("ACCL_SED_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "SED")
ACCLSED_BFP <- read_csv("ACCL_SED_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "SED")
ACCLSED_VAT <- read_csv("ACCL_SED_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "SED")

ACCLMVPA_LBM <- read_csv("ACCL_MVPA_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "MVPA")
ACCLMVPA_BMI <- read_csv("ACCL_MVPA_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "MVPA")
ACCLMVPA_WC <- read_csv("ACCL_MVPA_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "MVPA")
ACCLMVPA_WHR <- read_csv("ACCL_MVPA_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "MVPA")
ACCLMVPA_BFP <- read_csv("ACCL_MVPA_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "MVPA")
ACCLMVPA_VAT <- read_csv("ACCL_MVPA_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "MVPA")

ACCLVIGR_LBM <- read_csv("ACCL_VIGR_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "VIGR")
ACCLVIGR_BMI <- read_csv("ACCL_VIGR_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "VIGR")
ACCLVIGR_WC <- read_csv("ACCL_VIGR_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "VIGR")
ACCLVIGR_WHR <- read_csv("ACCL_VIGR_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "VIGR")
ACCLVIGR_BFP <- read_csv("ACCL_VIGR_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "VIGR")
ACCLVIGR_VAT <- read_csv("ACCL_VIGR_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "VIGR")

ACCLAVRG_LBM <- read_csv("ACCL_AVRG_LBM.csv") %>% mutate(., outcome = "LBM") %>% mutate(., exp = "AVRG")
ACCLAVRG_BMI <- read_csv("ACCL_AVRG_BMI.csv") %>% mutate(., outcome = "BMI") %>% mutate(., exp = "AVRG")
ACCLAVRG_WC <- read_csv("ACCL_AVRG_WC.csv") %>% mutate(., outcome = "WC") %>% mutate(., exp = "AVRG")
ACCLAVRG_WHR <- read_csv("ACCL_AVRG_WHR.csv") %>% mutate(., outcome = "WHR") %>% mutate(., exp = "AVRG")
ACCLAVRG_BFP <- read_csv("ACCL_AVRG_BFP.csv") %>% mutate(., outcome = "BFP") %>% mutate(., exp = "AVRG")
ACCLAVRG_VAT <- read_csv("ACCL_AVRG_VAT.csv") %>% mutate(., outcome = "VAT") %>% mutate(., exp = "AVRG")

df <- rbind(
  ACCLSED_LBM,
  ACCLSED_BMI,
  ACCLSED_WHR,
  ACCLSED_WC,
  ACCLSED_BFP,
  ACCLSED_VAT,
  
  ACCLMVPA_LBM,
  ACCLMVPA_BMI,
  ACCLMVPA_WHR,
  ACCLMVPA_WC,
  ACCLMVPA_BFP,
  ACCLMVPA_VAT,
  
  ACCLVIGR_LBM,
  ACCLVIGR_BMI,
  ACCLVIGR_WHR,
  ACCLVIGR_WC,
  ACCLVIGR_BFP,
  ACCLVIGR_VAT,
  
  ACCLAVRG_LBM,
  ACCLAVRG_BMI,
  ACCLAVRG_WHR,
  ACCLAVRG_WC,
  ACCLAVRG_BFP,
  ACCLAVRG_VAT
) %>%
  mutate(., p = nafill(p, type = "locf")) %>%
  filter(., model == "Causal") %>%
  mutate(b = word(gamma, 1, sep = "\\(")) %>%
  mutate(ci = word(gamma, 2, sep = "\\(")) %>%
  mutate(ci = str_sub(ci, 1, nchar(ci)-1) ) %>%
  mutate(sig = ifelse(p < 0.05, "*",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.001, "***", "")))) %>%
  mutate(ci = paste0(ci,sig)) %>%
  dplyr::select(., c("outcome","exp","b", "ci"))

write.csv(df, "C:\\Users\\Person\\OneDrive\\z1_Research\\THESIS\\manuscript\\final\\revision\\t2.csv")
