#CAUSE summary statistics files and names 

##################
# ACCCELEROMETER #
##################

setwd("C:/Users/Person/Desktop/MR/summary_statistics/")

###################
# Exposure: LIGHT #
###################

doherty_light <- data.table::fread("ACCL_Doherty_sedentary.txt", header=TRUE)
dim(doherty_light)
names(doherty_light)

exposure_data <- format_data(
  doherty_light,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "
  ",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  #min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)

##################
# Exposure: MVPA #
##################

load("ACCL_Ramadan_MVPA.rda")
ramadan_accl <- dat2
dim(ramadan_accl)
head(ramadan_accl)

exposure_data <- format_data(
  ramadan_accl,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "N",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)

##################
# Exposure: VIGR #
##################

klimentidis_vpa <- data.table::fread("ACCL_klimentidis_VPA.txt", header=TRUE)
dim(klimentidis_vpa)

#Use this to reformat data for MR Analyses
exposure_data <- format_data(
  klimentidis_vpa,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)

##################
# Exposure: AVRG #
##################

doherty_average <- data.table::fread("ACCL_Doherty_average.txt", header=TRUE)
dim(doherty_average)
names(doherty_average)

exposure_data <- format_data(
  doherty_average,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  #min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)


###############
# SELF REPORT #
###############

###################
# Exposure: LIGHT #
###################

exposure_data <- data.table::fread("vandev_SNPs.txt", header=TRUE)

exposure_data2 <- exposure_data[exposure_data$pval.exposure<5*10^-8,] 

##################
# Exposure: MODR #
##################

load("moderate_Ramadan.rda")
ramadan_modr <- dat2
dim(ramadan_modr)
names(ramadan_modr)

exposure_data <- format_data(
  ramadan_modr,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "N",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)

##################
# Exposure: MVPA #
##################

klimentidis_mvpa <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)
dim(klimentidis_mvpa)

exposure_data <- format_data(
  klimentidis_mvpa,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)  


###################
# Exposure: VIGR #
###################

klimentidis_vigr <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
dim(klimentidis_vigr)

exposure_data <- format_data(
  klimentidis_vigr,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)

##################
# Exposure: SSOE #
##################

klimentidis_ssoe <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
dim(klimentidis_ssoe)

exposure_data <- format_data(
  klimentidis_ssoe,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "CHR",
  pos_col = "BP",
  log_pval = FALSE
)

####################
# BODY COMPOSITION #
####################

#######
# LBM #
#######

gefos_lbm <- data.table::fread("gefos_LBM.txt", header=TRUE)
gefos_lbm2 <- gefos_lbm[gefos_lbm$`P-value`<5*10^-6,] 

exposure_lbm3 <- format_data(
  gefos_lbm2,
  type = "exposure",
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "TotalSampleSize",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  #chr_col = "hm_chrom",
  #pos_col = "hm_pos",
  log_pval = FALSE
)

#######
# BMI #
#######

locke_bmi <- data.table::fread("locke_bmi.txt", header=TRUE)
dim(locke_bmi)

outcome_data_bmi <- format_data(
  locke_bmi,
  type = "outcome",
  snps = exposure_data2$SNP,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  eaf_col = "eaf_ref",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  pval_col = "p_value",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  chr_col = "hm_chrom",
  pos_col = "hm_pos",
  log_pval = FALSE
)

######
# WC #
######

shungin_wc <- data.table::fread("shungin_wc.txt", header=TRUE)
dim(shungin_wc)

outcome_data_wc <- format_data(
  shungin_wc,
  type = "outcome",
  snps = exposure_data2$SNP,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "MarkerName",
  beta_col = "b",
  se_col = "se",
  eaf_col = "FreqAllele1HapMapCEU",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "N",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  #chr_col = "hm_chrom",
  #pos_col = "hm_pos",
  log_pval = FALSE
)

#######
# WHR #
#######  

shungin_whr <- data.table::fread("shungin_whr.txt", header=TRUE)
dim(shungin_whr)

outcome_data_whr <- format_data(
  shungin_whr,
  type = "outcome",
  snps = exposure_data2$SNP,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "MarkerName",
  beta_col = "b",
  se_col = "se",
  eaf_col = "FreqAllele1HapMapCEU",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "p",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "N",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  #chr_col = "hm_chrom",
  #pos_col = "hm_pos",
  log_pval = FALSE
)

#######
# BF% #
#######

lu_bfp <- data.table::fread("lu_bodyfatpercent.txt", header=TRUE)
dim(lu_bfp)

outcome_data_bfp <- format_data(
  lu_bfp,
  type = "outcome",
  snps = exposure_data2$SNP,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "hm_rsid",
  beta_col = "hm_beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "hm_effect_allele",
  other_allele_col = "hm_other_allele",
  pval_col = "p_value",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  samplesize_col = "n",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  #chr_col = "hm_chrom",
  #pos_col = "hm_pos",
  log_pval = FALSE
)

#######
# VAT #
#######

#Read in summary statistics 
load("VAT_karlsson.rda")
karlsson_vat <- dat2
dim(karlsson_vat)

outcome_data_vat <- format_data(
  karlsson_vat,
  type = "outcome",
  snps = exposure_data2$SNP,
  header = TRUE,
  #phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  #units_col = "units",
  #ncase_col = "ncase",
  #ncontrol_col = "ncontrol",
  #samplesize_col = "N",
  #gene_col = "gene",
  #id_col = "hm_variant_id",
  min_pval = 1e-200,
  #z_col = "z",
  #info_col = "info",
  #chr_col = "hm_chrom",
  #pos_col = "hm_pos",
  log_pval = FALSE
)



