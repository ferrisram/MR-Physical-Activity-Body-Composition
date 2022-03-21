#CAUSE Analysis for Accelerometer-measured PA and body composition exposures

library(TwoSampleMR)
library(readr)
library(dplyr)
#devtools::install_github("jean997/cause@v1.2.0")
library(cause)

bfile <- "C:\\Users\\Person\\Desktop\\gene\\ref\\EUR"
setwd("C:/Users/Person/Desktop/MR/summary_statistics/")

#############
# SEDENTARY #
#############

#SEDENTARY -> LBM
X1 <- data.table::fread("vandev_Television_LIGHT.txt", header=TRUE)
X2 <- data.table::fread("gefos_LBM.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "Effect"), 
                 se_cols = c("SE", "StdErr"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "P-value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SED_LBM.csv")



#SEDENTARY -> BMI
rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("vandev_Television_LIGHT.txt", header=TRUE)
X2 <- data.table::fread("locke_bmi.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SED_BMI.csv")


#SEDENTARY -> WHR
rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("vandev_Television_LIGHT.txt", header=TRUE)
X2 <- data.table::fread("shungin_whr.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SED_WHR.csv")




#SEDENTARY -> WC
rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("vandev_Television_LIGHT.txt", header=TRUE)
X2 <- data.table::fread("shungin_wc.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SED_WC.csv")




#SEDENTARY -> %TBF
rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("vandev_Television_LIGHT.txt", header=TRUE)
X2 <- data.table::fread("lu_bodyfatpercent.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SED_BFP.csv")




#SEDENTARY -> VAT
rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("vandev_Television_LIGHT.txt", header=TRUE)
load("VAT_karlsson.rda")
X2 <- dat2

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "SNP"), 
                 beta_hat_cols = c("BETA", "BETA"), 
                 se_cols = c("SE", "SE"), 
                 A1_cols = c("ALLELE1", "ALLELE1"), 
                 A2_cols = c("ALLELE0", "ALLELE0"),
                 pval_cols = c("P_BOLT_LMM_INF", "P_BOLT_LMM_INF")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SED_VAT.csv")






############
# MODERATE #
############

#MODERATE -> LBM
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

load("moderate_Ramadan.rda")
X1 <- dat2
X2 <- data.table::fread("gefos_LBM.txt", header=TRUE)

dim(X1)
dim(X2)

head(X1)
head(X2)

names(X1)
names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "Effect"), 
                 se_cols = c("SE", "StdErr"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "P-value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MODR_LBM.csv")



#MODERATE -> BMI
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

load("moderate_Ramadan.rda")
X1 <- dat2
X2 <- data.table::fread("locke_bmi.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MODR_BMI.csv")


#MODDERATE -> WHR
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

load("moderate_Ramadan.rda")
X1 <- dat2
X2 <- data.table::fread("shungin_whr.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MODR_WHR.csv")




#MODDERATE -> WC
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

load("moderate_Ramadan.rda")
X1 <- dat2
X2 <- data.table::fread("shungin_wc.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MODR_WC.csv")




#MODDERATE -> %TBF
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

load("moderate_Ramadan.rda")
X1 <- dat2
X2 <- data.table::fread("lu_bodyfatpercent.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MODR_BFP.csv")




#MODDERATE -> VAT
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

load("moderate_Ramadan.rda")
X1 <- dat2

load("VAT_karlsson.rda")
X2 <- dat2

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "SNP"), 
                 beta_hat_cols = c("BETA", "BETA"), 
                 se_cols = c("SE", "SE"), 
                 A1_cols = c("ALLELE1", "ALLELE1"), 
                 A2_cols = c("ALLELE0", "ALLELE0"),
                 pval_cols = c("P_BOLT_LMM_INF", "P_BOLT_LMM_INF")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MODR_VAT.csv")







########
# MVPA #
########

#MVPA -> LBM
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

X1 <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)
X2 <- data.table::fread("gefos_LBM.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "Effect"), 
                 se_cols = c("SE", "StdErr"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "P-value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)
res2

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MVPA_LBM.csv")



#MVPA -> BMI
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

X1 <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)
X2 <- data.table::fread("locke_bmi.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MVPA_BMI.csv")


#MVPA -> WHR
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

X1 <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)
X2 <- data.table::fread("shungin_whr.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MVPA_WHR.csv")




#MVPA -> WC
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

X1 <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)
X2 <- data.table::fread("shungin_wc.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MVPA_WC.csv")




#MVPA -> %TBF
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

X1 <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)
X2 <- data.table::fread("lu_bodyfatpercent.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MVPA_BFP.csv")




#MVPA -> VAT
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist); rm(dat2)

X1 <- data.table::fread("klimentidis_MVPA_Model1_BOLTLMM_500K.txt", header=TRUE)

load("VAT_karlsson.rda")
X2 <- dat2

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "SNP"), 
                 beta_hat_cols = c("BETA", "BETA"), 
                 se_cols = c("SE", "SE"), 
                 A1_cols = c("ALLELE1", "ALLELE1"), 
                 A2_cols = c("ALLELE0", "ALLELE0"),
                 pval_cols = c("P_BOLT_LMM_INF", "P_BOLT_LMM_INF")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_MVPA_VAT.csv")




########
# VIGR #
########

#VIGOROUS -> LBM
X1 <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
X2 <- data.table::fread("gefos_LBM.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "Effect"), 
                 se_cols = c("SE", "StdErr"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "P-value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)
res2

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_VIGR_LBM.csv")



#VIGOROUS -> BMI
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
X2 <- data.table::fread("locke_bmi.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_VIGR_BMI.csv")


#VIGOROUS -> WHR
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
X2 <- data.table::fread("shungin_whr.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_VIGR_WHR.csv")




#VIGOROUS -> WC
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
X2 <- data.table::fread("shungin_wc.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_VIGR_WC.csv")




#VIGOROUS -> %TBF
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
X2 <- data.table::fread("lu_bodyfatpercent.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_VIGR_BFP.csv")




#VIGOROUS -> VAT
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_VPA.txt", header=TRUE)
load("VAT_karlsson.rda")
X2 <- dat2

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "SNP"), 
                 beta_hat_cols = c("BETA", "BETA"), 
                 se_cols = c("SE", "SE"), 
                 A1_cols = c("ALLELE1", "ALLELE1"), 
                 A2_cols = c("ALLELE0", "ALLELE0"),
                 pval_cols = c("P_BOLT_LMM_INF", "P_BOLT_LMM_INF")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_VIGR_VAT.csv")




########
# SSOE #
########

#SSOE -> LBM
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
X2 <- data.table::fread("gefos_LBM.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "Effect"), 
                 se_cols = c("SE", "StdErr"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "P-value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)
res2

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SSOE_LBM.csv")




#SSOE -> BMI
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
X2 <- data.table::fread("locke_bmi.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SSOE_BMI.csv")




#SSOE -> WHR
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
X2 <- data.table::fread("shungin_whr.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SSOE_WHR.csv")




#SSOE -> WC
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
X2 <- data.table::fread("shungin_wc.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "MarkerName"), 
                 beta_hat_cols = c("BETA", "b"), 
                 se_cols = c("SE", "se"), 
                 A1_cols = c("ALLELE1", "Allele1"), 
                 A2_cols = c("ALLELE0", "Allele2"),
                 pval_cols = c("P_BOLT_LMM_INF", "p")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SSOE_WC.csv")




#SSOE -> %TBF
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
X2 <- data.table::fread("lu_bodyfatpercent.txt", header=TRUE)

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "hm_rsid"), 
                 beta_hat_cols = c("BETA", "hm_beta"), 
                 se_cols = c("SE", "standard_error"), 
                 A1_cols = c("ALLELE1", "hm_effect_allele"), 
                 A2_cols = c("ALLELE0", "hm_other_allele"),
                 pval_cols = c("P_BOLT_LMM_INF", "p_value")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SSOE_BFP.csv")




#SSOE -> VAT
rm(X1); rm(X2); rm(X3); rm(X_clump); rm(res); rm(res2); rm(params); rm(varlist)

X1 <- data.table::fread("klimentidis_SSOE.txt", header=TRUE)
load("VAT_karlsson.rda")
X2 <- dat2

dim(X1); dim(X2); head(X1); head(X2); names(X1); names(X2)

X3 <- gwas_merge(X1,
                 X2,
                 snp_name_cols = c( "SNP", "SNP"), 
                 beta_hat_cols = c("BETA", "BETA"), 
                 se_cols = c("SE", "SE"), 
                 A1_cols = c("ALLELE1", "ALLELE1"), 
                 A2_cols = c("ALLELE0", "ALLELE0"),
                 pval_cols = c("P_BOLT_LMM_INF", "P_BOLT_LMM_INF")
)

head(X3)

set.seed(100)
varlist <- with(X3, sample(snp, size = 100000, replace=FALSE))
params <- est_cause_params(X3, varlist)

class(params)
names(params)
head(params)
params$rho
head(params$mix_grid)

r2_thresh <- 1e-2
pval_thresh <- 1e-3

X_clump <- X3 %>% rename(rsid = snp,
                         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile
  )

keep_snps <- X_clump$rsid

res <- cause(X = X3, variants = keep_snps, param_ests = params, force = T)
res2 <- summary(res, ci_size=0.95)

res2 

write.csv(res2$tab, "C:\\Users\\Person\\Desktop\\MR\\cause\\SELF_SSOE_VAT.csv")




