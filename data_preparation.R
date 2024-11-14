library(tidyverse)
library(data.table)
library(broom)
library(here)
library(huxtable)
library(haven)
library(ivreg)

###
### Data construction part
###

### Get MoBa IDs
pregids <- fread(here("~/data/durable/data/moba/linkage/PDB2601_kobling_SSB_v12.csv"))
parents <- fread(here("~/data/durable/data/moba/Original files/csv/PDB2601_SV_INFO_v12.csv"))

### Get basic demography data
fasteoppl <- fread("~/data/durable/data/registers/original/csv/w19_0634_faste_oppl_ut.csv") %>% 
  mutate(cohort = as.numeric(str_sub(foedsels_aar_mnd, 1,4))) %>% 
  mutate(sex = 2-kjoenn) %>% 
  rename(cob = fodeland) %>% 
  select(w19_0634_lnr, cohort, cob, sex, mor_lnr, far_lnr)

### Fathers and mothers fertiliy history, select parity of focal child.
paternal_birthorder <- select(fasteoppl, cohort, w19_0634_lnr, far_lnr) %>% 
  arrange(far_lnr, cohort) %>% 
  filter(far_lnr!="") %>% 
  group_by(far_lnr) %>% 
  mutate(birthorder = row_number()) %>% 
  ungroup %>% 
  select(child_w19_0634_lnr = w19_0634_lnr, paternal_birthorder = birthorder)

maternal_birthorder <- select(fasteoppl, cohort, w19_0634_lnr, mor_lnr) %>% 
  arrange(mor_lnr, cohort) %>% 
  filter(mor_lnr!="") %>% 
  group_by(mor_lnr) %>% 
  mutate(birthorder = row_number()) %>% 
  ungroup %>% 
  select(child_w19_0634_lnr = w19_0634_lnr, maternal_birthorder = birthorder)

### Outcome data: National std. tests
###
nasjprov <- fread("~/data/durable/data/registers/original/csv/w19_0634_nasjonale_prover_ut.csv") 
outcomes <- c("NPREG05","NPREG08","NPREG09","NPLES05","NPLES08","NPLES09","NPENG05","NPENG08")
nprover <-   select(nasjprov, w19_0634_lnr, PROVE, AARGANG, POENG, DELTATTSTATUS) %>% 
  filter(PROVE %in% c("NPENG08","NPREG08","NPLES08","NPENG05","NPREG09","NPLES05","NPREG05","NPLES09")) %>% 
  filter(DELTATTSTATUS=="D", POENG>0) %>% 
  group_by(w19_0634_lnr, PROVE) %>% 
  arrange(desc(AARGANG)) %>% 
  slice(1) %>% 
  group_by(PROVE, AARGANG) %>% 
  mutate(zresult = scale.default(POENG)) %>% 
  ungroup %>% 
  gather(info, value, -w19_0634_lnr, -PROVE) %>% 
  mutate(info = if_else(info=="zresult", "",info)) %>% 
  unite(col = varname, PROVE, info, sep="", remove=T) %>% 
  spread(varname, value) %>% 
  mutate_at(vars(starts_with("NP")), as.numeric) %>% 
  rename(child_w19_0634_lnr = w19_0634_lnr) %>% 
  select(child_w19_0634_lnr, one_of(outcomes))
nprover$avg5test <- rowMeans(select(nprover, NPREG05, NPLES05, NPENG05), na.rm=T)
nprover$num5test <- rowSums(!is.na(select(nprover, NPREG05, NPLES05, NPENG05)), na.rm=T)

education <- fread("/tsd/p805/data/durable/data/registers/original/csv/w19_0634_utd_1970_2018_ut.csv")
gp_education <- education %>% 
  select(w19_0634_lnr,BU_2000) %>% 
  rename(mother_far_lnr = w19_0634_lnr) %>% 
  mutate(grandparent_educ = str_sub(BU_2000,1,1)) %>% 
  select(mother_far_lnr, grandparent_educ) %>% 
  group_by(mother_far_lnr) %>% 
  slice(1) %>% ungroup


### Create main data set of families, with ID numbers, child outcomes,
### demographic variables and links to biobank data
pars <- filter(pregids, rolle!="SU2PT_CHILD") %>% 
  select(-barn_nr) %>% 
  spread(rolle, w19_0634_lnr) %>% 
  rename(father_w19_0634_lnr = SU2PT_FATHER, mother_w19_0634_lnr = SU2PT_MOTHER)

### MoBa children except co-twins (barn_nr==1), with links to their parents
all_families <- filter(pregids, rolle=="SU2PT_CHILD", barn_nr==1) %>% 
  select(-rolle) %>% rename(child_w19_0634_lnr = w19_0634_lnr) %>% 
  left_join(pars) %>%
  left_join(prefixit(fasteoppl, "child_")) %>% 
  left_join(prefixit(fasteoppl, "father_")) %>% 
  left_join(prefixit(fasteoppl, "mother_")) %>% 
  left_join(parents, by=c("PREG_ID_2601")) %>% 
  rename(moba_child_cohort = FAAR) %>% 
  left_join(paternal_birthorder) %>% 
  left_join(maternal_birthorder) %>% 
  left_join(nprover) %>% 
  left_join(gp_education, by="mother_far_lnr") %>% 
  mutate(mothers_age = moba_child_cohort-mother_cohort) %>% 
  mutate(fathers_age = moba_child_cohort-father_cohort)

parents_ed <- all_families %>% select(mother_w19_0634_lnr, father_w19_0634_lnr, moba_child_cohort) %>% 
  gather(parent, w19_0634_lnr, -moba_child_cohort) %>% 
  select(-parent) %>% 
  filter(w19_0634_lnr!="") %>% 
  group_by(w19_0634_lnr) %>% 
  slice(1) %>% 
  ungroup %>% 
  left_join(education) %>% 
  select(-contains("igang")) %>% 
  gather(eduyr, edulevel, contains("BU_")) %>% 
  mutate(yr = as.numeric(str_sub(string = eduyr,start=4))) %>% 
  filter(yr==moba_child_cohort) %>% 
  mutate(edlevel = floor(edulevel/100000)) %>% 
  mutate(eduyears = case_when(edlevel==1 ~ 6.0,
                              edlevel==2 ~ 9.0,
                              edlevel==3 ~ 10.0,
                              edlevel==4 ~ 12.0,
                              edlevel==5 ~ 14.0,
                              edlevel==6 ~ 16.0,
                              edlevel==7 ~ 18.0,
                              edlevel==8 ~ 21.0,
                              TRUE ~ NA)) %>%
  select(w19_0634_lnr, edulevel, eduyears) %>% 
  group_by(w19_0634_lnr) %>% 
  slice(1) %>% 
  ungroup
  
firstborns <- all_families %>%  
  filter(paternal_birthorder==1, maternal_birthorder==1) %>% 
  left_join(prefixit(parents_ed, "mother_")) %>% 
  left_join(prefixit(parents_ed, "father_"))
  
### Add genetic data
###
### First, generate genetic IDs 
geninfo <- fread("~/data/durable/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt") %>% 
  select(SENTRIXID, FID, IID, num_range(prefix = "PC",range = 1:10))
child_sentrix <- fread(here("~/data/durable/data/moba/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Child_20220829.csv")) %>% 
  select(PREG_ID_2601, BARN_NR, child_SENTRIXID = SENTRIX_ID) %>% 
  inner_join(prefixit(geninfo,"child_"))
father_sentrix <- fread(here("~/data/durable/data/moba/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Father_20220829.csv"))     %>% 
  select(F_ID_2601, father_SENTRIXID = SENTRIX_ID) %>% 
  inner_join(prefixit(geninfo,"father_"))
mother_sentrix <- fread(here("~/data/durable/data/moba/MoBaGenetics/key_CENTRIX_ID_Aug_22/PDB2601_MoBaGeneticsTot_Mother_20220829.csv")) %>% 
  select(M_ID_2601, mother_SENTRIXID = SENTRIX_ID) %>% 
  inner_join(prefixit(geninfo,"mother_"))

### Add genetic IDs to data  
with_genetics <- firstborns %>% 
  left_join(child_sentrix, by=c("PREG_ID_2601",barn_nr = "BARN_NR")) %>% 
  left_join(mother_sentrix, by=c("M_ID_2601")) %>% 
  left_join(father_sentrix, by=c("F_ID_2601"))

### Add polygenic scores to data
#
#
# Get scores from cluster
score_names <- c("age1contraception","EA3","age1sex","A1B_NHB", "ADHD","smokeinit")
myscale <- function(x) return(as.vector(scale.default(x)))

# Read in prepared individual-PGIs
ordinary_pgis <- fread("/tsd/p805/data/durable/projects/openflux/data/prepared-genetic-data-individuals.csv")
relevant_pgis <- with_genetics %>% 
  select(child_IID, mother_IID, father_IID) %>% 
  gather(person, IID) %>% 
  select(-person) %>% 
  distinct %>% 
  left_join(ordinary_pgis) %>% 
  group_by(IID) %>% 
  slice(1) %>% 
  ungroup %>% 
  select(IID, 
         pgs_A1B = AFBPooled_Mills_2021_PRSice_PGS, 
         pgs_ADHD = ADHD_Demontis_et_al_2023_PGS, 
         pgs_smokeinit = SmokingInitiation_PRSice_PGS,
         pgs_contraception = age1contraception_PRSice_PGS,
         pgs_EA3 = EA4_excl23andMe_mobaref_PGS,
         pgs_age1sex = age1sex_PRSice_PGS) %>% 
  mutate_at(vars(contains("pgs_")), myscale)

with_scores <- with_genetics %>% 
  filter(child_IID!="") %>% 
  filter(mother_IID!="") %>% 
  filter(father_IID!="") %>% 
  left_join(prefixit(relevant_pgis,"child_")) %>% 
  left_join(prefixit(relevant_pgis,"father_")) %>% 
  left_join(prefixit(relevant_pgis,"mother_"))

final_df <- with_scores %>% 
  rename(mother_age = mothers_age) %>% 
  rename(father_age = fathers_age) %>% 
  rename(grand_educ = grandparent_educ)



# ### Sib analysis
# grandparents <- complete_trios %>% 
#   select(ends_with("_lnr")) %>% 
#   filter(father_mor_lnr!="") %>% 
#   filter(father_far_lnr!="") %>% 
#   filter(mother_mor_lnr!="") %>% 
#   filter(mother_far_lnr!="") %>% 
#   gather(var,lnr, -child_w19_0634_lnr, -father_w19_0634_lnr, -mother_w19_0634_lnr) %>%
#   filter(!str_detect(var,"child")) %>% 
#   mutate(branch = str_sub(var,start=1,end=6)) %>% 
#   mutate(grandp = str_sub(var,start=8,end=15)) %>% 
#   select(-var) %>% group_by(child_w19_0634_lnr, branch) %>% 
#   spread(grandp,lnr) %>% 
#   unite(col="grandparents", sep=":",far_lnr, mor_lnr, remove=T)
# 
# pgidevs <- fasteoppl %>% select(w19_0634_lnr, mor_lnr,far_lnr) %>% 
#   left_join(select(ordinary_pgis, w19_0634_lnr, pgi = AFBPooled_Mills_2021_PRSice_PGS)) %>% 
#   mutate(pgi = myscale(pgi)) %>% 
#   filter(!is.na(pgi)) %>% 
#   filter(mor_lnr!="", far_lnr!="") %>% 
#   group_by(w19_0634_lnr) %>% 
#   slice(1) %>% 
#   ungroup %>% 
#   unite(col="grandparents", sep=":",far_lnr, mor_lnr, remove=T) %>% 
#   group_by(grandparents) %>% 
#   mutate(pgimean = mean(pgi), pgidev = pgi-pgimean) %>% 
#   filter(pgidev!=0)
# 
# withdevs <- with_scores %>% 
#   left_join(prefixit(pgidevs, "mother_")) %>% 
#   left_join(prefixit(pgidevs, "father_"))
# 
# # II. Variables for analysis
# # a. Outcome: Children’s test score
# # b. Exposure: Parental AFB
# # c. IV: Parental PGI for AFB (including full child’s PGI for AFB)
# # d. Control Variables: 
# #   i. First 10 PC
# # ii. Grandparental education
# # iii. Partners full PGI for AFB
# # iv. Parental educational attainment
# # v. Birth year of child
# # vi. Sex of child
# # vii. Focal parent’s PGI for age at first sexual intercourse
# # viii. Focal parent’s PGI for age at smoking initiation
# # ix. Focal parent’s PGI for age at first use of oral contraceptives
# # x. Focal parent’s PGI for attention deficit hyperactivity disorder
# # xi. Focal parent’s PGI for educational attainment
# 
#   
### Save data set
fwrite(final_df, here("revised-data.csv"))
#   
