library(tidyverse)
library(data.table)
library(broom)
library(here)
library(huxtable)
library(haven)
library(AER)
library(xlsx)
library(strip)
#source(here("../my-useful-R-functions.R"))

complete_trios <- fread(here("data/revised-data.csv"))
nrow(complete_trios)
names(complete_trios)

###
### ANALYSIS PART STARTS HERE: Definitions
###
df <- complete_trios

### Mother models

### Estimate models and store results
outcome <- "avg5test"
endogenous_var <- "parent_age"
instrument_var <- "parent_pgs_A1B"
popstrat_vars <- str_c("parent_PC",1:10)
grandparent = "grand_educ"
child_sex <- "child_sex"
child_cohort <- "child_cohort"
partner_pgi <- "partner_pgs_A1B"
child_pgi <- "child_pgs_A1B"
indep_var = "parent_age"
pareduc <- "parent_eduyears"
parent_pgi_a1s <- "parent_pgs_age1sex"
parent_pgi_smokeinit <- "parent_pgs_smokeinit"
parent_pgi_contraception <- "parent_pgs_contraception"
parent_pgi_ea3 <- "parent_pgs_EA3"
parent_pgi_adhd <- "parent_pgs_ADHD"

variable_sets <- list(
  "popstrat_vars",
  "grandparent",
  "partner_pgi",
  "pareduc",
  "child_pgi",
  "child_cohort",
  "child_sex",
  "parent_pgi_a1s",
  "parent_pgi_smokeinit",
  "parent_pgi_contraception",
  "parent_pgi_adhd",
  "parent_pgi_ea3"
) %>% setNames(nm=.) %>% map(get)

nsets <- length(variable_sets)

model_list <- expand.grid(V1 = c(F,T), V2 = c(F,T), V3 = c(F,T),
            V4 = c(F,T), V5 = c(F,T), V6 = c(F,T),
            V7 = c(F,T), V8 = c(F,T), V9 = c(F,T),
            V10 = c(F,T), V11 = c(F,T), V12 = c(F,T)) %>% setNames(nm=names(variable_sets))

gen_varlist <- function(row) {
  myvars <- row %>% t %>% as.vector
  varlist <- unlist(variable_sets[myvars])
  
  return(varlist)
}
varlist <- apply(model_list, 1, gen_varlist) 
varsummary <- apply(model_list, 1, function(row) paste(ifelse(row,"T","F"), collapse=""))

model_list$varlist <- varlist
model_list$varsummary <- varsummary
head(model_list)
class(model_list)
names(model_list)

ols_mother_models <- model_list %>% mutate(parent="mother",model_type="ols")
ols_father_models <- model_list %>% mutate(parent="father",model_type="ols")
iv_mother_models <- model_list %>% mutate(parent="mother",model_type="iv")
iv_father_models <- model_list %>% mutate(parent="father",model_type="iv")
all_models <- bind_rows(ols_mother_models, ols_father_models, 
                        iv_mother_models, iv_father_models)

generate_model <- function(row) {
  varlist <- row$varlist
  parent <- row$parent
  model_type <- row$model_type
  partner <- if_else(parent=="mother","father","mother")
  
  ols_stage <- str_c(c(str_c(parent,"_age", collapse=""),varlist), collapse=" + ")
  iv_stage <-  str_c(c(" | parent_pgs_A1B",varlist), collapse=" + ")
  iv_stage <- if_else(model_type=="iv", iv_stage, "")
  lhs <- str_c(ols_stage, iv_stage)
  
  lhs <- gsub(pattern="parent",replacement=parent, lhs)
  lhs <- gsub(pattern="partner",replacement=partner, lhs)
  m <- str_c(outcome, " ~ ", lhs) 
  
  return(m)
}

all_models$model_def <- apply(all_models, 1, generate_model)


estimate_model <- function(row, ...) {
  model_def <- row$model_def
  m <- NULL
  if (row$model_type=="iv") {
    m <- ivreg(data=df, formula=as.formula(model_def))
  } else {
    m <- lm(data=df, formula=as.formula(model_def))
  }

  estimate = NA; std.error = NA; p.value = NA
  tm <- broom::tidy(m) %>% 
    filter(str_detect(term,".*_age$"))
  if (nrow(tm)==1) {
    estimate = tm$estimate
    std.error = tm$std.error
    p.value= tm$p.value
  }
  newrow = list(estimate=estimate, std.error=std.error, p.value=p.value)
      
  return(newrow)
}

all_models$model_results <- apply(all_models, 1, estimate_model)
all_models <- all_models %>% unnest_wider(model_results)

saveRDS(all_models, file = here("output/model_results.RDS"))

### Collected estimates
all_models <- readRDS(file=here("output/model_results.RDS"))


### Selected models
selected_model_ids <- list("FFFFFFFFFFFF","TFFFFFFFFFFF","TTTTTTTTTTTF","TTTTTTTTTTTT")
selected_models <- all_models %>% filter(varsummary %in% selected_model_ids)
nrow(selected_models)

get_estimated_model <- function(parent, model_type, model_def, ...) {
    m <- NULL
    if (model_type=="iv") {
      m <- ivreg(data=df, formula=as.formula(model_def))
    } else {
      m <- lm(data=df, formula=as.formula(model_def))
    }
    
    return(m)
}

estimated_selected_models <- selected_models %>% pmap(get_estimated_model)
selected_names <- selected_models %>% unite(parent,model_type,varsummary, col=model, sep=":") %>% pull(model)
names(estimated_selected_models) <- selected_names
coefdf <- map(estimated_selected_models, ~tidy(.x,conf.int=TRUE)) %>% bind_rows(.id="model")
statdf <- map(estimated_selected_models, performance::model_performance) %>% map(as.data.frame) %>% bind_rows(.id="model")
summary(estimated_selected_models$`mother:iv:TTTTTTTTTTTT`)


table1 <- coefdf %>% 
  filter(term=="mother_age"|term=="father_age") %>% 
  select(-p.value,-statistic,-std.error) %>% 
  left_join(select(statdf, model, BIC, F_test = weak_instruments)) %>% 
  separate(model, into = c("parent","model_type","spec"), remove=T) %>% 
  arrange(desc(parent),desc(model_type),spec) %>% 
  select(-term)

write.xlsx(file = here("output/table1.xlsx"), x = table1, sheetName = "Table 1")

    nest(results = estimate, conf.low, conf.high, BIC, F_test)

  gather(key,ar, estimate, conf.low, conf.high, BIC, F_test)


### First-stage regressions
generate_first_stage <- function(row) {
    varlist <- row$varlist
    parent <- row$parent
    model_type <- row$model_type
    partner <- if_else(parent=="mother","father","mother")
    
    first_stage <- str_c(c("parent_pgs_A1B",varlist), collapse=" + ")
    lhs <- str_c(str_c(parent,"_age", collapse="")," ~ ",first_stage)
    
    lhs <- gsub(pattern="parent",replacement=parent, lhs)
    lhs <- gsub(pattern="partner",replacement=partner, lhs)

    return(lhs)
}

first_stage_models <- selected_models %>% filter(model_type=="ols") %>% unite(col="model", varsummary, parent, sep=":")
first_stage_models <- apply(filter(selected_models, model_type=="ols"), 1, generate_first_stage) %>% setNames(first_stage_models$model)

estimate_first_stage <- function(row, ...) {
  model_def <- row
  m <- lm(data=df, formula=as.formula(model_def))
  return(m)
}

first_stages <- map(first_stage_models, estimate_first_stage)

first_coefdf <- map(first_stages, ~tidy(.x,conf.int=TRUE)) %>% bind_rows(.id="model")
first_statdf <- map(first_stages, performance::model_performance) %>% map(as.data.frame) %>% bind_rows(.id="model")


table2a <- first_coefdf %>% 
  separate(model, into = c("spec","parent"), remove=F) %>% 
  filter(term=="mother_pgs_A1B"|term=="father_pgs_A1B") %>% 
  filter(str_detect(term, parent)) %>% 
  select(-p.value,-statistic,-std.error) %>% 
  left_join(select(first_statdf, model, BIC)) %>% 
  arrange(desc(parent),spec) %>% 
  select(-term)

write.xlsx(file = here("output/appendix_tableA.xlsx"), x = table2a, sheetName = "Appendix Table A")

### 
generate_reduced <- function(row) {
  varlist <- row$varlist
  parent <- row$parent
  model_type <- row$model_type
  partner <- if_else(parent=="mother","father","mother")
  
  reduced <- str_c(c("parent_pgs_A1B",varlist), collapse=" + ")
  lhs <- str_c("avg5test ~ ",reduced)
  
  lhs <- gsub(pattern="parent",replacement=parent, lhs)
  lhs <- gsub(pattern="partner",replacement=partner, lhs)
  
  return(lhs)
}

reduced_models <- selected_models %>% filter(model_type=="ols") %>% unite(col="model", varsummary, parent, sep=":")
reduced_models <- apply(filter(selected_models, model_type=="ols"), 1, generate_reduced) %>% setNames(reduced_models$model)

estimate_reduced <- function(row, ...) {
  model_def <- row
  m <- lm(data=df, formula=as.formula(model_def))
  return(m)
}

reduceds <- map(reduced_models, estimate_reduced)

reduced_coefdf <- map(reduceds, ~tidy(.x,conf.int=TRUE)) %>% bind_rows(.id="model")
reduced_statdf <- map(reduceds, performance::model_performance) %>% map(as.data.frame) %>% bind_rows(.id="model")


table2b <- reduced_coefdf %>% 
  separate(model, into = c("spec","parent"), remove=F) %>% 
  filter(term=="mother_pgs_A1B"|term=="father_pgs_A1B") %>% 
  filter(str_detect(term, parent)) %>% 
  select(-p.value,-statistic,-std.error) %>% 
  left_join(select(reduced_statdf, model, BIC)) %>% 
  arrange(desc(parent),spec) %>% 
  select(-term)

write.xlsx(file = here("output/appendix_tableB.xlsx"), x = table2b, sheetName = "Appendix Table B.")


#### Corr matrix
cmatrix <- select(df, -contains("lnr"), -contains("_FID"), -contains("_IID"), -contains("_PC"), -
                    contains("SENTRIX"), -contains("2601"), -ends_with("_cob"), -
                    starts_with("NP"), -father_sex, -child_sex, -mother_sex, -contains("cohort"), 
                  -contains("birthorder"), -contains("edulevel"), -barn_nr,
                  -contains("eduyears"),-grand_educ, -avg5test, -num5test)
thematrix <- cor(cmatrix) %>% round(digits = 2)
write.xlsx2(thematrix, file = here("output/corr-matrix.xlsx"))

